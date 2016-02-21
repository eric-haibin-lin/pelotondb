//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// BWTree.h
//
// Identification: src/backend/index/BWTree.h
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once
#include "backend/common/types.h"
#include <map>
#include <vector>
#include <climits>
#include <string.h>

namespace peloton {
namespace index {

//===--------------------------------------------------------------------===//
// Types and Enums
//===--------------------------------------------------------------------===//
typedef uint64_t LPID;

enum BWTreeNodeType {
  TYPE_LPAGE = 0,
  TYPE_IPAGE = 1,
  TYPE_OTHER = 2,
};

#define IPAGE_ARITY 5
#define LPAGE_ARITY 5

#define INVALID_LPID ULLONG_MAX

#define DELTA_CHAIN_LIMIT 5

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode;

template <typename KeyType, typename ValueType, class KeyComparator>
class IPage;

template <typename KeyType, typename ValueType, class KeyComparator>
class LPage;

//===--------------------------------------------------------------------===//
// NodeStateBuilder
//===--------------------------------------------------------------------===//
// TODO add access methods for LNode scan
// TODO add methods for use by split deltas, etc
// TODO add constructors to LNode and INode to build node based on state
// TODO add methods on each node to build the NodeState
/**
 * Builder for a node state, to be used with Delta Compression and
 * scans on indexes with multiple keys
 */
template <typename KeyType, typename ValueType, class KeyComparator>
class NodeStateBuilder {
 protected:
  // number of k-v pairs. everything after size is ignored
  oid_t size;

  // pointer to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  bool is_separated = false;

  KeyType separator_key;
  LPID split_new_page_id = 0;

 public:
  NodeStateBuilder(oid_t size, BWTree<KeyType, ValueType, KeyComparator> *map)
      : size(size), map(map){};

  virtual BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage() = 0;

  inline bool IsSeparated() { return is_separated; }

  inline KeyType GetSeparatorKey() { return separator_key; }

  inline LPID GetSplitNewPageId() { return split_new_page_id; }

  virtual ~NodeStateBuilder() { ; };

  friend class LPage<KeyType, ValueType, KeyComparator>;
  friend class IPage<KeyType, ValueType, KeyComparator>;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class INodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // IPage children nodes
  std::pair<KeyType, LPID> *children_ = nullptr;
  IPage<KeyType, ValueType, KeyComparator> *new_page_ = nullptr;

 public:
  // IPage constructor
  INodeStateBuilder(std::pair<KeyType, LPID> *children, int children_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(children_len, map) {
    children_ = new std::pair<KeyType, LPID>[IPAGE_ARITY + DELTA_CHAIN_LIMIT]();
    for (int i = 0; i < children_len; i++) {
      children_[i] = children[i];
    }
  };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage();

  ~INodeStateBuilder(){
      // TODO delete arrays
  };

  //***************************************************
  // IPage Methods
  //***************************************************
  void AddChild(std::pair<KeyType, LPID> &new_pair);

  void RemoveChild(KeyType &key_to_remove);

  void SeparateFromKey(KeyType separator_key, LPID split_new_page_id);

  friend class LPage<KeyType, ValueType, KeyComparator>;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class LNodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // LPage members
  LPID left_sibling_ = 0;
  LPID right_sibling_ = 0;
  LPage<KeyType, ValueType, KeyComparator> *new_page_ = nullptr;
  ValueType separator_location_;

  // LPage members
  std::pair<KeyType, ValueType> *locations_ = nullptr;

 public:
  // LPage constructor
  LNodeStateBuilder(LPID left_sibling, LPID right_sibling,
                    std::pair<KeyType, ValueType> *locations, int location_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(location_len, map),
        left_sibling_(left_sibling),
        right_sibling_(right_sibling) {
    locations_ =
        new std::pair<KeyType, ValueType>[IPAGE_ARITY + DELTA_CHAIN_LIMIT]();
    for (int i = 0; i < location_len; i++) {
      locations_[i] = locations[i];
    }
  };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage();

  ~LNodeStateBuilder(){
      // TODO delete arrays
      // should the builder delete the page?
  };

  //***************************************************
  // LPage Methods
  //***************************************************

  void UpdateLeftSib(LPID new_left_sib) { left_sibling_ = new_left_sib; }

  void UpdateRightSib(LPID new_right_sib) { right_sibling_ = new_right_sib; }

  void AddLeafData(std::pair<KeyType, ValueType> &new_pair);

  // only called when keys are unique
  void RemoveLeafData(KeyType &key_to_remove);

  // only called when keys are non unique
  void RemoveLeafData(std::pair<KeyType, ValueType> &entry_to_remove);

  void SeparateFromKey(KeyType separator_key, ValueType location,
                       LPID split_new_page_id);

 private:
  bool ItemPointerEquals(ValueType v1, ValueType v2);

  friend class LPage<KeyType, ValueType, KeyComparator>;
};

//===--------------------------------------------------------------------===//
// BWTree
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree {
 private:
  // the logical page id for the root node
  LPID root_;

  KeyComparator comparator;

  // the mapping table
  unsigned int mapping_table_cap_ = 128;
  unsigned int mapping_table_size_ = 0;
  LPID *free_LPIDs = new LPID[mapping_table_cap_];
  int free_LPID_index = 0;
  //  BWTreeNode<KeyType, ValueType, KeyComparator> ** mapping_table_ =
  //		  (BWTreeNode<KeyType, ValueType, KeyComparator> **)
  //		  malloc(sizeof(BWTreeNode<KeyType, ValueType, KeyComparator>
  //*)*mapping_table_cap_);

  BWTreeNode<KeyType, ValueType, KeyComparator> **mapping_table_;
  LPID next_LPID_ = 0;
  int current_readers = 0;
  int current_writers = 0;

 public:
  /*
   * On construction, BWTree will create a IPage and an empty LPage.
   * The IPage has only one pointer, which points to the empty LPage.
   * The IPage serves as the root of all other nodes.
   */
  BWTree(bool unique_keys, KeyComparator comparator)
      : comparator(comparator), unique_keys(unique_keys) {
    // this->unique_keys = unique_keys;
    // BWTreeNode<KeyType, ValueType, KeyComparator>::comparator = comparator;
    mapping_table_ =
        new BWTreeNode<KeyType, ValueType, KeyComparator> *[mapping_table_cap_];

    // TODO @abj initialize the root IPage (and maybe a LPage?)

    IPage<KeyType, ValueType, KeyComparator> *root = new IPage<KeyType, ValueType, KeyComparator>(this);

    root_ = InstallPage(root);

    LPage<KeyType, ValueType, KeyComparator> *first_lpage = new LPage<KeyType, ValueType, KeyComparator>(this);

    std::pair<KeyType, LPID> first_lpage_pair;

    LPID first_lpage_lpid;

    first_lpage_lpid = InstallPage(first_lpage);
    first_lpage_pair.second = first_lpage_lpid;

    root->GetChildren()[0] = first_lpage_pair;

    // with the given comparator
  };

  ~BWTree() {
    delete[] free_LPIDs;
    delete[] mapping_table_;
  }

  // get the index of the first occurrence of the given key
  template <typename PairSecond>
  inline int BinarySearch(__attribute__((unused)) KeyType key,
                          __attribute__((unused))
                          std::pair<KeyType, PairSecond> *locations,
                          __attribute__((unused)) oid_t len);

  bool InsertEntry(KeyType key, ValueType location);

  bool DeleteEntry(KeyType key, ValueType location) {
    LPID child_lpid;
    child_lpid = root_;

    return GetNode(child_lpid)->DeleteEntry(key, location, child_lpid);
  };

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);
  std::vector<ValueType> ScanAllKeys();
  std::vector<ValueType> ScanKey(KeyType key);

 private:
  inline void AquireRead() {
    while (true) {
      while (current_writers == 1)
        ;
      __sync_add_and_fetch(&current_readers, 1);
      if (current_writers == 0)
        break;
      else
        __sync_add_and_fetch(&current_readers, -1);
    }
  }
  inline void ReleaseRead() { __sync_add_and_fetch(&current_readers, -1); }
  inline void AquireWrite() {
    while (__sync_bool_compare_and_swap(&current_writers, 0, 1))
      ;
    while (current_readers > 0)
      ;
  }
  inline void ReleaseWrite() {
    assert(__sync_bool_compare_and_swap(&current_writers, 1, 0));
  }

 public:
  // whether unique key is required
  bool unique_keys;

  // static bool unique_keys;
  // return 0 if the page install is not successful

  LPID InstallPage(BWTreeNode<KeyType, ValueType, KeyComparator> *node) {
    LPID newLPID = __sync_fetch_and_add(&next_LPID_, 1);
    // table grew too large, expand it
    while (newLPID >= mapping_table_cap_) {
      // only one thread should expand the table
      AquireWrite();
      if (newLPID < mapping_table_cap_) {
        ReleaseWrite();
        break;
      }
      int new_mapping_table_cap = mapping_table_cap_ * 2;
      auto new_mapping_table =
          new BWTreeNode<KeyType, ValueType,
                         KeyComparator> *[new_mapping_table_cap];
      memcpy(new_mapping_table, mapping_table_, mapping_table_cap_);
      delete[] mapping_table_;
      mapping_table_ = new_mapping_table;
      mapping_table_cap_ = new_mapping_table_cap;
      ReleaseWrite();
    }
    AquireRead();
    mapping_table_[newLPID] = node;
    ReleaseRead();
    return newLPID;
  }

  bool SwapNode(LPID id, BWTreeNode<KeyType, ValueType, KeyComparator> *oldNode,
                BWTreeNode<KeyType, ValueType, KeyComparator> *newNode) {
    AquireRead();
    bool ret =
        __sync_bool_compare_and_swap(mapping_table_ + id, oldNode, newNode);
    ReleaseRead();
    return ret;
  }

  // assumes that LPID is valid
  BWTreeNode<KeyType, ValueType, KeyComparator> *GetNode(LPID id) {
    AquireRead();
    auto ret = mapping_table_[id];
    ReleaseRead();
    return ret;
  }

  // return 0 if equal, -1 if left < right, 1 otherwise
  int CompareKey(KeyType left, KeyType right) {
    bool less_than_right = comparator(left, right);
    bool greater_than_right = comparator(right, left);
    if (!less_than_right && !greater_than_right) {
      return 0;
    } else if (less_than_right) {
      return -1;
    }
    return 1;
  }

  // compress delta chain
  bool CompressDeltaChain(LPID page_to_compress) {
    auto old_node_ptr = GetNode(page_to_compress);
    auto old_node_state = old_node_ptr->BuildNodeState();
    auto new_node_ptr = old_node_state->GetPage();

    // delete the temporary state
    delete old_node_state;

    bool completed = SwapNode(page_to_compress, old_node_ptr, new_node_ptr);
    // if we didn't failed to install we should clean up the page we created
    if (completed) {
      // TODO: add old node to epoch
    } else {
      // if we failed we should clean up the new page
      delete new_node_ptr;
    }

    return completed;
  }
};

//===--------------------------------------------------------------------===//
// BWTreeNode
// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode {
 public:
  BWTreeNode(BWTree<KeyType, ValueType, KeyComparator> *map,
             int delta_chain_len)
      : map(map), delta_chain_len_(delta_chain_len){};

  // These are virtual methods which child classes have to implement.
  // They also have to be redeclared in the child classes
  virtual bool InsertEntry(KeyType key, ValueType location, LPID self) = 0;

  virtual bool DeleteEntry(KeyType key, ValueType location, LPID self) = 0;

  virtual std::vector<ValueType> Scan(
      const std::vector<Value> &values,
      const std::vector<oid_t> &key_column_ids,
      const std::vector<ExpressionType> &expr_types,
      const ScanDirectionType &scan_direction) = 0;

  virtual std::vector<ValueType> ScanAllKeys() = 0;

  virtual std::vector<ValueType> ScanKey(KeyType key) = 0;

  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *
  BuildNodeState() = 0;

  // for scanKey, we have to build the node state as well. But we only care
  // about
  // the keys we want to scan
  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      __attribute__((unused)) KeyType key) {
    return nullptr;
  };

  //    virtual NodeStateBuilder<KeyType, ValueType, KeyComparator>
  //    *BuildScanState(
  //        __attribute__((unused)) const std::vector<Value> &values,
  //        __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
  //        __attribute__((unused)) const std::vector<ExpressionType>
  //        &expr_types,
  //        __attribute__((unused)) const ScanDirectionType &scan_direction) {
  //      return nullptr;
  //    };

  // Each sub-class will have to implement this function to return their type
  virtual BWTreeNodeType GetTreeNodeType() const = 0;

  virtual ~BWTreeNode(){

  };

  int GetDeltaChainLen() { return delta_chain_len_; }

 protected:
  // the handler to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  // length of delta chain
  int delta_chain_len_;
};

//===--------------------------------------------------------------------===//
// IPage
// The IPage hold pointers to all its children. An IPage with n + 1 keys
// k1, k2, .... kn has (n + 1) pointers to its children, where
// p1 for (-infinity, k1), p2 for [k1, k2), p3 for [k2, k3) ...
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  IPage(BWTree<KeyType, ValueType, KeyComparator> *map)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    size_ = 0;
    // children_ = new std::pair<KeyType, LPID>();
  };

  ~IPage(){};

  bool InsertEntry(KeyType key, ValueType location, LPID self);

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    LPID child_lpid = GetChild(key, children_, size_);
    return this->map->GetNode(child_lpid)
        ->DeleteEntry(key, location, child_lpid);
  };

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  // get the index of the child at next level, which contains the given key
  static int GetChild(__attribute__((unused)) KeyType key,
                      __attribute__((unused))
                      std::pair<KeyType, LPID> *children,
                      __attribute__((unused)) oid_t len) {
    return 0;
    // TODO @Matt implement this
  };

  std::pair<KeyType, LPID> *GetChildren() { return children_; }

 private:
  std::pair<KeyType, LPID> children_[IPAGE_ARITY];

  oid_t size_;
};

//===--------------------------------------------------------------------===//
// Delta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class Delta : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  Delta(BWTree<KeyType, ValueType, KeyComparator> *map,
        BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(
            map, modified_node->GetDeltaChainLen() + 1),
        modified_node(modified_node){};

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

 protected:
  // the modified node could either be a LPage or IPage or Delta
  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node;
};

//===--------------------------------------------------------------------===//
// IPageSplitDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageSplitDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  IPageSplitDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType key, LPID value)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node),
        modified_key_(key),
        modified_val_(value){};

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    // TODO implement this
    return false;
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location) {
    // TODO implement this
    return false;
  };

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

 private:
  // This key excluded in left child, included in the right child
  KeyType modified_key_;

  // The LPID of the new LPage
  LPID modified_val_;
};
//===--------------------------------------------------------------------===//
// LPageSplitDelta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageSplitDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  LPageSplitDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                  KeyType splitterKey, ValueType splitterVal,
                  LPID rightSplitPage)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node),
        modified_key_(splitterKey),
        modified_key_location_(splitterVal),
        modified_val_(rightSplitPage){};

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    // TODO implement this
    return false;
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location) {
    // TODO implement this
    return false;
  };

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

 private:
  // This key excluded in left child, included in the right child
  KeyType modified_key_;
  ValueType modified_key_location_;

  // The LPID of the new LPage
  LPID modified_val_;
};

//===--------------------------------------------------------------------===//
// LPageUpdateDelta
// LPageUpdateDelta either represents a insert delta or delete delta on LPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageUpdateDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  // TODO @abj initialize "is_delete_" to the desired value as well
  LPageUpdateDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                   BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                   KeyType key, ValueType value)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node),
        modified_key_(key),
        modified_val_(value){};

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    // TODO implement this
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta = new
    		LPageUpdateDelta<KeyType, ValueType, KeyComparator> (
        this->map, this, key, location);

    return this->map->SwapNode(self, this, new_delta);
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location, LPID self) {
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta = new
    		LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
        this->map, this, key, location);

    new_delta->SetDeleteFlag();
    return this->map->SwapNode(self, this, new_delta);
  };

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      __attribute__((unused)) KeyType key) {
    return nullptr;
  };

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      __attribute__((unused)) const std::vector<Value> &values,
      __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
      __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
      __attribute__((unused)) const ScanDirectionType &scan_direction) {
    return nullptr;
  };

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  void SetKey(KeyType modified_key) { modified_key_ = modified_key; }

  void SetValue(ValueType modified_val) { modified_val_ = modified_val; }

  void SetDeleteFlag() { is_delete_ = true; };

 private:
  // The key which is modified
  KeyType modified_key_;

  // The item pointer of the updated key
  ValueType modified_val_;

  // Whether it's a delete delta
  bool is_delete_ = false;
};

//===--------------------------------------------------------------------===//
// IPageUpdateDelta
//
// IPageUpdateDelta either represents a insert delta or delete delta on IPage
// The insert delta is essentially a separator in the case of a split SMO,
// where the new LPID has to be inserted. The delete delta is created in the
// case of a merge SMO, where an old LPID has to be removed.
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPageUpdateDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  // TODO @abj initialize "is_delete_" to the desired value as well

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    // TODO implement this
    return false;
  };
  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location) {
    // TODO implement this
    return false;
  };

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

 private:
  // The key which is modified
  KeyType modified_key_;

  // The logical page id of the updated key. Set to 0 for delete delta
  LPID modified_id_;

  // Whether it's a delete delta
  bool is_delete_ = false;
};

// TODO More delta classes such as
// merge, split, remove_page, separator

//===--------------------------------------------------------------------===//
// LPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  // get the index of the first occurrence of the given <k, v> pair
  // used when deleting a <k-v> entry from non-unique keys
  static int BinarySearch(__attribute__((unused))
                          std::pair<KeyType, ValueType> pair,
                          __attribute__((unused))
                          std::pair<KeyType, ValueType> *locations,
                          __attribute__((unused)) oid_t len);

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    // TODO initialize these with the proper values
    left_sib_ = 0;
    right_sib_ = 0;
    // locations_ = new std::pair<KeyType, ValueType>[LPAGE_ARITY]();
    size_ = 0;
  };

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map,
        NodeStateBuilder<KeyType, ValueType, KeyComparator> *state)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    size_ = state->size;
    LNodeStateBuilder<KeyType, ValueType, KeyComparator> *lstate =
        reinterpret_cast<
            LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(state);
    for (oid_t index = 0; index < size_; index++) {
      locations_[index] = lstate->locations_[index];
    }
    left_sib_ = 0;
    right_sib_ = 0;
    // TODO handle split case
  };

  ~LPage(){};

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
    		new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
        this->map, this, key, location);

    return this->map->SwapNode(self, this, new_delta);
    // return false;
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    // TODO implement this

    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
    		new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
        this->map, this, key, location);

    new_delta->SetDeleteFlag();
    return this->map->SwapNode(self, this, new_delta);
  };

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  void MergeNodes(LPID self, LPID parent) {
    LPID newLpageLPID;
    int j = 0;
    KeyType splitterKey;
    ValueType splitterVal;
    bool swapSuccess;

    LPage<KeyType, ValueType, KeyComparator> newLpage(this->map, 0);

    for (int i = size_ / 2; i < size_; i++) {
      newLpage.locations_[j++] = locations_[i];
    }

    // TODO Why am I able to access this size_ private variable?
    newLpage.size_ = j;

    // TODO left_sib is set to self
    newLpage.right_sib_ = right_sib_;

    splitterKey = locations_[size_ / 2 + 1].first;
    splitterVal = locations_[size_ / 2 + 1].second;

    newLpageLPID = this->map->InstallPage(&newLpage);

    LPageSplitDelta<KeyType, ValueType, KeyComparator> splitDelta(
        this->map, this, splitterKey, splitterVal, newLpageLPID);

    swapSuccess = this->map->SwapNode(self, this, &splitDelta);

    if (swapSuccess == false) {
      // What should we do on failure? This means that someone else succeeded in
      // doing the
      // atomic half split. Now if try and install our own Insert / Delete /
      // Update delta, it
      // will automatically be created on top of the LPageSplitDelta
      return;
    }

    // This completes the atomic half split
    // At this point no one else can succeed with the complete split because
    // this guy won in the half split
    // Now we still have to update the size field and the right_sib of this
    // node... how can we do it atomically?
    // Edit: No need to do that! Because the consolidation will do that lol you
    // dumbo abj

    // Now start with the second half
  };

 private:
  // return a vector of indices of the matched slots
  std::vector<oid_t> ScanKeyInternal(KeyType key);

  // return true if a given ItemPointer is invalid
  //  bool IsInvalidItemPointer(ValueType val);

  // left_sibling pointer is used for reverse scan
  LPID left_sib_;
  LPID right_sib_;

  // the key-value pairs
  std::pair<KeyType, ValueType> locations_[LPAGE_ARITY];

  // the size of stored locations
  oid_t size_;
};

}  // End index namespace
}  // End peloton namespace
