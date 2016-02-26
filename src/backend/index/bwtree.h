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
#include "backend/common/logger.h"
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

#define IPAGE_ARITY 256
#define LPAGE_ARITY 256

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
  std::pair<KeyType, LPID> children_[IPAGE_ARITY + DELTA_CHAIN_LIMIT];
  IPage<KeyType, ValueType, KeyComparator> *new_page_ = nullptr;

 public:
  // IPage constructor
  INodeStateBuilder(std::pair<KeyType, LPID> *children, int children_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(children_len, map) {
    memcpy(children_, children,
           sizeof(std::pair<KeyType, LPID>) * children_len);
  };

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage();

  ~INodeStateBuilder(){
      // TODO delete arrays
  };

  //***************************************************
  // IPage Methods
  //***************************************************

  /* AddChild adds a new <key, LPID> pair to its children if key is not
   * found. Otherwise overwrite the pair with the same key */
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
  LPID left_sibling_ = INVALID_LPID;
  LPID right_sibling_ = INVALID_LPID;
  LPage<KeyType, ValueType, KeyComparator> *new_page_ = nullptr;
  ValueType separator_location_;

  // LPage members
  std::pair<KeyType, ValueType> locations_[IPAGE_ARITY + DELTA_CHAIN_LIMIT];

 public:
  // LPage constructor
  LNodeStateBuilder(LPID left_sibling, LPID right_sibling,
                    std::pair<KeyType, ValueType> *locations, int location_len,
                    BWTree<KeyType, ValueType, KeyComparator> *map)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(location_len, map),
        left_sibling_(left_sibling),
        right_sibling_(right_sibling) {
    memcpy(locations_, locations,
           sizeof(std::pair<KeyType, ValueType>) * location_len);
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

  void RemoveLeafData(std::pair<KeyType, ValueType> &entry_to_remove);

  void SeparateFromKey(KeyType separator_key, ValueType location,
                       LPID split_new_page_id);

 private:
  bool ItemPointerEquals(ValueType v1, ValueType v2);

  friend class LPage<KeyType, ValueType, KeyComparator>;
};

//===--------------------------------------------------------------------===//
// ReadWriteLatch
//===--------------------------------------------------------------------===//
class ReadWriteLatch {
 public:
  inline void AquireRead() {
    LOG_INFO("Aquiring read Lock");
    while (true) {
      while (current_writers_ == 1)
        ;
      __sync_add_and_fetch(&current_readers_, 1);
      if (current_writers_ == 0)
        break;
      else
        __sync_add_and_fetch(&current_readers_, -1);
    }
  }
  inline void ReleaseRead() {
    LOG_INFO("Releasing read Lock");
    __sync_add_and_fetch(&current_readers_, -1);
  }
  inline void AquireWrite() {
    LOG_INFO("Aquiring write Lock");
    while (!__sync_bool_compare_and_swap(&current_writers_, 0, 1))
      ;
    while (current_readers_ > 0)
      ;
  }
  inline void ReleaseWrite() {
    LOG_INFO("Releasing write Lock");
    assert(__sync_bool_compare_and_swap(&current_writers_, 1, 0));
  }

 private:
  int current_readers_ = 0;
  int current_writers_ = 0;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class MappingTable {
 private:
  // the mapping table
  unsigned int mapping_table_cap_ = 128;
  unsigned int mapping_table_size_ = 0;

  BWTreeNode<KeyType, ValueType, KeyComparator> **mapping_table_;
  LPID next_LPID_ = 0;
  ReadWriteLatch latch_;

 public:
  MappingTable() {
    LOG_INFO("Constructing Mapping Table with initial capacity %d",
             mapping_table_cap_);
    mapping_table_ =
        new BWTreeNode<KeyType, ValueType, KeyComparator> *[mapping_table_cap_];
  };
  ~MappingTable() { delete[] mapping_table_; }
  LPID InstallPage(BWTreeNode<KeyType, ValueType, KeyComparator> *node) {
    LOG_INFO("Installing page in mapping table");
    LPID newLPID = __sync_fetch_and_add(&next_LPID_, 1);
    // table grew too large, expand it
    while (newLPID >= mapping_table_cap_) {
      LOG_INFO("mapping table has grown too large");
      // only one thread should expand the table
      latch_.AquireWrite();
      if (newLPID < mapping_table_cap_) {
        latch_.ReleaseWrite();
        break;
      }

      int new_mapping_table_cap = mapping_table_cap_ * 2;
      LOG_INFO("doubling size of mapping table capacity from %d to %d",
               mapping_table_cap_, new_mapping_table_cap);
      auto new_mapping_table =
          new BWTreeNode<KeyType, ValueType,
                         KeyComparator> *[new_mapping_table_cap];
      memcpy(new_mapping_table, mapping_table_,
             mapping_table_cap_ *
                 sizeof(new BWTreeNode<KeyType, ValueType, KeyComparator> *));
      delete[] mapping_table_;
      mapping_table_ = new_mapping_table;
      mapping_table_cap_ = new_mapping_table_cap;
      latch_.ReleaseWrite();
    }
    latch_.AquireRead();
    LOG_INFO("adding LPID: %lu to mapping table", newLPID);
    mapping_table_[newLPID] = node;
    latch_.ReleaseRead();
    return newLPID;
  }

  bool SwapNode(LPID id, BWTreeNode<KeyType, ValueType, KeyComparator> *oldNode,
                BWTreeNode<KeyType, ValueType, KeyComparator> *newNode) {
    LOG_INFO("swapping node for LPID: %lu into mapping table", id);
    latch_.AquireRead();
    bool ret =
        __sync_bool_compare_and_swap(mapping_table_ + id, oldNode, newNode);
    latch_.ReleaseRead();
    return ret;
  }

  // assumes that LPID is valid
  BWTreeNode<KeyType, ValueType, KeyComparator> *GetNode(LPID id) {
    LOG_INFO("getting node for LPID: %lu frommapping table", id);
    latch_.AquireRead();
    auto ret = mapping_table_[id];
    latch_.ReleaseRead();
    return ret;
  }
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
  MappingTable<KeyType, ValueType, KeyComparator> *mapping_table_;

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
    LOG_INFO("Inside BWTree Constructor");
    mapping_table_ = new MappingTable<KeyType, ValueType, KeyComparator>();

    // TODO @abj initialize the root IPage (and maybe a LPage?)

    IPage<KeyType, ValueType, KeyComparator> *root =
        new IPage<KeyType, ValueType, KeyComparator>(this);

    root_ = GetMappingTable()->InstallPage(root);

    LOG_INFO("Root got LPID: %lu", root_);

    LPage<KeyType, ValueType, KeyComparator> *first_lpage =
        new LPage<KeyType, ValueType, KeyComparator>(this);

    std::pair<KeyType, LPID> first_lpage_pair;

    LPID first_lpage_lpid;

    first_lpage_lpid = GetMappingTable()->InstallPage(first_lpage);

    LOG_INFO("The first LPage got LPID: %d", (int)first_lpage_lpid);
    first_lpage_pair.second = first_lpage_lpid;

    root->GetChildren()[0] = first_lpage_pair;

    LOG_INFO("Leaving BWTree constructor");

    // with the given comparator
  };

  ~BWTree() { delete mapping_table_; }

  // get the index of the first occurrence of the given key
  template <typename PairSecond>

  // positive index indicates found, negative indicates not found. 0 could be
  // either case
  inline int BinarySearch(KeyType key,

                          std::pair<KeyType, PairSecond> *locations, oid_t len);

  bool InsertEntry(KeyType key, ValueType location);

  bool DeleteEntry(KeyType key, ValueType location) {
    LPID child_lpid;
    child_lpid = root_;

    while (!GetMappingTable()
                ->GetNode(child_lpid)
                ->DeleteEntry(key, location, child_lpid))
      ;
    return true;
  };

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);
  std::vector<ValueType> ScanAllKeys();
  std::vector<ValueType> ScanKey(KeyType key);

 public:
  // whether unique key is required
  bool unique_keys;

  // static bool unique_keys;
  // return 0 if the page install is not successful

  inline MappingTable<KeyType, ValueType, KeyComparator> *GetMappingTable() {
    return mapping_table_;
  }

  // return 0 if equal, -1 if left < right, 1 otherwise
  inline int CompareKey(KeyType left, KeyType right) {
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
    LOG_INFO("Compressing delta chain for LPID: %lu", page_to_compress);
    auto old_node_ptr = GetMappingTable()->GetNode(page_to_compress);
    auto old_node_state = old_node_ptr->BuildNodeState();
    auto new_node_ptr = old_node_state->GetPage();

    // delete the temporary state
    delete old_node_state;

    bool completed = GetMappingTable()->SwapNode(page_to_compress, old_node_ptr,
                                                 new_node_ptr);
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

  // TODO to be implemented in the future
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

  inline int GetDeltaChainLen() {
    LOG_INFO("got delta chain len %d\n", delta_chain_len_);
    return delta_chain_len_;
  }

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
    int child_index = GetChild(key, children_, size_);
    LPID child_lpid = this->children_[child_index].second;
    auto child = this->map->GetMappingTable()->GetNode(child_lpid);
    while (child->GetDeltaChainLen() > DELTA_CHAIN_LIMIT) {
      this->map->CompressDeltaChain(child_lpid);
      child = this->map->GetMappingTable()->GetNode(child_lpid);
    }
    return child->DeleteEntry(key, location, child_lpid);
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
  int GetChild(__attribute__((unused)) KeyType key,
               __attribute__((unused)) std::pair<KeyType, LPID> *children,
               __attribute__((unused)) oid_t len);

  std::pair<KeyType, LPID> *GetChildren() { return children_; }

  void SplitNodes(LPID self, LPID parent);

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
        right_split_page_lpid_(rightSplitPage){};

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

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

 private:
  // This key included in left child, excluded in the right child
  KeyType modified_key_;
  ValueType modified_key_location_;

  // The LPID of the new LPage
  LPID right_split_page_lpid_;
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
        modified_val_(value) {
    LOG_INFO("Inside LPageUpdateDelta Constructor");
  };

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
        new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
                                                                key, location);
    bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
    if (!status) {
      delete new_delta;
    }
    return status;
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location, LPID self) {
    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
        new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
                                                                key, location);

    new_delta->SetDeleteFlag();
    bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
    if (!status) {
      delete new_delta;
    }
    return status;
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
  IPageUpdateDelta(BWTree<KeyType, ValueType, KeyComparator> *map,
                   BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node,
                   KeyType max_key_left_split_node,
                   KeyType max_key_right_split_node, LPID left_split_node_lpid,
                   LPID right_split_node_lpid, bool is_delete)
      : Delta<KeyType, ValueType, KeyComparator>(map, modified_node),
        max_key_left_split_node_(max_key_left_split_node),
        max_key_right_split_node_(max_key_right_split_node),
        left_split_node_lpid_(left_split_node_lpid),
        right_split_node_lpid_(right_split_node_lpid),
        is_delete_(is_delete) {
    LOG_INFO("Inside IPageUpdateDelta Constructor");
  };

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

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

 private:
  // The key which is modified
  KeyType max_key_left_split_node_, max_key_right_split_node_;

  // two members
  LPID left_split_node_lpid_, right_split_node_lpid_;

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
    left_sib_ = INVALID_LPID;
    right_sib_ = INVALID_LPID;
    size_ = 0;
  };

  LPage(BWTree<KeyType, ValueType, KeyComparator> *map,
        NodeStateBuilder<KeyType, ValueType, KeyComparator> *state)
      : BWTreeNode<KeyType, ValueType, KeyComparator>(map, 0) {
    LOG_INFO(" ");
    size_ = state->size;
    LNodeStateBuilder<KeyType, ValueType, KeyComparator> *lstate =
        reinterpret_cast<
            LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(state);

    memcpy(locations_, lstate->locations_,
           sizeof(std::pair<KeyType, ValueType>) * size_);
    left_sib_ = lstate->left_sibling_;
    if (state->IsSeparated()) {
      right_sib_ = state->split_new_page_id;
    } else {
      right_sib_ = lstate->right_sibling_;
    }
  };

  ~LPage(){};

  bool InsertEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    LOG_INFO("Inside LPage InsertEntry");

    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
        new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
                                                                key, location);

    bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
    if (!status) {
      delete new_delta;
    }
    return status;
  };

  bool DeleteEntry(__attribute__((unused)) KeyType key,
                   __attribute__((unused)) ValueType location,
                   __attribute__((unused)) LPID self) {
    // TODO implement this

    LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
        new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
                                                                key, location);

    new_delta->SetDeleteFlag();
    bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
    if (!status) {
      delete new_delta;
    }
    return status;
  };

  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  void SplitNodes(LPID self, LPID parent);

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

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
