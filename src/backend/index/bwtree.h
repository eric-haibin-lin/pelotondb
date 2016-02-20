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

namespace peloton {
namespace index {

//===--------------------------------------------------------------------===//
// Types and Enums
//===--------------------------------------------------------------------===//
typedef uint64_t LPID;

enum BWTreeNodeType {
  TYPE_BWTREE_NODE = 0,
  TYPE_LPAGE = 1,
  TYPE_IPAGE = 2,
  TYPE_DELTA = 3,
  TYPE_LPAGE_UPDATE_DELTA = 4,
  TYPE_IPAGE_UPDATE_DELTA = 5,
  // more types to add
};

#define IPAGE_ARITY 5
#define LPAGE_ARITY 5
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
  // number of k-v pairs
  oid_t size;

 public:
  NodeStateBuilder(oid_t size) : size(size){};

  virtual BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage() = 0;

  virtual ~NodeStateBuilder() = 0;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class INodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // IPage children nodes
  std::pair<KeyType, LPID> *children_ = nullptr;

 public:
  // IPage constructor
  INodeStateBuilder(std::pair<KeyType, LPID> *children, int children_len)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(children_len) {
    children_ = new std::pair<KeyType, LPID>[IPAGE_ARITY + DELTA_CHAIN_LIMIT]();
    for (int i = 0; i < children_len; i++) {
      children_[i] = children[i];
    }
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage() {
    // TODO: init an IPage based on collected state
    return nullptr;
  }

  //***************************************************
  // IPage Methods
  //***************************************************
  void AddChild(std::pair<KeyType, LPID> &new_pair) {
    assert(children_ != nullptr);
    KeyType key = new_pair.first;
    int index = IPage<KeyType, ValueType, KeyComparator>::GetChild(
        key, children_, this->size);
    // shift every element to the right
    for (int i = this->size; i > index; i--) {
      children_[i] = children_[i - 1];
    }
    // insert at the found index
    children_[index] = new_pair;
    // increase size
    this->size++;
  }

  void RemoveChild(KeyType key_to_remove) {
    assert(children_ != nullptr);
    int index = IPage<KeyType, ValueType, KeyComparator>::GetChild(
        key_to_remove, children_, this->size);
    // delete at the found index
    for (int i = index; i < this->size - 1; i++) {
      children_[i] = children_[i + 1];
    }
    // decrement size
    this->size--;
  }
};

template <typename KeyType, typename ValueType, class KeyComparator>
class LNodeStateBuilder
    : public NodeStateBuilder<KeyType, ValueType, KeyComparator> {
 private:
  // LPage members
  LPID left_sibling_ = 0;
  LPID right_sibling_ = 0;
  bool unique_keys_ = true;

  // LPage members
  std::pair<KeyType, ValueType> *locations_ = nullptr;

 public:
  // LPage constructor
  LNodeStateBuilder(LPID left_sibling, LPID right_sibling,
                    std::pair<KeyType, ValueType> *locations, int location_len,
                    bool unique_keys)
      : NodeStateBuilder<KeyType, ValueType, KeyComparator>(location_len),
        left_sibling_(left_sibling),
        right_sibling_(right_sibling),
        unique_keys_(unique_keys) {
    locations_ =
        new std::pair<KeyType, ValueType>[IPAGE_ARITY + DELTA_CHAIN_LIMIT]();
    for (int i = 0; i < location_len; i++) {
      locations_[i] = locations[i];
    }
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetPage() {
    // TODO: init an LPage based on collected state
    return nullptr;
  }

  //***************************************************
  // LPage Methods
  //***************************************************

  void UpdateLeftSib(LPID new_left_sib) { left_sibling_ = new_left_sib; }

  void UpdateRightSib(LPID new_right_sib) { right_sibling_ = new_right_sib; }

  void AddLeafData(std::pair<KeyType, ValueType> &new_pair) {
    assert(locations_ != nullptr);
    KeyType key = new_pair.first;
    int index = LPage<KeyType, ValueType, KeyComparator>::BinarySearch(
        locations_, this->size);
    // shift every element to the right
    for (int i = this->size; i > index; i--) {
      locations_[i] = locations_[i - 1];
    }
    // insert at the found index
    locations_[index] = new_pair;
    // increase size
    this->size++;
  }

  void RemoveLeafData(std::pair<KeyType, ValueType> &entry_to_remove) {
    assert(locations_ != nullptr);
    KeyType key = entry_to_remove.first;
    int index = LPage<KeyType, ValueType, KeyComparator>::BinarySearch(
        key, locations_, this->size);
    // keys are unique
    if (unique_keys_) {
      // delete at the found index
      for (int i = index; i < this->size - 1; i++) {
        locations_[i] = locations_[i + 1];
      }
      // decrement size
      this->size--;
      return;
    }
    // keys are not unique
    // we have the first appearance of the given key, do linear scan to see
    // which one matches exactly
    bool found_exact_key = false;
    for (; index < this->size; index++) {
      std::pair<KeyType, ValueType> pair = locations_[index];
      if (BWTreeNode<KeyType, ValueType, KeyComparator>::comparator(pair.first,
                                                                    key)) {
        if (pair.second == entry_to_remove.second) {
          found_exact_key = true;
        }
      } else {
        // not found
        break;
      }
    }
    if (found_exact_key) {
      for (int i = index; i < this->size - 1; i++) {
        locations_[i] = locations_[i + 1];
      }
      // decrement size
      this->size--;
    }
  }
};
//===--------------------------------------------------------------------===//
// BWTree
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree {
 private:
  // the default LPID for the root node is 2^64 - 1
  static const LPID DEFAULT_ROOT_LPID = ULLONG_MAX;

  // the logical page id for the root node
  LPID root_;

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
  BWTree() : root_(DEFAULT_ROOT_LPID) {
    mapping_table_ =
        new BWTreeNode<KeyType, ValueType, KeyComparator> *[mapping_table_cap_];
    // TODO @abj initialize the root IPage (and maybe a LPage?) ---
    // done!
    // with the given comparator
    // TODO Comparator too?

    /* First create an object of the IPAGE class, which then
     * acts as the root.
     */
    //@abj please fix the compiler warnings :P
    //    IPage<KeyType, ValueType, KeyComparator> *root =
    //        new IPage<KeyType, ValueType, KeyComparator>(this, true);
    //
    //    // TODO: do we need the () at the end? Also, this should be
    //    moved to
    //    // the constructor of IPAGE
    //    root->children_map_ = new std::pair<KeyType,
    //    LPID>[IPAGE_ARITY]();
    //
    //    /* Install the root in the mapping table */
    //    root_ = InstallPage(root);
    //
    //    /* Now create the first LPAGE */
    //    LPage<KeyType, ValueType, KeyComparator> *first_lpage =
    //        new LPage<KeyType, ValueType, KeyComparator>(this, true);
    //
    //    LPID first_lpage_lpid;
    //
    //    first_lpage_lpid = InstallPage(first_lpage);
    //
    //    /* Now grow a pair (:P) to store the first LPAGE pointer */
    //    std::pair<KeyType, LPID> *first_lpage_pair = new
    //    std::pair<KeyType, LPID>;
    //
    //    // TODO: some way to denote the max KeyType value
    //    // first_lpage_pair->first =
    //    first_lpage_pair->second = first_lpage_lpid;
    //
    //    root->children_map_[0] = *first_lpage_pair;
  };

  /*
   * On construction, BWTree will create a IPage and an empty LPage.
   * The IPage has only one pointer, which points to the empty LPage.
   * The IPage serves as the root of all other nodes.
   */
  BWTree(bool unique_keys, KeyComparator comparator)
      : root_(DEFAULT_ROOT_LPID) {
    mapping_table_ =
        new BWTreeNode<KeyType, ValueType, KeyComparator> *[mapping_table_cap_];
    // TODO @abj initialize the root IPage (and maybe a LPage?)
    // with the given comparator
  };

  ~BWTree() {
    delete[] free_LPIDs;
    delete[] mapping_table_;
  }
  // instead of explicitly use ValueType as the type for location, we should
  // use the template type ValueType instead (although it's ValueType is always
  // instantiated as ValueType class
  // see index/index_factory.cpp IndexFactory::GetInstance()
  bool InsertEntry(KeyType key, ValueType location);

  bool DeleteEntry(KeyType key, ValueType location);

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
  // return 0 if the page install is not successful

  LPID InstallPage(BWTreeNode<KeyType, ValueType, KeyComparator> *node) {
    LPID newLPID = __sync_fetch_and_add(&next_LPID_, 1);
    // table grew too large, expand it
    while (newLPID >= mapping_table_cap_) {
      // only one thread should expand the table
      AquireWrite();
      int new_mapping_table_cap = mapping_table_cap_ * 2;
      auto new_mapping_table =
          new (BWTreeNode<KeyType, ValueType, KeyComparator> *
               [new_mapping_table_cap]);
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
};

//===--------------------------------------------------------------------===//
// BWTreeNode
// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode {
 public:
  // the comparator for key
  static KeyComparator comparator;

  BWTreeNode(BWTree<KeyType, ValueType, KeyComparator> *map, bool unique_keys)
      : map(map), unique_keys(unique_keys){};

  // These are virtual methods which child classes have to implement.
  // They also have to be redeclared in the child classes
  virtual bool InsertEntry(KeyType key, ValueType location) = 0;

  virtual bool DeleteEntry(KeyType key, ValueType location);

  virtual std::vector<ValueType> Scan(
      const std::vector<Value> &values,
      const std::vector<oid_t> &key_column_ids,
      const std::vector<ExpressionType> &expr_types,
      const ScanDirectionType &scan_direction) = 0;

  virtual std::vector<ValueType> ScanAllKeys() = 0;

  virtual std::vector<ValueType> ScanKey(KeyType key) = 0;

  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *
  BuildNodeState() = 0;

  // for scan, we have to build the node state as well. But we only care about
  // the keys we want to scan
  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      KeyType key) = 0;

  virtual NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      const std::vector<Value> &values,
      const std::vector<oid_t> &key_column_ids,
      const std::vector<ExpressionType> &expr_types,
      const ScanDirectionType &scan_direction) = 0;

  // Each sub-class will have to implement this function to return their type
  virtual BWTreeNodeType GetTreeNodeType() const = 0;

  virtual ~BWTreeNode() = 0;

 protected:
  // the handler to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  // whether unique key is required
  bool unique_keys;
};

//===--------------------------------------------------------------------===//
// IPage
// The IPage hold pointers to all its children. An IPage with n keys
// k1, k2, .... kn actually has (n + 1) pointers to its children, where
// p1 for (-infinity, k1], p2 for (k1, k2], p3 for (k2, k3] ...
// pn for (k_n-1, kn], p_n+1 for (kn, +infinity). Variable right_most_child_
// stores the value of p_n+1
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class IPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  bool InsertEntry(KeyType key, ValueType location);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

  // get the index of the child at next level, which contains the given key
  static int GetChild(KeyType key, std::pair<KeyType, LPID> *children,
                      oid_t len);

 private:
  std::pair<KeyType, LPID> *children_;

  oid_t size_;
};

//===--------------------------------------------------------------------===//
// Delta
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class Delta : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_DELTA; };

 protected:
  // the modified node could either be a LPage or IPage or Delta
  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node;
};

//===--------------------------------------------------------------------===//
// LPageUpdateDelta
// LPageUpdateDelta either represents a insert delta or delete delta on LPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPageUpdateDelta : public Delta<KeyType, ValueType, KeyComparator> {
 public:
  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const {
    return TYPE_LPAGE_UPDATE_DELTA;
  };

 private:
  // The key which is modified
  KeyType modified_key_;

  // The location of the updated key
  // Set to INVALID_ITEMPOINTER for delete delta
  ValueType modified_val_;
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
  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  inline BWTreeNodeType GetTreeNodeType() const {
    return TYPE_IPAGE_UPDATE_DELTA;
  };

 private:
  // The key which is modified
  KeyType modified_key_;

  // The logical page id of the updated key. Set to 0 for delete delta
  LPID modified_id_;
};

// TODO More delta classes such as
// merge, split, remove_page, separator

//===--------------------------------------------------------------------===//
// LPage
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
class LPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  std::vector<ValueType> Scan(const std::vector<Value> &values,
                              const std::vector<oid_t> &key_column_ids,
                              const std::vector<ExpressionType> &expr_types,
                              const ScanDirectionType &scan_direction);

  std::vector<ValueType> ScanAllKeys();

  std::vector<ValueType> ScanKey(KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildNodeState();

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      KeyType key);

  NodeStateBuilder<KeyType, ValueType, KeyComparator> *BuildScanState(
      const std::vector<Value> &values,
      const std::vector<oid_t> &key_column_ids,
      const std::vector<ExpressionType> &expr_types,
      const ScanDirectionType &scan_direction);

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

  // get the index of the first occurrence of the given key
  static int BinarySearch(KeyType key, std::pair<KeyType, ValueType> *locations,
                          oid_t len);

 private:
  // return a vector of indices of the matched slots
  // or return an iterator?
  std::vector<int> ScanKeyInternal(KeyType key);

  // left_sibling pointer is used to do reverse iterate
  LPID left_sib_;
  LPID right_sib_;

  std::pair<KeyType, ValueType> *locations_;
  // the size of stored locations
  oid_t size_;
};

}  // End index namespace
}  // End peloton namespace
