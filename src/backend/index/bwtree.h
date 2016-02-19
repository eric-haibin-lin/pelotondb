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
 private:
  // type of node, should be IPage or LPage
  BWTreeNodeType output_type_;
  // LPage members
  LPID left_sibling_ = 0;
  LPID right_sibling_ = 0;
  // TODO: change from maps to data structures of fixed length (arrays)
  std::map<KeyType, LPID> children_map_;

  // IPage members
  std::map<KeyType, ValueType> leaf_data_;

 public:
  // LPage constructor
  NodeStateBuilder(LPID left_sibling, LPID right_sibling,
                   std::pair<KeyType, ValueType> *leaf_data, int leaf_data_len)
      : output_type_(TYPE_LPAGE),
        left_sibling_(left_sibling),
        right_sibling_(right_sibling) {
    for (int i = 0; i < leaf_data_len; i++) {
      leaf_data_[leaf_data[i]->first] = leaf_data[i]->second;
    }
  }

  // IPage constructor
  NodeStateBuilder(std::pair<KeyType, LPID> *children_map, int children_map_len)
      : output_type_(TYPE_IPAGE) {
    for (int i = 0; i < children_map_len; i++) {
      children_map_[children_map[i]->first] = children_map[i]->second;
    }
  }

  //***************************************************
  // IPage Methods
  //***************************************************
  void AddChild(std::pair<KeyType, LPID> &new_pair) {
    children_map_[new_pair->first] = new_pair->second;
  }

  void RemoveChild(KeyType key_to_remove) {
    children_map_.erase(key_to_remove);
  }

  //***************************************************
  // LPage Methods
  //***************************************************

  void UpdateLeftSib(LPID new_left_sib) { left_sibling_ = new_left_sib; }

  void UpdateRightSib(LPID new_right_sib) { right_sibling_ = new_right_sib; }

  void AddLeafData(std::pair<KeyType, ValueType> &new_entry) {
    leaf_data_[new_entry->first] = new_entry->second;
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
  void AquireRead() {
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
  void ReleaseRead() { __sync_add_and_fetch(&current_readers, -1); }
  void AquireWrite() {
    while (__sync_bool_compare_and_swap(&current_writers, 0, 1))
      ;
    while (current_readers > 0)
      ;
  }
  void ReleaseWrite() {
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

  // Each sub-class will have to implement this function to return their type
  virtual BWTreeNodeType GetTreeNodeType() const = 0;

  virtual ~BWTreeNode() = 0;

 protected:
  // the handler to the mapping table
  BWTree<KeyType, ValueType, KeyComparator> *map;

  // whether unique key is required
  bool unique_keys;

  // the comparator for key
  KeyComparator comparator;
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

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_IPAGE; };

 private:
  std::pair<KeyType, LPID> *children_map_;
  LPID right_most_child_;

  // get the LPID of the child at next level, which contains the given key
  // TODO implement this
  BWTreeNode<KeyType, ValueType, KeyComparator> GetChild(KeyType key);
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

  inline BWTreeNodeType GetTreeNodeType() const { return TYPE_LPAGE; };

 private:
  // left_sibling pointer is used to do reverse iterate
  LPID left_sib_;
  LPID right_sib_;

  std::pair<KeyType, ValueType> *locations_;
  // the size of stored locations
  oid_t size_;
};

}  // End index namespace
}  // End peloton namespace
