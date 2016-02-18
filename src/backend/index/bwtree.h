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

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode;

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
  std::map<LPID, BWTreeNode<KeyType, ValueType, KeyComparator> **>
      mapping_table_;

 public:
  BWTree()
      : root_(DEFAULT_ROOT_LPID){
            // TODO @abj initialize the root IPage (and maybe a LPage?)
        };

  /*
   * On construction, BWTree will create a IPage and an empty LPage.
   * The IPage has only one pointer, which points to the empty LPage.
   * The IPage serves as the root of all other nodes.
   */
  BWTree(bool unique_keys)
      : root_(DEFAULT_ROOT_LPID){
            // TODO @abj initialize the root IPage (and maybe a LPage?)
        };

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

  // return 0 if the page install is not successful
  LPID InstallPage(BWTreeNode<KeyType, ValueType, KeyComparator> *node);

  bool SwapNode(LPID id, BWTreeNode<KeyType, ValueType, KeyComparator> *node);

  BWTreeNode<KeyType, ValueType, KeyComparator> *GetNode(LPID id);
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
      : map_(map), unique_keys_(unique_keys){};

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
  BWTree<KeyType, ValueType, KeyComparator> *map_;
  // whether unique key is required
  bool unique_keys_;
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
  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node_;
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

  // The location of the updated key. Set to INVALID for delete delta
  ValueType modified_val_;
};

//===--------------------------------------------------------------------===//
// IPageUpdateDelta
// IPageUpdateDelta either represents a insert delta or delete delta on IPage
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

  // scan with information of deleted/inserted keys/pages
  // std::vector<ValueType> ScanKey(const std::vector<Delta *> deltas);

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
