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
  TYPE_LPAGE_UPDATE_DELTA = 3,
  TYPE_IPAGE_UPDATE_DELTA = 4,
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
            // TODO initialize the root IPage (and maybe a LPage?)
        };

  BWTree(bool unique_keys)
      : root_(DEFAULT_ROOT_LPID){
            // TODO initialize the root IPage (and maybe a LPage?)
        };

  // instead of explicitly use ValueType as the type for location, we should
  // use the template type ValueType instead (although it's ValueType is always
  // templated as ValueType class
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
  // This is better than having to store redundant types in all the objects
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

 private:
  std::pair<KeyType, LPID> *children_map_;

  // get the LPID of the child at next level, which contains the given key
  // TODO implement this!
  BWTreeNode<KeyType, ValueType, KeyComparator> GetChild(KeyType key);
};

template <typename KeyType, typename ValueType, class KeyComparator>
class Delta : public BWTreeNode<KeyType, ValueType, KeyComparator> {
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
  KeyType modified_key_;
  LPID modified_id_;  // set to 0 for delete delta
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
  KeyType modified_key_;
  ValueType modified_val_;  // set to INVALID for delete delta
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
