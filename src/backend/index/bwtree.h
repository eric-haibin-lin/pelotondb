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
#include "backend/storage/tuple.h"
#include <map>
#include <vector>

namespace peloton {
namespace index {

typedef uint32_t LPID;

// declaration of all the classes
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree {
 private:
  LPID root;
  std::map<LPID, BWTree<KeyType, ValueType, KeyComparator> **> mapping_table_;

 public:
  bool InsertEntry(const storage::Tuple *key, const ItemPointer location);

  bool DeleteEntry(const storage::Tuple *key, const ItemPointer location);

  std::vector<ItemPointer> Scan(const std::vector<Value> &values,
                                const std::vector<oid_t> &key_column_ids,
                                const std::vector<ExpressionType> &expr_types,
                                const ScanDirectionType &scan_direction);
  std::vector<ItemPointer> ScanAllKeys();
  std::vector<ItemPointer> ScanKey(const storage::Tuple *key);

  // TODO pick one value for failure status
  LPID InstallPage(BWTreeNode<KeyType, ValueType, KeyComparator> *node);
  bool SwapNode(LPID id, BWTreeNode<KeyType, ValueType, KeyComparator> *node);
  BWTreeNode<KeyType, ValueType, KeyComparator> *Lookup(LPID id);
};

// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode {
 protected:
  BWTree<KeyType, ValueType, KeyComparator> *map_;

 public:
  BWTreeNode(BWTree<KeyType, ValueType, KeyComparator> *map) : map_(map){};

  // These are virtual methods which child classes have to implement.
  // They also have to be redeclared in the child classes
  virtual bool InsertEntry(const storage::Tuple *key,
                           const ItemPointer location) = 0;

  virtual bool DeleteEntry(const storage::Tuple *key,
                           const ItemPointer location);

  virtual std::vector<ItemPointer> Scan(
      const std::vector<Value> &values,
      const std::vector<oid_t> &key_column_ids,
      const std::vector<ExpressionType> &expr_types,
      const ScanDirectionType &scan_direction) = 0;

  virtual std::vector<ItemPointer> ScanAllKeys() = 0;

  virtual std::vector<ItemPointer> ScanKey(const storage::Tuple *key) = 0;

  virtual ~BWTreeNode() = 0;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class IPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  std::vector<ItemPointer> Scan(const std::vector<Value> &values,
                                const std::vector<oid_t> &key_column_ids,
                                const std::vector<ExpressionType> &expr_types,
                                const ScanDirectionType &scan_direction);

  std::vector<ItemPointer> ScanAllKeys();

  std::vector<ItemPointer> ScanKey(const storage::Tuple *key);

 private:
  std::pair<KeyType, LPID> *children_map_;
};

template <typename KeyType, typename ValueType, class KeyComparator>
class Delta : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  std::vector<ItemPointer> Scan(const std::vector<Value> &values,
                                const std::vector<oid_t> &key_column_ids,
                                const std::vector<ExpressionType> &expr_types,
                                const ScanDirectionType &scan_direction);

  std::vector<ItemPointer> ScanAllKeys();

  std::vector<ItemPointer> ScanKey(const storage::Tuple *key);

 protected:
  BWTreeNode<KeyType, ValueType, KeyComparator> *modified_node_;
};

// TODO More delta classes such as
// delete, insert, merge, split, remove_page, separator

template <typename KeyType, typename ValueType, class KeyComparator>
class LPage : public BWTreeNode<KeyType, ValueType, KeyComparator> {
 public:
  std::vector<ItemPointer> Scan(const std::vector<Value> &values,
                                const std::vector<oid_t> &key_column_ids,
                                const std::vector<ExpressionType> &expr_types,
                                const ScanDirectionType &scan_direction);

  std::vector<ItemPointer> ScanAllKeys();

  std::vector<ItemPointer> ScanKey(const storage::Tuple *key);

 private:
  // we don't support left_sibling pointer for the moment
  LPID left_sib_;
  LPID right_sib_;
  std::pair<KeyType, ItemPointer> *locations_;
};

}  // End index namespace
}  // End peloton namespace
