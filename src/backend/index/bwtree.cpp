//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// bwtree.cpp
//
// Identification: src/backend/index/bwtree.cpp
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "backend/index/bwtree.h"

namespace peloton {
namespace index {

/*
 * methods implementation for BWTree
 */

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> BWTree<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  assert(mapping_table_.size() > 0);
  assert(unique_keys_ == true);

  // TODO define a macro for easier access to the pages' & addresses
  result = (*(mapping_table_[root_]))
               ->Scan(values, key_column_ids, expr_types, scan_direction);
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  assert(mapping_table_.size() > 0);
  assert(unique_keys_ == true);
  result = (*(mapping_table_[root_]))->ScanAllKeys();
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  assert(mapping_table_.size() > 0);
  assert(unique_keys_ == true);
  result = GetNode(root_)->ScanKey(key);
  return result;
};

/*
 * methods implementation for IPage
 */
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> IPage<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
IPage<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> IPage<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  assert(key != nullptr);
  BWTreeNode<KeyType, ValueType, KeyComparator> *child = GetChild(key);
  if (child == nullptr) {
    // Key is not included in the tree
  } else {
    result = child->ScanKey(key);
  }
  return result;
};

/*
 * methods implementation for DeleteDelta
 */
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> DeleteDelta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
DeleteDelta<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
DeleteDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  // TODO get the result from child, delete the key
  return result;
};

/*
 * methods implementation for InsertDelta
 */
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> InsertDelta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
InsertDelta<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
InsertDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

/*
 * method implementation for LPage
 */
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> LPage<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
LPage<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> LPage<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  // TODO implement this
  return result;
};
}  // End index namespace
}  // End peloton namespace
