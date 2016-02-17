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

//===--------------------------------------------------------------------===//
// BWTree Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  assert(mapping_table_.size() > 0);

  // TODO define a macro for easier access to the pages' & addresses
  result = (*(mapping_table_[root_]))
               ->Scan(values, key_column_ids, expr_types, scan_direction);
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;
  assert(mapping_table_.size() > 0);
  result = (*(mapping_table_[root_]))->ScanAllKeys();
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key) {
  std::vector<ValueType> result;
  assert(mapping_table_.size() > 0);
  result = GetNode(root_)->ScanKey(key);
  return result;
};

//===--------------------------------------------------------------------===//
// IPage Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> IPage<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> IPage<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> IPage<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) KeyType key) {
  std::vector<ValueType> result;
  assert(key != nullptr);
  BWTreeNode<KeyType, ValueType, KeyComparator> *child = GetChild(key);
  if (child == nullptr) {
    // Key is not included in the tree
  } else {
    result = child->ScanKey(key);
  }
  return result;
};

//===--------------------------------------------------------------------===//
// IPageUpdateDelta Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
IPageUpdateDelta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
IPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
IPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
  std::vector<ValueType> result;
  assert(key != nullptr);
  assert(this->unique_keys_ == true);
  return result;
};

//===--------------------------------------------------------------------===//
// LPageUpdateDelta Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) KeyType key) {
  std::vector<ValueType> result;
  if (this->unique_keys_) {
    // check the underlying type
  }
  // TODO get the result from child, delete the key
  return result;
};

//===--------------------------------------------------------------------===//
// LPage Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> LPage<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> LPage<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> LPage<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) KeyType key) {
  std::vector<ValueType> result;
  // in the first version, the LPage has no content at all
  assert(size_ == 0);
  return result;
};
}  // End index namespace
}  // End peloton namespace
