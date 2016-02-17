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
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  //TODO implement this
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
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
IPage<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> IPage<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

/*
 * methods implementation for Delta
 */

/*
 * methods implementation for LPage
 */

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> Delta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
Delta<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> Delta<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  //TODO implement this
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
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer>
LPage<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ItemPointer> LPage<KeyType, ValueType, KeyComparator>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {
  std::vector<ItemPointer> result;
  //TODO implement this
  return result;
};
}  // End index namespace
}  // End peloton namespace
