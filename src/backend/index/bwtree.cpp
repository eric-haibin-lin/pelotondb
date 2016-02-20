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
// NodeStateBuilder
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// INodeStateBuilder
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LNodeStateBuilder
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
BWTreeNode<KeyType, ValueType, KeyComparator> *
LNodeStateBuilder<KeyType, ValueType, KeyComparator>::GetPage() {
  if (new_page == nullptr) {
    // TODO call LPage constructor
    // new_page = new LPage<KeyType, ValueType, KeyComparator>();
  }
  return new_page;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveLeafData(
    KeyType &key_to_remove) {
  // TODO remove entry based on key
  assert(locations_ != nullptr);
  // keys are unique
  assert(this->map->unique_keys);

  int index = LPage<KeyType, ValueType, KeyComparator>::BinarySearch(
      key_to_remove, locations_, this->size);
  // delete at the found index
  for (int i = index; i < this->size - 1; i++) {
    locations_[i] = locations_[i + 1];
  }
  // decrement size
  this->size--;
  return;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveLeafData(
    std::pair<KeyType, ValueType> &entry_to_remove) {
  assert(locations_ != nullptr);
  // keys are not unique
  assert(this->map->unique_keys == false);

  KeyType key = entry_to_remove.first;
  int index = LPage<KeyType, ValueType, KeyComparator>::BinarySearch(
      key, locations_, this->size);
  // we have the first appearance of the given key, do linear scan to see
  // which one matches exactly
  bool found_exact_key = false;
  for (; index < this->size; index++) {
    std::pair<KeyType, ValueType> pair = locations_[index];
    if (this->map->comparator(pair.first, key)) {
      if (ItemPointerEquals(pair.second, entry_to_remove.second)) {
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

template <typename KeyType, typename ValueType, class KeyComparator>
bool LNodeStateBuilder<KeyType, ValueType, KeyComparator>::ItemPointerEquals(
    ValueType v1, ValueType v2) {
  return v1.block == v2.block && v1.offset == v2.offset;
};

//===--------------------------------------------------------------------===//
// BWTree Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result =
      GetNode(root_)->Scan(values, key_column_ids, expr_types, scan_direction);

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result = GetNode(root_)->ScanAllKeys();

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key) {
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result = GetNode(root_)->ScanKey(key);

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool BWTree<KeyType, ValueType, KeyComparator>::InsertEntry(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) ValueType location) {
  // just call InsertEntry on root
  // return false;
  return GetNode(root_)->InsertEntry(key, location);
};

//===--------------------------------------------------------------------===//
// IPage Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> IPage<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
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
    KeyType key) {
  std::vector<ValueType> result;

  int child_idx = GetChild(key, this->children_, this->size_);
  LPID child_id = this->children_[child_idx].second;

  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetNode(child_id);
  if (child == nullptr) {
    // Key is not included in the tree, do nothing
  } else {
    result = child->ScanKey(key);
  }
  return result;
};
//@abj please fix the warnings :P
template <typename KeyType, typename ValueType, class KeyComparator>
bool IPage<KeyType, ValueType, KeyComparator>::InsertEntry(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) ValueType location) {
  /* int i;
   bool last_level_page;
   LPID target_child_lpid;  // TODO: write code to get this -- partially done*/

  return this->map->GetNode(GetChild(key, children_, size_))
      ->InsertEntry(key, location);
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

  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
IPage<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder;
  // build node state for IPage
  builder = new INodeStateBuilder<KeyType, ValueType, KeyComparator>(
      children_, size_, this->map);
  return builder;
};

//===--------------------------------------------------------------------===//
// LPageUpdateDelta Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
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
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
  std::vector<ValueType> result;
  if (this->map->unique_keys) {
    // the modified key matches the scanKey
    if (this->map->comparator(modified_key_, key) == true) {
      if (!is_delete_) {
        // the modified key is inserted, add to result vector
        result.push_back(modified_val_);
      }
    }
    return result;
  }
  // non unique key. we have to build the state
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->modified_node->BuildScanState(key);
  BWTreeNode<KeyType, ValueType, KeyComparator> *page = builder->GetPage();
  assert(page != nullptr);
  // do scan on the new state
  result = page->ScanKey(key);
  //release builder
  delete(builder);
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  // Children of LPageDelta always return a LNodeStateBuilder
  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
  // delete delta
  if (is_delete_) {
    if (this->map->unique_keys) {
      builder->RemoveLeafData(this->modified_key_);
    } else {
      std::pair<KeyType, ValueType> pair(modified_key_, modified_val_);
      builder->RemoveLeafData(pair);
    }

  } else {
    // insert delta
    std::pair<KeyType, ValueType> pair(modified_key_, modified_val_);
    builder->AddLeafData(pair);
  }
  return builder;
};

//===--------------------------------------------------------------------===//
// LPage Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
int LPage<KeyType, ValueType, KeyComparator>::BinarySearch(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) std::pair<KeyType, ValueType> *locations,
    __attribute__((unused)) oid_t len) {
  // TODO @Matt implement this
  return 0;
}

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> LPage<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
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
    KeyType key) {
  std::vector<ValueType> result;
  std::vector<oid_t> indices = ScanKeyInternal(key);
  // we only need the values
  oid_t index;
  for (index = 0; index < indices.size(); index++) {
    result.push_back((locations_)[index].second);
  }
  // reach the end of current LPage, go to next Lpage for more results
  if (index == size_) {
    std::vector<ValueType> sib_result =
        this->map->GetNode(right_sib_)->ScanKey(key);
    result.insert(result.end(), sib_result.begin(), sib_result.end());
  }
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<oid_t> LPage<KeyType, ValueType, KeyComparator>::ScanKeyInternal(
    KeyType key) {
  assert(locations_ != nullptr);
  std::vector<oid_t> result;
  // empty LPage
  if (size_ == 0) {
    return result;
  }
  assert(size_ > 0);
  // do a binary search on locations to get the key
  int index = BinarySearch(key, locations_, size_);
  if (index == -1) {
    // key not found, return empty result
    return result;
  }

  // try to collect all matching keys. If unique_keys, only one key matches
  while (index < size_) {
    std::pair<KeyType, ValueType> location = (locations_)[index];
    if (this->map->comparator(location.first, key) == true) {
      // found a matching key
      result.push_back(index++);
    } else {
      // key not found, return result
      return result;
    }
  }
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPage<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder;
  // build node state for LPage
  builder = new LNodeStateBuilder<KeyType, ValueType, KeyComparator>(
      left_sib_, right_sib_, locations_, size_, this->map);
  return builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPage<KeyType, ValueType, KeyComparator>::BuildScanState(__attribute__((unused))
                                                         KeyType key) {
  // TODO call ScanKeyInternal(key);
  return nullptr;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPage<KeyType, ValueType, KeyComparator>::BuildScanState(
    __attribute__((unused)) const std::vector<Value> &values,
    __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  // TODO call ScanKeyInternal(values, ...);
  return nullptr;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPage<KeyType, ValueType, KeyComparator>::IsInvalidItemPointer(
    ValueType val) {
  return val.block == INVALID_OID || val.offset == INVALID_OID;
}

}  // End index namespace
}  // End peloton namespace
