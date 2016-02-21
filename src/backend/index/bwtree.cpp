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
#include "backend/common/exception.h"

namespace peloton {
namespace index {

//===--------------------------------------------------------------------===//
// NodeStateBuilder
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// INodeStateBuilder
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
BWTreeNode<KeyType, ValueType, KeyComparator> *
INodeStateBuilder<KeyType, ValueType, KeyComparator>::GetPage() {
  if (new_page_ == nullptr) {
    // TODO call IPage constructor
    // new_page = new IPage<KeyType, ValueType, KeyComparator>();
  }
  return new_page_;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::AddChild(
    std::pair<KeyType, LPID> &new_pair) {
  assert(children_ != nullptr);
  KeyType key = new_pair.first;
  int index = IPage<KeyType, ValueType, KeyComparator>::GetChild(key, children_,
                                                                 this->size);
  assert(index < IPAGE_ARITY + DELTA_CHAIN_LIMIT);
  // shift every element to the right
  for (int i = this->size; i > index; i--) {
    children_[i] = children_[i - 1];
  }
  // insert at the found index
  children_[index] = new_pair;
  // increase size
  this->size++;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveChild(
    KeyType &key_to_remove) {
  assert(children_ != nullptr);
  int index = IPage<KeyType, ValueType, KeyComparator>::GetChild(
      key_to_remove, children_, this->size);
  // if key found
  if (index < this->size && index >= 0) {
    for (int i = index; i < this->size - 1; i++) {
      children_[i] = children_[i + 1];
    }
    // decrement size
    this->size--;
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::SeparateFromKey(
    KeyType separator_key, LPID split_new_page_id) {
  assert(children_ != nullptr);
  int index = IPage<KeyType, ValueType, KeyComparator>::GetChild(
      separator_key, children_, this->size);
  assert(index < this->size && index >= 0);
  // assume we exclude the key at the split page
  // decrement size
  this->size = index;
  // update separator info
  this->is_separated = true;
  this->separator_key = separator_key;
  this->split_new_page_id = split_new_page_id;
}

//===--------------------------------------------------------------------===//
// LNodeStateBuilder
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
BWTreeNode<KeyType, ValueType, KeyComparator> *
LNodeStateBuilder<KeyType, ValueType, KeyComparator>::GetPage() {
  if (new_page_ == nullptr) {
    new_page_ = new LPage<KeyType, ValueType, KeyComparator>(this->map, this);
  }
  return new_page_;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::AddLeafData(
    std::pair<KeyType, ValueType> &new_pair) {
  assert(locations_ != nullptr);
  KeyType key = new_pair.first;
  int index = this->map->BinarySearch(key, locations_, this->size);
  assert(index < IPAGE_ARITY + DELTA_CHAIN_LIMIT);
  // not found. shift every element to the right
  if (index < 0) {
    index = -1 * index;
    for (int i = this->size; i > index; i--) {
      locations_[i] = locations_[i - 1];
    }
  } else {
    // found
    if (!this->map->unique_keys) {
      // shift every element to the right if we allow non-unique keys
      for (int i = this->size; i > index; i--) {
        locations_[i] = locations_[i - 1];
      }
    }
  }
  // insert at the found index
  locations_[index] = new_pair;
  // increase size
  this->size++;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveLeafData(
    KeyType &key_to_remove) {
  // TODO remove entry based on key
  assert(locations_ != nullptr);
  // keys are unique
  assert(this->map->unique_keys);
  // TODO implement non-unique-key case

  int index = this->map->BinarySearch(key_to_remove, locations_, this->size);
  // if key found
  if (index < this->size && index >= 0) {
    for (int i = index; i < this->size - 1; i++) {
      locations_[i] = locations_[i + 1];
    }
    // decrement size
    this->size--;
  }
  return;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveLeafData(
    std::pair<KeyType, ValueType> &entry_to_remove) {
  assert(locations_ != nullptr);
  // keys are not unique
  assert(this->map->unique_keys == false);

  KeyType key = entry_to_remove.first;
  static int index = this->map->BinarySearch(key, locations_, this->size);
  // we have the first appearance of the given key, do linear scan to see
  // which one matches exactly
  bool found_exact_key = false;
  for (; index < this->size; index++) {
    std::pair<KeyType, ValueType> pair = locations_[index];
    if (this->map->CompareKey(key, pair.first) == 0) {
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
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::SeparateFromKey(
    KeyType separator_key, ValueType location, LPID split_new_page_id) {
  assert(locations_ != nullptr);

  int index = this->map->BinarySearch(
      std::pair<KeyType, ValueType>(separator_key, location), locations_,
      this->size);
  assert(index < this->size && index >= 0);
  // assume we exclude the key at the split page
  // decrement size
  this->size = index;
  // update separator info
  this->is_separated = true;
  this->separator_key = separator_key;
  separator_location_ = location;
  this->split_new_page_id = split_new_page_id;
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
bool BWTree<KeyType, ValueType, KeyComparator>::InsertEntry(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) ValueType location) {
  // just call InsertEntry on root
  // return false;
  LPID child_lpid;
  child_lpid = root_;

  return GetNode(child_lpid)->InsertEntry(key, location, child_lpid);
};

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

//===--------------------------------------------------------------------===//
// IPage Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
bool IPage<KeyType, ValueType, KeyComparator>::InsertEntry(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) ValueType location,
    __attribute__((unused)) LPID self) {
  /* int i;
   bool last_level_page;
   LPID target_child_lpid;  // TODO: write code to get this -- partially done*/
	LOG_INFO("Inside IPage InsertEntry");
	int child_lpid_index = GetChild(key, children_, size_);
	LOG_INFO("Got child_lpid_index as %d", child_lpid_index);
  //LPID child_lpid = GetChild(key, children_, size_);
  return this->map->GetNode(children_[child_lpid].second)->InsertEntry(key, location, child_lpid);
};

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
  // locate the child who covers the key
  int child_idx = GetChild(key, this->children_, this->size_);
  LPID child_id = this->children_[child_idx].second;

  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetNode(child_id);
  assert(child != nullptr);
  result = child->ScanKey(key);
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

template <typename KeyType, typename ValueType, class KeyComparator>
// get the index of the child at next level, which contains the given key
int IPage<KeyType, ValueType, KeyComparator>::GetChild(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) std::pair<KeyType, LPID> *children,
    __attribute__((unused)) oid_t len) {
  int index = this->map->BinarySearch(key, children, len);
  return index >= 0 ? index : -index;
};

//===--------------------------------------------------------------------===//
// Delta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> Delta<KeyType, ValueType, KeyComparator>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> Delta<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;
  // TODO implement this
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> Delta<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key) {
  std::vector<ValueType> result;
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->BuildNodeState();
  BWTreeNode<KeyType, ValueType, KeyComparator> *page = builder->GetPage();
  assert(page != nullptr);
  if (!builder->IsSeparated()) {
    // the Page doesn't split, scan left page
    result = page->ScanKey(key);
  } else {
    // the Page splits
    KeyType separator_key = builder->GetSeparatorKey();
    // the scan key is in the left page, scan left page
    // if the scan == separator key, we still scan from the left page since non
    // duplicate key and go across boarder
    if (this->map->CompareKey(separator_key, key) <= 0) {
      result = page->ScanKey(key);
    } else {
      // the desired key is in the right page
      LPID right_page_id = builder->GetSplitNewPageId();
      // scan right page
      result = this->map->GetNode(right_page_id)->ScanKey(key);
    }
  }
  // release builder
  delete (builder);
  return result;
};
//===--------------------------------------------------------------------===//
// Delta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// IPageSplitDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
IPageSplitDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  // Children of IPageDelta always return a INodeStateBuilder
  INodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<INodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
  assert(builder != nullptr);
  std::pair<KeyType, LPID> separator(modified_key_, modified_val_);
  builder->SeparateFromKey(separator);

  return builder;
}
//===--------------------------------------------------------------------===//
// IPageSplitDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPageSplitDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPageSplitDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  // Children of IPageDelta always return a INodeStateBuilder
  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
  assert(builder != nullptr);
  builder->SeparateFromKey(modified_key_, modified_key_location_,
                           modified_val_);

  return builder;
}
//===--------------------------------------------------------------------===//
// LPageSplitDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// IPageUpdateDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
IPageUpdateDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  // Children of IPageDelta always return a INodeStateBuilder
  INodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<INodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
  assert(builder != nullptr);
  // delete delta
  if (is_delete_) {
    builder->RemoveChild(this->modified_key_);
  } else {
    // insert delta
    std::pair<KeyType, LPID> pair(modified_key_, modified_id_);
    builder->AddChild(pair);
  }
  return builder;
};
//===--------------------------------------------------------------------===//
// IPageUpdateDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPageUpdateDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
  std::vector<ValueType> result;
  LOG_INFO(" ");
  if (this->map->unique_keys) {
    // the modified key matches the scanKey
    if (this->map->CompareKey(modified_key_, key) == 0) {
      LOG_INFO("Found a matching key");
      if (!is_delete_) {
        // the modified key is inserted, add to result vector
        LOG_INFO("Inserts the matching key");
        result.push_back(modified_val_);
      }
      return result;
    }
  }
  LOG_INFO("Building NodeState of all children nodes");

  // non unique key. we have to build the state
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->modified_node->BuildScanState(key);

  BWTreeNode<KeyType, ValueType, KeyComparator> *page = builder->GetPage();
  assert(page != nullptr);
  // TODO check split status of builder for shortcut
  // do scan on the new state
  LOG_INFO("Scan on the new NodeState");
  result = page->ScanKey(key);
  // release builder
  delete (builder);
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  LOG_INFO(" ");
  // Children of LPageDelta always return a LNodeStateBuilder
  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
  assert(builder != nullptr);
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
// LPageUpdateDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPage Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
template <typename PairSecond>
int BWTree<KeyType, ValueType, KeyComparator>::BinarySearch(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) std::pair<KeyType, PairSecond> *locations,
    __attribute__((unused)) oid_t len) {
  int low = 0, high = len - 1, first = -1;
  int mid;
  while (low <= high) {
    mid = low + ((high - low) >> 1);
    switch (CompareKey(locations[mid].first, key)) {
      case -1:
        low = mid + 1;
        break;
      case 0:
        first = mid;
      // no break
      default:
        // positive case
        high = mid - 1;
        break;
    }
  }
  if (first < 0) {
    return first;
  } else {
    return -low;
  }
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
  LOG_INFO(" ");
  std::vector<ValueType> result;
  std::vector<oid_t> indices = ScanKeyInternal(key);
  // we only need the values
  oid_t index;
  for (index = 0; index < indices.size(); index++) {
    result.push_back((locations_)[index].second);
  }
  // reach the end of current LPage, go to next LPage for more results
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
  LOG_INFO(" ");
  assert(locations_ != nullptr);
  std::vector<oid_t> result;
  // empty LPage
  if (size_ == 0) {
    return result;
  }
  assert(size_ > 0);
  // do a binary search on locations to get the key
  int index = this->map->BinarySearch(key, locations_, size_);
  if (index == -1) {
    // key not found, return empty result
    return result;
  }

  // try to collect all matching keys. If unique_keys, only one key matches
  while (index < size_) {
    std::pair<KeyType, ValueType> location = (locations_)[index];
    if (this->map->CompareKey(location.first, key) == 0) {
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

// DO NOT DELETE
/*
 template <typename KeyType, typename ValueType, class KeyComparator>
 bool LPage<KeyType, ValueType, KeyComparator>::IsInvalidItemPointer(
    ValueType val) {
  return val.block == INVALID_OID || val.offset == INVALID_OID;
}*/

//===--------------------------------------------------------------------===//
// LPage Methods End
//===--------------------------------------------------------------------===//

}  // End index namespace
}  // End peloton namespace
