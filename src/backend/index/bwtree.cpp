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
#include <iostream>

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
  return new IPage<KeyType, ValueType, KeyComparator>(this->map, this);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::AddChild(
    std::pair<KeyType, LPID> &new_pair) {
  assert(children_ != nullptr);
  KeyType key = new_pair.first;
  int index = this->map->BinarySearch(key, children_, this->size - 1);
  assert(index < IPAGE_ARITY + IPAGE_DELTA_CHAIN_LIMIT);
  // Key not found. shift every element to the right
  if (index < 0 || this->size == 0 ||
      (index == 0 && this->map->CompareKey(children_[0].first, key) != 0)) {
    index = -index;
    for (int i = this->size; i > index; i--) {
      children_[i] = children_[i - 1];
    }
    // increase size
    this->size++;
  }
  // insert at the found index
  children_[index] = new_pair;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::ReplaceLastChild(
    LPID &inf_LPID) {
  assert(children_ != nullptr);
  children_[this->size - 1].second = inf_LPID;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveLastChild() {
  assert(children_ != nullptr);
  this->size--;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveChild(
    KeyType &key_to_remove) {
  assert(children_ != nullptr);
  int index = this->map->BinarySearch(key_to_remove, children_, this->size - 1);
  // if key found
  if ((index == 0 &&
       this->map->CompareKey(children_[0].first, key_to_remove) == 0) ||
      (index < this->size && index > 0)) {
    for (int i = index; i < this->size - 1; i++) {
      children_[i] = children_[i + 1];
    }
    // decrement size
    this->size--;
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::ScanAllKeys(
    std::vector<ValueType> &result) {
  LPID child_id = this->children_[0].second;
  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetMappingTable()->GetNode(child_id);
  child->ScanAllKeys(result);
}

//===--------------------------------------------------------------------===//
// LNodeStateBuilder
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
BWTreeNode<KeyType, ValueType, KeyComparator> *
LNodeStateBuilder<KeyType, ValueType, KeyComparator>::GetPage() {
  return new LPage<KeyType, ValueType, KeyComparator>(this->map, this);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::AddLeafData(
    std::pair<KeyType, ValueType> &new_pair) {
  assert(locations_ != nullptr);
  KeyType key = new_pair.first;
  int index = this->map->BinarySearch(key, locations_, this->size);

  LOG_INFO(
      "LNodeStateBuilder::AddLeafData - BinarySearch (%lu, %lu) result %d, "
      "original size: %lu",
      new_pair.second.block, new_pair.second.offset, index, this->size);

  assert(index < LPAGE_ARITY + LPAGE_DELTA_CHAIN_LIMIT);
  // not found. shift every element to the right
  if (index < 0 || this->size == 0 ||
      (index == 0 && this->map->CompareKey(locations_[0].first, key) != 0)) {
    // LOG_INFO("LNodeStateBuilder::AddLeafData - key not found");
    index = -1 * index;
    for (int i = this->size; i > index; i--) {
      locations_[i] = locations_[i - 1];
    }
    this->size++;
  } else {
    // LOG_INFO("LNodeStateBuilder::AddLeafData - key found");
    if (!this->map->unique_keys) {
      // shift every element to the right if we allow non-unique keys
      for (int i = this->size; i > index; i--) {
        locations_[i] = locations_[i - 1];
      }
      this->size++;
    }
  }
  // insert at the found index
  locations_[index] = new_pair;
  LOG_INFO("LNodeStateBuilder::AddLeafData - size after addition %lu",
           this->size);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveLeafData(
    std::pair<KeyType, ValueType> &entry_to_remove) {
  assert(locations_ != nullptr);

  // key doesn't have to be unique
  KeyType key = entry_to_remove.first;
  int index = this->map->BinarySearch(key, locations_, this->size);
  LOG_INFO(
      "LNodeStateBuilder::RemoveLeafData - BinarySearch result %d, original "
      "size: %lu",
      index, this->size);
  // we have the first appearance of the given key, do linear scan to see
  // which one matches exactly
  int found_exact_entry_count = 0;
  int dest = index;
  bool key_matches = true;
  for (int src = index; src >= 0 && src < this->size; src++) {
    std::pair<KeyType, ValueType> pair = locations_[src];
    if (key_matches && this->map->CompareKey(key, pair.first) == 0) {
      if (ItemPointerEquals(pair.second, entry_to_remove.second)) {
        // LOG_INFO(
        //    "LNodeStateBuilder::RemoveLeafData - found exact entry at index
        //    %d",
        //    src);
        found_exact_entry_count++;
      } else {
        // value doesn't match
        // LOG_INFO("LNodeStateBuilder::RemoveLeafData - not exact entry");
        locations_[dest++] = pair;
      }
    } else {
      key_matches = false;
      // LOG_INFO("LNodeStateBuilder::RemoveLeafData - key not match");
      locations_[dest++] = pair;
    }
  }

  // decrement size
  this->size -= found_exact_entry_count;
  LOG_INFO("LNodeStateBuilder::RemoveLeafData - size after removal %lu",
           this->size);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::SeparateFromKey(
    KeyType separator_key, int sep_index, LPID split_new_page_id) {
  assert(locations_ != nullptr);

  int index = this->map->BinarySearch(separator_key, locations_, this->size);
  assert(index < this->size && index >= 0);
  // assume we include the key at the split page
  // decrement size
  this->size = index + 1;
  // update separator info
  this->is_separated = true;
  this->separator_key = separator_key;
  separator_index_ = sep_index;
  this->split_new_page_id = split_new_page_id;
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool LNodeStateBuilder<KeyType, ValueType, KeyComparator>::ItemPointerEquals(
    ValueType v1, ValueType v2) {
  return v1.block == v2.block && v1.offset == v2.offset;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction, std::vector<ValueType> &result,
    const KeyType *index_key) {
  LOG_INFO("LNodeStateBuilder::Scan, size: %lu", this->size);
  this->map->ScanHelper(values, key_column_ids, expr_types, scan_direction,
                        index_key, result, locations_, this->size,
                        right_sibling_);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::ScanAllKeys(
    std::vector<ValueType> &result) {
  LOG_INFO("LNodeStateBuilder::ScanAllKeys, size: %lu", this->size);
  this->map->ScanAllKeysHelper(this->size, locations_, right_sibling_, result);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LNodeStateBuilder<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  LOG_INFO("LNodeStateBuilder::ScanKey, size: %lu", this->size);
  this->map->ScanKeyHelper(key, this->size, locations_, right_sibling_, result);
}
//===--------------------------------------------------------------------===//
// BWTree Methods
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
bool BWTree<KeyType, ValueType, KeyComparator>::InsertEntry(
    KeyType key, ValueType location) {
  // return false;
  LPID child_lpid;
  LOG_INFO("\nBWTree::InsertEntry, %s", this->ToString(key).c_str());
  bool complete = false;
  // just call InsertEntry on root
  while (!complete) {
    auto epochNum = epoch_manager_.GetCurrentEpoch();
    child_lpid = root_;
    complete = GetMappingTable()
                   ->GetNode(child_lpid)
                   ->InsertEntry(key, location, child_lpid, child_lpid);
    epoch_manager_.ReleaseEpoch(epochNum);
  }
  return true;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction, const KeyType *index_key) {
  std::vector<ValueType> result;
  LOG_INFO("Enter BWTree::Scan");
  // recursive call scan from the root of BWTree
  auto curr_epoch = epoch_manager_.GetCurrentEpoch();
  GetMappingTable()->GetNode(root_)->Scan(values, key_column_ids, expr_types,
                                          scan_direction, result, index_key);
  epoch_manager_.ReleaseEpoch(curr_epoch);

  // reverse the result if scan in backward direction. inefficient
  // implementation
  if (scan_direction == SCAN_DIRECTION_TYPE_BACKWARD) {
    std::reverse(result.begin(), result.end());
  }

  LOG_INFO("Leave BWTree::Scan");
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  auto curr_epoch = epoch_manager_.GetCurrentEpoch();
  GetMappingTable()->GetNode(root_)->ScanAllKeys(result);
  epoch_manager_.ReleaseEpoch(curr_epoch);
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key) {
  // LOG_INFO
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  LOG_INFO("Inside ScanKey of BWTree");
  auto curr_epoch = epoch_manager_.GetCurrentEpoch();
  GetMappingTable()->GetNode(root_)->ScanKey(key, result);
  epoch_manager_.ReleaseEpoch(curr_epoch);
  LOG_INFO("Leaving ScanKey of BWTree");

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void BWTree<KeyType, ValueType, KeyComparator>::Debug() {
  // recursive call Debug from the root of BWTree
  std::string info = GetMappingTable()->GetNode(root_)->Debug(1, root_);
  LOG_INFO("\n%s", info.c_str());
}

template <typename KeyType, typename ValueType, class KeyComparator>
void BWTree<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  // recursive call Check from the root of BWTree
  LOG_INFO("BWTree::BWTreeCheck");
  GetMappingTable()->GetNode(root_)->BWTreeCheck();
}

template <typename KeyType, typename ValueType, class KeyComparator>
template <typename PairSecond>
int BWTree<KeyType, ValueType, KeyComparator>::BinarySearch(
    KeyType key, std::pair<KeyType, PairSecond> *locations, oid_t len) {
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
  if (first != -1) {
    return first;
  } else {
    return -low;
  }
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<oid_t> BWTree<KeyType, ValueType, KeyComparator>::ScanKeyInternal(
    KeyType key, std::pair<KeyType, ValueType> *locations, oid_t size) {
  LOG_INFO("BWTree::ScanKeyInternal");
  assert(locations != nullptr);
  std::vector<oid_t> result;
  // empty LPage
  if (size == 0) {
    return result;
  }
  assert(size > 0);
  // do a binary search on locations to get the key
  int index = BinarySearch(key, locations, size);
  if (index < 0 || size == 0 ||
      (index == 0 && CompareKey(locations[0].first, key) != 0)) {
    // key not found, return empty result
    return result;
  }

  // try to collect all matching keys. If unique_keys, only one key matches
  while (index < size) {
    std::pair<KeyType, ValueType> location = (locations)[index];
    if (CompareKey(location.first, key) == 0) {
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
void BWTree<KeyType, ValueType, KeyComparator>::ScanHelper(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction, const KeyType *index_key,
    std::vector<ValueType> &result, std::pair<KeyType, ValueType> *locations,
    oid_t size, LPID right_sibling) {
  LOG_INFO("BWTree::ScanHelper");
  int index = 0;
  // equality constraint on index_key
  if (index_key != nullptr) {
    index = BinarySearch(*index_key, locations, size);
    LOG_INFO("Scan from position %d, total len: %lu", index, size);
    // key not found.
    if (index < 0) {
      index = -index;
    }

    for (; index < size; index++) {
      std::pair<KeyType, ValueType> pair = locations[index];
      KeyType key = pair.first;
      auto tuple = key.GetTupleForComparison(GetKeySchema());

      // key doesn't match index_key
      if (!MatchLeadingColumn(tuple, key_column_ids, values)) {
        LOG_INFO("Scan key's leading column doesn't match");
        return;
      }
      if (Index::Compare(tuple, key_column_ids, expr_types, values) == true) {
        ItemPointer location = pair.second;
        result.push_back(location);
      }
    }
    // reach the end of current LPage, go to next LPage for more results
    if (index == size && right_sibling != INVALID_LPID) {
      GetMappingTable()
          ->GetNode(right_sibling)
          ->Scan(values, key_column_ids, expr_types, scan_direction, result,
                 index_key);
    }
    return;
  }

  // no constraint for index_key equality check
  for (; index < size; index++) {
    std::pair<KeyType, ValueType> pair = locations[index];
    KeyType key = pair.first;
    auto tuple = key.GetTupleForComparison(GetKeySchema());
    // Compare the current key in the scan with "values" based on "expression
    // types" For instance, "5" EXPR_GREATER_THAN "2" is true
    if (Index::Compare(tuple, key_column_ids, expr_types, values) == true) {
      ItemPointer location = pair.second;
      result.push_back(location);
    }
  }
  // reach the end of current LPage, go to next LPage for more results
  if (index == size && right_sibling != INVALID_LPID) {
    GetMappingTable()
        ->GetNode(right_sibling)
        ->Scan(values, key_column_ids, expr_types, scan_direction, result,
               index_key);
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
void BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeysHelper(
    oid_t size, std::pair<KeyType, ValueType> *locations, oid_t right_sibling,
    std::vector<ValueType> &result) {
  oid_t index;
  for (index = 0; index < size; index++) {
    std::pair<KeyType, ValueType> result_pair = locations[index];
    result.push_back(result_pair.second);
  }
  // reach the end of current LPage, go to next LPage for more results
  if (index == size && right_sibling != INVALID_LPID) {
    GetMappingTable()->GetNode(right_sibling)->ScanAllKeys(result);
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
void BWTree<KeyType, ValueType, KeyComparator>::ScanKeyHelper(
    KeyType key, oid_t size, std::pair<KeyType, ValueType> *locations,
    oid_t right_sibling, std::vector<ValueType> &result) {
  std::vector<oid_t> indices = ScanKeyInternal(key, locations, size);
  // we only need the values
  oid_t index;
  for (index = 0; index < indices.size(); index++) {
    std::pair<KeyType, ValueType> result_pair = (locations[indices[index]]);
    result.push_back(result_pair.second);
  }
  LOG_INFO("BWTree::ScanKeyHelper found %lu results", indices.size());
  // reach the end of current LPage, go to next LPage for more results
  if (indices.size() > 0 && indices[indices.size() - 1] == size - 1 &&
      right_sibling != INVALID_LPID) {
    LOG_INFO("BWTree::ScanKeyHelper go to the right sibling for more result");
    GetMappingTable()->GetNode(right_sibling)->ScanKey(key, result);
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool BWTree<KeyType, ValueType, KeyComparator>::MatchLeadingColumn(
    const AbstractTuple &index_key, const std::vector<oid_t> &key_column_ids,
    const std::vector<Value> &values) {
  int diff;
  oid_t leading_column_id = 0;
  auto key_column_itr = std::find(key_column_ids.begin(), key_column_ids.end(),
                                  leading_column_id) -
                        key_column_ids.begin();
  assert(key_column_itr < key_column_ids.size());

  const Value &rhs = values[key_column_itr];
  const Value &lhs = index_key.GetValue(0);

  diff = lhs.Compare(rhs);
  if (diff == VALUE_COMPARE_EQUAL) {
    return true;
  }
  return false;
}

//===--------------------------------------------------------------------===//
// BWTree Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// IPage Methods Begin
//===--------------------------------------------------------------------===//

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPage<KeyType, ValueType, KeyComparator>::AddINodeEntry(
    LPID self, KeyType max_key_left_split_node,
    KeyType max_key_right_split_node, bool right_node_is_infinity,
    LPID left_split_node_lpid, LPID right_split_node_lpid, bool is_delete) {
  IPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new IPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, max_key_left_split_node, max_key_right_split_node,
          right_node_is_infinity, left_split_node_lpid, right_split_node_lpid,
          is_delete, this->GetRightMostKey(), this->IsInifinity());

  bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);

  if (!status) {
    LOG_INFO("IPage AddInodeEntry failed");
    delete new_delta;
  }
  return status;
}
//
// template <typename KeyType, typename ValueType, class KeyComparator>
// bool IPage<KeyType, ValueType, KeyComparator>::AddINodeSplit(KeyType key,
// LPID value, int modified_index){
//	IPageSplitDelta<KeyType, ValueType, KeyComparator> *new_delta =
//		      new IPageSplitDelta<KeyType, ValueType, KeyComparator>(
//		          this->map, this, key, value, modified_index,
// this->GetRightMostKey(),
//	              this->IsInifinity());
//
//		  bool status = this->map->GetMappingTable()->SwapNode(self, this,
// new_delta);
//
//		  if (!status) {
//		    LOG_INFO("LPage InsertEntry failed");
//		    delete new_delta;
//		  }
//		  return status;
//}

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPage<KeyType, ValueType, KeyComparator>::InsertEntry(
    __attribute__((unused)) KeyType key,
    __attribute__((unused)) ValueType location,
    __attribute__((unused)) LPID self, __attribute__((unused)) LPID parent) {
  /* int i;
   bool last_level_page;
   LPID target_child_lpid;  // TODO: write code to get this -- partially done*/
  if (should_split_) {
    SplitNodes(self, parent);
    return false;
  }
  LOG_INFO("Inside IPage InsertEntry");
  assert(size_ <= IPAGE_ARITY);
  int child_lpid_index = GetChild(key, children_, size_);
  LOG_INFO("Got child_lpid_index as %d", child_lpid_index);
  // LPID child_lpid = GetChild(key, children_, size_);
  LPID child_lpid = children_[child_lpid_index].second;
  return this->map->GetMappingTable()
      ->GetNode(child_lpid)
      ->InsertEntry(key, location, child_lpid, self);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void IPage<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction, std::vector<ValueType> &result,
    const KeyType *index_key) {
  LOG_INFO("Enter IPage::Scan");
  LPID child_id = this->children_[0].second;
  if (index_key != nullptr) {
    // special case
    int child_idx = GetChild(*index_key, children_, size_);
    LOG_INFO("Got child_idx as %d", child_idx);
    child_id = this->children_[child_idx].second;
  }
  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetMappingTable()->GetNode(child_id);
  assert(child != nullptr);
  child->Scan(values, key_column_ids, expr_types, scan_direction, result,
              index_key);
  LOG_INFO("Leave IPage::Scan");
};

template <typename KeyType, typename ValueType, class KeyComparator>
void IPage<KeyType, ValueType, KeyComparator>::ScanAllKeys(
    std::vector<ValueType> &result) {
  LPID child_id = this->children_[0].second;
  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetMappingTable()->GetNode(child_id);
  child->ScanAllKeys(result);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void IPage<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  LOG_INFO("Enter IPage::ScanKey");
  // locate the child who covers the key
  int child_idx = GetChild(key, this->children_, this->size_);
  LOG_INFO("Got child_idx as %d", child_idx);
  LPID child_id = this->children_[child_idx].second;

  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetMappingTable()->GetNode(child_id);
  assert(child != nullptr);
  child->ScanKey(key, result);
  LOG_INFO("Leave IPage::ScanKey");
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
IPage<KeyType, ValueType, KeyComparator>::BuildNodeState(int max_index) {
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder;
  // build node state for IPage
  builder = new INodeStateBuilder<KeyType, ValueType, KeyComparator>(
      children_, max_index == -1 ? size_ : max_index + 1, this->map,
      this->right_most_key, this->infinity);
  return builder;
};
template <typename KeyType, typename ValueType, class KeyComparator>
std::string IPage<KeyType, ValueType, KeyComparator>::Debug(int depth,
                                                            LPID self) {
  std::string blank = this->map->GetStringPrefix("  ", depth) + " ";
  std::string info;
  if (self != INVALID_LPID) {
    std::string prefix = this->map->GetStringPrefix("==", depth);
    info += (prefix + " IPage - ");
    info += "LPID: " + std::to_string(self) + " ";
  } else {
    info += blank + "IPage - ";
  }
  info += " should_split: " + std::to_string(this->should_split_) + " inf: " +
          std::to_string(this->infinity) + " right_most_key: " +
          this->map->ToString(this->right_most_key) + " size: " +
          std::to_string(size_) + "\n" + blank;
  for (oid_t i = 0; i < size_; i++) {
    std::pair<KeyType, LPID> pair = this->children_[i];
    if (i + 1 == size_) {
      info += "?," + std::to_string(pair.second) + "\t";
    } else {
      info += this->map->ToString(pair.first) + "," +
              std::to_string(pair.second) + "\t";
    }
    if (i % 5 == 4) {
      info += "\n" + blank;
    }
  }
  info += "\n";
  for (int i = 0; i < size_; i++) {
    LPID child_id = children_[i].second;
    BWTreeNode<KeyType, ValueType, KeyComparator> *child =
        this->map->GetMappingTable()->GetNode(child_id);
    info += child->Debug(depth + 1, child_id);
  }
  return info;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void IPage<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  LOG_INFO("IPage::BWTreeCheck");

  for (oid_t i = 0; i < size_; i++) {
    std::pair<KeyType, LPID> pair = this->children_[i];
    KeyType key = pair.first;
    // 1. check responsibility
    if (i + 1 < size_) {
      assert(this->map->KeyNotGreaterThan(key, this->right_most_key,
                                          this->infinity));
    }
    // 2. check increasing
    if (i + 2 < size_) {
      assert(this->map->CompareKey(key, this->children_[i + 1].first) < 0);
    }
    // 3. check child's responsibility
    LPID child_id = pair.second;
    BWTreeNode<TEMPLATE_TYPE> *child =
        this->map->GetMappingTable()->GetNode(child_id);
    BWTreeNode<TEMPLATE_TYPE> *page = child->BuildNodeState(-1)->GetPage();
    auto child_right_most_key = page->GetRightMostKey();
    auto child_inf = page->IsInifinity();
    if (i < size_ - 1) {
      assert(this->map->CompareKey(child_right_most_key, key) == 0);
    } else {
      assert(this->map->KeyNotGreaterThan(
          child_right_most_key, this->right_most_key, this->infinity));
    }
  }

  for (oid_t i = 0; i < size_; i++) {
    std::pair<KeyType, LPID> pair = this->children_[i];
    // 4. check children
    LPID child_id = pair.second;
    BWTreeNode<TEMPLATE_TYPE> *child =
        this->map->GetMappingTable()->GetNode(child_id);
    child->BWTreeCheck();
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
// get the index of the child at next level, which contains the given key
int IPage<KeyType, ValueType, KeyComparator>::GetChild(
    KeyType key, std::pair<KeyType, LPID> *children, oid_t len) {
  int index = this->map->BinarySearch(key, children, len - 1);
  return index >= 0 ? index : -index;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void IPage<KeyType, ValueType, KeyComparator>::SplitNodes(LPID self,
                                                          LPID parent) {
  auto currNode = this->map->GetMappingTable()->GetNode(self);
  if (currNode != this) {
    this->map->CompressDeltaChain(self, currNode, currNode);
    return false;
  }

  LPID newIpageLPID;
  int newPageIndex = 0;
  KeyType maxLeftSplitNodeKey;  //, maxRightSplitNodeKey;
  // ValueType leftSplitNodeVal;
  bool swapSuccess;
  LOG_INFO("Splitting Node with LPID: %lu, whose parent is %lu", self, parent);

  LOG_INFO("The size of this node (LPID %lu) is %d", self, (int)size_);
  IPage<KeyType, ValueType, KeyComparator> *newIpage =
      new IPage<KeyType, ValueType, KeyComparator>(
          this->map, this->right_most_key, this->infinity);

  for (int i = size_ / 2 + 1; i < size_; i++) {
    newIpage->children_[newPageIndex++] = children_[i];
  }

  newIpage->size_ = newPageIndex;
  LOG_INFO("The size of the new right split node is %d", (int)newIpage->size_);

  // Assuming we have ( .. ] ranges
  maxLeftSplitNodeKey = children_[size_ / 2].first;
  int modifiedIndex = size_ / 2;

  newIpageLPID = this->map->GetMappingTable()->InstallPage(newIpage);

  LOG_INFO("This newly created right split node got LPID: %d",
           (int)newIpageLPID);

  IPageSplitDelta<KeyType, ValueType, KeyComparator> *splitDelta =
      new IPageSplitDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, maxLeftSplitNodeKey, newIpageLPID, modifiedIndex,
          this->right_most_key, this->infinity);

  swapSuccess = this->map->GetMappingTable()->SwapNode(self, this, splitDelta);

  if (swapSuccess == false) {
    LOG_INFO("This SwapNode attempt for split failed");
    delete splitDelta;
    // delete newIpage;
    return;
  }

  should_split_ = false;

  LOG_INFO("The SwapNode attempt for split succeeded.");
  // This completes the atomic half split

  LOG_INFO("Split page %lu, into new page %lu", self, newIpageLPID);

  if (self != parent)  // Internal IPage
  {
    // Now start with the second half
    LOG_INFO("Now try to create a new IPageUpdateDelta");

    do {
      swapSuccess =
          this->map->GetMappingTable()->GetNode(parent)->AddINodeEntry(
              parent, maxLeftSplitNodeKey, this->GetRightMostKey(),
              this->IsInifinity(), self, newIpageLPID, false);

      // This SwapNode has to be successful. No one else can do this for now.
      // TODO if we allow any node to complete a partial SMO, then this will
      // change

    } while (!swapSuccess);
    LOG_INFO("Returning from the split.");
    splitDelta->SetSplitCompleted();
    return;
  }

  LOG_INFO("This is a split on the root node. Must handle it separately");

  IPage<KeyType, ValueType, KeyComparator> *new_root =
      new IPage<KeyType, ValueType, KeyComparator>(
          this->map, this->right_most_key, this->infinity);

  LPID new_root_LPID;

  new_root_LPID = this->map->GetMappingTable()->InstallPage(new_root);
  LOG_INFO("The new root's LPID is %d", (int)new_root_LPID);

  new_root->size_ = 2;

  new_root->children_[0].first = maxLeftSplitNodeKey;  // key
  new_root->children_[0].second = self;                // LPID

  //  new_root->children_[1].first = maxRightSplitNodeKey;  // key
  new_root->children_[1].second = newIpageLPID;  // LPID

  // First, we must install this new IPage as the new root in this->map
  bool successFlag;
  //  LPID *rootLPIDAddress = this->map->GetRootLPIDAddress();
  //  LPID oldRootLPID = *rootLPIDAddress;

  successFlag = __sync_bool_compare_and_swap(this->map->GetRootLPIDAddress(),
                                             self, new_root_LPID);
  if (!successFlag)  // someone else succeeded
  {
    LOG_INFO("The CAS to change root LPID in the BWTree failed. Return");
    return;
  }
  splitDelta->SetSplitCompleted();
  LOG_INFO(
      "The CAS to change root LPID in the BWTree succeeded. We now have a new "
      "root.");
  return;
};

//===--------------------------------------------------------------------===//
// Delta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
void Delta<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction, std::vector<ValueType> &result,
    const KeyType *index_key) {
  LOG_INFO("Delta::Scan");
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->BuildNodeState();
  auto page = builder->GetPage();
  page->Scan(values, key_column_ids, expr_types, scan_direction, result,
             index_key);
  // release builder
  delete page;
  delete builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void Delta<KeyType, ValueType, KeyComparator>::ScanAllKeys(
    std::vector<ValueType> &result) {
  LOG_INFO("Delta::ScanAllKeys");
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->BuildNodeState();
  builder->ScanAllKeys(result);
  // release builder
  delete (builder);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void Delta<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  LOG_INFO("Delta::ScanKey");
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->BuildNodeState(-1);
  auto page = builder->GetPage();
  assert(builder != nullptr);
  if (!builder->IsSeparated()) {
    // the Page doesn't split, scan left page
    // TODO support ScanKey on IBuilder
    page->ScanKey(key, result);
  } else {
    // the Page splits
    KeyType separator_key = builder->GetSeparatorKey();
    // the scan key is in the left page, scan left page
    // if the scan == separator key, we still scan from the left page since non
    // duplicate key and go across boarder
    if (this->map->CompareKey(separator_key, key) <= 0) {
      page->ScanKey(key, result);
    } else {
      // the desired key is in the right page
      LPID right_page_id = builder->GetSplitNewPageId();
      // scan right page
      this->map->GetMappingTable()
          ->GetNode(right_page_id)
          ->ScanKey(key, result);
    }
  }
  // release builder
  delete page;
  delete builder;
};
//===--------------------------------------------------------------------===//
// Delta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// IPageSplitDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageSplitDelta<KeyType, ValueType, KeyComparator>::AddINodeEntry(
    LPID self, KeyType max_key_left_split_node,
    KeyType max_key_right_split_node, bool right_node_is_infinity,
    LPID left_split_node_lpid, LPID right_split_node_lpid, bool is_delete) {
  // TODO Is this flag useful??
  if (!this->split_completed_) {
    LOG_INFO("Returning early because split in progress");
    return false;
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *old_child_node_hard_ptr =
      this->modified_node;

  if (old_child_node_hard_ptr->GetDeltaChainLen() + 1 >
      IPAGE_DELTA_CHAIN_LIMIT) {
    this->map->CompressDeltaChain(self, this, this);
    return false;
  }

  // This Delta must now be inserted BELOW this delta
  IPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new IPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this->modified_node, max_key_left_split_node,
          max_key_right_split_node, right_node_is_infinity,
          left_split_node_lpid, right_split_node_lpid, is_delete,
          this->GetRightMostKey(), this->IsInifinity());

  bool status = __sync_bool_compare_and_swap(
      &(this->modified_node), old_child_node_hard_ptr, new_delta);
  if (!status) {
    delete new_delta;
  }
  return status;
}

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
IPageSplitDelta<KeyType, ValueType, KeyComparator>::BuildNodeState(
    __attribute__((unused)) int max_index) {
  // Children of IPageDelta always return a INodeStateBuilder
  INodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<INodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState(modified_index_));
  int index =
      this->map->BinarySearch(modified_key_, builder->children_, builder->size);
  assert(index > 0);
  builder->size = index + 1;
  builder->SetInfinity(false);
  builder->SetRightMostKey(modified_key_);
  assert(builder != nullptr);
  return builder;
}

template <typename KeyType, typename ValueType, class KeyComparator>
std::string IPageSplitDelta<KeyType, ValueType, KeyComparator>::Debug(
    int depth, LPID self) {
  std::string blank = this->map->GetStringPrefix("  ", depth) + " ";
  std::string info;
  if (self != INVALID_LPID) {
    std::string prefix = this->map->GetStringPrefix("==", depth);
    info += (prefix + " IPageSplitDelta - ");
    info += "LPID: " + std::to_string(self) + " ";
  } else {
    info += blank + "IPageSplitDelta - ";
  }

  info += "split_completed_: " + std::to_string(split_completed_) +
          " right_most_key: " + this->map->ToString(this->right_most_key) +
          " inf: " + std::to_string(this->infinity) + " index: " +
          std::to_string(modified_index_) + " key: " +
          this->map->ToString(modified_key_) + ", " +
          std::to_string(modified_val_) + " right_lpid: " +
          std::to_string(this->modified_val_) + "\n";
  return info + this->modified_node->Debug(depth, INVALID_LPID);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void IPageSplitDelta<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  LOG_INFO("IPageSplitDelta::BWTreeCheck");

  assert(this->map->KeyNotGreaterThan(modified_key_, this->right_most_key,
                                      this->infinity));
  assert(modified_index_ < IPAGE_ARITY + IPAGE_DELTA_CHAIN_LIMIT);
  assert(modified_index_ > 0);

  // check left page
  BWTreeNode<TEMPLATE_TYPE> *page = BuildNodeState(-1)->GetPage();
  page->BWTreeCheck();
  delete page;

  // check right page
  BWTreeNode<TEMPLATE_TYPE> *right_node =
      this->map->GetMappingTable()->GetNode(modified_val_);
  right_node->BWTreeCheck();
}

template <typename KeyType, typename ValueType, class KeyComparator>
void IPageSplitDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  // Children of IPageDelta always return a INodeStateBuilder
  // this key is greater than split key
  if (this->map->CompareKey(key, modified_key_) > 0) {
    return this->map->GetMappingTable()
        ->GetNode(modified_val_)
        ->ScanKey(key, result);
  }

  return this->modified_node->ScanKey(key, result);
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageSplitDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(
    KeyType key, ValueType location, LPID self, LPID parent) {
  // This function will just call delete entry on the appropriate child

  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }
  if (this->map->CompareKey(key, modified_key_) >
      0)  // this key is greater than split key
  {
    return this->map->GetMappingTable()
        ->GetNode(modified_val_)
        ->DeleteEntry(key, location, modified_val_, parent);
  }

  return this->modified_node->DeleteEntry(key, location, self, parent);
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageSplitDelta<KeyType, ValueType, KeyComparator>::InsertEntry(
    KeyType key, ValueType location, LPID self, LPID parent) {
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }
  // This function will just call delete entry on the appropriate child
  if (this->map->CompareKey(key, modified_key_) >
      0)  // this key is greater than split key
  {
    return this->map->GetMappingTable()
        ->GetNode(modified_val_)
        ->InsertEntry(key, location, modified_val_, parent);
  }

  return this->modified_node->InsertEntry(key, location, self, parent);
};

//===--------------------------------------------------------------------===//
// IPageSplitDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPageSplitDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
void LPageSplitDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  assert(this->modified_node != nullptr);
  LOG_INFO("LPageSplitDelta::ScanKey");

  bool greater_than_left_key = this->map->CompareKey(key, modified_key_) > 0;

  if (greater_than_left_key) {
    LOG_INFO(
        "LPageSplitDelta::ScanKey Found a matching key for right split page");
    this->map->GetMappingTable()
        ->GetNode(right_split_page_lpid_)
        ->ScanKey(key, result);
  } else {
    // Scan the modified node
    auto builder = this->BuildNodeState(-1);
    builder->ScanKey(key, result);
  }
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPageSplitDelta<KeyType, ValueType, KeyComparator>::BuildNodeState(int) {
  // Children of IPageDelta always return a INodeStateBuilder
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      this->modified_node->BuildNodeState(modified_key_index_);
  builder->SetInfinity(false);
  builder->SetRightMostKey(modified_key_);
  assert(builder != nullptr);

  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *lbuilder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          builder);
  lbuilder->UpdateRightSib(right_split_page_lpid_);

  return lbuilder;
}
template <typename KeyType, typename ValueType, class KeyComparator>
std::string LPageSplitDelta<KeyType, ValueType, KeyComparator>::Debug(
    int depth, LPID self) {
  std::string blank = this->map->GetStringPrefix("  ", depth) + " ";
  std::string info;
  if (self != INVALID_LPID) {
    std::string prefix = this->map->GetStringPrefix("==", depth);
    info += (prefix + " LPageSplitDelta - ");
    info += "LPID: " + std::to_string(self) + " ";
  } else {
    info += blank + "LPageSplitDelta - ";
  }

  info += " right_most_key: " + this->map->ToString(this->right_most_key) +
          " inf: " + std::to_string(this->infinity) + "split_completed_: " +
          std::to_string(split_completed_) + " index: " +
          std::to_string(modified_key_index_) + " key: " +
          this->map->ToString(modified_key_) + ", " +
          std::to_string(right_split_page_lpid_) + " right_lpid: " +
          std::to_string(this->right_split_page_lpid_) + "\n";
  return info + this->modified_node->Debug(depth, INVALID_LPID);
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LPageSplitDelta<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  LOG_INFO("LPageSplitDelta::BWTreeCheck");

  assert(this->map->KeyNotGreaterThan(modified_key_, this->right_most_key,
                                      this->infinity));
  assert(modified_key_index_ < LPAGE_ARITY + LPAGE_DELTA_CHAIN_LIMIT);
  assert(modified_key_index_ > 0);

  // check left page
  BWTreeNode<TEMPLATE_TYPE> *page = BuildNodeState(-1)->GetPage();
  page->BWTreeCheck();
  delete page;

  // check right page
  BWTreeNode<TEMPLATE_TYPE> *right_node =
      this->map->GetMappingTable()->GetNode(right_split_page_lpid_);
  right_node->BWTreeCheck();
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPageSplitDelta<KeyType, ValueType, KeyComparator>::InsertEntry(
    KeyType key, ValueType location, LPID self, LPID parent) {
  // if the key falls out of responsible range, retry
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }

  // this key is greater than modified_key_
  if (this->map->CompareKey(key, modified_key_) == 1) {
    return this->map->GetMappingTable()
        ->GetNode(right_split_page_lpid_)
        ->InsertEntry(key, location, right_split_page_lpid_, parent);
  }

  if (!this->split_completed_) {
    LOG_INFO("Returning early because split in progress");
    return false;
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *old_child_node_hard_ptr =
      this->modified_node;

  if (old_child_node_hard_ptr->GetDeltaChainLen() + 1 >
      LPAGE_DELTA_CHAIN_LIMIT) {
    this->map->CompressDeltaChain(self, this, this);
    return false;
  }

  // This Delta must now be inserted BELOW this delta
  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, old_child_node_hard_ptr, key, location,
          this->right_most_key, this->infinity, right_split_page_lpid_);

  bool status = __sync_bool_compare_and_swap(
      &(this->modified_node), old_child_node_hard_ptr, new_delta);
  if (!status) {
    delete new_delta;
  }
  return status;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPageSplitDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(
    KeyType key, ValueType location, LPID self,
    __attribute__((unused)) LPID parent) {
  // if the key falls out of responsible range, retry
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }

  // this key is greater than modified_key_
  if (this->map->CompareKey(key, modified_key_) == 1) {
    return this->map->GetMappingTable()
        ->GetNode(right_split_page_lpid_)
        ->DeleteEntry(key, location, right_split_page_lpid_, parent);
  }

  BWTreeNode<KeyType, ValueType, KeyComparator> *old_child_node_hard_ptr =
      this->modified_node;

  if (!this->split_completed_) {
    LOG_INFO("Returning early because split in progress");
    return false;
  }

  if (old_child_node_hard_ptr->GetDeltaChainLen() + 1 >
      LPAGE_DELTA_CHAIN_LIMIT) {
    this->map->CompressDeltaChain(self, this, this);
    return false;
  }

  LOG_INFO("LPageSplitDelta inserting delete delta below itself");
  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, old_child_node_hard_ptr, key, location,
          this->right_most_key, this->infinity, right_split_page_lpid_);
  new_delta->SetDeleteFlag();

  bool status = __sync_bool_compare_and_swap(
      &(this->modified_node), old_child_node_hard_ptr, new_delta);
  if (!status) {
    delete new_delta;
    return false;
  }

  if (this->right_sibling != INVALID_LPID &&
      this->map->CompareKey(key, this->modified_key_) == 0) {
    LOG_INFO("Cascading deleteEntry to right sib..");
    do {
      BWTreeNode<TEMPLATE_TYPE> *right_node =
          this->map->GetMappingTable()->GetNode(this->right_sibling);
      status = right_node->DeleteEntry(key, location, this->right_sibling,
                                       INVALID_LPID);
    } while (!status);
  }
  return true;
};
//===--------------------------------------------------------------------===//
// LPageSplitDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// IPageUpdateDelta Methods Begin
//===--------------------------------------------------------------------===//

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageUpdateDelta<KeyType, ValueType, KeyComparator>::AddINodeEntry(
    LPID self, KeyType max_key_left_split_node,
    KeyType max_key_right_split_node, bool right_node_is_infinity,
    LPID left_split_node_lpid, LPID right_split_node_lpid, bool is_delete) {
  IPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new IPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, max_key_left_split_node, max_key_right_split_node,
          right_node_is_infinity, left_split_node_lpid, right_split_node_lpid,
          is_delete, this->GetRightMostKey(), this->IsInifinity());
  bool status = this->PerformDeltaInsert(self, new_delta);
  if (!status) {
    delete new_delta;
  }
  return status;
}

// template <typename KeyType, typename ValueType, class KeyComparator>
// bool IPageUpdateDelta<KeyType, ValueType,
// KeyComparator>::AddINodeSplit(KeyType key, LPID value, int modified_index){
//	IPageSplitDelta<KeyType, ValueType, KeyComparator> *new_delta =
//		      new IPageSplitDelta<KeyType, ValueType, KeyComparator>(
//		          this->map, this, key, value, modified_index,
// this->GetRightMostKey(),
//	              this->IsInifinity());
//
//	bool status = this->PerformDeltaInsert(my_lpid, new_delta);
//		    if (!status) {
//		      delete new_delta;
//		    }
//		    return status;
//}

template <typename KeyType, typename ValueType, class KeyComparator>
void IPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  assert(this->modified_node != nullptr);
  LOG_INFO("IPageUpdateDelta::ScanKey");

  bool greater_than_left_key =
      this->map->CompareKey(key, max_key_left_split_node_) > 0;
  bool not_greater_than_right_key =
      this->IsInifinity() ||
      this->map->CompareKey(key, max_key_right_split_node_) <= 0;

  if (!is_delete_) {
    // the scanKey is in the range of the right split page
    if (greater_than_left_key && not_greater_than_right_key) {
      LOG_INFO("IPageUpdateDelta::ScanKey Found a matching key in range");
      this->map->GetMappingTable()
          ->GetNode(right_split_node_lpid_)
          ->ScanKey(key, result);
    } else {
      // ScanKey on modified node
      this->modified_node->ScanKey(key, result);
    }

  } else {
    // TODO handle the merge case
  }
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *IPageUpdateDelta<
    KeyType, ValueType, KeyComparator>::BuildNodeState(int max_index) {
  // Children of IPageDelta always return a INodeStateBuilder
  INodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<INodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState(max_index));
  assert(builder != nullptr);
  // delete delta
  if (is_delete_) {
    if (right_node_is_infinity_ ||
        this->map->CompareKey(this->GetRightMostKey(),
                              max_key_right_split_node_) == 0) {
      builder->RemoveLastChild();
    } else {
      builder->RemoveChild(max_key_right_split_node_);
    }
  } else {
    // insert delta
    if (right_node_is_infinity_ ||
        this->map->CompareKey(this->GetRightMostKey(),
                              max_key_right_split_node_) == 0) {
      builder->ReplaceLastChild(right_split_node_lpid_);
    } else {
      std::pair<KeyType, LPID> right_pair(max_key_right_split_node_,
                                          right_split_node_lpid_);
      builder->AddChild(right_pair);
    }
    if (this->map->CompareKey(this->GetRightMostKey(),
                              max_key_left_split_node_) == 0) {
      builder->ReplaceLastChild(right_split_node_lpid_);
    } else {
      std::pair<KeyType, LPID> left_pair(max_key_left_split_node_,
                                         left_split_node_lpid_);
      builder->AddChild(left_pair);
    }
  }

  return builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::string IPageUpdateDelta<KeyType, ValueType, KeyComparator>::Debug(
    int depth, LPID self) {
  std::string blank = this->map->GetStringPrefix("  ", depth) + " ";
  std::string info;
  if (self != INVALID_LPID) {
    std::string prefix = this->map->GetStringPrefix("==", depth);
    info += (prefix + " IPageUpdateDelta - ");
    info += "LPID: " + std::to_string(self) + " ";
  } else {
    info += blank + "IPageUpdateDelta - ";
  }
  info += " inf: " + std::to_string(this->infinity) + " right_most_key: " +
          this->map->ToString(this->right_most_key) + "is_delete: " +
          std::to_string(is_delete_) + " left: " +
          this->map->ToString(max_key_left_split_node_) + ", " +
          std::to_string(left_split_node_lpid_) + "\t" + "right: " +
          this->map->ToString(max_key_right_split_node_) + ", " +
          std::to_string(right_split_node_lpid_) + "\n";
  info += this->modified_node->Debug(depth, INVALID_LPID);
  info += this->map->GetMappingTable()
              ->GetNode(right_split_node_lpid_)
              ->Debug(depth + 1, right_split_node_lpid_);
  return info;
}
template <typename KeyType, typename ValueType, class KeyComparator>
void IPageUpdateDelta<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  LOG_INFO("IPageUpdateDelta::BWTreeCheck");
  // check the left split key is not greater than the right one
  assert(this->map->KeyNotGreaterThan(max_key_left_split_node_,
                                      this->max_key_right_split_node_,
                                      this->right_node_is_infinity_));
  assert(this->map->KeyNotGreaterThan(max_key_left_split_node_,
                                      this->right_most_key, this->infinity));
  assert(this->map->KeyNotGreaterThan(max_key_right_split_node_,
                                      this->right_most_key, this->infinity));

  BWTreeNode<TEMPLATE_TYPE> *page = BuildNodeState(-1)->GetPage();
  page->BWTreeCheck();
  delete page;
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageUpdateDelta<KeyType, ValueType, KeyComparator>::InsertEntry(
    KeyType key, ValueType location, LPID self, LPID parent) {
  assert(this->GetDeltaChainLen() <= this->GetDeltaChainLimit());
  LOG_INFO("Delta chain len: %d", this->GetDeltaChainLen());
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }
  LOG_INFO("Inside IPageUpdateDelta InsertEntry");
  if (this->map->CompareKey(key, max_key_left_split_node_) <=
      0)  // should go to the underlying Ipage
  {
    return this->modified_node->InsertEntry(key, location, self, parent);
  }
  if (right_node_is_infinity_ ||
      this->map->CompareKey(key, max_key_right_split_node_) <=
          0)  // should go to the right node
  {
    return this->map->GetMappingTable()
        ->GetNode(right_split_node_lpid_)
        ->InsertEntry(key, location, right_split_node_lpid_, self);
  }
  // Less than or equal to max_key_left_split_node_, go directly down
  return this->modified_node->InsertEntry(key, location, self, parent);
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageUpdateDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(
    KeyType key, ValueType location, LPID self, LPID parent) {
  //	if (!this->infinity && this->map->CompareKey(key, this->right_most_key)
  //> 0) {
  //	    return false;
  //	  }
  if (this->map->CompareKey(key, max_key_left_split_node_) <=
      0)  // should go to the right split child
  {
    return this->modified_node->DeleteEntry(key, location, self, parent);
  }
  if (right_node_is_infinity_ ||
      this->map->CompareKey(key, max_key_right_split_node_) <=
          0)  // should go down to the lower level IPage
  {
    return this->map->GetMappingTable()
        ->GetNode(right_split_node_lpid_)
        ->DeleteEntry(key, location, right_split_node_lpid_, self);
  }
  // Less than or equal to max_key_left_split_node_, go directly down
  return this->modified_node->DeleteEntry(key, location, self, parent);
}

//===--------------------------------------------------------------------===//
// IPageUpdateDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPageUpdateDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
void LPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  assert(this->modified_node != nullptr);
  LOG_INFO("LPageUpdateDelta::ScanKey");
  // if unique keys, just walk down the delta chain without building node state
  if (this->map->unique_keys) {
    // the modified key matches the scanKey
    if (this->map->CompareKey(modified_key_, key) == 0) {
      LOG_INFO("LPageUpdateDelta::ScanKey Found a matching key");
      if (!is_delete_) {
        // the modified key is insert, add to result vector
        LOG_INFO("LPageUpdateDelta::ScanKey Inserts the matching item pointer");
        result.push_back(modified_val_);
        return;
      }
    }
  }
  LOG_INFO("Building NodeState of all children nodes");

  // we have to build the state
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder;
  builder = this->BuildNodeState(-1);

  assert(builder != nullptr);
  LOG_INFO("Scan on the new NodeState");
  builder->ScanKey(key, result);
  // release builder
  delete (builder);
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *LPageUpdateDelta<
    KeyType, ValueType, KeyComparator>::BuildNodeState(int max_index) {
  LOG_INFO("LPageUpdateDelta::BuildNodeState");
  // Children of LPageDelta always return a LNodeStateBuilder
  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState(max_index));
  assert(builder != nullptr);
  // delete delta
  if (is_delete_) {
    std::pair<KeyType, ValueType> pair(modified_key_, modified_val_);
    builder->RemoveLeafData(pair);
  } else {
    // insert delta
    std::pair<KeyType, ValueType> pair(modified_key_, modified_val_);
    builder->AddLeafData(pair);
  }
  return builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::string LPageUpdateDelta<KeyType, ValueType, KeyComparator>::Debug(
    int depth, LPID self) {
  std::string blank = this->map->GetStringPrefix("  ", depth) + " ";
  std::string info;
  if (self != INVALID_LPID) {
    std::string prefix = this->map->GetStringPrefix("==", depth);
    info += (prefix + " LPageUpdateDelta - ");
    info += "LPID: " + std::to_string(self) + " ";
  } else {
    info += blank + "LPageUpdateDelta - ";
  }
  info += " inf: " + std::to_string(this->infinity) + " right_most_key: " +
          this->map->ToString(this->right_most_key) + "is_delete: " +
          std::to_string(is_delete_) + " " +
          this->map->ToString(modified_key_) + ", " +
          this->map->ToString(modified_val_) + "\n";
  return info + this->modified_node->Debug(depth, INVALID_LPID);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LPageUpdateDelta<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  // check the modified key is not greater than the right most one
  LOG_INFO("LPageUpdateDelta::BWTreeCheck");
  assert(this->map->KeyNotGreaterThan(modified_key_, this->right_most_key,
                                      this->infinity));
  BWTreeNode<TEMPLATE_TYPE> *page = BuildNodeState(-1)->GetPage();
  page->BWTreeCheck();
  delete page;
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPageUpdateDelta<KeyType, ValueType, KeyComparator>::InsertEntry(
    KeyType key, ValueType location, LPID my_lpid,
    __attribute__((unused)) LPID parent) {
  // if the key falls out of responsible range, retry
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }

  if (should_split_) {
    this->map->CompressDeltaChain(my_lpid, this, this);
    return false;
  }

  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, key, location, this->right_most_key, this->infinity,
          this->right_sibling);
  bool status = this->PerformDeltaInsert(my_lpid, new_delta);
  if (!status) {
    delete new_delta;
  }
  return status;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPageUpdateDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(
    KeyType key, ValueType location, LPID my_lpid,
    __attribute__((unused)) LPID parent) {
  // if the key falls out of responsible range, retry
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }
  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, key, location, this->right_most_key, this->infinity,
          this->right_sibling);

  new_delta->SetDeleteFlag();
  new_delta->SetShouldSplit(should_split_);
  bool status = this->PerformDeltaInsert(my_lpid, new_delta);
  if (!status) {
    delete new_delta;
    return false;
  }
  // cascade delete
  if (this->right_sibling != INVALID_LPID &&
      this->map->CompareKey(key, this->right_most_key) == 0) {
    LOG_INFO("Cascading deleteEntry to right sib..");
    do {
      BWTreeNode<TEMPLATE_TYPE> *right_node =
          this->map->GetMappingTable()->GetNode(this->right_sibling);
      status = right_node->DeleteEntry(key, location, this->right_sibling,
                                       INVALID_LPID);
    } while (!status);
  }
  return true;
};

//===--------------------------------------------------------------------===//
// LPageUpdateDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPage Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
void LPage<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction, std::vector<ValueType> &result,
    const KeyType *index_key) {
  LOG_INFO("LPage::Scan");
  this->map->ScanHelper(values, key_column_ids, expr_types, scan_direction,
                        index_key, result, locations_, size_, right_sib_);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LPage<KeyType, ValueType, KeyComparator>::ScanAllKeys(
    std::vector<ValueType> &result) {
  LOG_INFO("LPage::ScanAllKeys");
  this->map->ScanAllKeysHelper(size_, locations_, right_sib_, result);
};

template <typename KeyType, typename ValueType, class KeyComparator>
void LPage<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key, std::vector<ValueType> &result) {
  LOG_INFO("LPage::ScanKey");
  this->map->ScanKeyHelper(key, size_, locations_, right_sib_, result);
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::string LPage<KeyType, ValueType, KeyComparator>::Debug(int depth,
                                                            LPID self) {
  std::string blank = this->map->GetStringPrefix("  ", depth) + " ";
  std::string info;
  if (self != INVALID_LPID) {
    std::string prefix = this->map->GetStringPrefix("==", depth);
    info += (prefix + " LPage - ");
    info += "LPID: " + std::to_string(self) + " ";
  } else {
    info += blank + "LPage - ";
  }
  info += " inf: " + std::to_string(this->infinity) + " right_most_key: " +
          this->map->ToString(this->right_most_key) + " size: " +
          std::to_string(size_) + " right_sib: " + std::to_string(right_sib_) +
          "\n" + blank;
  for (oid_t i = 0; i < size_; i++) {
    std::pair<KeyType, ValueType> pair = this->locations_[i];
    info += this->map->ToString(pair.first) + "," +
            this->map->ToString(pair.second) + "\t";
    if (i % 5 == 4) {
      info += "\n" + blank;
    }
  }
  info += "\n";
  return info;
}

template <typename KeyType, typename ValueType, class KeyComparator>
void LPage<KeyType, ValueType, KeyComparator>::BWTreeCheck() {
  LOG_INFO("LPage::BWTreeCheck");
  for (oid_t i = 0; i < size_; i++) {
    KeyType key = this->locations_[i].first;
    // 1. check responsibility
    if (i + 1 < size_) {
      assert(this->map->KeyNotGreaterThan(key, this->right_most_key,
                                          this->infinity));
    }
    // 2. check non decreasing
    if (i + 2 < size_) {
      assert(this->map->CompareKey(key, this->locations_[i + 1].first) <= 0);
    }
  }
}

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPage<KeyType, ValueType, KeyComparator>::BuildNodeState(int max_index) {
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder;
  // build node state for LPage
  builder = new LNodeStateBuilder<KeyType, ValueType, KeyComparator>(
      left_sib_, right_sib_, locations_,
      max_index == -1 ? size_ : max_index + 1, this->map, this->right_most_key,
      this->infinity);
  return builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPage<KeyType, ValueType, KeyComparator>::InsertEntry(KeyType key,
                                                           ValueType location,
                                                           LPID self,
                                                           LPID parent) {
  if (this->size_ > LPAGE_SPLIT_THRESHOLD) {
    this->SplitNodes(self, parent);
    return false;
  }
  LOG_INFO("Inside LPage InsertEntry. It's size is %d", (int)this->size_);
  // if the key falls out of responsible range, retry
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }
  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, key, location, this->right_most_key, this->infinity,
          right_sib_);

  bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);

  if (!status) {
    LOG_INFO("LPage InsertEntry failed");
    delete new_delta;
  }
  return status;
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPage<KeyType, ValueType, KeyComparator>::DeleteEntry(KeyType key,
                                                           ValueType location,
                                                           LPID self,
                                                           LPID parent) {
  LOG_INFO("LPage::DeleteEntry, %lu, %lu", self, parent);
  bool should_split = this->size_ > LPAGE_SPLIT_THRESHOLD;
  if (should_split) {
    if (parent != INVALID_LPID) {
      this->SplitNodes(self, parent);
      return false;
    }
  }
  // if the key falls out of responsible range, retry
  if (!this->infinity && this->map->CompareKey(key, this->right_most_key) > 0) {
    return false;
  }
  LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
      new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, key, location, this->right_most_key, this->infinity,
          right_sib_);

  new_delta->SetDeleteFlag();
  new_delta->SetShouldSplit(should_split);
  bool status = this->map->GetMappingTable()->SwapNode(self, this, new_delta);
  if (!status) {
    delete new_delta;
    return false;
  }
  if (right_sib_ != INVALID_LPID &&
      this->map->CompareKey(key, this->right_most_key) == 0) {
    LOG_INFO("Cascading deleteEntry to right sib..");
    do {
      BWTreeNode<TEMPLATE_TYPE> *right_node =
          this->map->GetMappingTable()->GetNode(right_sib_);
      status = right_node->DeleteEntry(key, location, right_sib_, INVALID_LPID);
    } while (!status);
  }
  return true;
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPage<KeyType, ValueType, KeyComparator>::SplitNodes(LPID self,
                                                          LPID parent) {
  LPID newLpageLPID, left_page_lpid = self;
  int newPageIndex = 0;
  KeyType maxLeftSplitNodeKey;  //, maxRightSplitNodeKey;
  int leftSplitNodeIndex;
  bool swapSuccess;

  LOG_INFO("Splitting Node with LPID: %lu, whose parent is %lu", self, parent);

  LOG_INFO("The size of this node (LPID %lu) is %lu", self, size_);
  LPage<KeyType, ValueType, KeyComparator> *newLpage =
      new LPage<KeyType, ValueType, KeyComparator>(
          this->map, this->right_most_key, this->infinity);

  for (int i = size_ / 2 + 1; i < size_; i++) {
    newLpage->locations_[newPageIndex++] = locations_[i];
  }

  newLpage->size_ = newPageIndex;
  LOG_INFO("The size of the new right split node is %lu", newLpage->size_);

  // TODO left_sib is set to self
  newLpage->right_sib_ = right_sib_;

  // Assuming we have ( .. ] ranges
  maxLeftSplitNodeKey = locations_[size_ / 2].first;
  // this ensures that the split is inclusive
  leftSplitNodeIndex = size_ / 2;

  //  maxRightSplitNodeKey = locations_[size_ - 1].first;

  //  assert(this->map->CompareKey(maxLeftSplitNodeKey, maxRightSplitNodeKey) !=
  //  1);

  for (int i = 0; i <= size_ / 2; i++) {
    assert(this->map->CompareKey(locations_[i].first, maxLeftSplitNodeKey) <=
           0);
    //    assert(this->map->CompareKey(locations_[i].first,
    //    maxRightSplitNodeKey) <=
    //           0);
  }

  for (int i = 0; i < newLpage->GetSize(); i++) {
    assert(this->map->CompareKey(newLpage->GetLocationsArray()[i].first,
                                 maxLeftSplitNodeKey) != -1);
    //    assert(this->map->CompareKey(newLpage->GetLocationsArray()[i].first,
    //                                 maxRightSplitNodeKey) != 1);
  }

  newLpageLPID = this->map->GetMappingTable()->InstallPage(newLpage);

  LOG_INFO("This newly created right split node got LPID: %lu", newLpageLPID);

  LPageSplitDelta<KeyType, ValueType, KeyComparator> *splitDelta =
      new LPageSplitDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, maxLeftSplitNodeKey, leftSplitNodeIndex,
          newLpageLPID, this->right_most_key, this->infinity, newLpageLPID);

  swapSuccess = this->map->GetMappingTable()->SwapNode(self, this, splitDelta);

  if (swapSuccess == false) {
    LOG_INFO("This SwapNode attempt for split failed");
    delete splitDelta;
    // delete newLpage;
    // What should we do on failure? This means that someone else succeeded in
    // doing the
    // atomic half split. Now if try and install our own Insert / Delete /
    // Update delta, it
    // will automatically be created on top of the LPageSplitDelta
    return false;
  }

  LOG_INFO("The SwapNode attempt for split succeeded.");

  // This completes the atomic half split
  // At this point no one else can succeed with the complete split because
  // this guy won in the half split
  // Now we still have to update the size field and the right_sib of this
  // node... how can we do it atomically?
  // Edit: No need to do that! Because the consolidation will do that.

  LOG_INFO("Split page %lu, into new page %lu", self, newLpageLPID);

  // Now start with the second half
  LOG_INFO("Now try to create a new IPageUpdateDelta");

  do {
    swapSuccess = this->map->GetMappingTable()->GetNode(parent)->AddINodeEntry(
        parent, maxLeftSplitNodeKey, this->GetRightMostKey(),
        this->IsInifinity(), left_page_lpid, newLpageLPID, false);

    // This SwapNode has to be successful. No one else can do this for now.
    // TODO if we allow any node to complete a partial SMO, then this will
    // change

  } while (!swapSuccess);
  assert(swapSuccess == true);

  splitDelta->SetSplitCompleted();
  LOG_INFO("Split finished");
  return true;
}

//===--------------------------------------------------------------------===//
// LPage Methods End
//===--------------------------------------------------------------------===//

}  // End index namespace
}  // End peloton namespace
