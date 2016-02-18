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
  assert(mapping_table_.size() > 0);

  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result =
      GetNode(root_)->Scan(values, key_column_ids, expr_types, scan_direction);

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  assert(mapping_table_.size() > 0);
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result = GetNode(root_)->ScanAllKeys();

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key) {
  assert(mapping_table_.size() > 0);
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result = GetNode(root_)->ScanKey(key);

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool BWTree<KeyType, ValueType, KeyComparator>::InsertEntry(KeyType key, ValueType location) {
	//just call InsertEntry on root

	return GetNode(root_)->InsertEntry(key, location);
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
    KeyType key) {
  assert(key != nullptr);

  std::vector<ValueType> result;

  BWTreeNode<KeyType, ValueType, KeyComparator> *child = GetChild(key);
  if (child == nullptr) {
    // Key is not included in the tree, do nothing
  } else {
    result = child->ScanKey(key);
  }
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool
IPage<KeyType, ValueType, KeyComparator>::InsertEntry(
		const storage::Tuple *key, const ItemPointer location) {
	//TODO Now call InsertEntry on the appropriate child (also note that the LPAGE
	//has no InsertEntry)
	//Shouldn't this node know that it is the last level IPAGE, and perform the
	//delta insert?

	int i;
	bool last_level_page;
	LPID target_child_lpid; //TODO: write code to get this -- partially done

	for (i = 0; i < children_map.size(); i++)
		if (KeyComparator(key, children_map[i].first))
			continue;
		else
		{
			target_child_lpid = children_map[i].second;
			break;
		}

	// Ideally, i should never be equal to children_map.size(), because that would mean somewhere
	// the ranges are not correct
	assert(i != children_map.size());

	if (last_level_ipage)
	{
		LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
				new LPageUpdateDelta<KeyType, ValueType, KeyComparator>();

		new_delta->modified_node = GetNode(target_child_lpid);
		new_delta->modified_key_ = key;
		new_delta->modified_val_ = location;
		//TODO: handle failed SwapNode
		return SwapNode(target_child_lpid, new_delta->modified_node, new_delta);
	}
	else
		return GetNode(target_child_lpid).InsertEntry(key, location);

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
  assert(this->unique_keys == true);
  // TODO implement this
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
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
  std::vector<ValueType> result;
  if (this->unique_keys) {
    // the modified key matches the scanKey
    if (this->comparator(modified_key_, key) == true) {
      // the modified key is deleted, return empty result
      if (modified_val_ == INVALID_ITEMPOINTER) {
        // do nothing
      } else {
        result.push_back(modified_val_);
      }
    }
  } else {
    // TODO pass down the delta info to delete the key at bottom LPage
  }

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
  // in the first version, the LPage has no content at all
  assert(size_ == 0);
  std::vector<ValueType> result;
  return result;
};
}  // End index namespace
}  // End peloton namespace
