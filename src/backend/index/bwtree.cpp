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
  int index = this->map->BinarySearch(key, children_, this->size);
  assert(index < IPAGE_ARITY + IPAGE_DELTA_CHAIN_LIMIT);
  // Key not found. shift every element to the right
  if (index < 0 || this->size == 0 ||
      (index == 0 && CompareKey(children_[0].first, key) != 0)) {
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
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::RemoveChild(
    KeyType &key_to_remove) {
  assert(children_ != nullptr);
  int index = this->map->BinarySearch(key_to_remove, children_, this->size);
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
void INodeStateBuilder<KeyType, ValueType, KeyComparator>::SeparateFromKey(
    KeyType separator_key, LPID split_new_page_id) {
  assert(children_ != nullptr);
  int index = this->map->BinarySearch(separator_key, children_, this->size);
  assert(index < this->size && index >= 0);
  // assume we include the key at the split page
  // decrement size
  this->size = index + 1;
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

  LOG_INFO(
      "LNodeStateBuilder::AddLeafData - BinarySearch (%lu, %lu) result %d, "
      "original size: %lu",
      new_pair.second.block, new_pair.second.offset, index, this->size);

  assert(index < LPAGE_ARITY + LPAGE_DELTA_CHAIN_LIMIT);
  // not found. shift every element to the right
  if (index < 0 || this->size == 0 ||
      (index == 0 && this->map->CompareKey(locations_[0].first, key) != 0)) {
    LOG_INFO("LNodeStateBuilder::AddLeafData - key not found");
    index = -1 * index;
    for (int i = this->size; i > index; i--) {
      locations_[i] = locations_[i - 1];
    }
    this->size++;
  } else {
    LOG_INFO("LNodeStateBuilder::AddLeafData - key found");
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
        LOG_INFO(
            "LNodeStateBuilder::RemoveLeafData - found exact entry at index %d",
            src);
        found_exact_entry_count++;
      } else {
        // value doesn't match
        LOG_INFO("LNodeStateBuilder::RemoveLeafData - not exact entry");
        locations_[dest++] = pair;
      }
    } else {
      key_matches = false;
      LOG_INFO("LNodeStateBuilder::RemoveLeafData - key not match");
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
    KeyType separator_key, ValueType location, LPID split_new_page_id) {
  assert(locations_ != nullptr);

  int index = this->map->BinarySearch(
      std::pair<KeyType, ValueType>(separator_key, location), locations_,
      this->size);
  assert(index < this->size && index >= 0);
  // assume we include the key at the split page
  // decrement size
  this->size = index + 1;
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

  while (!GetMappingTable()
              ->GetNode(child_lpid)
              ->InsertEntry(key, location, child_lpid))
    ;
  return true;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::Scan(
    const std::vector<Value> &values, const std::vector<oid_t> &key_column_ids,
    const std::vector<ExpressionType> &expr_types,
    const ScanDirectionType &scan_direction) {
  std::vector<ValueType> result;
  LOG_INFO("Enter BWTree::Scan");
  // recursive call scan from the root of BWTree
  result = GetMappingTable()->GetNode(root_)->Scan(values, key_column_ids,
                                                   expr_types, scan_direction);

  LOG_INFO("Leave BWTree::Scan");
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
BWTree<KeyType, ValueType, KeyComparator>::ScanAllKeys() {
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  result = GetMappingTable()->GetNode(root_)->ScanAllKeys();

  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType> BWTree<KeyType, ValueType, KeyComparator>::ScanKey(
    KeyType key) {
  // LOG_INFO
  std::vector<ValueType> result;

  // recursive call scan from the root of BWTree
  LOG_INFO("Inside ScanKey of BWTree");
  result = GetMappingTable()->GetNode(root_)->ScanKey(key);
  LOG_INFO("Leaving ScanKey of BWTree");

  return result;
};

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
  // LPID child_lpid = GetChild(key, children_, size_);
  LPID child_lpid = children_[child_lpid_index].second;
  return this->map->GetMappingTable()
      ->GetNode(child_lpid)
      ->InsertEntry(key, location, child_lpid);
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
  LOG_INFO("Enter IPage::ScanKey");
  std::vector<ValueType> result;
  // locate the child who covers the key
  int child_idx = GetChild(key, this->children_, this->size_);
  LOG_INFO("Got child_idx as %d", child_idx);
  LPID child_id = this->children_[child_idx].second;

  BWTreeNode<KeyType, ValueType, KeyComparator> *child =
      this->map->GetMappingTable()->GetNode(child_id);
  assert(child != nullptr);
  result = child->ScanKey(key);
  LOG_INFO("Leave IPage::ScanKey");
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
    KeyType key, std::pair<KeyType, LPID> *children, oid_t len) {
  int index = this->map->BinarySearch(key, children, len - 1);
  return index >= 0 ? index : -index;
};

template <typename KeyType, typename ValueType, class KeyComparator>
void IPage<KeyType, ValueType, KeyComparator>::SplitNodes(LPID self,
                                                          LPID parent) {
  LPID newIpageLPID;
  int newPageIndex = 0;
  KeyType maxLeftSplitNodeKey, maxRightSplitNodeKey;
  // ValueType leftSplitNodeVal;
  bool swapSuccess;

  LOG_INFO("Splitting Node with LPID: %lu, whose parent is %lu", self, parent);

  LOG_INFO("The size of this node (LPID %lu) is %d", self, size_);
  IPage<KeyType, ValueType, KeyComparator> *newIpage =
      new IPage<KeyType, ValueType, KeyComparator>(this->map);

  for (int i = size_ / 2 + 1; i < size_; i++) {
    newIpage->children_[newPageIndex++] = children_[i];
  }

  newIpage->size_ = newPageIndex;
  LOG_INFO("The size of the new right split node is %d", newIpage->size_);

  // Assuming we have ( .. ] ranges
  maxLeftSplitNodeKey = children_[size_ / 2].first;
  maxRightSplitNodeKey = children_[size_ - 1].first;

  newIpageLPID = this->map->GetMappingTable()->InstallPage(newIpage);

  LOG_INFO("This newly created right split node got LPID: %d", newIpageLPID);

  IPageSplitDelta<KeyType, ValueType, KeyComparator> *splitDelta =
      new IPageSplitDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, maxLeftSplitNodeKey, newIpageLPID);

  swapSuccess = this->map->GetMappingTable()->SwapNode(self, this, splitDelta);

  if (swapSuccess == false) {
    LOG_INFO("This SwapNode attempt for split failed");
    delete splitDelta;
    delete newIpage;
    // What should we do on failure? This means that someone else succeeded in
    // doing the
    // atomic half split. Now if try and install our own Insert / Delete /
    // Update delta, it
    // will automatically be created on top of the LPageSplitDelta
    return;
  }

  LOG_INFO("The SwapNode attempt for split succeeded.");
  // This completes the atomic half split
  // At this point no one else can succeed with the complete split because
  // this guy won in the half split
  // Now we still have to update the size field and the right_sib of this
  // node... how can we do it atomically?
  // Edit: No need to do that! Because the consolidation will do that.

  LOG_INFO("Split page %lu, into new page %lu", self, newIpageLPID);

  if (self != parent) //Internal IPage
  {
	  // Now start with the second half
	  LOG_INFO("Now try to create a new IPageUpdateDelta");

	  BWTreeNode<KeyType, ValueType, KeyComparator> *parentHardPtr;

  parentHardPtr = this->map->GetMappingTable()->GetNode(parent);

	  IPageUpdateDelta<KeyType, ValueType, KeyComparator> *parentUpdateDelta =
	      new IPageUpdateDelta<KeyType, ValueType, KeyComparator>(
	          this->map, parentHardPtr, maxLeftSplitNodeKey, maxRightSplitNodeKey,
	          newIpageLPID);

	  LOG_INFO("Now Doing a SwapNode to install this IPageUpdateDelta");

  // This SwapNode has to be successful. No one else can do this for now.
  // TODO if we allow any node to complete a partial SMO, then this will change
  swapSuccess = this->map->GetMappingTable()->SwapNode(parent, parentHardPtr,
                                                       parentUpdateDelta);
  }


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
  LOG_INFO(" ");
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
      result =
          this->map->GetMappingTable()->GetNode(right_page_id)->ScanKey(key);
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
  std::pair<KeyType, LPID> splitter(modified_key_, modified_val_);
  builder->SeparateFromKey(splitter);

  return builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageSplitDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(KeyType key,
		   ValueType location, LPID self)
	 {
		  //This function will just call delete entry on the appropriate child
		  if (this->map->CompareKey(key, modified_key_) > 1) //this key is greater than split key
		  {
			  return this->map->GetNode(modified_val_)->DeleteEntry(key, location, modified_val_);
		  }

		  return this->modified_node->DeleteEntry(key, location, self);
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageSplitDelta<KeyType, ValueType, KeyComparator>::InsertEntry(KeyType key,
		   ValueType location, LPID self)
	 {
		  //This function will just call delete entry on the appropriate child
		  if (this->map->CompareKey(key, modified_key_) > 1) //this key is greater than split key
		  {
			  return this->map->GetNode(modified_val_)->InsertEntry(key, location, modified_val_);
		  }

		  return this->modified_node->InsertEntry(key, location, self);
};

//===--------------------------------------------------------------------===//
// IPageSplitDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPageSplitDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
LPageSplitDelta<KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
  std::vector<ValueType> result;
  assert(this->modified_node != nullptr);
  LOG_INFO("LPageSplitDelta::ScanKey");

  bool greater_than_left_key = this->map->CompareKey(key, modified_key_) > 0;

  if (greater_than_left_key) {
    LOG_INFO(
        "LPageSplitDelta::ScanKey Found a matching key for right split page");
    result = this->map->GetMappingTable()
                 ->GetNode(right_split_page_lpid_)
                 .ScanKey(key);
  } else {
    // Scan the modified node
    result = this->modified_node->ScanKey(key);
  }
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPageSplitDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  // Children of IPageDelta always return a INodeStateBuilder
  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
  assert(builder != nullptr);
  builder->SeparateFromKey(modified_key_, modified_key_location_,
                           right_split_page_lpid_);

  return builder;
}

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPageSplitDelta<KeyType, ValueType, KeyComparator>::InsertEntry(KeyType key,
		ValueType location, LPID self) {

	if (this->map->CompareKey(key, modified_key_) == 1) //this key is greater than modified_key_
	{
		return this->map->GetNode(right_split_page_lpid_)->InsertEntry(key, location, right_split_page_lpid_);
	}

	LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
			new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
					key, location);
	bool status = this->map->SwapNode(self, this, new_delta);
	if (!status) {
		delete new_delta;
	}
	return status;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool LPageSplitDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(KeyType key,
		ValueType location, LPID self) {

	if (this->map->CompareKey(key, modified_key_) == 1) //this key is greater than modified_key_
	{
		return this->map->GetNode(right_split_page_lpid_)->InsertEntry(key, location, right_split_page_lpid_);
	}

	LPageUpdateDelta<KeyType, ValueType, KeyComparator> *new_delta =
			new LPageUpdateDelta<KeyType, ValueType, KeyComparator>(this->map, this,
					key, location);
	new_delta->SetDeleteFlag();
	bool status = this->map->SwapNode(self, this, new_delta);
	if (!status) {
		delete new_delta;
	}
	return status;
};
//===--------------------------------------------------------------------===//
// LPageSplitDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// IPageUpdateDelta Methods Begin
//===--------------------------------------------------------------------===//
template <typename KeyType, typename ValueType, class KeyComparator>
std::vector<ValueType>
IPageUpdateDelta<KeyType, ValueType, KeyComparator>::ScanKey(KeyType key) {
  std::vector<ValueType> result;
  assert(this->modified_node != nullptr);
  LOG_INFO("IPageUpdateDelta::ScanKey");

  bool greater_than_left_key =
      this->map->CompareKey(key, max_key_left_split_node_) > 0;
  bool not_greater_than_right_key =
      this->map->CompareKey(key, max_key_right_split_node_) <= 0;

  if (!is_delete_) {
    // the scanKey is in the range of the right split page
    if (greater_than_left_key && not_greater_than_right_key) {
      LOG_INFO("IPageUpdateDelta::ScanKey Found a matching key in range");
      result = this->map->GetMappingTable()
                   ->GetNode(right_split_node_lpid_)
                   .ScanKey(key);
    } else {
      // ScanKey on modified node
      result = this->modified_node->ScanKey(key);
    }

  } else {
    // TODO handle the merge case
  }
  return result;
};

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
    builder->RemoveChild(max_key_right_split_node_);
  } else {
    // insert delta
    std::pair<KeyType, LPID> right_pair(max_key_right_split_node_,
                                        right_split_node_lpid_);
    builder->AddChild(right_pair);

    std::pair<KeyType, LPID> left_pair(max_key_left_split_node_,
                                       right_split_node_lpid_);
    builder->AddChild(left_pair);
  }
  return builder;
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageUpdateDelta<KeyType, ValueType, KeyComparator>::InsertEntry(KeyType key,
		ValueType location, LPID self) {
	  if (this->map->CompareKey(key, max_key_right_split_node_) == 1) //should go down to the lower level IPage
	  {
		  return this->modified_node->InsertEntry(key, location, self);
	  }
	  if (this->map->CompareKey(key, max_key_left_split_node_) == 1) //should go to the right split child
	  {
		  return this->map->GetNode(right_split_node_lpid_)->InsertEntry(key, location, right_split_node_lpid_);
	  }
	  // Less than or equal to max_key_left_split_node_, go directly down
	  return this->modified_node->InsertEntry(key, location, self);
};

template <typename KeyType, typename ValueType, class KeyComparator>
bool IPageUpdateDelta<KeyType, ValueType, KeyComparator>::DeleteEntry(KeyType key,
		 ValueType location, LPID self) {
	  if (this->map->CompareKey(key, max_key_right_split_node_) == 1) //should go down to the lower level IPage
	  {
		  return this->modified_node->DeleteEntry(key, location, self);
	  }
	  if (this->map->CompareKey(key, max_key_left_split_node_) == 1) //should go to the right split child
	  {
		  return this->map->GetNode(right_split_node_lpid_)->DeleteEntry(key, location, right_split_node_lpid_);
	  }
	  // Less than or equal to max_key_left_split_node_, go directly down
	  return this->modified_node->DeleteEntry(key, location, self);
 }

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
      }
      // if the modified key is delete, then we have to build the state to see
      // if the delete<k, v> is valid.
      return result;
    }
  }
  LOG_INFO("Building NodeState of all children nodes");

  // non unique key. we have to build the state
  NodeStateBuilder<KeyType, ValueType, KeyComparator> *builder;

  if (this->map->unique_keys) {
    builder = this->modified_node->BuildNodeState();
  } else {
    builder = this->BuildNodeState();
  }

  assert(builder != nullptr);
  BWTreeNode<KeyType, ValueType, KeyComparator> *page = builder->GetPage();
  assert(page != nullptr);
  // TODO check split status of builder for shortcut
  // do scan on the new state
  LOG_INFO("Scan on the new NodeState");
  result = page->ScanKey(key);
  // release builder
  delete page;
  delete (builder);
  return result;
};

template <typename KeyType, typename ValueType, class KeyComparator>
NodeStateBuilder<KeyType, ValueType, KeyComparator> *
LPageUpdateDelta<KeyType, ValueType, KeyComparator>::BuildNodeState() {
  LOG_INFO("LPageUpdateDelta::BuildNodeState");
  // Children of LPageDelta always return a LNodeStateBuilder
  LNodeStateBuilder<KeyType, ValueType, KeyComparator> *builder =
      reinterpret_cast<LNodeStateBuilder<KeyType, ValueType, KeyComparator> *>(
          this->modified_node->BuildNodeState());
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
//===--------------------------------------------------------------------===//
// LPageUpdateDelta Methods End
//===--------------------------------------------------------------------===//

//===--------------------------------------------------------------------===//
// LPage Methods Begin
//===--------------------------------------------------------------------===//

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
    std::pair<KeyType, ValueType> result_pair = (locations_[indices[index]]);
    result.push_back(result_pair.second);
  }
  // reach the end of current LPage, go to next LPage for more results
  if (index == size_ && right_sib_ != INVALID_LPID) {
    std::vector<ValueType> sib_result =
        this->map->GetMappingTable()->GetNode(right_sib_)->ScanKey(key);
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
  if (index < 0) {
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

template <typename KeyType, typename ValueType, class KeyComparator>
void LPage<KeyType, ValueType, KeyComparator>::SplitNodes(LPID self,
                                                          LPID parent) {
  LPID newLpageLPID, left_page_lpid = self;
  int newPageIndex = 0;
  KeyType maxLeftSplitNodeKey, maxRightSplitNodeKey;
  ValueType leftSplitNodeVal;
  bool swapSuccess;

  LOG_INFO("Splitting Node with LPID: %lu, whose parent is %lu", self, parent);

  LOG_INFO("The size of this node (LPID %lu) is %d", self, size_);
  LPage<KeyType, ValueType, KeyComparator> *newLpage =
      new LPage<KeyType, ValueType, KeyComparator>(this->map);

  for (int i = size_ / 2 + 1; i < size_; i++) {
    newLpage->locations_[newPageIndex++] = locations_[i];
  }

  newLpage->size_ = newPageIndex;
  LOG_INFO("The size of the new right split node is %d", newLpage->size_);

  // TODO left_sib is set to self
  newLpage->right_sib_ = right_sib_;

  // Assuming we have ( .. ] ranges
  maxLeftSplitNodeKey = locations_[size_ / 2].first;
  leftSplitNodeVal = locations_[size_ / 2].second;

  maxRightSplitNodeKey = locations_[size_ - 1].first;

  newLpageLPID = this->map->GetMappingTable()->InstallPage(newLpage);

  LOG_INFO("This newly created right split node got LPID: %d", newLpageLPID);

  LPageSplitDelta<KeyType, ValueType, KeyComparator> *splitDelta =
      new LPageSplitDelta<KeyType, ValueType, KeyComparator>(
          this->map, this, maxLeftSplitNodeKey, leftSplitNodeVal, newLpageLPID);

  swapSuccess = this->map->GetMappingTable()->SwapNode(self, this, splitDelta);

  if (swapSuccess == false) {
    LOG_INFO("This SwapNode attempt for split failed");
    delete splitDelta;
    delete newLpage;
    // What should we do on failure? This means that someone else succeeded in
    // doing the
    // atomic half split. Now if try and install our own Insert / Delete /
    // Update delta, it
    // will automatically be created on top of the LPageSplitDelta
    return;
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

  BWTreeNode<KeyType, ValueType, KeyComparator> *parentHardPtr;

  parentHardPtr = this->map->GetMappingTable()->GetNode(parent);
  IPageUpdateDelta<KeyType, ValueType, KeyComparator> *parentUpdateDelta =
      new IPageUpdateDelta<KeyType, ValueType, KeyComparator>(
          this->map, parentHardPtr, maxLeftSplitNodeKey, maxRightSplitNodeKey,
          left_page_lpid, newLpageLPID, false);

  LOG_INFO("Now Doing a SwapNode to install this IPageUpdateDelta");

  // This SwapNode has to be successful. No one else can do this for now.
  // TODO if we allow any node to complete a partial SMO, then this will
  // change
  swapSuccess = this->map->GetMappingTable()->SwapNode(parent, parentHardPtr,
                                                       parentUpdateDelta);

  assert(swapSuccess == true);

  LOG_INFO("Split finished");
}

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
