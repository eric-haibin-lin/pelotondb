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


namespace peloton {
namespace index {


typedef uint32_t LPID;

template <typename KeyType, typename ValueType, class KeyComparator>
class BWTree {
private:
  LPID root;
  std::map<uint32_t, BWTree*> mapping_table_;
  bool InsertEntry(const storage::Tuple *key,
      const ItemPointer location);

  bool DeleteEntry(const storage::Tuple *key, const ItemPointer location);
  std::vector<ItemPointer> Scan(const std::vector<Value> &values,
                                const std::vector<oid_t> &key_column_ids,
                                const std::vector<ExpressionType> &expr_types,
                                const ScanDirectionType& scan_direction);

  std::vector<ItemPointer> ScanAllKeys();

  std::vector<ItemPointer> ScanKey(const storage::Tuple *key);
};


// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
template <typename KeyType, typename ValueType, class KeyComparator>
class BWTreeNode {

protected:
  std::map<uint32_t, BWTree*> *mapping_table_;

public:

  BWTreeNode(std::map<uint32_t, BWTree*> *mapping_table)
    : mapping_table_ (mapping_table){};
  virtual bool InsertEntry(const storage::Tuple *key,
      const ItemPointer location)=0;

  virtual bool DeleteEntry(const storage::Tuple *key, const ItemPointer location);
  std::vector<ItemPointer> Scan(const std::vector<Value> &values,
                                const std::vector<oid_t> &key_column_ids,
                                const std::vector<ExpressionType> &expr_types,
                                const ScanDirectionType& scan_direction)=0;

  virtual std::vector<ItemPointer> ScanAllKeys()=0;

  virtual std::vector<ItemPointer> ScanKey(const storage::Tuple *key)=0;

  virtual ~BWTreeNode()=0;

};




}  // End index namespace
}  // End peloton namespace
