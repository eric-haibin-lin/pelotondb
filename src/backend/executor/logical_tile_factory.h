//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// logical_tile_factory.h
//
// Identification: src/backend/executor/logical_tile_factory.h
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <vector>
#include <memory>

#include "backend/common/types.h"

namespace peloton {

namespace storage {
class Tile;
class TileGroup;
class AbstractTable;
}

//===--------------------------------------------------------------------===//
// Logical Tile Factory
//===--------------------------------------------------------------------===//

namespace executor {
class LogicalTile;

class LogicalTileFactory {
 public:
  static LogicalTile *GetTile();

  static LogicalTile *WrapTiles(
      const std::vector<std::shared_ptr<storage::Tile>> &base_tile_refs);

  static LogicalTile *WrapTileGroup(
      const std::shared_ptr<storage::TileGroup> &tile_group,
      txn_id_t txn_id);

  static std::vector<LogicalTile *> WrapTileGroups(
      const std::vector<ItemPointer> tuple_locations,
      const std::vector<oid_t> column_ids, txn_id_t txn_id, cid_t commit_id);
};

}  // namespace executor
}  // namespace peloton
