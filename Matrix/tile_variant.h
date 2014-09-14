#ifndef TCC_MATRIX_TILEVARIANT_H
#define TCC_MATRIX_TILEVARIANT_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"

template <typename T>
class TileVariant {
  public:
    TileVariant() : tag_(FullRank), ftile_() {}
    ~TileVariant() {
        switch (tag_) {
        case TileVariant::LowRank:
            lrtile_.~LowRankTile<T>();
            break;
        case TileVariant::FullRank:
            ftile_.~FullRankTile<T>();
            break;
        }
    }

    TileVariant(TileVariant const &t) : tag_(t.tag_) { copyTileVariant(t); }
    // TODO figure out what to do about move constructor.

    void copyTileVariant(TileVariant const &t) {
        switch (t.tag_) {
        case TileVariant::LowRank:
            new (&lrtile_) LowRankTile<T>(t.lrtile_);
            break;
        case TileVariant::FullRank:
            new (&ftile_) FullRankTile<T>(t.ftile_);
            break;
        }
    }

  private:
    enum {
        LowRank,
        FullRank
    } tag_;

    union {
        LowRankTile<T> lrtile_;
        FullRankTile<T> ftile_;
    };
};

#endif // TTC_MATRIX_TILEVARIANT_H
