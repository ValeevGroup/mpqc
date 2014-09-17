#ifndef TCC_MATRIX_TILEVARIANT_H
#define TCC_MATRIX_TILEVARIANT_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include "../include/tcc_all_same.h"
#include <cstdint>

template <typename T>
class TileVariant {
  public:
    enum TileType : std::int8_t {
        LowRank = 0,
        FullRank = 1
    };

    TileVariant() : tag_(FullRank), ftile_() {}
    ~TileVariant() {
        switch (tag_) {
        case LowRank:
            lrtile_.~LowRankTile<T>();
            break;
        case FullRank:
            ftile_.~FullRankTile<T>();
            break;
        }
    }


    TileVariant(TileVariant const &t) : tag_(t.tag_) { copyTileVariant(t); }
    void copyTileVariant(TileVariant const &t) {
        switch (t.tag_) {
        case LowRank:
            new (&lrtile_) LowRankTile<T>(t.lrtile_);
            break;
        case FullRank:
            new (&ftile_) FullRankTile<T>(t.ftile_);
            break;
        }
    }

    // TODO figure out what to do about move constructor.

    TileVariant(const LowRankTile<T> &l) : tag_(LowRank), lrtile_(l) {}
    TileVariant(LowRankTile<T> &&l) : tag_(LowRank), lrtile_(std::move(l)) {}

    TileVariant(const FullRankTile<T> &f) : tag_(FullRank), ftile_(f) {}
    TileVariant(FullRankTile<T> &&f) : tag_(FullRank), ftile_(std::move(f)) {}

    LowRankTile<T> const &lrtile() const {
        assert(tag_ == LowRank);
        return lrtile_;
    }

    FullRankTile<T> const &ftile() const {
        assert(tag_ == FullRank);
        return ftile_;
    }

    template <typename Func>
    void apply_binary_op_to(const TileVariant &left, const TileVariant &right,
                            Func op) {
        switch ((tag() << 2) | (left.tag() << 1) | right.tag()) {
        case 0: // Low Low Low
            op(lrtile_, left.lrtile(), right.lrtile());
            break;
        case 1: // Low Low Full
            op(lrtile_, left.lrtile(), right.ftile());
            break;
        case 2: // Low Full Low
            op(lrtile_, left.ftile(), right.lrtile());
            break;
        case 3: // Full Low Low
            op(ftile_, left.lrtile(), right.lrtile());
            break;
        case 4: // Low Full Full
            op(lrtile_, left.ftile(), right.ftile());
            break;
        case 5: // Full Low Full
            op(ftile_, left.lrtile(), right.ftile());
            break;
        case 6: // Full Full Low
            op(ftile_, left.ftile(), right.lrtile());
            break;
        case 7: // Full Full Full
            op(ftile_, left.ftile(), right.ftile());
            break;
        }
    }

    template <typename Func>
    TileVariant apply_binary_op(const TileVariant &right, Func op) const {
        switch ((tag() << 1) | right.tag()) {
        case 0: // Low Low
            return TileVariant{op(lrtile(), right.lrtile())};
        case 1: // Low Full
            return TileVariant{op(lrtile(), right.ftile())};
        case 2: // Full Low
            return TileVariant{op(ftile(), right.lrtile())};
        case 3: // Full Full
            return TileVariant{op(ftile(), right.ftile())};
        }
    }

    template <typename Func>
    auto apply_binary_transform_op(const TileVariant &right, Func op)
        -> decltype(op(lrtile(), lrtile())) const {

        static_assert(
            tcc::all_same
            <decltype(op(lrtile(), lrtile())), decltype(op(lrtile(), ftile())),
             decltype(op(ftile(), lrtile())),
             decltype(op(ftile(), ftile()))>::value,
            "All return types of functor for binary transform op must have the "
            "same type.");

        switch ((tag() << 1) | right.tag()) {
        case 0: // Low Low
            return op(lrtile(), right.lrtile());
        case 1: // Low Full
            return op(lrtile(), right.ftile());
        case 2: // Full Low
            return op(ftile(), right.lrtile());
        case 3: // Full Full
            return op(ftile(), right.ftile());
        }
    }

    template <typename Func>
    void apply_unary_op_to(Func op) {
        switch (tag()) {
        case LowRank:
            return TileVariant{op(lrtile_)};
        case FullRank:
            return TileVariant{op(ftile_)};
        }
    }

    template <typename Func>
    TileVariant apply_unary_op(Func op) const {
        switch (tag()) {
        case LowRank:
            return TileVariant{op(lrtile())};
        case FullRank:
            return TileVariant{op(ftile())};
        }
    }

    // This function is for operations which don't modify the tile, but need
    // to read the tile's data, an example would be 2norm.
    template <typename Func>
    auto apply_unary_transform_op(Func op) -> decltype(op(lrtile())) const {
        static_assert(tcc::all_same
                      <decltype(op(lrtile())), decltype(op(ftile()))>::value,
                      "Unary Transform op must return the same type for every "
                      "tile type.");
        switch (tag()) {
        case LowRank:
            return op(lrtile());
        case FullRank:
            return op(ftile());
        }
    }


    TileType tag() const { return tag_; }

  private:
    TileType tag_;

    union {
        LowRankTile<T> lrtile_;
        FullRankTile<T> ftile_;
    };
};

#endif // TTC_MATRIX_TILEVARIANT_H
