#ifndef TCC_MATRIX_TILE_VARIANT_H
#define TCC_MATRIX_TILE_VARIANT_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include "../include/tcc_all_same.h"
#include <cstdint>

template <typename T>
class TileVariant {
  public:
    enum TileType : std::uint8_t { LowRank = 0, FullRank = 1 };

    using scaler_type = T;

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

    TileVariant &operator=(TileVariant const &t) {
        if (tag_ == t.tag()) {
            if (tag_ == LowRank) {
                lrtile_ = t.lrtile_;
            } else {
                ftile_ = t.ftile_;
            }
        } else {
            if (tag_ == LowRank) {
                lrtile_.~LowRankTile<T>();
            } else {
                ftile_.~FullRankTile<T>();
            }
            tag_ = t.tag_;
            copyTileVariant(t);
        }
        return *this;
    }


    TileVariant(TileVariant &&t) : tag_(std::move(t.tag_)) {
        moveTileVariant(std::move(t));
    }

    TileVariant &operator=(TileVariant &&t) {
        if (tag_ == t.tag()) {
            if (tag_ == LowRank) {
                lrtile_ = std::move(t.lrtile_);
            } else {
                ftile_ = std::move(t.ftile_);
            }
        } else {
            if (tag_ == LowRank) {
                lrtile_.~LowRankTile<T>();
            } else {
                ftile_.~FullRankTile<T>();
            }
            tag_ = std::move(t.tag_);
            moveTileVariant(std::move(t));
        }
        return *this;
    }


    explicit TileVariant(const LowRankTile<T> &l) : tag_(LowRank), lrtile_(l) {}
    explicit TileVariant(LowRankTile<T> &&l)
        : tag_(LowRank), lrtile_(std::move(l)) {}

    explicit TileVariant(const FullRankTile<T> &f)
        : tag_(FullRank), ftile_(f) {}
    explicit TileVariant(FullRankTile<T> &&f)
        : tag_(FullRank), ftile_(std::move(f)) {}

    LowRankTile<T> const &lrtile() const {
        assert(tag_ == LowRank);
        return lrtile_;
    }

    FullRankTile<T> const &ftile() const {
        assert(tag_ == FullRank);
        return ftile_;
    }

    template <typename Func>
    TileVariant &apply_binary_mutation(const TileVariant &left,
                                       const TileVariant &right, Func op) {

        switch ((tag() << 2) | (left.tag() << 1) | right.tag()) {
        case 0: // Low Low Low
            *this = op(std::move(lrtile_), left.lrtile(), right.lrtile());
            return *this;
        case 1: // Low Low Full
            *this = op(std::move(lrtile_), left.lrtile(), right.ftile());
            return *this;
        case 2: // Low Full Low
            *this = op(std::move(lrtile_), left.ftile(), right.lrtile());
            return *this;
        case 3: // Full Low Low
            *this = op(std::move(ftile_), left.lrtile(), right.lrtile());
            return *this;
        case 4: // Low Full Full
            *this = op(std::move(lrtile_), left.ftile(), right.ftile());
            return *this;
        case 5: // Full Low Full
            *this = op(std::move(ftile_), left.lrtile(), right.ftile());
            return *this;
        case 6: // Full Full Low
            *this = op(std::move(ftile_), left.ftile(), right.lrtile());
            return *this;
        case 7: // Full Full Full
            *this = op(std::move(ftile_), left.ftile(), right.ftile());
            return *this;
        default: // Should never be reached
            assert(false);
            return *this;
        }
    }

    template <typename Func>
    auto apply_binary_op(const TileVariant &right,
                         Func op) const -> decltype(op(lrtile(), lrtile())) {

        static_assert(
            tcc::all_same<decltype(op(lrtile(), lrtile())),
                          decltype(op(lrtile(), ftile())),
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
        default: // Should never be reached
            assert(false);
            return decltype(op(ftile(), ftile())){};
        }
    }

    template <typename Func>
    TileVariant &apply_unary_mutation(Func op) {
        static_assert(
            tcc::all_same<decltype(op(lrtile())), decltype(op(ftile())),
                          decltype(*this)>::value,
            "Return types for a mutation must be TileVariants.");


        switch (tag()) {
        case LowRank:
            *this = op(std::move(lrtile_));
            return *this;
        case FullRank:
            *this = op(std::move(ftile_));
            return *this;
        }
    }

    template <typename Func>
    auto apply_unary_op(Func op) const -> decltype(op(lrtile())) {
        static_assert(
            tcc::all_same<decltype(op(lrtile())), decltype(op(ftile()))>::value,
            "Unary Transform op must return the same type for every "
            "tile type.");

        switch (tag()) {
        case LowRank:
            return op(lrtile());
        case FullRank:
            return op(ftile());
        }
    }

    unsigned long rank() const {
        return apply_unary_op([](auto const &t) { return t.rank(); });
    }

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix() const {
        return apply_unary_op([](auto const &t) { return t.matrix(); });
    }


    TileType tag() const { return tag_; }

  private:
    TileType tag_;

    union {
        LowRankTile<T> lrtile_;
        FullRankTile<T> ftile_;
    };

    /*
     * Utililty Functions
     */
  private:
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

    void moveTileVariant(TileVariant &&t) {
        switch (t.tag_) {
        case LowRank:
            new (&lrtile_) LowRankTile<T>(std::move(t.lrtile_));
            break;
        case FullRank:
            new (&ftile_) FullRankTile<T>(std::move(t.ftile_));
            break;
        }
    }
};

#endif // TTC_MATRIX_TILE_VARIANT_H
