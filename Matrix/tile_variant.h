#ifndef TCC_MATRIX_TILE_VARIANT_H
#define TCC_MATRIX_TILE_VARIANT_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include <cstdint>

template <typename T>
class TileVariant {
  public:
    enum TileType : std::uint8_t { LowRank = 0, FullRank = 1 };
    using scaler_type = T;

    TileVariant() : tag_(FullRank), ftile_() {}
    TileVariant(TileVariant const &t) : tag_(t.tag_) { copyTile(t); }
    TileVariant(TileVariant &&t) noexcept : tag_(t.tag_) {
        moveTile(std::move(t));
    }
    ~TileVariant() { destroy_tile(); }

    explicit TileVariant(LowRankTile<T> const &l) : tag_(LowRank), lrtile_(l) {}
    explicit TileVariant(LowRankTile<T> &&l) noexcept : tag_(LowRank),
                                                        lrtile_(std::move(l)) {}

    explicit TileVariant(FullRankTile<T> const &f)
        : tag_(FullRank), ftile_(f) {}
    explicit TileVariant(FullRankTile<T> &&f) noexcept : tag_(FullRank),
                                                         ftile_(std::move(f)) {}

    TileVariant &operator=(TileVariant const &t) {
        if (tag_ == t.tag()) {
            if (tag_ == LowRank) {
                lrtile_ = t.lrtile_;
            } else {
                ftile_ = t.ftile_;
            }
        } else {
            destroy_tile();
            tag_ = t.tag_;
            copyTile(t);
        }
        return *this;
    }

    TileVariant &operator=(TileVariant &&t) noexcept {
        if (tag_ == t.tag()) {
            if (tag_ == LowRank) {
                lrtile_ = std::move(t.lrtile_);
            } else {
                ftile_ = std::move(t.ftile_);
            }
        } else {
            destroy_tile();
            tag_ = t.tag_;
            moveTile(std::move(t));
        }
        return *this;
    }


    LowRankTile<T> const &lrtile() const {
        assert(tag_ == LowRank);
        return lrtile_;
    }

    FullRankTile<T> const &ftile() const {
        assert(tag_ == FullRank);
        return ftile_;
    }

    template <typename Func>
    TileVariant &apply_ternary_mutation(const TileVariant &left,
                                        const TileVariant &right, Func op) {

        switch ((tag() << 2) | (left.tag() << 1) | right.tag()) {
        case LowLowLow:
            *this = op(std::move(lrtile_), left.lrtile(), right.lrtile());
            return *this;
        case LowLowFull:
            *this = op(std::move(lrtile_), left.lrtile(), right.ftile());
            return *this;
        case LowFullLow:
            *this = op(std::move(lrtile_), left.ftile(), right.lrtile());
            return *this;
        case LowFullFull:
            *this = op(std::move(lrtile_), left.ftile(), right.ftile());
            return *this;
        case FullLowLow:
            *this = op(std::move(ftile_), left.lrtile(), right.lrtile());
            return *this;
        case FullLowFull:
            *this = op(std::move(ftile_), left.lrtile(), right.ftile());
            return *this;
        case FullFullLow:
            *this = op(std::move(ftile_), left.ftile(), right.lrtile());
            return *this;
        default: // Full Full Full
            *this = op(std::move(ftile_), left.ftile(), right.ftile());
            return *this;
        }
    }

    template <typename Func>
    TileVariant &apply_binary_mutation(TileVariant const &right, Func op) {
        switch ((tag() << 1) | right.tag()) {
        case LowLow:
            *this = op(std::move(lrtile_), right.lrtile());
            return *this;
        case LowFull:
            *this = op(std::move(lrtile_), right.ftile());
            return *this;
        case FullLow:
            *this = op(std::move(ftile_), right.lrtile());
            return *this;
        default: // Full Full
            *this = op(std::move(ftile_), right.ftile());
            return *this;
        }
    }

    template <typename Func>
    TileVariant &apply_unary_mutation(Func op) {
        if (tag() == LowRank) {
            *this = op(std::move(lrtile_));
        } else {
            *this = op(std::move(ftile_));
        }
        return *this;
    }

    template <typename Func>
    auto apply_binary_op(const TileVariant &right,
                         Func op) const -> decltype(op(lrtile(),
                                                       lrtile())) const {
        switch ((tag() << 1) | right.tag()) {
        case LowLow:
            return op(lrtile(), right.lrtile());
        case LowFull: // Low Full
            return op(lrtile(), right.ftile());
        case FullLow:
            return op(ftile(), right.lrtile());
        default:
            return op(ftile(), right.ftile());
        }
    }

    template <typename Func>
    auto apply_unary_op(Func op) const -> decltype(op(lrtile())) const {
        if (tag() == LowRank) {
            return op(lrtile());
        } else {
            return op(ftile());
        }
    }

    TileType tag() const { return tag_; }

    /*
     * Some Predefined operations implemented using the interface the op
     * interface
     */
    unsigned long rank() const { return apply_unary_op(low_rank_functor); }

    unsigned long full_rank() const {
        return apply_unary_op(full_rank_functor);
    }

    typename FullRankTile<T>::template Matrix<T> matrix() const {
        return apply_unary_op(matrix_functor);
    }


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
    /*
     * Functions that use the unary op interface.
     */
    struct {
        unsigned long operator()(FullRankTile<T> const &t) {
            return std::min(t.Rows(), t.Cols());
        }

        unsigned long operator()(LowRankTile<T> const &t) {
            return std::min(t.Rows(), t.Cols());
        }
    } full_rank_functor;

    struct {
        unsigned long operator()(FullRankTile<T> const &t) { return t.rank(); }

        unsigned long operator()(LowRankTile<T> const &t) { return t.rank(); }
    } low_rank_functor;

    struct {
        typename FullRankTile<T>::template Matrix<T>
        operator()(FullRankTile<T> const &t) {
            return t.matrix();
        }

        // Using FullRankTile because tile.matrix() is always a full matrix.
        typename FullRankTile<T>::template Matrix<T>
        operator()(LowRankTile<T> const &t) {
            return t.matrix();
        }
    } matrix_functor;


    /*
     * Functions to help with copying and moving tiles
     */
    void copyTile(TileVariant const &t) {
        if (t.tag_ == LowRank) {
            new (&lrtile_) LowRankTile<T>{t.lrtile_};
        } else {
            new (&ftile_) FullRankTile<T>{t.ftile_};
        }
    }

    void moveTile(TileVariant &&t) noexcept {
        if (t.tag_ == LowRank) {
            new (&lrtile_) LowRankTile<T>{std::move(t.lrtile_)};
        } else {
            new (&ftile_) FullRankTile<T>{std::move(t.ftile_)};
        }
    }

    // Calls destructors of union members
    void destroy_tile() noexcept {
        if (tag_ == LowRank) {
            lrtile_.~LowRankTile<T>();
        } else {
            ftile_.~FullRankTile<T>();
        }
    }

    enum VariantSwitchId : std::uint8_t {
        LowLow = 0,
        LowFull = 1,
        FullLow = 2,
        LowLowLow = 0,
        LowLowFull = 1,
        LowFullLow = 2,
        LowFullFull = 3,
        FullLowLow = 4,
        FullLowFull = 5,
        FullFullLow = 6
    };

    // Check binary switch
    static_assert((LowRank << 1 | LowRank) == LowLow,
                  "Low Low switch is incorrect");
    static_assert((LowRank << 1 | FullRank) == LowFull,
                  "Low Full switch is incorrect");
    static_assert((FullRank << 1 | LowRank) == FullLow,
                  "Full Low switch is incorrect");

    // Check ternary switch
    static_assert((LowRank << 2 | LowRank << 1 | LowRank) == LowLowLow,
                  "Low Low Low switch is incorrect");
    static_assert((LowRank << 2 | LowRank << 1 | FullRank) == LowLowFull,
                  "Low Low Full switch is incorrect");
    static_assert((LowRank << 2 | FullRank << 1 | LowRank) == LowFullLow,
                  "Low Full Low switch is incorrect");
    static_assert((LowRank << 2 | FullRank << 1 | FullRank) == LowFullFull,
                  "Low Full Full switch is incorrect");
    static_assert((FullRank << 2 | LowRank << 1 | LowRank) == FullLowLow,
                  "Full Low Low switch is incorrect");
    static_assert((FullRank << 2 | LowRank << 1 | FullRank) == FullLowFull,
                  "Full Low Full switch is incorrect");
    static_assert((FullRank << 2 | FullRank << 1 | LowRank) == FullFullLow,
                  "Full Full Low switch is incorrect");
};

#endif // TTC_MATRIX_TILE_VARIANT_H
