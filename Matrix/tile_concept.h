#ifndef TILE_CONCEPT_H
#define TILE_CONCEPT_H

#include "../include/tiledarray.h"
#include <memory>

class TileConcept {
  public:
    virtual ~TileConcept() = default;
    virtual TileConcept *clone() const = 0;

    virtual void const *data() const = 0;
    virtual void *data() = 0;
    virtual std::type_info const &type() const noexcept = 0;

};

template <typename T>
class TileModel : public TileConcept {
  public:
    TileModel(T t) : tile_(std::move(t)) {}
    TileModel(const TileModel &t) = default;
    TileModel(TileModel &&t) = default;
    TileModel &operator=(TileModel &&t) = default;
    TileModel &operator=(TileModel t) {
        tile_ = std::move(t.tile_);
        return *this;
    }

    virtual TileConcept *clone() const override final {
        return new TileModel(*this);
    }

    virtual void const *data() const override final { return &tile_; }
    virtual void *data() override final { return &tile_; }

    virtual std::type_info const &type() const noexcept override final {
        return typeid(T);
    }

  private:
    T tile_;
};

class TileAble {
  public:
    TileAble(TiledArray::Range range) : range_(range), tile_ptr_() {}

    template <typename T>
    TileAble(TiledArray::Range range, T t)
        : range_(std::move(range)),
          tile_ptr_(std::make_shared<TileModel<T>>(std::move(t))) {}

    template <typename T>
    TileAble(T t)
        : range_(TiledArray::Range::Range()),
          tile_ptr_(std::make_shared<TileModel<T>>(std::move(t))) {}

    TileAble(const TileAble &t)
        : range_(t.range_), tile_ptr_(std::move(t.tile_ptr_->clone())) {}
    TileAble(TileAble &&t) = default;
    TileAble &operator=(TileAble &&t) = default;
    TileAble &operator=(TileAble t) {
        range_ = std::move(t.range_);
        tile_ptr_ = std::move(t.tile_ptr_);
        return *this;
    }

    inline TileAble clone() const {
        return TileAble(range_, std::move(tile_ptr_->clone()));
    }

    inline std::type_info const &type() const { return tile_ptr_->type(); }

    template <typename T>
    T *cast_data() {
        assert(tile_ptr_->type() == typeid(T));
        return static_cast<T *>(tile_ptr_->data());
    }

    template <typename T>
    T const *cast_data() const {
        assert(tile_ptr_->type() == typeid(T));
        return static_cast<T const *>(tile_ptr_->data());
    }

  private:
    TiledArray::Range range_;
    std::shared_ptr<TileConcept> tile_ptr_;
};

#endif // TILE_CONCEPT_H
