#ifndef TILE_CONCEPT_H
#define TILE_CONCEPT_H

#include <include/tiledarray.h>
#include <memory>

class TileConcept {
  public:
    virtual ~TileConcept() = default;
    virtual TileConcept clone() const = 0;
};

template <typename T>
class TileModel : public TileConcept {
  public:
    TileModel(T t) : tile_type_(std::move(t)) {}
    TileModel(const TileModel &t) = default;
    TileModel(TileModel &&t) = default;
    TileModel &operator=(TileModel &&t) = default;
    TileModel &operator=(TileModel t) {
        tile_ = std::move(t.tile_);
        return *this;
    }

    virtual TileConcept clone() const override final {
        return new TileModel(*this);
    }

  private:
    T tile_;
};

class TileAble {
  public:
    TileAble(TiledArray::Range range) : range_(range), tile_ptr_() {}

    template <typename T>
    TileAble(TiledArray::Range range, T t)
        : range_(range_),
          tile_ptr(std::make_shared<TileModel<T>>(std::move(t))) {
    }

    TileAble(const TileAble &t)
        : range_(range), tile_ptr(std::move(t.tile_ptr->clone())) {}
    TileAble(TileAble &&t) = default;
    TileAble &operator=(TileAble &&t) = default;
    TileAble &operator=(TileAble t) {
        *this = std::move(t);
        return *this;
    }

    inline TileAble clone() const {
        return TileAble(range_, std::move(tile_ptr->clone()));
    }

  private:
    TiledArray::Range range_;
    std::shared_ptr<TileConcept> tile_ptr;
};

#endif // TILE_CONCEPT_H
