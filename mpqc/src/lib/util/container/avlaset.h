
#include <util/container/pixintRAVLMap.h>
#include <util/container/intpixRAVLMap.h>

template <class Type>
class Arrayset : public AVLSet<Type>
{
  private:
    int nelement;
    PixintRAVLMap pixtoint;
    intPixRAVLMap inttopix;
  public:
    Arrayset(): pixtoint(-1), inttopix(0), nelement(0) {}
    Arrayset(const Arrayset<Type>&a):
      AVLSet<Type>(a),pixtoint(-1),inttopix(0),nelement(0) {}
    ~Arrayset() {}
    Type& operator[](int i)
    {
      return operator()(inttopix[i]);
    }
    const Type& operator[](int i) const
    {
      // must cast away the constness
      Arrayset<Type> *ncthis = (Arrayset<Type>*) this;
      return ncthis->operator()(ncthis->inttopix[i]);
    }
    int iseek(const Type &item)
    {
      return pixtoint[seek(item)];
    }
    Pix add(Type& item) {
        Pix r = AVLSet<Type>::add(item);
        pixtoint[r] = nelement;
        inttopix[nelement] = r;
        nelement++;
      }
    void del(Type& item) {
        Pix p = seek(item);
        int i = pixtoint[p];
        AVLSet<Type>::del(item);
        pixtoint.del(p);
        inttopix.del(i);
      }
};
