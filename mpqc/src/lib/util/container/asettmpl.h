
template <class Type>
class Arrayset : public Set<Type>
{
 public:
    Arrayset() {}
    Arrayset(const Arrayset<Type>&a): Set<Type>(a) {}
    ~Arrayset() {}
    Type& operator[](int i)
    {
      if (i<0 || i>=nelement) {
          fprintf(stderr,"Arrayset::operator[] out of range: "
                  "%d (nelement = %d)\n", i, nelement);
          abort();
        };
      return element[i];
    }
    const Type& operator[](int i) const
    {
      if (i<0 || i>=nelement) {
          fprintf(stderr,"Arrayset::operator[] out of range: "
                  "%d (nelement = %d)\n", i, nelement);
          abort();
        };
      return element[i];
    }
    int iseek(const Type &item)
    {
      for (int i=0; i<nelement; i++) {
          if (item == element[i]) return i;
        }
      return -1;
    }
};
