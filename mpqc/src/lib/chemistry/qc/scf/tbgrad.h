
#ifndef _chemistry_qc_scf_tbgrad_h
#define _chemistry_qc_scf_tbgrad_h

#ifdef __GNUC__
#pragma interface
#endif

template<class T>
class TBGrad {
  protected:
    T& contribution;

  public:
    TBGrad(T&t) : contribution(t) {}
    virtual ~TBGrad() {}

    virtual void build_tbgrad(const RefSCVector&, double, double accuracy) =0;

    inline void set_scale(double& coulombscale, double& exchangescale,
                          int i, int j, int k, int l) const
    {
      double scale = 1.0;

      if ((i!=k)||(j!=l))
        scale *= 2.0;

      if (i!=j)
        scale *= 2.0;

      coulombscale = 0.5*scale;
      exchangescale = -0.25*scale;

      if (k!=l)
        coulombscale *= 2.0;

      if ((k!=l)&&(i==j))
        exchangescale *= 2.0;
    }
};

#endif
