
#ifndef _chemistry_qc_scf_gbuild_h
#define _chemistry_qc_scf_gbuild_h

#ifdef __GNUC__
#pragma interface
#endif

template<class T>
class GBuild {
  protected:
    T& contribution;

  public:
    GBuild(T&t) : contribution(t) {}
    virtual ~GBuild() {}

    virtual void build_gmat(double accuracy) =0;
};

#endif
