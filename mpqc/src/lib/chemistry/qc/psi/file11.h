#ifndef _chemistry_qc_psi_file11_h
#define _chemistry_qc_psi_file11_h

class FILE11 {
  private:
    char label_[132];
    char theory_[132];
    char dertype_[132];
    int nat;
    double energy_;
    double *coord_[3];
    int *charges_;
    double *grad_[3];
    int ngrad_;
  public:
    FILE11(int=0);
    ~FILE11();
    int read(int=0);
    void print();
    const char *label(){return label_;}
    const char *theory(){return theory_;}
    const char *dertype(){return dertype_;}
    int natoms(){return nat;}
    int ngrad(){return ngrad_;}
    double energy(){return energy_;}
    double coordinate(int i, int j){return coord_[i][j];}
    double charge(int i){return charges_[i];}
    double gradient(int i, int j){return grad_[i][j];}
    };


#endif
