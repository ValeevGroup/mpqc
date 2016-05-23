#include <memory>
#include <array>
#include <vector>

class Atom {
 private:
  int isotope_;
  std::array<double, 3> r_;
  int Z_;

 public:
  Atom(int iso, int charge, double x, double y, double z)
      : isotope_(iso), r_{{x, y, z}}, Z_(charge) {}

  int isotope() const { return isotope_; }
  int Z() const { return Z_; }
  std::array<double, 3> const &r() const { return r_; }
};

class Molecule {
 private:
  std::vector<Atom> atoms_;

 public:
  Molecule(std::vector<Atom> atoms) : atoms_(atoms) {}
  std::vector<Atom> const &atoms() const { return atoms_; }
};

class Shell {
 private:
  std::array<double, 3> center_;

 public:
  Shell(double x, double y, double z) : center_{{x, y, z}} {}
  std::array<double, 3> const &center() const { return center_; }
};

class Basis {
 private:
  std::vector<Shell> shells_;

 public:
  Basis(Molecule const &mol) {
    for (auto atom : mol.atoms()) {
      auto const &r = atom.r();
      shells_.push_back(Shell(r[0], r[1], r[2]));
    }
  }

  std::vector<Shell> const &shells() const { return shells_; }
};

class Wavefunction {
 public:
  virtual Molecule const &molecule() const = 0;
  virtual Basis const &basis() const = 0;

  virtual void update(Molecule const &) = 0;
};

class MinimalWFN : public Wavefunction {
 private:
  std::shared_ptr<Molecule> mol_;
  std::shared_ptr<Basis> basis_;

 public:
  MinimalWFN(Molecule const &mol)
      : mol_(std::make_shared<Molecule>(mol)),
        basis_(std::make_shared<Basis>(*mol_)) {}

  Molecule const &molecule() const override { return *mol_; }
  Basis const &basis() const override { return *basis_; }
  void update(Molecule const &mol) override {
    *mol_ = mol;
    *basis_ = Basis(*mol_);
  }
};

class Integrals {
 public:
  Integrals(Basis const &) {}
};

class OneBodyWFN : public Wavefunction {
 private:
  std::shared_ptr<Wavefunction> ref_;
  std::shared_ptr<Integrals> ints_;

 public:
  Molecule const &molecule() const override { return ref_->molecule(); }
  Basis const &basis() const override { return ref_->basis(); }

  void update(Molecule const &mol) override {
    ref_->update(mol);
    *ints_ = Integrals(basis());
  }

  OneBodyWFN(std::shared_ptr<Wavefunction> ref)
      : ref_(std::move(ref)),
        ints_(std::make_shared<Integrals>(ref_->basis())) {}
};

int main() {
  std::vector<Atom> atoms = {
      Atom(1, 1, 0.0, 0.0, 0.0), Atom(1, 1, 1.0, 0.0, 0.0),
      Atom(1, 1, 2.0, 0.0, 0.0),
  };

  auto min_wfn = std::make_shared<MinimalWFN>(MinimalWFN(atoms));
  auto ob_wfn = std::make_shared<OneBodyWFN>(OneBodyWFN(std::move(min_wfn)));

  atoms[1] = Atom(1, 1, 2.0, 1.0, 4.0);

  ob_wfn->update(Molecule(atoms));

  auto &mol = ob_wfn->molecule();
  auto &basie = ob_wfn->basis();

  return 0;
};
