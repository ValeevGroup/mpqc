//
// findisp.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _chemistry_molecule_findisp_h
#define _chemistry_molecule_findisp_h

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include <util/misc/scexception.h>
#include <chemistry/molecule/deriv.h>
#include <chemistry/molecule/energy.h>

namespace sc {

  /// @addtogroup ChemistryMolecule
  /// @{

  /// Maps displacements in terms of symmetrized coordinates to property values
  /// @tparam Value the property
  template <typename Value>
    class Displacements {
      public:
      Displacements() {}
      ~Displacements() {}

      void save(StateOut& s) const {
        s.put(disps_);
      }
      void restore(StateIn& s) {
        s.get(disps_);
      }

      typedef int Coord;
      typedef double DisplacementSize;
      typedef std::pair<Coord,DisplacementSize> KeyElem;
      typedef std::vector<KeyElem> Key;
      typedef typename std::map< Key, Value>::iterator iter;
      typedef typename std::map< Key, Value>::const_iterator citer;

      citer begin() const { return disps_.begin(); }
      citer end() const { return disps_.end(); }
      std::size_t size() const { return disps_.size(); }

      std::pair<Key,Value> find(Coord c1, DisplacementSize d1) const {
        Key key;
        key.push_back(std::make_pair(c1,d1));
        return find(key);
      }
      std::pair<Key,Value> find(Coord c1, DisplacementSize d1,
                                Coord c2, DisplacementSize d2) const {
        Key key;
        key.push_back(std::make_pair(c1,d1));
        key.push_back(std::make_pair(c2,d2));
        std::stable_sort(key.begin(), key.end());
        return find(key);
      }
      /// assumes that key is sorted
      bool has(const Key& key) const {
        citer v = disps_.find(key);
        return v != disps_.end();
      }
      /// assumes that key is sorted
      std::pair<Key,Value> find(const Key& key) const {
        citer v = disps_.find(key);
        if (v == disps_.end())
          throw sc::ProgrammingError("Displacement::find -- displacement not found",__FILE__,__LINE__);
        return *v;
      }
      /// finds the first instance of value (implemented as dumb O(N) search)
      citer find(const Value& value) const {
        for(citer v = disps_.begin(); v != disps_.end(); ++v) {
          if ((*v).second == value) return v;
        }
        throw sc::ProgrammingError("Displacement::find -- value not found",__FILE__,__LINE__);
      }

      void push(Coord c1, DisplacementSize d1,
                const Value& v) {
        Key key;
        key.push_back(std::make_pair(c1,d1));
        return push(key,v);
      }
      void push(Coord c1, DisplacementSize d1,
                Coord c2, DisplacementSize d2,
                const Value& v) {
        Key key;
        key.push_back(std::make_pair(c1,d1));
        key.push_back(std::make_pair(c2,d2));
        std::stable_sort(key.begin(), key.end());
        return push(key,v);
      }
      /// assumes that key is sorted
      void push(const Key& key, const Value& v) {
        iter i = disps_.find(key);
        if (i != disps_.end())
          throw sc::ProgrammingError("Displacement::push -- displacement already exists",__FILE__,__LINE__);
        disps_[key] = v;
      }

      private:
      /** displacements are represented as a key -> index map
       *  keys are formatted as comma-separated list of integer,double pairs
       *  each pair contains the index of SALC (0-based) + displacement size
       *  SALC indices are nondecreasing
       */
      std::map< Key, Value> disps_;
  };

#if 0
  /** FinDispDerivative computes derivatives of functions using finite-difference formulas
   * @tparam TargetOrder desired order of the derivative
   * @tparam Function class that describes the function (see traits class)
  */
  template <unsigned int TargetOrder, typename Function>
  class FinDispDerivative : virtual public SavableState {
    public:
      /**
       * Constructor
       * @param function The function (since it must have state, it is a functor) whose derivative is to be evaluated. Its parameters must be set to the "reference" values
       * @return none
       */
      FinDispDerivative(const Ref<Function>& function);
      FinDispDerivative(StateIn&);
      ~FinDispDerivative();
      void save_data_state(StateOut&);

      typedef typename Function::Value FunctionValue;
      typedef typename Function::Parameter FunctionParameter;
      static const unsigned int TensorTraits<FunctionParameter>::rank param_rank;
      const ParameterDerivative<FunctionValue,IntegerProduct<TargetOrder,param_rank>::value>::result Result;

      const Result& result();
    private:
      static ClassDesc class_desc_;
      Ref<Function> function_;
      Result result_;

      void compute();
  };
#endif

  /// energy + gradient + hessian
  struct EGH {
    public:
      EGH();
      EGH(double e, const RefSCVector& g, const RefSymmSCMatrix& h);
      ~EGH();

      double energy() const { return energy_; }
      const RefSCVector& gradient() const { return gradient_; }
      const RefSymmSCMatrix& hessian() const { return hessian_; }

    private:
      double energy_;
      RefSCVector gradient_;
      RefSymmSCMatrix hessian_;
  };
  //@{ specify how to serialize/deserialize EGH
  template <> void FromStateIn<EGH>(EGH& v, StateIn& s, int& count);
  template <> void ToStateOut<EGH>(const EGH& v, StateOut& so, int& count);
  //@}


/** Computes the molecular hessian by finite displacements of gradients (or, if not available, energies).
    This will use the minimum number of displacements, each in the
    highest possible point group. */
class FinDispMolecularHessian: public MolecularHessian {

  /// Params encpasulates parameters of the finite-difference molecular hessian evaluator
  class Params : virtual public SavableState {
    public:
      Params();
      Params(const Ref<KeyVal>& kv);
      Params(StateIn&);
      ~Params();
      void save_data_state(StateOut&);

      const Ref<PointGroup>& disp_pg() const { return disp_pg_; }
      double disp_size() const { return disp_; }
      bool only_totally_symmetric() const { return only_totally_symmetric_; }
      bool eliminate_quadratic_terms() const { return eliminate_quadratic_terms_; }
      bool do_null_displacement() const { return do_null_displacement_; }
      int debug() const { return debug_; }
      bool checkpoint() const { return checkpoint_; }
      const std::string& checkpoint_file() const { return checkpoint_file_; }
      bool restart() const { return restart_; }
      const std::string& restart_file() const { return restart_file_; }
      bool use_energies() const { return use_energies_; }
      double gradient_accuracy() const { return gradient_accuracy_; }
      double energy_accuracy() const { return energy_accuracy_; }
      int nirrep() const { return disp_pg_->char_table().nirrep(); }

      void set_eliminate_quadratic_terms(bool e) { eliminate_quadratic_terms_ = e; }
      void set_disp_size(double s);
      void set_disp_pg(const Ref<PointGroup>& pg) { disp_pg_ = pg; }
      void set_restart(bool r = true) { restart_ = r; }
      void set_checkpoint(bool c = true) { checkpoint_ = c; }
      void set_desired_accuracy(double acc);

    private:
      static ClassDesc class_desc_;

      // In case molecule must be given in lower symmetry, its actual
      // symmetry and the symmetry used to compute displacements is this
      Ref<PointGroup> disp_pg_;
      // the cartesian displacement size in bohr
      double disp_;
      // only do the totally symmetric displacements (makes sense if the Hessian is to be used in geometry optimization)
      bool only_totally_symmetric_;
      // eliminate the quadratic terms in the hessian by doing an extra displacement for
      // each of the totally-symmetric coordinates
      bool eliminate_quadratic_terms_;
      // use the gradient at the initial geometry to remove first order terms
      // (important if not at equilibrium geometry)
      bool do_null_displacement_;
      // print flag
      int debug_;
      // whether or not to checkpoint
      bool checkpoint_;
      // the name of the checkpoint file
      std::string checkpoint_file_;
      // whether or not to attempt a restart
      bool restart_;
      // the name of the restart file
      std::string restart_file_;
      // force computation from energies
      bool use_energies_;
      // the accuracy for energy calculations
      double energy_accuracy_;
      // the accuracy for gradient calculations
      double gradient_accuracy_;
  };


    // Implements FinDispMolecularHessian
    class Impl : virtual public SavableState {
      public:
      typedef Displacements<EGH>::Key Displacement;

      Impl(const Ref<MolecularEnergy>& e,
           const Ref<Params>& params);
      Impl(StateIn& s);
      ~Impl();
      void save_data_state(StateOut& s);

      void checkpoint_displacements(StateOut&);
      void restore_displacements(StateIn&);

      const Ref<Params>& params() const { return params_; }
      const Ref<MolecularEnergy>& mole() const { return mole_; }
      void set_mole(const Ref<MolecularEnergy>& mole) { mole_ = mole; }
      virtual int ndisplace() const =0;
      // returns a matrix whose columns are SALC displacements
      RefSCMatrix displacements(int irrep) const;
      void displace(const Displacement& d);
      void original_geometry();

      virtual void init();
      virtual void restart();

      /** This returns the cartesian hessian.  If it has not yet been
          computed, it will be computed by finite displacements. */
      RefSymmSCMatrix cartesian_hessian();

      // transforms hessian in symmetrized coordinates to cartesian coordinates
      void do_hess_for_irrep(int irrep,
                             const RefSymmSCMatrix &dhessian,
                             const RefSymmSCMatrix &xhessian);

      protected:
      Ref<Params> params_;

      Ref<MolecularEnergy> mole_;
      // The molecule's original point group for restoration at the end.
      Ref<PointGroup> original_point_group_;
      // The molecule's original geometry for restoration at the end and
      //computing displacements.
      RefSCVector original_geometry_;
      // a basis for the symmetrized cartesian coordinates (rows)
      RefSCMatrix symbasis_;

      // values, gradients, and hessians at various (displaced) geometries
      Displacements<EGH> values_;

      /** computes a cartesian hessian from energies or gradients at displaced geometries */
      virtual RefSymmSCMatrix compute_hessian() =0;
      /// computes MolecularEnergy object
      virtual void compute_mole(const Displacement& d) =0;

      Ref<SCMatrixKit> matrixkit() const { return mole_->matrixkit(); }
      RefSCDimension d3natom() const { return mole_->moldim(); }
      /// given symmetrized coordinate index return its irrep
      unsigned int coor_to_irrep(unsigned int symm_coord) const;

      static ClassDesc class_desc_;
    };

    // Implementation of finite-difference hessians from gradients
    class GradientsImpl : public Impl {
      public:
        GradientsImpl(const Ref<MolecularEnergy>& e,
                      const Ref<Params>& params);
        GradientsImpl(StateIn& s);
        ~GradientsImpl();
        void save_data_state(StateOut& s);

      private:

        // will throw if mole lacks gradients capability; do nothing if mole is null
        static void validate_mole(const Ref<MolecularEnergy>& e);
        int ndisplace() const;
        void set_gradient(const Displacement& d, double energy, const RefSCVector &grad);
        RefSymmSCMatrix compute_hessian();
        void compute_mole(const Displacement& d);

        static ClassDesc class_desc_;
    };

    // Implementation of finite-difference hessians from energies
    class EnergiesImpl : public Impl {
      public:
        // will throw if mole lacks energies capability; do nothing if mole is null
        EnergiesImpl(const Ref<MolecularEnergy>& e,
                     const Ref<Params>& params);
        EnergiesImpl(StateIn& s);
        ~EnergiesImpl();
        void save_data_state(StateOut& s);

      private:

        static void validate_mole(const Ref<MolecularEnergy>& e);
        int ndisplace() const;
        RefSymmSCMatrix compute_hessian();
        void compute_mole(const Displacement& d);
        const Displacements<EGH>& values() const { return values_; }

        struct Eij {
            Eij(int i, int j, EnergiesImpl& eimpl) :
              i_(i), j_(j), eimpl_(eimpl)
            {
            }

            // return energy computed at geometry displaced by di (in units of disp_)
            // along i and dj along j
            double operator()(int di, int dj) {
              Displacement disp;
              if (i_ != j_) {
                if (di != 0) disp.push_back(std::make_pair(i_,double(di)));
                if (dj != 0) disp.push_back(std::make_pair(j_,double(dj)));
              }
              else { // i_ == j_
                if (di+dj != 0) disp.push_back(std::make_pair(i_,double(di+dj)));
              }
              eimpl_.compute_mole(disp);
              const double energy = eimpl_.values().find(disp).second.energy();
              return energy;
            }

            int i_, j_;
            EnergiesImpl& eimpl_;
        };

        static ClassDesc class_desc_;
    };

  private:
    Ref<Impl> pimpl_; //<< initliazed lazily
    Ref<MolecularEnergy> mole_init_;   //< pimpl_ is initalized lazily, so this is used to hold the MolecularEnergy object used to initialize the pimpl_
    Ref<Params> params_;
    bool user_provided_eliminate_quadratic_terms_; // default for eliminate_quadratic_terms depends on type of Impl ...
                                               // must know whether the default needs to vary

    void override_default_params(); // override some defaults based on the properties of MolecularEnergy

  protected:

    /// initializes pimpl_, it should not be called until e is fully initalized, hence use this lazily
    void init_pimpl(const Ref<MolecularEnergy>& e);
    void restart();

  public:
    FinDispMolecularHessian(const Ref<MolecularEnergy>&);
    /** The FinDispMolecularHessian KeyVal constructor is used to generate a
        FinDispMolecularHessian object from the input.  It reads the keywords
        below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>energy</tt><td>MolecularEnergy<td>none<td>This gives an
        object which will be used to compute the gradients (or energies) needed to form
        the hessian.  If this is not specified, the object using
        FinDispMolecularHessian will, in some cases, fill it in
        appropriately.  However, even in these cases, it may be desirable
        to specify this keyword.  For example, this could be used in an
        optimization to compute frequencies using a lower level of theory.

        <tr><td><tt>debug</tt><td>boolean<td>false<td>If true,
        print out debugging information.

        <tr><td><tt>point_group</tt><td>PointGroup<td>none<td>
        The point group to use for generating the displacements.

        <tr><td><tt>restart</tt><td>boolean<td>true<td>If true, and a
        checkpoint file exists, restart from that file.

        <tr><td><tt>restart_file</tt><td>string
        <td><em>basename</em><tt>.ckpt.hess</tt><td>The name of
        the file where checkpoint information is written to or read from.

        <tr><td><tt>checkpoint</tt><td>boolean<td>false<td>If true,
        checkpoint intermediate data.

        <tr><td><tt>only_totally_symmetric</tt><td>boolean<td>false
        <td>If true, only follow totally symmetric displacments.  The
        hessian will not be complete, but it has enough information
        to use it in a geometry optimization.

        <tr><td><tt>eliminate_quadratic_terms</tt><td>boolean<td>see notes<td>
        If <tt>true</tt>, then contributions to the hessian quadratic in the displacement
        will be eliminated (i.e. the leading-order errors will be <em>quartic</em> in the displacement).
        If implemented in terms of gradients,
        this requires that two displacements are done for each totally symmetric
        coordinate, rather than one. If implemented in terms of energies,
        this requires twice as many displacements as the standard algorithm for each force constant, regardless
        of its symmetry. If using gradients, the default setting is <tt>true</tt>
        (in such case setting this keyword to <tt>false</tt> will produces lower
        accuracy which will only be sufficient for geometry optimizations).
        If using energies, the default setting is <tt>false</tt>.
        Benchmark calculations should always set this to <tt>true</tt>.

        <tr><td><tt>do_null_displacement</tt><td>boolean<td>true<td>Run
        the calculation at the given geometry as well.

        <tr><td><tt>displacement</tt><td>double<td>1.0e-2<td>The size of
        the displacement in Bohr.

        <tr><td><tt>gradient_accuracy</tt><td>double<td><tt>accuracy</tt>
        * <tt>displacement</tt><td>The accuracy to which the gradients will be computed.

        <tr><td><tt>energy_accuracy</tt><td>double<td><tt>accuracy</tt>
        * <tt>displacement</tt>^2<td>The accuracy to which the energies will be computed.

        <tr><td><tt>use_energies</tt><td>boolean<td>false<td>Setting to true will
        force computation from energies.

        </table>
    */
    FinDispMolecularHessian(const Ref<KeyVal>&);
    FinDispMolecularHessian(StateIn&);
    ~FinDispMolecularHessian();
    void save_data_state(StateOut&);

    /** This returns the cartesian hessian.  If it has not yet been
        computed, it will be computed by finite displacements. */
    RefSymmSCMatrix cartesian_hessian();

    void set_energy(const Ref<MolecularEnergy> &energy);
    MolecularEnergy* energy() const { return pimpl_ ? pimpl_->mole().pointer() : 0; }

    const Ref<Params>& params() const { return params_; }

    void set_desired_accuracy(double acc);
};

/** Computes the molecular gradient by finite differences of energies.
    This will use the minimum number of displacements, each in the
    same point group as the reference geometry. */
class FinDispMolecularGradient: public MolecularGradient {
//
//  public:
//
//    class Displacements {
//      public:
//
//      private:
//        /** displacements are represented as a key -> index map
//         *  keys are formatted as comma-separated list of integer,double pairs
//         *  each pair contains the index of SALC (0-based) + displacement size
//         *  SALC indices are nondecreasing
//         */
//        std::map< std::list< std::pair<int,double> >, > disps;
//    };

  protected:
    Ref<MolecularEnergy> mole_;
    // In case molecule must be given in lower symmetry, its actual
    // symmetry and the symmetry used to compute displacements is this
    Ref<PointGroup> displacement_point_group_;
    // The molecule's original geometry for restoration at the end and
    //computing displacements.
    RefSCVector original_geometry_;
    // the cartesian displacement size in bohr
    double disp_;
    // the accuracy of the energy computations
    double energy_accuracy_;
    // whether or not to attempt a restart
    int restart_;
    // the name of the restart file
    std::string restart_file_;
    // whether or not to checkpoint
    int checkpoint_;
    // the name of the checkpoint file
    std::string checkpoint_file_;
    // 2-pt formula for gradient is accurate to O(h^2) (contaminated by 3rd derivatives)
    // make it O(h^4), i.e. eliminate the h^2 terms, by using a 4-pt formula
    int eliminate_quadratic_terms_;
    // print flag
    int debug_;
    // a basis for the symmetrized cartesian coordinates
    RefSCMatrix symbasis_;
    // the energies at each of the completed displacements
    std::vector<double> energies_;

    // given displacement # disp returns internal coordinate index and displacement size (in units of disp_)
    void get_disp(int disp, int &index, double &dispsize);
//    void do_grad_for_irrep(int irrep,
//                           const RefSymmSCMatrix &dhessian,
//                           const RefSymmSCMatrix &xhessian);
    void init();
    void restart();

    /** These members are used to compute a cartesian gradient from
        gradients at finite displacements. */
    RefSCVector compute_gradient();
    int ndisplace() const;
    int ndisplacements_done() const { return energies_.size(); }
    RefSCMatrix displacements(int irrep) const;
    void displace(int disp);
    void original_geometry();
    //void set_gradient(int disp, const RefSCVector &grad);
    void checkpoint_displacements(StateOut&);
    void restore_displacements(StateIn&);

  public:
    FinDispMolecularGradient(const Ref<MolecularEnergy>&);
    /** The FinDispMolecularGradient KeyVal constructor is used to generate a
        FinDispMolecularGradient object from the input.  It reads the keywords
        below.

        <table border="1">

        <tr><td>%Keyword<td>Type<td>Default<td>Description

        <tr><td><tt>energy</tt><td>MolecularEnergy<td>none<td>This gives an
        object which will be used to compute the energies needed to form
        the gradient.  If this is not specified, the object using
        FinDispMolecularGradient will, in some cases, fill it in
        appropriately.

        <tr><td><tt>debug</tt><td>boolean<td>false<td>If true,
        print out debugging information.

        <tr><td><tt>point_group</tt><td>PointGroup<td>none<td>
        The point group to use for generating the displacements.

        <tr><td><tt>restart</tt><td>boolean<td>true<td>If true, and a
        checkpoint file exists, restart from that file.

        <tr><td><tt>restart_file</tt><td>string
        <td><em>basename</em><tt>.ckpt.grad</tt><td>The name of
        the file where checkpoint information is written to or read from.

        <tr><td><tt>checkpoint</tt><td>boolean<td>false<td>If true,
        checkpoint intermediate data.

        <tr><td><tt>eliminate_quadratic_terms</tt><td>boolean<td>false<td>
        If true, then the terms quadratic in the displacement will be eliminated (i.e. the error
        in the gradient will be quartic in the displacement size).  This requires
        that four displacements are done for each (totally symmetric)
        coordinate, rather than two.  Setting this to true only makes sense
        for benchmark computations of high precision and should not
        be necessary in routine computations.

        <tr><td><tt>displacement</tt><td>double<td>1.0e-2<td>The size of
        the displacement in Bohr.

        <tr><td><tt>energy_accuracy</tt><td>double<td><tt>accuracy</tt> * <tt>displacement</tt>
        <td>The accuracy to which the energies will be computed.

        </table>
    */
    FinDispMolecularGradient(const Ref<KeyVal>&);
    FinDispMolecularGradient(StateIn&);
    ~FinDispMolecularGradient();
    void save_data_state(StateOut&);


    /** This returns the cartesian gradient.  If it has not yet been
        computed, it will be computed by finite displacements. */
    RefSCVector cartesian_gradient();

    /// Set checkpoint option.
    void set_checkpoint(int c) { checkpoint_ = c; }
    /// Return the current value of the checkpoint option.
    int checkpoint() const { return checkpoint_; }

    void set_energy(const Ref<MolecularEnergy> &energy);
    MolecularEnergy* energy() const;

    Ref<SCMatrixKit> matrixkit() const { return mole_->matrixkit(); }
    RefSCDimension d3natom() const { return mole_->moldim(); }

    void set_desired_accuracy(double acc);

    void set_eliminate_quadratic_terms(bool e) { eliminate_quadratic_terms_ = e; }
    void set_disp_size(double s);
};

/// @}
// end of addtogroup ChemistryMolecule

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
