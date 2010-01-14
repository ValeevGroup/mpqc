
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_extrap_h
#define _chemistry_qc_lmp2_extrap_h

#include <math/optimize/scextrap.h>
#include <chemistry/qc/lmp2/sma.h>

namespace sc {

namespace sma2 {

/** \brief This permits an Array<2> and an Array<4> to be used with
    SelfConsistentExtrapolation derivatives. */
class Array24SCExtrapData: public sc::SCExtrapData {
  private:
    Array<2> array2_;
    Array<4> array4_;
    double bound_;
    bool distrib2_;
    bool distrib4_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array24SCExtrapData(Array<2>&, bool distrib2,
                        Array<4>&, bool distrib4,
                        const sc::Ref<sc::MessageGrp>&msg);
    ~Array24SCExtrapData();
    void save_data_state(sc::StateOut&);

    sc::SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const sc::Ref<sc::SCExtrapData>&);
    void update(Array<2> &array2, Array<4> &array4);
};

/** \brief This permits an Array<2> and an Array<4> to be used with
    SelfConsistentExtrapolation derivatives. */
class Array24SCExtrapError: public sc::SCExtrapError {
  private:
    Array<2> array2_;
    Array<4> array4_;
    double bound_;
    bool distrib2_;
    bool distrib4_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array24SCExtrapError(Array<2>&, bool distrib2,
                         Array<4>&, bool distrib4,
                         const sc::Ref<sc::MessageGrp>&msg);
    ~Array24SCExtrapError();
    void save_data_state(sc::StateOut&);
    
    double error();
    double scalar_product(const sc::Ref<sc::SCExtrapError>&);
};

/** \brief This permits Array<2>'s to be used with
    SelfConsistentExtrapolation derivatives. */
class Array2SCExtrapData: public sc::SCExtrapData {
  private:
    Array<2> array_;
    double bound_;
    bool distrib_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array2SCExtrapData(Array<2>&, bool distrib,
                       const sc::Ref<sc::MessageGrp>&msg);
    ~Array2SCExtrapData();
    void save_data_state(sc::StateOut&);

    sc::SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const sc::Ref<sc::SCExtrapData>&);
    void update(Array<2> &array);
};

/** \brief This permits Array<2>'s to be used with
    SelfConsistentExtrapolation derivatives. */
class Array2SCExtrapError: public sc::SCExtrapError {
  private:
    Array<2> array_;
    double bound_;
    bool distrib_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array2SCExtrapError(Array<2>&,bool distrib,
                        const sc::Ref<sc::MessageGrp>&msg);
    ~Array2SCExtrapError();
    void save_data_state(sc::StateOut&);
    
    double error();
    double scalar_product(const sc::Ref<sc::SCExtrapError>&);
};

/** \brief This permits Array<4>'s to be used with SelfConsistentExtrapolation
    derivatives. */
class Array4SCExtrapData: public sc::SCExtrapData {
  private:
    Array<4> array_;
    double bound_;
    bool distrib_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array4SCExtrapData(Array<4>&,bool distrib,
                       const sc::Ref<sc::MessageGrp>&msg);
    ~Array4SCExtrapData();
    void save_data_state(sc::StateOut&);

    sc::SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const sc::Ref<sc::SCExtrapData>&);
    void update(Array<4> &array);
};

/** \brief This permits Array<4>'s to be used with
    SelfConsistentExtrapolation derivatives. */
class Array4SCExtrapError: public sc::SCExtrapError {
  private:
    Array<4> array_;
    double bound_;
    bool distrib_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array4SCExtrapError(Array<4>&,bool distrib,
                        const sc::Ref<sc::MessageGrp>&msg);
    ~Array4SCExtrapError();
    void save_data_state(sc::StateOut&);
    
    double error();
    double scalar_product(const sc::Ref<sc::SCExtrapError>&);
};

/** \brief This permits Array<6>'s to be used with
    SelfConsistentExtrapolation derivatives. */
class Array6SCExtrapData: public sc::SCExtrapData {
  private:
    Array<6> array_;
    double bound_;
    bool distrib_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array6SCExtrapData(Array<6>&,bool distrib,
                       const sc::Ref<sc::MessageGrp>&msg);
    ~Array6SCExtrapData();
    void save_data_state(sc::StateOut&);

    sc::SCExtrapData* copy();
    void zero();
    void accumulate_scaled(double, const sc::Ref<sc::SCExtrapData>&);
    void update(Array<6> &array);
};

/** \brief This permits Array<6>'s to be used with
    SelfConsistentExtrapolation derivatives. */
class Array6SCExtrapError: public sc::SCExtrapError {
  private:
    Array<6> array_;
    double bound_;
    bool distrib_;
    sc::Ref<sc::MessageGrp> msg_;
  public:
    Array6SCExtrapError(Array<6>&,bool distrib,
                        const sc::Ref<sc::MessageGrp>&msg);
    ~Array6SCExtrapError();
    void save_data_state(sc::StateOut&);
    
    double error();
    double scalar_product(const sc::Ref<sc::SCExtrapError>&);
};

}

}

#endif
