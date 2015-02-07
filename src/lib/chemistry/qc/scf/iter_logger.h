//
// iter_logger.h
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: Jun 13, 2014
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


#ifndef _lib_chemistry_qc_scf_h
#define _lib_chemistry_qc_scf_h

#include <functional>
#include <utility>
#include <chemistry/qc/scf/scf.h>

#include <util/misc/xml.h>

namespace sc {

/////////////////////////////////////////////////////////////////////////////

class SCFIterationData;

class SCFIterationLogger : public XMLWritable, public DescribedClass {
  public:

    typedef std::function<void(boost::property_tree::ptree&, const XMLWriter&)> element_write_function;

  private:

    bool log_evals_;
    bool log_density_;
    bool log_coeffs_;

    std::vector<element_write_function> other_details_;
    std::vector<std::vector<element_write_function>> other_iter_details_;

    std::vector<SCFIterationData> iterations_;

  public:

    SCFIterationLogger(const Ref<KeyVal>& keyval);

    boost::property_tree::ptree& write_xml(boost::property_tree::ptree& parent, const XMLWriter& writer);

    // TODO handle alpha and beta for evals and density

    void new_iteration();

    void log_density(RefSymmSCMatrix density, SpinCase1 spin_case=AnySpinCase1);

    void log_evals(RefDiagSCMatrix evals, SpinCase1 spin_case=AnySpinCase1);

    void log_coeffs(RefSCMatrix evals, SpinCase1 spin_case=AnySpinCase1);

    bool log_coeffs_enabled(){ return log_coeffs_; }

    template <typename Func, typename... Args>
    void log_global_misc(Func&& func, Args&&... args) {
      other_details_.emplace_back(std::bind(
          func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...
      ));
    }

    template <typename Func, typename... Args>
    void log_iter_misc(Func&& func, Args&&... args) {
      other_iter_details_.back().emplace_back(std::bind(
          func, std::placeholders::_1, std::placeholders::_2, std::forward<Args>(args)...
      ));
    }

};

/////////////////////////////////////////////////////////////////////////////

class SCFIterationData : public XMLWritable {

  public:

    SCFIterationLogger* parent;
    int number;
    RefDiagSCMatrix evals = 0;
    RefDiagSCMatrix alpha_evals = 0;
    RefDiagSCMatrix beta_evals = 0;
    RefSymmSCMatrix density = 0;
    RefSymmSCMatrix alpha_density = 0;
    RefSymmSCMatrix beta_density = 0;
    RefSCMatrix coeffs = 0;
    RefSCMatrix alpha_coeffs = 0;
    RefSCMatrix beta_coeffs = 0;


    boost::property_tree::ptree& write_xml(boost::property_tree::ptree& parent, const XMLWriter& writer);

};


}


#endif /* _lib_chemistry_qc_scf_h */
