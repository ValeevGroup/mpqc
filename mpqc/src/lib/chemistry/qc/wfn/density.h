//
// density.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_wfn_density_h
#define _chemistry_qc_wfn_density_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/volume.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/molecule/molrender.h>

class ElectronDensity: public Volume {
  protected:
    Ref<Wavefunction> wfn_;
    virtual void compute();
  public:
    ElectronDensity(const Ref<KeyVal>&);
    ElectronDensity(const Ref<Wavefunction>&);
    ~ElectronDensity();
    virtual void boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2);
};

class DensityColorizer: public MoleculeColorizer {
  protected:
    Ref<Wavefunction> wfn_;
    double scale_;
    double reference_;
    int have_scale_;
    int have_reference_;
  public:
    DensityColorizer(const Ref<KeyVal>&);
    ~DensityColorizer();

    void colorize(const Ref<RenderedPolygons> &);
};

class GradDensityColorizer: public MoleculeColorizer {
  protected:
    Ref<Wavefunction> wfn_;
    double scale_;
    double reference_;
    int have_scale_;
    int have_reference_;
  public:
    GradDensityColorizer(const Ref<KeyVal>&);
    ~GradDensityColorizer();

    void colorize(const Ref<RenderedPolygons> &);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
