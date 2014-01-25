//
// molrender.cc
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

#include <math.h>

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/render/sphere.h>
#include <util/render/polygons.h>
#include <util/render/polylines.h>
#include <util/render/color.h>
#include <chemistry/molecule/molrender.h>
#include <math/scmat/vector3.h>

using namespace sc;

////////////////////////////////////////////////////////////////
// RenderedMolecule

static ClassDesc RenderedMolecule_cd(
  typeid(RenderedMolecule),"RenderedMolecule",1,"public RenderedObject",
  0, 0, 0);

RenderedMolecule::RenderedMolecule(const Ref<KeyVal>& keyval):
  RenderedObject(keyval)
{
  mol_ << keyval->describedclassvalue("molecule");
  atominfo_ << keyval->describedclassvalue("atominfo");
  if (atominfo_.null()) {
      atominfo_ = new AtomInfo();
    }

  if (mol_.null()) {
      throw InputError("missing required input of type Molecule",
                       __FILE__, __LINE__, "molecule", 0,
                       class_desc());
    }
}

RenderedMolecule::~RenderedMolecule()
{
}

void
RenderedMolecule::render(const Ref<Render>& render)
{
  object_->render(render);
}

////////////////////////////////////////////////////////////////
// RenderedBallMolecule

static ClassDesc RenderedBallMolecule_cd(
  typeid(RenderedBallMolecule),"RenderedBallMolecule",1,"public RenderedMolecule",
  0, create<RenderedBallMolecule>, 0);

RenderedBallMolecule::RenderedBallMolecule(const Ref<KeyVal>& keyval):
  RenderedMolecule(keyval)
{
  init();
}

RenderedBallMolecule::~RenderedBallMolecule()
{
}

void
RenderedBallMolecule::init()
{
  Ref<RenderedObjectSet> set = new RenderedObjectSet;

  for (int i=0; i<mol_->natom(); i++) {
      Ref<RenderedObject> atom = new RenderedSphere;

      int Z = mol_->Z(i);

      Ref<Material> material = new Material;
      Color color(atominfo_->red(Z),
                  atominfo_->green(Z),
                  atominfo_->blue(Z));
      material->diffuse().set(color);
      material->ambient().set(color);

      Ref<Transform> transform = new Transform;
      transform->scale(atominfo_->vdw_radius(Z));
      transform->translate(mol_->r(i,0), mol_->r(i,1), mol_->r(i,2));

      atom->material(material);
      atom->transform(transform);

      set->add(atom);
    }

  object_ = set.pointer();
}

////////////////////////////////////////////////////////////////
// RenderedStickMolecule

static ClassDesc RenderedStickMolecule_cd(
  typeid(RenderedStickMolecule),"RenderedStickMolecule",1,"public RenderedMolecule",
  0, create<RenderedStickMolecule>, 0);

RenderedStickMolecule::RenderedStickMolecule(const Ref<KeyVal>& keyval):
  RenderedMolecule(keyval)
{
  use_color_ = keyval->booleanvalue("color");
  if (keyval->error() != KeyVal::OK) use_color_ = 1;
  init();
}

RenderedStickMolecule::~RenderedStickMolecule()
{
}

static int
bonding(const Ref<Molecule>& m, const Ref<AtomInfo>& a, int i, int j)
{
  SCVector3 ri(m->r(i));
  SCVector3 rj(m->r(j));
  double maxbonddist = 1.1*(m->atominfo()->atomic_radius(m->Z(i))
                            +m->atominfo()->atomic_radius(m->Z(j)));
  SCVector3 r(ri-rj);
  if (r.dot(r) <= maxbonddist*maxbonddist) return 1;
  return 0;
}

void
RenderedStickMolecule::init()
{
  int i,j;
  int nbonds;
  int natoms = mol_->natom();
  
  Ref<RenderedPolylines> o = new RenderedPolylines;

  // count the number of bonds
  nbonds = 0;
  for (i=0; i<natoms; i++) {
      for (j=0; j<i; j++) {
          if (bonding(mol_, atominfo_, i, j)) nbonds++;
        }
    }

  int nvertex = natoms;
  if (use_color_) nvertex += 2*nbonds;

  // initialize the polylines
  o->initialize(nvertex, nbonds, RenderedPolylines::Vertex);

  // put the atoms in the vertex list
  for (i=0; i<natoms; i++) {
      o->set_vertex(i,
                    mol_->r(i,0),
                    mol_->r(i,1),
                    mol_->r(i,2));
      if (use_color_) {
          int Z = mol_->Z(i);
          o->set_vertex_rgb(i,
                            atominfo_->red(Z),
                            atominfo_->green(Z),
                            atominfo_->blue(Z));
        }
      else {
          o->set_vertex_rgb(i, 0.0, 0.0, 0.0);
        }
    }

  // put the bonds in the line list
  nbonds = 0;
  int ibonds2 = natoms;
  for (i=0; i<natoms; i++) {
      SCVector3 ri(mol_->r(i));
      int Zi = mol_->Z(i);
      for (j=0; j<i; j++) {
          if (bonding(mol_, atominfo_, i, j)) {
              if (use_color_) {
                  SCVector3 rj(mol_->r(j));
                  int Zj = mol_->Z(j);
                  SCVector3 v = 0.5*(ri+rj);
                  o->set_vertex(ibonds2, v.x(), v.y(), v.z());
                  o->set_vertex_rgb(ibonds2,
                                    atominfo_->red(Zi),
                                    atominfo_->green(Zi),
                                    atominfo_->blue(Zi));
                  o->set_vertex(ibonds2+1, v.x(), v.y(), v.z());
                  o->set_vertex_rgb(ibonds2+1,
                                    atominfo_->red(Zj),
                                    atominfo_->green(Zj),
                                    atominfo_->blue(Zj));
                  o->set_polyline(nbonds, i, ibonds2, ibonds2+1, j);
                  ibonds2 += 2;
                }
              else {
                  o->set_polyline(nbonds, i, j);
                }
              nbonds++;
            }
        }
    }

  object_ = o.pointer();
}

////////////////////////////////////////////////////////////////
// RenderedMolecularSurface

static ClassDesc RenderedMolecularSurface_cd(
  typeid(RenderedMolecularSurface),"RenderedMolecularSurface",1,"public RenderedMolecule",
  0, create<RenderedMolecularSurface>, 0);

RenderedMolecularSurface::RenderedMolecularSurface(const Ref<KeyVal>& keyval):
  RenderedMolecule(keyval)
{
  surf_ << keyval->describedclassvalue("surface");
  colorizer_ << keyval->describedclassvalue("colorizer");
  if (colorizer_.null())
      colorizer_ = new AtomProximityColorizer(mol_,atominfo_);
  init(0);
}

RenderedMolecularSurface::~RenderedMolecularSurface()
{
}

void
RenderedMolecularSurface::init()
{
  init(1);
}

void
RenderedMolecularSurface::init(int reinit_surf)
{
  int i, ij, j;
  if (reinit_surf || !surf_->inited()) surf_->init();
  int nvertex = surf_->nvertex();
  int ntriangle = surf_->ntriangle();
  int natom = mol_->natom();

  Ref<RenderedPolygons> o = new RenderedPolygons;

  o->initialize(nvertex, ntriangle, RenderedPolygons::Vertex);

  // extract the atomic positions and colors into an array for rapid access
  double *axyz = new double[3*natom];
  double *argb = new double[3*natom];
  double *arad = new double[natom];
  ij = 0;
  for (i=0; i<natom; i++) {
      int Z = mol_->Z(i);
      arad[i] = atominfo_->vdw_radius(Z);
      for (j=0; j<3; j++,ij++) {
          axyz[ij] = mol_->r(i,j);
          argb[ij] = atominfo_->rgb(Z, j);
        }
    }

  for (i=0; i<nvertex; i++) {
      const SCVector3& v = surf_->vertex(i)->point();
      double x = v[0];
      double y = v[1];
      double z = v[2];
      o->set_vertex(i, x, y, z);
    }
  colorizer_->colorize(o);

  delete[] axyz;
  delete[] argb;
  delete[] arad;

  for (i=0; i<ntriangle; i++) {
      o->set_face(i,
                  surf_->triangle_vertex(i,0),
                  surf_->triangle_vertex(i,1),
                  surf_->triangle_vertex(i,2));
    }

  object_ = o.pointer();
}

/////////////////////////////////////////////////////////////////////////////
// MoleculeColorizer

static ClassDesc MoleculeColorizer_cd(
  typeid(MoleculeColorizer),"MoleculeColorizer",1,"public DescribedClass",
  0, 0, 0);

MoleculeColorizer::MoleculeColorizer(const Ref<Molecule>&mol)
{
  mol_ = mol;
}

MoleculeColorizer::MoleculeColorizer(const Ref<KeyVal>&keyval)
{
  mol_ << keyval->describedclassvalue("molecule");
}

MoleculeColorizer::~MoleculeColorizer()
{
}

/////////////////////////////////////////////////////////////////////////////
// AtomProximityColorizer

static ClassDesc AtomProximityColorizer_cd(
  typeid(AtomProximityColorizer),"AtomProximityColorizer",1,"public MoleculeColorizer",
  0, create<AtomProximityColorizer>, 0);

AtomProximityColorizer::AtomProximityColorizer(const Ref<Molecule> &mol,
                                               const Ref<AtomInfo> &ai):
  MoleculeColorizer(mol)
{
  atominfo_ = ai;
}

AtomProximityColorizer::AtomProximityColorizer(const Ref<KeyVal>&keyval):
  MoleculeColorizer(keyval)
{
  atominfo_ << keyval->describedclassvalue("atominfo");
  if (atominfo_.null()) {
      atominfo_ = new AtomInfo();
    }
}

AtomProximityColorizer::~AtomProximityColorizer()
{
}

static void
compute_color(int n, double* axyz, double* argb, double* arad,
              double x, double y, double z, Color& c)
{
  int i, j;
  const int maxclosest = 10;
  int closest[maxclosest];
  double distance2[maxclosest]; // the distance squared - radius squared
  int nclosest;

  if (n == 0) {
      c.set_rgb(1.0, 1.0, 1.0);
      return;
    }

  // find the closest atoms
  nclosest = 0;
  for (i=0; i<n; i++) {
      SCVector3 r(axyz[0] - x, axyz[1] - y, axyz[2] - z);
      double tmpdist2 = r.dot(r) - arad[i]*arad[i];
//       if (tmpdist2 < 1.e-6) {
//           c.set_rgb(argb[3*i], argb[3*i+1], argb[3*i+2]);
//           return;
//         }
      if (tmpdist2 < 0.0) tmpdist2 = 0.0;
      for (j=nclosest-1; j>=0; j--) {
          if (distance2[j] <= tmpdist2) break;
          if (j+1 < maxclosest) {
              distance2[j+1] = distance2[j];
              closest[j+1] = closest[j];
            }
        }
      if (j+1 < maxclosest) {
          distance2[j+1] = tmpdist2;
          closest[j+1] = i;
          if (maxclosest > nclosest) nclosest++;
        }
      axyz += 3;
    }

  if (nclosest == 1) {
      c.set_rgb(argb[3*closest[0]],argb[3*closest[0]]+1,argb[3*closest[0]]+2);
      return;
    }

  // average the colors of the closest atoms
  for (i=0; i<nclosest; i++) distance2[i] = sqrt(distance2[i]);
  for (i=1; i<nclosest; i++) distance2[i] -= distance2[0];
  distance2[0] = 0.0;
  for (i=0; i<nclosest; i++) distance2[i] = 2.0*exp(-distance2[i]);
  double sum = 0.0;
  for (i=0; i<nclosest; i++) sum += distance2[i];
  sum = 1.0/sum;
  for (i=0; i<nclosest; i++) distance2[i] *= sum;

  double rgb[3] = {0.0, 0.0, 0.0};
  for (i=0; i<nclosest; i++) {
      for (j=0; j<3; j++) {
          rgb[j] += distance2[i]*argb[3*closest[i] + j];
        }
    }
  c.set_rgb(rgb[0], rgb[1], rgb[2]);
}

void
AtomProximityColorizer::colorize(const Ref<RenderedPolygons> &poly)
{
  int natom = mol_->natom();
  int nvertex = poly->nvertex();

  int i,j,ij;
  // extract the atomic positions and colors into an array for rapid access
  double *axyz = new double[3*natom];
  double *argb = new double[3*natom];
  double *arad = new double[natom];
  ij = 0;
  for (i=0; i<natom; i++) {
      int Z = mol_->Z(i);
      arad[i] = atominfo_->vdw_radius(Z);
      for (j=0; j<3; j++,ij++) {
          axyz[ij] = mol_->r(i,j);
          argb[ij] = atominfo_->rgb(Z, j);
        }
    }

  for (i=0; i<nvertex; i++) {
      const double *v = poly->vertex(i);
      double x = v[0];
      double y = v[1];
      double z = v[2];
      Color c;
      compute_color(natom, axyz, argb, arad, x, y, z, c);
      poly->set_vertex_color(i, c);
    }

  delete[] axyz;
  delete[] argb;
  delete[] arad;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
