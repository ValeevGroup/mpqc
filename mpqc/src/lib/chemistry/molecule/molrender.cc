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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <util/render/sphere.h>
#include <util/render/polygons.h>
#include <util/render/polylines.h>
#include <util/render/color.h>
#include <chemistry/molecule/molrender.h>
#include <math/scmat/vector3.h>

////////////////////////////////////////////////////////////////
// RenderedMolecule

#define CLASSNAME RenderedMolecule
#define PARENTS public RenderedObject
#include <util/class/classia.h>
void *
RenderedMolecule::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedObject::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedMolecule::RenderedMolecule(const RefKeyVal& keyval):
  RenderedObject(keyval),
  atominfo_(keyval->describedclassvalue("atominfo")),
  mol_(keyval->describedclassvalue("molecule"))
{
  if (atominfo_.null()) {
      atominfo_ = new AtomInfo();
    }

  if (mol_.null()) {
      cerr << node0 << indent
           << "RenderedMolecule: no \"molecule\" in keyval\n";
      abort();
    }
}

RenderedMolecule::~RenderedMolecule()
{
}

void
RenderedMolecule::render(const RefRender& render)
{
  object_->render(render);
}

////////////////////////////////////////////////////////////////
// RenderedBallMolecule

#define CLASSNAME RenderedBallMolecule
#define PARENTS public RenderedMolecule
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
RenderedBallMolecule::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedMolecule::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedBallMolecule::RenderedBallMolecule(const RefKeyVal& keyval):
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
  RefRenderedObjectSet set = new RenderedObjectSet;

  for (int i=0; i<mol_->natom(); i++) {
      RefRenderedObject atom = new RenderedSphere;

      AtomicCenter& atomcent = mol_->operator[](i);
      ChemicalElement& chemelem = atomcent.element();

      RefMaterial material = new Material;
      Color color(atominfo_->red(chemelem),
                  atominfo_->green(chemelem),
                  atominfo_->blue(chemelem));
      material->diffuse().set(color);
      material->ambient().set(color);

      RefTransform transform = new Transform;
      transform->scale(atominfo_->radius(chemelem));
      transform->translate(atomcent[0], atomcent[1], atomcent[2]);

      atom->material(material);
      atom->transform(transform);

      set->add(atom);
    }

  object_ = set.pointer();
}

////////////////////////////////////////////////////////////////
// RenderedStickMolecule

#define CLASSNAME RenderedStickMolecule
#define PARENTS public RenderedMolecule
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
RenderedStickMolecule::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedMolecule::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedStickMolecule::RenderedStickMolecule(const RefKeyVal& keyval):
  RenderedMolecule(keyval)
{
  init();
}

RenderedStickMolecule::~RenderedStickMolecule()
{
}

static int
bonding(const RefMolecule& m, const RefAtomInfo& a, int i, int j)
{
  SCVector3 ri(m->atom(i)[0], m->atom(i)[1], m->atom(i)[2]);
  SCVector3 rj(m->atom(j)[0], m->atom(j)[1], m->atom(j)[2]);
  //double maxbonddist = 1.1*(a->radius(m->atom(i).element())
  //                          +a->radius(m->atom(j).element()));
  double maxbonddist = 1.1*(m->atom(i).element().atomic_radius()
                            +m->atom(j).element().atomic_radius());
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
  
  RefRenderedPolylines o = new RenderedPolylines;

  // count the number of bonds
  nbonds = 0;
  for (i=0; i<natoms; i++) {
      for (j=0; j<i; j++) {
          if (bonding(mol_, atominfo_, i, j)) nbonds++;
        }
    }

  // initialize the polylines
  o->initialize(natoms+2*nbonds, nbonds, RenderedPolylines::Vertex);

  // put the atoms in the vertex list
  for (i=0; i<natoms; i++) {
      o->set_vertex(i,
                    mol_->atom(i)[0],
                    mol_->atom(i)[1],
                    mol_->atom(i)[2]);
      o->set_vertex_rgb(i,
                        atominfo_->red(mol_->atom(i).element()),
                        atominfo_->green(mol_->atom(i).element()),
                        atominfo_->blue(mol_->atom(i).element()));
    }

  // put the bonds in the line list
  nbonds = 0;
  int ibonds2 = natoms;
  for (i=0; i<natoms; i++) {
      SCVector3 ri(mol_->atom(i)[0],
                   mol_->atom(i)[1],
                   mol_->atom(i)[2]);
      for (j=0; j<i; j++) {
          if (bonding(mol_, atominfo_, i, j)) {
              SCVector3 rj(mol_->atom(j)[0],
                           mol_->atom(j)[1],
                           mol_->atom(j)[2]);
              SCVector3 v = 0.5*(ri+rj);
              o->set_vertex(ibonds2, v.x(), v.y(), v.z());
              o->set_vertex_rgb(ibonds2,
                                atominfo_->red(mol_->atom(i).element()),
                                atominfo_->green(mol_->atom(i).element()),
                                atominfo_->blue(mol_->atom(i).element()));
              o->set_vertex(ibonds2+1, v.x(), v.y(), v.z());
              o->set_vertex_rgb(ibonds2+1,
                                atominfo_->red(mol_->atom(j).element()),
                                atominfo_->green(mol_->atom(j).element()),
                                atominfo_->blue(mol_->atom(j).element()));
              o->set_polyline(nbonds, i, ibonds2, ibonds2+1, j);
              ibonds2 += 2;
              nbonds++;
            }
        }
    }

  object_ = o.pointer();
}

////////////////////////////////////////////////////////////////
// RenderedMolecularSurface

#define CLASSNAME RenderedMolecularSurface
#define PARENTS public RenderedMolecule
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
RenderedMolecularSurface::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RenderedMolecule::_castdown(cd);
  return do_castdowns(casts,cd);
}

RenderedMolecularSurface::RenderedMolecularSurface(const RefKeyVal& keyval):
  RenderedMolecule(keyval)
{
  surf_ = keyval->describedclassvalue("surface");
  colorizer_ = keyval->describedclassvalue("colorizer");
  if (colorizer_.null())
      colorizer_ = new AtomProximityColorizer(mol_,atominfo_);
  init();
}

RenderedMolecularSurface::~RenderedMolecularSurface()
{
}

void
RenderedMolecularSurface::init()
{
  int i, ij, j;
  surf_->init();
  int nvertex = surf_->nvertex();
  int ntriangle = surf_->ntriangle();
  int natom = mol_->natom();

  RefRenderedPolygons o = new RenderedPolygons;

  o->initialize(nvertex, ntriangle, RenderedPolygons::Vertex);

  // extract the atomic positions and colors into an array for rapid access
  double *axyz = new double[3*natom];
  double *argb = new double[3*natom];
  double *arad = new double[natom];
  ij = 0;
  for (i=0; i<natom; i++) {
      ChemicalElement& element = mol_->atom(i).element();
      arad[i] = atominfo_->radius(element);
      for (j=0; j<3; j++,ij++) {
          axyz[ij] = mol_->atom(i)[j];
          argb[ij] = atominfo_->rgb(element, j);
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

#define CLASSNAME MoleculeColorizer
#define PARENTS public DescribedClass
#include <util/class/classia.h>
void *
MoleculeColorizer::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

MoleculeColorizer::MoleculeColorizer(const RefMolecule&mol)
{
  mol_ = mol;
}

MoleculeColorizer::MoleculeColorizer(const RefKeyVal&keyval)
{
  mol_ = keyval->describedclassvalue("molecule");
}

MoleculeColorizer::~MoleculeColorizer()
{
}

/////////////////////////////////////////////////////////////////////////////
// AtomProximityColorizer

#define CLASSNAME AtomProximityColorizer
#define PARENTS public MoleculeColorizer
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
AtomProximityColorizer::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = MoleculeColorizer::_castdown(cd);
  return do_castdowns(casts,cd);
}

AtomProximityColorizer::AtomProximityColorizer(const RefMolecule &mol,
                                               const RefAtomInfo &ai):
  MoleculeColorizer(mol)
{
  atominfo_ = ai;
}

AtomProximityColorizer::AtomProximityColorizer(const RefKeyVal&keyval):
  MoleculeColorizer(keyval)
{
  atominfo_ = keyval->describedclassvalue("atominfo");
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
AtomProximityColorizer::colorize(const RefRenderedPolygons &poly)
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
      ChemicalElement& element = mol_->atom(i).element();
      arad[i] = atominfo_->radius(element);
      for (j=0; j<3; j++,ij++) {
          axyz[ij] = mol_->atom(i)[j];
          argb[ij] = atominfo_->rgb(element, j);
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
