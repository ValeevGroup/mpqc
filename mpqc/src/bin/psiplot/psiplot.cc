
extern "C" {
#include <math.h>
#include <stdio.h>
#include <cx/DataTypes.h>
#include <cx/DataAccess.h>
#include <cx/DataOps.h>
#include <cx/PortAccess.h>
  void psiplot();
};

#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/shape.h>
#include <chemistry/qc/mpqc/mpqc.h>

static int inited=0;

extern "C" {void __do_global_ctors();}

// orbital values:
//   0 compute total density
//   1 compute natural orbital orbnum's value
//   2 compute natural orbital orbnum's absolute value
//   3 compute natural orbital orbnum's density
//   4 compute orbital orbnum's value
//   5 compute orbital orbnum's absolute value
//   6 compute orbital orbnum's density
//   7 use the density's rms gradient
//   8 the distance to the VDW surface
//   9 the distance to the Connolly surface
void
psiplot(double margin, double incr, int orbital, int orbnum,
        cxLattice**outlat)
{

  int i,j,k,l;

  static MPQC* mpqc;
  static Molecule* mol;
  static Shape* con;
  static Shape* vdw;

  if (!inited) {
      // HACK FOLLOWS (DON'T LOOK)
      // call the gcc specific function that calls the global ctors
      __do_global_ctors();

      inited = 1;
      const char * infile = "mpqc++.in";

      // set up keyval
      printf("Reading from file \"%s\"\n",infile);
      static ParsedKeyVal keyval(infile);
      static PrefixKeyVal pkeyval("mpqc",keyval);

      // set up the molecule
      static Molecule newmol(pkeyval);
      mol = &newmol;
      mol->print();

      // set up the basis set
      static GaussianBasisSet gbs(pkeyval,*mol);
      gbs.print();

      // set up the wavefunction
      static MPQC newmpqc(pkeyval,*mol,gbs);
      mpqc = &newmpqc;
      mpqc->print();

      con = new ConnollyShape(*mol);
      vdw = new VDWShape(*mol);

      //printf("The eigenvectors:\n");
      //Print(mpqc->eigenvectors());
   }

  printf("energy = %12.8f\n",
        mpqc->energy());

  int norb = mpqc->eigenvectors().Ncols();

  if (orbnum < 0) {
      orbnum = 0;
      printf("resetting orbnum to %d\n",orbnum);
    }
  else if (orbnum >= norb) {
      orbnum = norb-1;
      printf("resetting orbnum to %d\n",orbnum);
    }

  // find a box that the molecule will fit in
  double startx;
  double lastx;
  double starty;
  double lasty;
  double startz;
  double lastz;
  startx = (*mol)[0][0];
  lastx  = (*mol)[0][0];
  starty = (*mol)[0][1];
  lasty  = (*mol)[0][1];
  startz = (*mol)[0][2];
  lastz  = (*mol)[0][2];
  for (i=1; i<mol->natom(); i++) {
      if (startx > (*mol)[i][0]) startx = (*mol)[i][0];
      else if (lastx  < (*mol)[i][0]) lastx  = (*mol)[i][0];
      if (starty > (*mol)[i][1]) starty = (*mol)[i][1];
      else if (lasty  < (*mol)[i][1]) lasty  = (*mol)[i][1];
      if (startz > (*mol)[i][2]) startz = (*mol)[i][2];
      else if (lastz  < (*mol)[i][2]) lastz  = (*mol)[i][2];
    }
  startx -= margin;
  lastx  += margin;
  starty -= margin;
  lasty  += margin;
  startz -= margin;
  lastz  += margin;
  //const double incr = 0.5;
  printf("the x box size is % 12.8f % 12.8f\n",startx,lastx);
  printf("the y box size is % 12.8f % 12.8f\n",starty,lasty);
  printf("the z box size is % 12.8f % 12.8f\n",startz,lastz);

  int nx,ny,nz;
  nx = (int)((lastx - startx)/incr);
  ny = (int)((lasty - starty)/incr);
  nz = (int)((lastz - startz)/incr);

  cxLattice *lattice;

  // create the lattice
  long int dims[3];
  dims[0] = nx; dims[1] = ny; dims[2] = nz;
  int ndata;
  ndata = 1;
  lattice = cxLatNew(3,dims,ndata,cx_prim_float,3,cx_coord_uniform);
  if (cxDataAllocErrorGet()) return;
  float *data;
  data = lattice->data->d.cx_prim_float.values;
  lattice->coord->c.cx_coord_uniform.bBox[0] = startx;
  lattice->coord->c.cx_coord_uniform.bBox[1] = lastx;
  lattice->coord->c.cx_coord_uniform.bBox[2] = starty;
  lattice->coord->c.cx_coord_uniform.bBox[3] = lasty;
  lattice->coord->c.cx_coord_uniform.bBox[4] = startz;
  lattice->coord->c.cx_coord_uniform.bBox[5] = lastz;
  //float *coord;
  //cxLatPtrGet(lattice,0,(void**)&data,0,(void**)&coord);

  cart_point r;
  //double nelec = 0.0;
  const double volume = incr*incr*incr;
  //int icoord=0;
  //int idata=0;
  for (i=0,r[0]=startx; i<nx; i++,r[0]+=incr) {
      for (j=0,r[1]=starty; j<ny; j++,r[1]+=incr) {
          for (k=0,r[2]=startz; k<nz; k++,r[2]+=incr) {
              double value;
              double orbval;
              double rmsgrad;
              //nelec +=  dens * volume;

              if (orbital == 1) {
                  value = mpqc->natural_orbital(r,orbnum);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
              else if (orbital == 2) {
                  value = mpqc->natural_orbital(r,orbnum);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = fabs(value);
                }
              else if (orbital == 3) {
                  value = mpqc->natural_orbital_density(r,orbnum,&orbval);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
              else if (orbital == 4) {
                  value = mpqc->orbital(r,orbnum);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
              else if (orbital == 5) {
                  value = mpqc->orbital(r,orbnum);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = fabs(value);
                }
              else if (orbital == 6) {
                  value = mpqc->orbital_density(r,orbnum,&orbval);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
              else if (orbital == 7) {
                  double grad[3];
                  value = mpqc->density_gradient(r,grad);
                  rmsgrad = sqrt(grad[0]*grad[0]
                                 +grad[1]*grad[1]
                                 +grad[2]*grad[2]);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = rmsgrad;
                }
              else if (orbital == 8) {
                  value = vdw->distance_to_surface(r);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
              else if (orbital == 9) {
                  value = con->distance_to_surface(r);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
              else {
                  value = mpqc->density(r);
                  data[k*nx*ny*ndata + j*nx*ndata + i*ndata] = value;
                }
            }
        }
    }
  //printf("nelec = %12.8f volume element = %12.8f\n",nelec,volume);

  *outlat = lattice;

}

