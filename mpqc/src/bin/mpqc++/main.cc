
#include <util/options/GetLongOpt.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/molecule/coor.h>
#include <chemistry/qc/mpqc/mpqc.h>
#include <chemistry/molecule/shape.h>

main(int argc,char** argv)
{
  // process the options
  GetLongOpt option;
  option.enroll("help",GetLongOpt::NoValue,
		"print this option summary",0);
  option.enroll("input", GetLongOpt::MandatoryValue,
		"read input from file $val", "mpqc++.in");
  int optind = option.parse(argc,argv);
  if (optind < 1) return -1;
  const char *infile = option.retrieve("input");
  if ( option.retrieve("help") ) {
      option.usage();
      return 0;
    }
  FILE*outfp = stdout;

  // set up keyval
  fprintf(outfp,"Reading from file \"%s\"\n",infile);
  ParsedKeyVal keyval(infile);
  PrefixKeyVal pkeyval("mpqc",keyval);

  // set up the molecule
  Molecule mol(pkeyval);

  // set up the basis set
  GaussianBasisSet gbs(pkeyval,mol);
  //gbs.print();

  // set up the molecular coordinate transformation (if any)
  //ProjectedCartesian trans(mol);

  // set up the wavefunction
  //MPQC mpqc(pkeyval,mol,gbs,trans);
  MPQC mpqc(pkeyval,mol,gbs);

  mpqc.print();
  mpqc.do_exchange_energy(1);

  int do_opt = pkeyval.booleanvalue("opt");
  if (pkeyval.error() != KeyVal::OK) do_opt = 0;
  int do_dens = pkeyval.booleanvalue("dens");
  if (pkeyval.error() != KeyVal::OK) do_dens = 0;
  int do_dens_grad = pkeyval.booleanvalue("dens_grad");
  if (pkeyval.error() != KeyVal::OK) do_dens_grad = 0;
  int do_shape = pkeyval.booleanvalue("shape");
  if (pkeyval.error() != KeyVal::OK) do_shape = 0;

  if (do_shape) {
      RefShape con;
      RefShape vdw;
      double radius = pkeyval.doublevalue("radius");
      if (pkeyval.error() != KeyVal::OK) radius = 2.0;
      vdw = new VDWShape(mol);
      con = new ConnollyShape(mol,radius);

      int singlepoint = pkeyval.booleanvalue("singlepoint");
      if (pkeyval.error() != KeyVal::OK) singlepoint = 0;
      if (singlepoint) {
          Point3 p;
          p[0] = pkeyval.doublevalue("point",0);
          p[1] = pkeyval.doublevalue("point",1);
          p[2] = pkeyval.doublevalue("point",2);
          con->is_outside(p);
          con->distance_to_surface(p);
        }

      int singlepoints = pkeyval.booleanvalue("singlepoints");
      if (pkeyval.error() != KeyVal::OK) singlepoints = 0;
      if (singlepoints) {
          while(1) {
              Point3 p;
              printf("x y z: "); fflush(stdout);
              scanf("%lf",&p[0]);
              scanf("%lf",&p[1]);
              scanf("%lf",&p[2]);
              printf("distance = %f\n",con->distance_to_surface(p));
            }
        }

      double incr = pkeyval.doublevalue("incr");
      if (pkeyval.error() != KeyVal::OK) incr = 1.0;
      int dovdw = pkeyval.booleanvalue("vdw");
      if (pkeyval.error() != KeyVal::OK) dovdw = 0;
      int docon = pkeyval.booleanvalue("con");
      if (pkeyval.error() != KeyVal::OK) docon = 0;

      double x0[3];
      if (pkeyval.exists("x")) {
          for (int i=0; i<3; i++) x0[i] = pkeyval.doublevalue("x",i);
        }

    next_x:
      cart_point r;
      printf("x: "); fflush(stdout);
      scanf("%lf",&r.x());
      for (int i=-20; i<20; i++) {
          for (int j=-40; j<39; j++) {
              //r.x() = 0.0;
              r.y() = x0[1] + incr*i;
              r.z() = x0[2] + 0.5*incr*j;
              if (docon||dovdw) {
                  RefShape shape;
                  if (docon) shape = con;
                  else shape = vdw;
                  double d = shape->distance_to_surface(r);
                  if (d>=1.0) printf(" ");
                  else if (d<0.0) printf(".");
                  else printf("%c",'0'+((int)(d*10)));
                }
              else {
                  if (con->is_outside(r) &&
                      vdw->is_outside(r)) printf(" ");
                  else if (con->is_outside(r) &&
                           !vdw->is_outside(r)) printf(".");
                  else if (!con->is_outside(r) &&
                           vdw->is_outside(r)) printf("'");
                  else printf(":");
                }
            }
          printf("\n");
        }
      goto next_x;
    }

  if (do_opt) {
      mpqc.do_gradient(1);  
    }

  //printf("The eigenvectors:\n");
  //Print(mpqc.eigenvectors());

  printf("energy = %12.8f, exchange_energy = %12.8f\n",
         mpqc.energy(),
         mpqc.exchange_energy());

  if (do_dens) {
      cart_point r;
      r[0] = r[1] = r[2] = 0.0;
      printf("density (0.0,0.0,0.0) = %10.8f\n",mpqc.density(r));
      double nelec = 0.0;
      const double range = 4.0;
      double incr = pkeyval.doublevalue("incr");
      if (pkeyval.error() != KeyVal::OK) incr = 0.2;
      printf("incr = %8.5f\n",incr);
      double volume = incr*incr*incr;
      for (r[0] = -range; r[0]<range; r[0] += incr) {
          for (r[1] = -range; r[1]<range; r[1] += incr) {
              for (r[2] = -range; r[2]<range; r[2] += incr) {
                  nelec += mpqc.density(r) * volume;
                }
            }
        }
      printf("nelec = %12.8f volume element = %12.8f\n",nelec,volume);
    }

  if (do_dens_grad) {
      cart_point r;
      double grad[3];
      r[0] = r[1] = r[2] = 0.0;
      printf("density (0.0,0.0,0.0) = %10.8f\n",mpqc.density(r));
      double nelec = 0.0;
      const double range = 4.0;
      double incr = pkeyval.doublevalue("incr");
      if (pkeyval.error() != KeyVal::OK) incr = 0.2;
      double disp = pkeyval.doublevalue("disp");
      if (pkeyval.error() != KeyVal::OK) disp = 0.01;
      printf("incr = %8.5f\n",incr);
      double volume = incr*incr*incr;
      for (r[0] = -range; r[0]<range; r[0] += incr) {
          for (r[1] = -range; r[1]<range; r[1] += incr) {
              for (r[2] = -range; r[2]<range; r[2] += incr) {
                  double dens = mpqc.density_gradient(r,grad);
                  nelec += dens * volume;
                  double testgrad[3];
                  for (int l=0; l<3; l++) {
                      cart_point rplus = r;
                      rplus[l] += disp;
                      double dplus = mpqc.density(rplus);
                      cart_point rminus = r;
                      rminus[l] -= disp;
                      double dminus = mpqc.density(rminus);
                      testgrad[l] = (dplus - dminus)/(2.0*disp);
                    }
                  printf(" g=(% 8.4f % 8.4f % 8.4f)\n",
                         grad[0],grad[1],grad[2]);
                  printf("tg=(% 8.4f % 8.4f % 8.4f)\n",
                         testgrad[0],testgrad[1],testgrad[2]);
                }
            }
        }
      printf("nelec = %12.8f volume element = %12.8f\n",nelec,volume);
    }

  if (do_opt) {
      // set up tolerances
      TOLS tol;
      tol.SetDefaultTol();
      tol.SetAndTol(1);
      tol.PrintTol();

      // do the newton steps using a generalized inverse
      GeneralizedNewtonStep step(3 * mol.natom() - 6);
      // use bfgs to update the hessian
      BFGSupdate bfgs;
      // optimize with a quasinewton method
      OptQNewton opt(&mpqc,&tol,&bfgs,&step);

      // optimize with a conjugate gradient method
      //OptCG opt(&mpqc,&tol);

      opt.optimize();
      opt.PrintStatus();

      // print results
      fprintf(outfp,"Optimized molecule:\n");
      mol.print();
    }
}
