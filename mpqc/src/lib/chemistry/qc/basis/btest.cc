
#include <iostream.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/integral/integralv2.h>

#include <util/misc/timer.h>
#include <util/group/picl.h>

int
main(int, char *argv[])
{
  char *filename = (argv[1]) ? argv[1] : SRCDIR "/btest.kv";
  
  RefKeyVal keyval = new ParsedKeyVal(filename);

  int nproc, me, host, top, ord, dir;
  RefMessageGrp grp = MessageGrp::get_default_messagegrp();
  open0_messagegrp(&nproc, &me, &host, grp);
  setarc0(&nproc,&top,&ord,&dir);

  for (int i=0; i<keyval->count("test"); i++) {
      RefGaussianBasisSet gbs = keyval->describedclassvalue("test", i);

      if (SymmGaussianBasisSet::castdown(gbs)) {
        SymmGaussianBasisSet *sgbs = SymmGaussianBasisSet::castdown(gbs);
        
        CharacterTable ct = gbs->molecule()->point_group().char_table();
        PetiteList& pl = sgbs->petite_list();

        pl.print(stdout,0);

        tim_enter("aotoso");
        RefSCMatrix aotoso(pl.AO_basisdim(), pl.SO_basisdim());
        RefSCElementOp sotrans = new AOSO_Transformation(gbs);
        aotoso.element_op(sotrans);
        sotrans=0;
        tim_exit("aotoso");

        tim_enter("overlap");
        RefSymmSCMatrix s(pl.AO_basisdim());
        s.assign(0.0);
        RefSCElementOp op =
          new GaussianOverlapIntv2(gbs, new SymmOneBodyIntIter(gbs));
        s.element_op(op);
        op=0;
        tim_exit("overlap");

        //s.print("overlap");
        tim_print(0);

        tim_enter("symmetrize");
        RefSymmSCMatrix ss(pl.SO_basisdim());
        pl.symmetrize(s,ss);
        tim_exit("symmetrize");
        //ss.print("blocked ss");
        tim_print(0);

        RefSCMatrix sotoao = aotoso.i();
        RefSymmSCMatrix us = s.clone();
        us.assign(0.0);
        us.accumulate_transform(sotoao.t(),ss);
        //us.print("us");

        tim_enter("full overlap");
        op = new GaussianOverlapIntv2(gbs);
        s.assign(0.0);
        s.element_op(op);
        op=0;
        tim_exit("full overlap");
        //s.print("real overlap");
        tim_print(0);

        RefSymmSCMatrix uss = ss.clone();
        uss.assign(0.0);
        tim_enter("transform");
        uss.accumulate_transform(aotoso.t(),s);
        tim_exit("transform");
        //uss.print("uss");
        tim_print(0);
        
        //uss = uss-ss;
        //uss.print("uss - ss");
      
        //us = us-s;
        //us.print("us - s");
      } else {
        SymmGaussianBasisSet *sgbs = new SymmGaussianBasisSet(*gbs.pointer());

        CharacterTable ct = gbs->molecule()->point_group().char_table();
        PetiteList& pl = sgbs->petite_list();

        pl.print(stdout,0);
      }
      
      //gbs->print();

      fflush(stdout);
      cout.flush();

      StateOutText out("btest.out");
      gbs.save_state(out);
      StateInText in("btest.out");
      gbs.restore_state(in);
      //gbs->print();
    }

  return 0;
}
