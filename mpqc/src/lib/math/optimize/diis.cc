
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/optimize/diis.h>

#define CLASSNAME DIIS
#define PARENTS public SelfConsistentExtrapolation
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>
void *
DIIS::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SelfConsistentExtrapolation::_castdown(cd);
  return do_castdowns(casts,cd);
}

void
DIIS::init()
{
  int i;
  int dim = ndiis+1;

  iter = 0;

  if ((allocbn_double_matrix(&bmat,"n1 n2",dim,dim) != 0) ||
      (allocbn_double_matrix(&bold,"n1 n2",ndiis,ndiis) != 0) ||
      (allocbn_double_vector(&btemp,"n",dim) != 0)) {
    fprintf(stderr,"DIIS::init:  alloc of bmat, bold, and btemp failed\n");
    abort();
  }

  diism_data = new RefSCExtrapData[ndiis];
  diism_error = new RefSCExtrapError[ndiis];
}

DIIS::DIIS()
{
  ndiis = 5;
  start = 1;

  init();
}

DIIS::DIIS(StateIn& s) :
  SelfConsistentExtrapolation(s)
{
  s.get(start);
  s.get(ndiis);
  s.get(iter);
  s.get(damping_factor);

  // hack for old sgen stuff...this will go away soon
  s.get(btemp.n);
  allocbn_double_vector(&btemp,"n",btemp.n);
  for (int i=0; i < btemp.n; i++)
    s.get(btemp.d[i]);

  s.get(bold.n1);
  s.get(bold.n2);
  allocbn_double_matrix(&bold,"n1 n2", bold.n1, bold.n2);
  for (int i=0; i < bold.n1; i++)
    for (int j=0; j < bold.n2; j++)
      s.get(bold.d[i][j]);
  
  s.get(bmat.n1);
  s.get(bmat.n2);
  allocbn_double_matrix(&bmat,"n1 n2", bmat.n1, bmat.n2);
  for (int i=0; i < bmat.n1; i++)
    for (int j=0; j < bmat.n2; j++)
      s.get(bmat.d[i][j]);
  
  diism_data = new RefSCExtrapData[ndiis];
  diism_error = new RefSCExtrapError[ndiis];
  for (int i=0; i < ndiis; i++) {
    diism_data[i].restore_state(s);
    diism_error[i].restore_state(s);
  }
}

DIIS::DIIS(const RefKeyVal& keyval):
  SelfConsistentExtrapolation(keyval)
{
  ndiis = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) ndiis = 5;

  start = keyval->intvalue("start");
  if (keyval->error() != KeyVal::OK) start = 1;

  if (ndiis <= 0) {
      fprintf(stderr, "DIIS::DIIS(const RefKeyVal& keyval): got ndiis = 0\n");
      abort();
    }

  init();
}

DIIS::~DIIS()
{
  free_double_matrix(&bmat);
  free_double_matrix(&bold);
  free_double_vector(&btemp);

  delete[] diism_data;
  delete[] diism_error;
}

void
DIIS::save_data_state(StateOut& s)
{
  SelfConsistentExtrapolation::save_data_state(s);
  s.put(start);
  s.put(ndiis);
  s.put(iter);
  s.put(damping_factor);

  // hack for old sgen stuff...this will go away soon
  s.put(btemp.n);
  for (int i=0; i < btemp.n; i++)
    s.put(btemp.d[i]);

  s.put(bold.n1);
  s.put(bold.n2);
  for (int i=0; i < bold.n1; i++)
    for (int j=0; j < bold.n2; j++)
      s.put(bold.d[i][j]);
  
  s.put(bmat.n1);
  s.put(bmat.n2);
  for (int i=0; i < bmat.n1; i++)
    for (int j=0; j < bmat.n2; j++)
      s.put(bmat.d[i][j]);
  
  for (int i=0; i < ndiis; i++) {
    diism_data[i].save_state(s);
    diism_error[i].save_state(s);
  }
}

int
DIIS::extrapolate(const RefSCExtrapData& data,
                  const RefSCExtrapError& error)
{
  int i, j, k, kk;
  int last = iter;
  int trial = 0;
  int col = iter + 2;
  double norm, determ;
  double scale;

  iter++;

  scale = 1.0 + damping_factor;

  if (iter > ndiis) {
      last = ndiis-1;
      col = ndiis+1;
      dtemp_data = diism_data[0];
      dtemp_error = diism_error[0];
      for (i=0; i < last ; i++) {
          diism_data[i] = diism_data[i+1];
          diism_error[i] = diism_error[i+1];
        }
      diism_data[last] = dtemp_data;
      diism_error[last] = dtemp_error;
    }

  diism_data[last] = data->copy();
  diism_error[last] = error;

  set_error(error->error());
               
  // then set up B matrix, where B(i,j) = <ei|ej>

  // move bold(i+1,j+1) to bold(i,j)
  if (iter > ndiis) {
      for (i=0; i < last ; i++) {
          for (j=0; j <= i ; j++) {
              bold.d[i][j]=bold.d[j][i]=bold.d[i+1][j+1];
            }
        }
    }

  // and set the current rows of bold
  for (i=0; i <= last ; i++)
      bold.d[i][last]=bold.d[last][i] = 
                      diism_error[i]->scalar_product(diism_error[last]);

  bmat.d[0][0] = 0.0;
  btemp.d[0] = -1.0;

  if (bold.d[0][0] > 1.e-10) {
      norm = 1.0/bold.d[0][0];
    }
  else {
      norm = 1.0;
    }

  for (i=1; i <= last+1 ; i++) {
      bmat.d[i][0]=bmat.d[0][i] = -1.0;
      btemp.d[i] = 0.0;
      for (j=1; j <= i ; j++) {
          bmat.d[i][j]=bmat.d[j][i] = bold.d[i-1][j-1]*norm;
          if (i==j) bmat.d[i][j] *= scale;
        }
    }

  // finally, solve the set of linear equations, obtain the coefficients,
  // and form the new fock matrix F= sum(i=1,n) ci*Fi

  if (iter-1) {
      determ = math_lin(&bmat,&btemp,col,1);

      // test for poorly conditioned equations */
      while (fabs(determ) < 1.0e-19 && trial < last) {

          trial++;
          col--;

          bmat.d[0][0] = 0.0;
          btemp.d[0] = -1.0;

          if (bold.d[trial][trial] > 1.e-10) {
              norm=1.0/bold.d[trial][trial];
            }
          else {
              norm = 1.0;
            }

          for (i=1; i <= ndiis-trial ; i++) {
              bmat.d[i][0]=bmat.d[0][i] = -1.0;
              for (j=1; j <= i ; j++) {
                  bmat.d[i][j]=bmat.d[j][i]=bold.d[i+trial-1][j+trial-1]*norm;
                  if (i==j) bmat.d[i][j] *= scale;
                }
              btemp.d[i] = 0.0;
            }

          determ = math_lin(&bmat,&btemp,col,1);
        }

      if (fabs(determ) < 10.0e-20) {
          fprintf(stderr,"DIIS::extrapolate:  trial %d no good\n",trial);
          return -1;
        }

      if (iter >= start) {
          int kk=1;

          data->zero();

          for (k=trial; k < last+1 ; k++) {
              data->accumulate_scaled(btemp.d[kk], diism_data[k]);
              kk++;
            }
        }
    }

  return 0;
}
