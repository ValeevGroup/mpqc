
#include <math.h>
#include <float.h>

#include <chemistry/qc/basis/fjt.h>
#include <util/misc/formio.h>

using namespace std;
using namespace sc;

// this should only be used for x < a + 1
static double
mgamma_small_x(double a, double x)
{
  const int nterm = 1000;
  // n = 0
  double xn = 1.0;
  double an = a;
  double c0 = xn/an;
  int n = 1;
  double sum = c0;
  double contrib;
  double c0eps = c0 * DBL_EPSILON;
  do {
      an *= a+n;
      xn *= x;
      contrib = xn/an;
      n++;
      sum += contrib;
    } while(contrib > c0eps);
  return sum * exp(-x);
}

static double
Gamma_m(int m, double x)
{
  double a = m+0.5;
  const double tiny = DBL_EPSILON * DBL_EPSILON;
  // iteration 0
  double f = tiny;
  double c = f;
  double d = 0.0;
  double delta;
  // iteration 1 +
  int j=1;
  double ac = 1.0;
  double bc = x + 1.0 - a;
  do {
      d = bc + ac * d;
      if (d == 0.0) d = tiny;
      c = bc + ac/c;
      if (c == 0.0) c = tiny;
      d = 1.0/d;
      delta = c*d;
      f = f*delta;
      ac = - j * (j - a);
      bc += 2.0;
      j++;
    } while (fabs(delta-1.0) > DBL_EPSILON);
  return exp(-x)*pow(x,a)*f;
}

static double
Gamma_m(int m)
{
  double a = m+0.5;
  const double c[6] = {  76.18009172947146,
                        -86.50532032941677,
                         24.01409824083091,
                         -1.231739572450155,
                          0.1208650973866179e-2,
                         -0.5395239384953e-5};
  double x,y,tmp,ser;
  int j;
  y = x = a;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) ser += c[j]/++y;
  return exp(-tmp+log(2.5066282746310005*ser/x));
}

static double
gamma_m_large_x(int m, double x)
{
  return Gamma_m(m) - Gamma_m(m,x);
}

static double
mgamma_m(int m, double x)
{
  double a = m+0.5;
  if (x < a + 1.0) return mgamma_small_x(a,x);
  return gamma_m_large_x(m,x)/pow(x,a);
}

static double
_Fjt(int j, double t)
{
  //return gamma(j+0.5,sqrt(t))/(2.0*pow(t,j+0.5));
  return 0.5 * mgamma_m(j,t);
}

int
main(int,char**)
{
  int maxj = 18;
  Ref<FJT> fjt = new FJT(maxj+1);

  double tinc = 0.1;
  for (double T=0.0; T<1000.0; T+=tinc) {
      double *values = fjt->values(maxj,T);
      for (int j=0; j<=maxj; j+=1) {
          double v1 = values[j];
          double v2 = _Fjt(j,T);
          double error = fabs((v1-v2)/v1);
          if (error > DBL_EPSILON*10.0) {
              cout << scprintf("F(%2d,%5.2f) = %15.12f %15.12f e = %18.15f",
                               j,T,v1,v2, error) << endl;
            }
          else {
              cout << scprintf("F(%2d,%5.2f) = %15.12f %15.12f",
                               j,T,v1,v2) << endl;
            }
        }
      if (T > 10.0) tinc = 10.0;
      if (T > 100.0) tinc = 100.0;
    }

  return 0;
}
