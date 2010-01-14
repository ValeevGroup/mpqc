
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

#include <map>
#include <iostream>
#include <iomanip>

#include <chemistry/qc/lmp2/dgemminfo.h>

namespace sc {

class dgemminfo: public std::pair<std::pair<int,int>,int> {
  public:
    dgemminfo(int m, int l, int n) {
      first.first = m;
      first.second = l;
      second = n;
    }
    int m() const { return first.first; }
    int l() const { return first.second; }
    int n() const { return second; }
};

std::ostream &
operator << (std::ostream &o, const dgemminfo &d) {
  o << d.m();
  o << "," << d.l();
  o << "," << d.n();
  return o;
}

class count {
    unsigned int i;
  public:
    count(): i(0) {}
    unsigned int n() const { return i; }
    unsigned int operator++() {
      i = i+1;
      return i;
    }
    unsigned int operator++(int) {
      unsigned int tmp = i;
      i = i+1;
      return tmp;
    }
};

std::ostream &
operator << (std::ostream &o, const count &c) {
  o << c.n();
  return o;
}

void print_dgemm();

class Double {
    double d;
  public:
    Double(): d(0.0) {}
    double val() const { return d; }
    void operator += (double v) { d += v; }
};


class dgemmmap_t: public std::map<dgemminfo,count> {
  public:
    ~dgemmmap_t() { print_dgemm(); }
    std::map<dgemminfo,Double> dgemmtime;
    void add_time(const dgemminfo &di, double t) {
      dgemmtime[di] += t;
    }
};

static dgemmmap_t dgemmmap;

void
count_dgemm(int m, int l, int n, double t)
{
  dgemminfo di(m<n?m:n,l,m<n?n:m);
  dgemmmap[di]++;
  dgemmmap.add_time(di,t);
}

void
print_dgemm()
{
  double sum_nflop = 0.0;
  double total_seconds = 0.0;
  std::multimap<double,dgemminfo> sorted;
  for (std::map<dgemminfo,count>::iterator i = dgemmmap.begin();
       i != dgemmmap.end();
       i++) {
      const dgemminfo &info = i->first;
      int m = info.m();
      int l = info.l();
      int n = info.n();
      int nflop = m*l*n*2;
      int ncall = dgemmmap[info].n();
      double seconds = dgemmmap.dgemmtime[info].val();
      double nflop_total = (double)nflop * ncall;
//       sorted.insert(std::make_pair(nflop_total,i->first));
      sorted.insert(std::make_pair(seconds,i->first));
      total_seconds += seconds;
      sum_nflop += nflop_total;
    }

  if (total_seconds == 0.0) return;

  double sum_percent_nflop = 0.0;
  double sum_t = 0.0;
  std::cout << " " << std::setw(4) << "m"
            << " " << std::setw(4) << "l"
            << " " << std::setw(4) << "n"
            << " " << std::setw(8) << "ncall"
            << " " << std::setw(6) << "%nflop"
//             << " " << std::setw(7) << "cum%nfl"
            << " " << std::setw(7) << "t(msec)"
            << " " << std::setw(7) << "%t"
            << " " << std::setw(7) << "cum%t"
            << " " << std::setw(7) << "MFLOP/s"
            << std::endl;
  for (std::multimap<double,dgemminfo>::reverse_iterator
           i = sorted.rbegin(); i != sorted.rend(); i++) {
      const dgemminfo &info = i->second;
      int m = info.m();
      int l = info.l();
      int n = info.n();
      int ncall = dgemmmap[info].n();
      double nflop_total = 2.0*ncall*m*l*n;
      sum_percent_nflop += nflop_total/sum_nflop;
      double seconds = dgemmmap.dgemmtime[info].val();
      sum_t += seconds;
      std::cout << " " << std::setw(4) << m
                << " " << std::setw(4) << l
                << " " << std::setw(4) << n
                << " " << std::setw(8) << ncall
                << " " << std::setw(6) << nflop_total/sum_nflop*100.0
//                 << " " << std::setw(7) << sum_percent_nflop*100.0
                << " " << std::setw(7) << std::setprecision(1)
                << seconds*1e3
                << " " << std::setw(7) << seconds/total_seconds*100.0
                << " " << std::setw(7) << std::setprecision(1)
                << sum_t/total_seconds*100.0
                << " " << std::setw(7) << std::setprecision(3)
                << 1e-6*nflop_total/seconds
                << std::endl;
    }

  std::cout << "The total number of flops was "
            << sum_nflop << std::endl;
  std::cout << "The total number of seconds was "
            << std::setprecision(3)
            << total_seconds << std::endl;
  std::cout << "Computation rate was "
            << 1e-6*sum_nflop/total_seconds
            << " MFLOP/s" << std::endl;
}

}
