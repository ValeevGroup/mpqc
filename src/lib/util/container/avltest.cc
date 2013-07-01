//
// avltest.h --- test program for avl maps and sets
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

#include <iostream>
#include <math.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif
#include <util/container/eavlmmap.h>
#include <util/container/avlmap.h>
#include <util/container/avlset.h>

using namespace std;
using sc::AVLMap;
using sc::AVLSet;
using sc::EAVLMMap;
using sc::EAVLMMapNode;
using sc::compare;

class Data {
  public:
    EAVLMMapNode<int,Data> map1;
    EAVLMMapNode<int,Data> map2;
  public:
    Data(int k1, int k2 = 0): map1(k1), map2(k2) {};
    void print(int indent = 0);
    void change1(int val) { map1.key = val; }
};

void
Data::print(int indent)
{
  for (int i=0; i<indent; i++) cout << " ";
  cout << map1.key;
}

#define TEST1 1
#define TEST2 1
#define TEST3 1
#define TEST4 1
#define TEST5 1
#define TEST6 1
#define TEST7 1
#define TEST8 0
#define TEST9 0

static int Ni = 0;
static int Nr = 0;
static int Nf = 0;

void
testmap(EAVLMMap<int, Data>& map, Data** data, int n)
{
  for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
          map.insert(data[j]);
          Ni++;
        }
#if 0
      if (i==0) {
          cout << "--------------------------------------------" << endl;
          map.check();
          map.print2();
        }
      cout << "............................. removing ";
      data[i]->print();
      cout << endl;
#endif
      map.remove(data[i]);
      Nr++;
#if 0
      map.print2();
#endif
      map.check();
      map.clear_without_delete();
    }
}

void
rantest(EAVLMMap<int, Data>&map1, Data** data, int n)
{
  int i;
  for (i=0; i<n; i++) {
      Data* d = data[i];
      d->change1(random());
      map1.insert(d);
      Ni++;
    }
  map1.check();
  for (i=0; i<n; i++) {
      map1.find(i);
      Nf++;
    }
  for (i=0; i<n; i++) {
      Data* d = data[i];
      map1.remove(d);
      Nr++;
    }
  map1.check();
  map1.clear_without_delete();
}

int
main(int argc, char* argv[])
{
  int i;
  const int maxkey = 9;
  EAVLMMap<int,Data> map1(&Data::map1);
  Data* data[maxkey][maxkey];
  Data* currentdata[maxkey];
  for (i=0; i<maxkey; i++) {
      for (int j=0; j<maxkey; j++) {
          data[i][j] = new Data(j);
        }
    }
  int max;

  const int unique = 1;

  cout << "emmap:" << endl;
  EAVLMMap<int,Data> emap(&Data::map1);
  for (i=0; i<maxkey; i++) {
      emap.insert(new Data(i/2,i));
      emap.check();
    }
  for (EAVLMMap<int,Data>::iterator im=emap.begin(); im!=emap.end(); im++) {
      cout << " " << im.key() << " " << im->map2.key << endl;
    }

  cout << "map:" << endl;
  AVLMap<int, char> icmap;
  for (i=0; i<maxkey; i++) {
      int d = random();
      icmap.insert(d, i+'a');
      icmap.insert(d, i+'a');
      icmap.check();
    }
  icmap.print();
  for (AVLMap<int,char>::iterator ic=icmap.begin(); ic!=icmap.end(); ic++) {
      cout << " " << ic.key() << " " << ic.data() << endl;
    }

  cout << "set:" << endl;
  AVLSet<int> iset;
  for (i=0; i<10; i++) {
      iset.insert(i);
    }
  if (iset.length() != 10) abort();
  iset.remove(100);
  if (iset.length() != 10) abort();
  if (!iset.contains(3)) abort();
  iset.remove(3);
  if (iset.contains(3)) abort();
  if (iset.length() != 9) abort();
  for (i=0; i<10; i++) iset.remove(i);
  if (iset.length() != 0) abort();
  if (iset.contains(3)) abort();
  if (iset.length() != 0) abort();
  iset.clear();
  for (i=0; i<maxkey; i++) {
      int d = random();
      iset.insert(d);
      iset.insert(d);
      iset.check();
    }
  iset.print();
  for (AVLSet<int>::iterator is=iset.begin(); is!=iset.end(); is++) {
      cout << " " << is.key() << endl;
    }

#if TEST1
  cout << "=================================================" << endl;
  max = 1;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];

      testmap(map1, currentdata, max);
    }
#endif

#if TEST2
  cout << "=================================================" << endl;
  max = 2;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          testmap(map1, currentdata, max);
        }
    }
#endif

#if TEST3
  cout << "=================================================" << endl;
  max = 3;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];

              testmap(map1, currentdata, max);
            }
        }
    }
#endif

#if TEST4
  cout << "=================================================" << endl;
  max = 4;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
                  currentdata[3] = data[0][l];

                  testmap(map1, currentdata, max);
                }
            }
        }
    }
#endif

#if TEST5
  cout << "=================================================" << endl;
  max = 5;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
                  currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
                  currentdata[4] = data[0][m];

                  testmap(map1, currentdata, max);
                  }
                }
            }
        }
    }
#endif

#if TEST6
  cout << "=================================================" << endl;
  max = 6;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
          currentdata[1] = data[0][j];
          cout << "6: i = " << i << " j = " << j << endl;
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
      currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];

                  testmap(map1, currentdata, max);

                }
                }
                }
            }
        }
    }
#endif

#if TEST7
  cout << "=================================================" << endl;
  max = 7;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
      currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
              currentdata[2] = data[0][k];
              cout << "7: i = " << i << " j = " << j << " k = " << k << endl;
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];
                  for (int o=0; o<max; o++) {
                  if (unique && o==i || o==j || o==k || o==l || o==m || o==n) continue;
      currentdata[6] = data[0][o];

                  testmap(map1, currentdata, max);

                }
                }
                }
                }
            }
        }
    }
#endif

#if TEST8
  cout << "=================================================" << endl;
  max = 8;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
      currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
      currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  cout << "7: i = " << i << " j = " << j << " k = " << k
                       << " l = " << l << endl;
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];
                  for (int o=0; o<max; o++) {
                  if (unique && o==i || o==j || o==k || o==l || o==m || o==n) continue;
      currentdata[6] = data[0][o];
                  for (int p=0; p<max; p++) {
                  if (unique && p==i||p==j||p==k||p==l||p==m||p==n||p==o) continue;
      currentdata[7] = data[0][p];

                  testmap(map1, currentdata, max);

                }
                }
                }
                }
                }
            }
        }
    }
#endif

#if TEST9
  cout << "=================================================" << endl;
  max = 9;
  for (i=0; i<max; i++) {
      currentdata[0] = data[0][i];
      for (int j=0; j<max; j++) {
          if (unique && i == j) continue;
      currentdata[1] = data[0][j];
          for (int k=0; k<max; k++) {
              if (unique && k==i || k==j) continue;
      currentdata[2] = data[0][k];
              for (int l=0; l<max; l++) {
                  if (unique && l==i || l==j || l==k) continue;
      currentdata[3] = data[0][l];
                  for (int m=0; m<max; m++) {
                  if (unique && m==i || m==j || m==k || m==l) continue;
      currentdata[4] = data[0][m];
                  cout << "7: i = " << i << " j = " << j << " k = " << k
                       << " l = " << l << " m = " << m << endl;
                  for (int n=0; n<max; n++) {
                  if (unique && n==i || n==j || n==k || n==l || n==m) continue;
      currentdata[5] = data[0][n];
                  for (int o=0; o<max; o++) {
                  if (unique && o==i || o==j || o==k || o==l || o==m || o==n) continue;
      currentdata[6] = data[0][o];
                  for (int p=0; p<max; p++) {
                  if (unique && p==i||p==j||p==k||p==l||p==m||p==n||p==o) continue;
      currentdata[7] = data[0][p];
                  for (int q=0; q<max; q++) {
                  if (unique && q==i||q==j||q==k||q==l||q==m||q==n||q==o||q==p) continue;
      currentdata[8] = data[0][q];

                  testmap(map1, currentdata, max);

                }
                }
                }
                }
                }
                }
            }
        }
    }
#endif

  cout << "Ni = " << Ni << ", Nr = " << Nr << ", N = " << Ni+Nr << endl;

  const int maxdat2 = 2000;
  Data * data2[maxdat2];
  for (i=0; i<maxdat2; i++) {
      data2[i] = new Data(i);
    }
  for (i=0; i<maxdat2; i++) {
      if (i%100 == 0) cout << "-";
    }
  cout << endl;
  for (i=0; i<maxdat2; i++) {
      if (i%100 == 0) {
          cout << ".";
        }
      rantest(map1, data2, i);
    }
  cout << endl;

  cout << "Ni = " << Ni << ", Nr = " << Nr << ", Nf = " << Nf
       << ", N = " << Ni+Nr << endl;

  return 0;
}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template class AVLMapNode<int, char>;
template class EAVLMMap<int,AVLMapNode<int, char> >;
template class AVLMapNode<int, int>;
template class EAVLMMap<int,AVLMapNode<int, int> >;
template class EAVLMMap<int, Data >;
#endif

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
