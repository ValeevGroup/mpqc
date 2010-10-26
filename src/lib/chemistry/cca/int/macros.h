//
// macros.h
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

/* True if the integral is nonzero. */
#define INT_NONZERO(x) (((x)< -1.0e-15)||((x)> 1.0e-15))

/* Computes an index to a Cartesian function within a shell given
 * am = total angular momentum
 * i = the exponent of x (i is used twice in the macro--beware side effects)
 * j = the exponent of y
 * formula: (am - i + 1)*(am - i)/2 + am - i - j unless i==am, then 0
 * The following loop will generate indices in the proper order:
 *  cartindex = 0;
 *  for (i=am; i>=0; i--) {
 *    for (j=am-i; j>=0; j--) {
 *      do_it_with(cartindex);
 *      cartindex++;
 *      }
 *    }
 */
#define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))

/* This sets up the above loop over cartesian exponents as follows
 * FOR_CART(i,j,k,am)
 *   Stuff using i,j,k.
 *   END_FOR_CART
 */
#define FOR_CART(i,j,k,am) for((i)=(am);(i)>=0;(i)--) {\
                           for((j)=(am)-(i);(j)>=0;(j)--) \
                           { (k) = (am) - (i) - (j);
#define END_FOR_CART }}

/* This sets up a loop over all of the generalized contractions
 * and all of the cartesian exponents.
 * gc is the number of the gen con
 * index is the index within the current gen con.
 * i,j,k are the angular momentum for x,y,z
 * sh is the shell pointer
 */
#define FOR_GCCART(gc,index,i,j,k,sh)\
    for ((gc)=0; (gc)<(sh)->ncon; (gc)++) {\
    (index)=0;\
    FOR_CART(i,j,k,(sh)->type[gc].am)

#define FOR_GCCART_GS(gc,index,i,j,k,sh)\
    for ((gc)=0; (gc)<(sh)->ncontraction(); (gc)++) {\
    (index)=0;\
    FOR_CART(i,j,k,(sh)->am(gc))

#define END_FOR_GCCART(index)\
    (index)++;\
    END_FOR_CART\
    }

#define END_FOR_GCCART_GS(index)\
    (index)++;\
    END_FOR_CART\
    }

/* These are like the above except no index is kept track of. */
#define FOR_GCCART2(gc,i,j,k,sh)\
    for ((gc)=0; (gc)<(sh)->ncon; (gc)++) {\
    FOR_CART(i,j,k,(sh)->type[gc].am)

#define END_FOR_GCCART2\
    END_FOR_CART\
    }

/* These are used to loop over shells, given the centers structure
 * and the center index, and shell index. */
#define FOR_SHELLS(c,i,j) for((i)=0;(i)<(c)->n;i++) {\
                          for((j)=0;(j)<(c)->center[(i)].basis.n;j++) {
#define END_FOR_SHELLS }}

/* Computes the number of Cartesian function in a shell given
 * am = total angular momentum
 * formula: (am*(am+1))/2 + am+1;
 */
#define INT_NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

/* Like INT_NCART, but only for nonnegative arguments. */
#define INT_NCART_NN(am) ((((am)+2)*((am)+1))>>1)

/* For a given ang. mom., am, with n cartesian functions, compute the
 * number of cartesian functions for am+1 or am-1
 */
#define INT_NCART_DEC(am,n) ((n)-(am)-1)
#define INT_NCART_INC(am,n) ((n)+(am)+2)

/* Computes the number of pure angular momentum functions in a shell
 * given am = total angular momentum
 */
#define INT_NPURE(am) (2*(am)+1)

/* Computes the number of functions in a shell given
 * pu = pure angular momentum boolean
 * am = total angular momentum
 */
#define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))

/* Given a centers pointer and a shell number, this evaluates the
 * pointer to that shell. */
#define INT_SH(c,s) ((c)->center[(c)->center_num[s]].basis.shell[(c)->shell_num[s]])

/* Given a centers pointer and a shell number, get the angular momentum
 * of that shell. */
#define INT_SH_AM(c,s) ((c)->center[(c)->center_num[s]].basis.shell[(c)->shell_num[s]].type.am)

/* Given a centers pointer and a shell number, get pure angular momentum
 * boolean for that shell. */
#define INT_SH_PU(c,s) ((c)->center[(c)->center_num[s]].basis.shell[(c)->shell_num[s]].type.puream)

/* Given a centers pointer, a center number, and a shell number,
 * get the angular momentum of that shell. */
#define INT_CE_SH_AM(c,a,s) ((c)->center[(a)].basis.shell[(s)].type.am)

/* Given a centers pointer, a center number, and a shell number,
 * get pure angular momentum boolean for that shell. */
#define INT_CE_SH_PU(c,a,s) ((c)->center[(a)].basis.shell[(s)].type.puream)

/* Given a centers pointer and a shell number, compute the number
 * of functions in that shell. */
/* #define INT_SH_NFUNC(c,s) INT_NFUNC(INT_SH_PU(c,s),INT_SH_AM(c,s)) */
#define INT_SH_NFUNC(c,s) ((c)->center[(c)->center_num[s]].basis.shell[(c)->shell_num[s]].nfunc)

/* These macros assist in looping over the unique integrals
 * in a shell quartet.  The exy variables are booleans giving
 * information about the equivalence between shells x and y.  The nx
 * variables give the number of functions in each shell, x. The
 * i,j,k are the current values of the looping indices for shells 1, 2, and 3.
 * The macros return the maximum index to be included in a summation
 * over indices 1, 2, 3, and 4.
 * These macros require canonical integrals.  This requirement comes
 * from the need that integrals of the shells (1 2|2 1) are not
 * used.  The integrals (1 2|1 2) must be used with these macros to
 * get the right nonredundant integrals.
 */
#define INT_MAX1(n1) ((n1)-1)
#define INT_MAX2(e12,i,n2) ((e12)?(i):((n2)-1))
#define INT_MAX3(e13e24,i,n3) ((e13e24)?(i):((n3)-1))
#define INT_MAX4(e13e24,e34,i,j,k,n4) \
  ((e34)?(((e13e24)&&((k)==(i)))?(j):(k)) \
        :((e13e24)&&((k)==(i)))?(j):(n4)-1)
/* A note on integral symmetries:
 *  There are 15 ways of having equivalent indices.
 *  There are 8 of these which are important for determining the
 *  nonredundant integrals (that is there are only 8 ways of counting
 *  the number of nonredundant integrals in a shell quartet)
 * Integral type   Integral    Counting Type
 *     1           (1 2|3 4)      1
 *     2           (1 1|3 4)      2
 *     3           (1 2|1 4)       ->1
 *     4           (1 2|3 1)       ->1
 *     5           (1 1|1 4)      3
 *     6           (1 1|3 1)       ->2
 *     7           (1 2|1 1)       ->5
 *     8           (1 1|1 1)      4
 *     9           (1 2|2 4)       ->1
 *    10           (1 2|3 2)       ->1
 *    11           (1 2|3 3)      5
 *    12           (1 1|3 3)      6
 *    13           (1 2|1 2)      7
 *    14           (1 2|2 1)      8    reduces to 7 thru canonicalization
 *    15           (1 2|2 2)       ->5
 */
