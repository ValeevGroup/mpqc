
#ifndef _mpqc_disclaimer_h
#define _mpqc_disclaimer_h

#include <util/misc/formio.h>

static char *disclaimer[] = {
  "!!!!!!! WARNING !!!!!!!!!! WARNING !!!!!!!!!!! WARNING !!!!!!!\n",
  "This is a quantum chemistry program which utilizes libraries",
  "developed at the Univerity of Georgia, Sandia National",
  "Laboratories, Stanford University, and the National Institutes",
  "of Health.  It is free.  You need pay a license fee to no one",
  "to use this program.  You may modify it at will.  We don\'t",
  "care.  Just understand one thing...\n",
  "This program comes with NO WARRANTY.  The answers it produces are",
  "NOT guaranteed to be better than random numbers.  If you lose a",
  "million because MPQC messes up, it\'s you that\'s out the million,",
  "not us.  If you don\'t like this disclaimer:  tough.  We reserve",
  "the right to do the absolute minimum provided by law, up to and",
  "including nothing.  Also, if you modify this code and publish",
  "results obtained with it, please mention in your paper that you",
  "modified the code so we won\'t be blamed if your results are garbage.\n",
  "!!!!!!! WARNING !!!!!!!!!! WARNING !!!!!!!!!!! WARNING !!!!!!!\n",
  0
};

static void
print_disclaimer(ostream& out)
{
  out << node0 << endl;
  
  for (int i=0; disclaimer[i]; i++)
    out << node0 << indent << disclaimer[i] << endl;
}

#endif
