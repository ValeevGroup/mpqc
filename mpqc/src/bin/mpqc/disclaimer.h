

static char *disclaimer[] = {
"   !!!!!!! WARNING !!!!!!!!!! WARNING !!!!!!!!!!! WARNING !!!!!!!\n\n",
"   This is a direct scf program which utilizes quantum chemistry\n",
"   libraries developed at the Univerity of Georgia, Sandia National\n",
"   Laboratories, and the National Institutes of Health.  It is free.\n",
"   You need pay a license to no one to use this program.  You may\n",
"   modify it at will.  We don\'t care.  Just understand one thing...\n",
"   This program comes with NO WARRANTY.  The answers it produces are\n",
"   NOT guaranteed to be better than random numbers.  If you lose a\n",
"   million because MPQC messes up, it\'s you that\'s out the million,\n",
"   not us.  If you don\'t like this disclaimer:  tough.  We reserve\n",
"   the right to do the absolute minimum provided by law, up to and\n",
"   including nothing.  Also, if you modify this code and publish\n",
"   results obtained with it, please mention in your paper that you\n",
"   modified the code so we won\'t be blamed if your results are garbage.\n\n",
"   !!!!!!! WARNING !!!!!!!!!! WARNING !!!!!!!!!!! WARNING !!!!!!!\n\n"};

print_disclaimer(FILE *out)
{
  int i;

  for(i=0; i < 15 ; i++)
    fprintf(out,"%s",disclaimer[i]);
  }
