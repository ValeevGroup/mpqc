
#include <stdio.h>
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>

int
main(int argc, char** argv)
{
  if (argc != 3) {
      fprintf(stderr,"Usage: %s inputfile keyword\n", argv[0]);
      exit(1);
    }
  char* inputfile = argv[1];
  char* keyword = argv[2];

  FILE*fp = fopen(inputfile,"r");
  if (!fp) {
      fprintf(stderr,"%s: error opening input file \"%s\"\n",
              argv[0], inputfile);
      perror("fopen");
      exit(1);
    }

  RefParsedKeyVal keyval = new ParsedKeyVal(fp);

  if (!keyval->exists(keyword)) {
      fprintf(stderr,"%s: keyword \"%s\" doesn't exist in file \"%s\"\n",
              argv[0], keyword, inputfile);
      exit(1);
    }

  if (!keyval->classname(keyword)) {
      fprintf(stderr,"%s: keyword \"%s\" in file \"%s\" is not an object\n",
              argv[0], keyword, inputfile);
      exit(1);
    }

  RefDescribedClass dc = keyval->describedclassvalue(keyword);

  if (dc.null()) {
      fprintf(stderr,"%s: keyword \"%s\" in file \"%s\" could not be"
              " converted into a DescribedClass object\n",
              argv[0], keyword, inputfile);
      exit(1);
    }

  RefOptimize opt = dc;

  if (opt.null()) {
      fprintf(stderr,"%s: keyword \"%s\" in file \"%s\" could not be"
              " converted into an Optimize object\n",
              argv[0], keyword, inputfile);
      exit(1);
    }

  opt->optimize();

  return 0;
}

