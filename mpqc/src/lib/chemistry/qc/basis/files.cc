
#include <ctype.h>
#include <string.h>
#include <chemistry/qc/basis/files.h>
#include <stdio.h>

BasisFileSet::BasisFileSet(const RefKeyVal& keyval)
{
  nbasissets_ = 0;
  basissets_ = 0;

  dir_[0] = keyval->pcharvalue("basisdir");
  dir_[1] = getenv("SCLIBDIR");
  if (dir_[1]) {
      char *tmp = strchr(dir_[1],'=');
      if (!tmp) tmp = dir_[1];
      else tmp = &tmp[1];
      dir_[1] = strcpy(new char[strlen(tmp)+6+1], tmp);
      strcat(dir_[1], "/basis");
    }
  else {
      dir_[1] = strcpy(new char[strlen(BASISDIR)+1], BASISDIR);
    }
}

BasisFileSet::~BasisFileSet()
{
  int i;
  for (i=0; i<nbasissets_; i++) delete[] basissets_[i];
  delete[] basissets_;
  delete[] dir_[0];
  delete[] dir_[1];
}

RefKeyVal
BasisFileSet::keyval(const RefKeyVal &keyval, const char *basisname)
{
  int i;

  // check if the basisname has already been located
  for (i=0; i<nbasissets_; i++) {
      if (!strcmp(basisname, basissets_[i])) return keyval;
    }

  char *filename = new char[strlen(basisname)+1];

  for (i=0; basisname[i] != '\0'; i++) {
      if (basisname[i] >= 'A' && basisname[i] <= 'Z') {
          filename[i] = tolower(basisname[i]);
        }
      else if (basisname[i] == ','
               || basisname[i] == ' ') {
          filename[i] = '_';
        }
      else if (basisname[i] == '+') {
          filename[i] = 'P';
        }
      else if (basisname[i] == '*') {
          filename[i] = 'S';
        }
      else if (basisname[i] == '(') {
          filename[i] = 'L';
        }
      else if (basisname[i] == ')') {
          filename[i] = 'R';
        }
      else {
          filename[i] = basisname[i];
        }
    }
  filename[i] = '\0';

  // find the basis file
  RefKeyVal newkeyval(keyval);
  for (i=0; i<2; i++) {
      if (!dir_[i]) continue;
      char *path = new char[strlen(dir_[i]) + strlen(filename) + 5];
      strcpy(path, dir_[i]);
      strcat(path, "/");
      strcat(path, filename);
      strcat(path, ".kv");

      // test to see if the file can be opened read only.
      FILE *fp = fopen(path, "r");
      if (fp) {
          fclose(fp);
          RefKeyVal libkeyval = new ParsedKeyVal(path);
          newkeyval = new AggregateKeyVal(keyval,libkeyval);
          delete[] path;
          break;
        }
      delete[] path;
    }

  // add the current basis set to basissets_
  char **newbasissets = new char*[nbasissets_+1];
  for (i=0; i<nbasissets_; i++) newbasissets[i] = basissets_[i];
  newbasissets[nbasissets_] = strcpy(new char[strlen(basisname)+1],
                                     basisname);
  nbasissets_++;
  delete[] basissets_;
  basissets_ = newbasissets;

  delete[] filename;

  return newkeyval;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
