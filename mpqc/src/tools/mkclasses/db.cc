
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <errno.h>

#include "mkclasses.h"

void
MkClasses::writedb()
{
  if (options_.dbname().empty() || options_.libname().empty()) return;

  string classes;

  vector<Class>::iterator ci;
  for (ci=classes_.begin(); ci != classes_.end(); ci++) {
      if (!classes.empty()) classes += " ";
      classes += (*ci).name();
    }

  // open the old database
  FILE* db;
  db = fopen(options_.dbname().c_str(), "r");
  if (!db && errno != ENOENT) {
      fprintf(stderr,"writedb: couldn't open dbfile \"%s\"\n",
              options_.dbname().c_str());
      perror("fopen");
      exit(1);
    }

  // open the new database (a temporary file)
  char* dbnewname = new char[strlen(options_.dbname().c_str())+4+1];
  strcpy(dbnewname, options_.dbname().c_str());
  strcat(dbnewname, ".tmp");
  FILE* dbnew = fopen(dbnewname,"w");
  if (!dbnew) {
      fprintf(stderr,"writedb: couldn't open dbfile tmp file \"%s\"\n",
              dbnewname);
      perror("fopen");
      exit(1);
    }

  // copy the old database to the new data base, making necessary
  // changes.
  const int bufsize = 10000;
  char buf[bufsize];
  int ncharlibname = strlen(options_.libname().c_str());
  fprintf(dbnew, "%s %s\n", options_.libname().c_str(), classes.c_str());
  if (db) {
      while(fgets(buf, bufsize, db)) {
          if (buf[0] != '\0' && buf[strlen(buf)-1] == '\n') {
              buf[strlen(buf)-1] = '\0';
            }
          if (strncmp(options_.libname().c_str(),buf,ncharlibname) == 0
              && buf[ncharlibname] == ' ') {
              // skip output, this library's classes have already been written
            }
          else if (buf[0] != '\0' && buf[0] != '\n') {
              fprintf(dbnew, "%s\n", buf);
            }
        }
    }

  // update the file
  if (rename(dbnewname, options_.dbname().c_str())) {
      fprintf(stderr,"writedb: couldn't rename dbfile\n");
      perror("rename");
      exit(1);
    }
}
