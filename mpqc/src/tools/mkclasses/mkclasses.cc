
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <errno.h>

int
main(int argc, char** argv)
{
  int i;

  // process the command line
  char* dbname;
  char* libname;
  char* classes;
  int classes_length;

  if (argc < 4) {
      fprintf(stderr,"Usage: %s dbfile libname class1 ...\n",argv[0]);
      return 1;
    }
  dbname = argv[1];
  libname = argv[2];

  classes_length = 1;
  for (i=3; i<argc; i++) {
      classes_length += strlen(argv[i]) + 1;
    }
  classes = new char[classes_length];
  classes[0] = '\0';
  for (i=3; i<argc; i++) {
      strcat(classes,argv[i]);
      strcat(classes," ");
    }
  classes[classes_length-1] = '\0';

  // open the old database
  FILE* db;
  db = fopen(dbname, "r");
  if (!db && errno != ENOENT) {
      fprintf(stderr,"%s: couldn't open dbfile \"%s\"\n", argv[0], dbname);
      perror("fopen");
      exit(1);
    }

  // open the new database (a temporary file)
  char* dbnewname = new char[strlen(dbname)+4+1];
  strcpy(dbnewname, dbname);
  strcat(dbnewname, ".tmp");
  FILE* dbnew = fopen(dbnewname,"w");
  if (!dbnew) {
      fprintf(stderr,"%s: couldn't open dbfile tmp file \"%s\"\n",
              argv[0], dbnewname);
      perror("fopen");
      exit(1);
    }

  // copy the old database to the new data base, making necessary
  // changes.
  const int bufsize = 10000;
  char buf[bufsize];
  int ncharlibname = strlen(libname);
  fprintf(dbnew, "%s %s\n", libname, classes);
  if (db) {
      while(fgets(buf, bufsize, db)) {
          if (buf[0] != '\0' && buf[strlen(buf)-1] == '\n') {
              buf[strlen(buf)-1] = '\0';
            }
          if (strncmp(libname,buf,ncharlibname) == 0
              && buf[ncharlibname] == ' ') {
              // skip output, this library's classes have already been written
            }
          else if (buf[0] != '\0' && buf[0] != '\n') {
              fprintf(dbnew, "%s\n", buf);
            }
        }
    }

  // update the file
  if (rename(dbnewname, dbname)) {
      fprintf(stderr,"%s: couldn't rename dbfile\n", argv[0]);
      perror("rename");
      exit(1);
    }

  return 0;
}
