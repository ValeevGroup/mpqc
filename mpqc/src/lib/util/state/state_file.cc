
#include <stdio.h>
#include <util/class/class.h>
#include "state.h"

#include "statenumSet.h"
#include "classdImplMap.h"

StateOutFile::StateOutFile() :
  opened_(0), fp_(stdout)
{
}

StateOutFile::StateOutFile(FILE* fp) :
  opened_(0), fp_(fp)
{
}

StateOutFile::StateOutFile(const char * path, const char * mode) :
  opened_(1)
{
  fp_ = fopen(path,mode);
}

StateOutFile::~StateOutFile()
{
  if (opened_) close();
}

void StateOutFile::flush() { fflush(fp_); }
void StateOutFile::close()
{
  if(opened_) fclose(fp_);
  opened_=0; fp_=0;

  _classidmap->clear(); _nextclassid=0;
  forget_references();
}

void StateOutFile::rewind() { if(fp_) fseek(fp_,0,0); }

int StateOutFile::open(const char *path, const char * mode)
{
  if (opened_) close();

  if ((fp_ = fopen(path,mode))==0) {
    fprintf(stderr,"StateOutFile::open(%s,%s) failed\n",path,mode);
    return -1;
    }

  opened_ = 1;
  return 0;
}

////////////////////////////////////

StateInFile::StateInFile() :
  opened_(0), fp_(stdin)
{
}

StateInFile::StateInFile(FILE* fp) :
  opened_(0), fp_(fp)
{
}

StateInFile::StateInFile(const char * path, const char * mode) :
  opened_(1)
{
  fp_ = fopen(path,mode);
}

StateInFile::~StateInFile()
{
  if (opened_) close();
}

void StateInFile::flush() { fflush(fp_); }
void StateInFile::close()
{
  if(opened_) fclose(fp_);
  opened_=0; fp_=0;

  _cd.clear();
  forget_references();
}
void StateInFile::rewind() { if(fp_) fseek(fp_,0,0); }

int StateInFile::open(const char *path, const char * mode)
{
  if (opened_) close();

  if ((fp_ = fopen(path,mode))==0) {
    fprintf(stderr,"StateInFile::open(%s,%s) failed\n",path,mode);
    return -1;
    }

  opened_ = 1;
  return 0;
}

