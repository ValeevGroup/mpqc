
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <util/class/scexception.h>
#include <chemistry/qc/ccr12/parenthesis2tnum.h>

using namespace sc;

static ClassDesc Parenthesis2tNum_cd(
  typeid(Parenthesis2tNum),"Parenthesis2tNum",1,"virtual public RefCount"
  ,0,0,0);

Parenthesis2tNum::Parenthesis2tNum(CCR12_Info* info): z(info){
}
    
Parenthesis2tNum::~Parenthesis2tNum(){
}

void Parenthesis2tNum::compute_amp(double*,const long,const long,const long,
                                   const long,const long,const long,const long){
} 

void Parenthesis2tNum::compute_amp(double*,const long,const long,const long,const long,
                                   const long,const long,const long,const long,const long){
} 
