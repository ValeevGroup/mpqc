//$$except.cxx                        Exception handler



#define WANT_STREAM                  // include.h will get stream fns


#include "include.h"                 // include standard files
#include "boolean.h"


#include "except.h"                  // for exception handling

//#define REG_DEREG                    // for print out uses of new/delete
//#define CLEAN_LIST                   // to print entries being added to
                                       // or deleted from cleanup list


void Throw()
{
   for (Janitor* jan = JumpBase::jl->janitor; jan; jan = jan->NextJanitor)
      jan->CleanUp();
   JumpBase::jl = JumpBase::jl->ji;
   if ( ! JumpBase::jl ) Terminate();
   Exception::last = JumpBase::jl->trace;
   longjmp(JumpBase::jl->env, 1);
}

void Throw(const Exception& exc) { JumpBase::type = exc.type(); Throw(); }


void Exception::PrintTrace(Boolean)
{
   cout << "\n";
   {
      for (Tracer* et = last; et; et=et->previous)
	 cout << "  * " << et->entry << "\n";
   }
}

Exception::Exception(int action)
{
   if (action)
   {
      cout << "\nAn exception has occurred: call trace follows.";
      PrintTrace();
      if (action < 0) exit(1);
   }
}

Janitor::Janitor()
{
   if (do_not_link)
   {
      do_not_link = FALSE; NextJanitor = 0; OnStack = FALSE;
#ifdef CLEAN_LIST
      cout << "Not added to clean-list " << (unsigned long)this << "\n";
#endif
   }
   else
   {
      OnStack = TRUE;
#ifdef CLEAN_LIST
      cout << "Add to       clean-list " << (unsigned long)this << "\n";
#endif
      NextJanitor = JumpBase::jl->janitor; JumpBase::jl->janitor=this;
   }
}

Janitor::~Janitor()
{
   // expect the item to be deleted to be first on list
   // but must be prepared to search list
   if (OnStack)
   {
#ifdef CLEAN_LIST
      cout << "Delete from  clean-list " << (unsigned long)this << "\n";
#endif
      Janitor* lastjan = JumpBase::jl->janitor;
      if (this == lastjan) JumpBase::jl->janitor = NextJanitor;
      else
      {
	 for (Janitor* jan = lastjan->NextJanitor; jan;
	    jan = lastjan->NextJanitor)
	 {
	    if (jan==this)
	       { lastjan->NextJanitor = jan->NextJanitor; return; }
	    lastjan=jan;
	 }

         cout << "\nCannot resolve memory linked list\n";
         cout << "See notes in except.cxx for details\n";
	 Throw(Exception(-1));
/*
This message occurs when a call to ~Janitor() occurs, apparently without a
corresponding call to Janitor(). This could happen if my way of deciding
whether a constructor is being called by new fails. Possibly also if
delete is applied an object on the stack (ie not called by new). Otherwise,
it is a bug in Newmat or your compiler. If you don't #define
TEMPS_DESTROYED_QUICKLY you will get this error with Microsoft C 7.0. There
are probably situations where you will get this when you do define
TEMPS_DESTROYED_QUICKLY. This is a bug in MSC. Beware of "operator" statements
for defining conversions; particularly for converting from a Base class to
a Derived class. 

You may get away with simply deleting this error message and Throw statement
if you can't find a better way of overcoming the problem. In any case please
tell me if you get this error message, particularly for compilers apart from
Microsoft C.
*/
      }
   }
}

JumpItem* JumpBase::jl = 0;
long JumpBase::type;
Tracer* Exception::last = 0;

Boolean Janitor::do_not_link = FALSE;

static JumpItem JI;                  // need JumpItem at head of list



void Terminate()
{
   cout << "\nThere has been an exception with no handler - exiting\n";
   exit(1);
}






#ifdef DO_FREE_CHECK
// Routines for tracing whether new and delete calls are balanced

FreeCheckLink::FreeCheckLink() : next(FreeCheck::next)
   { FreeCheck::next = this; }

FCLClass::FCLClass(void* t, char* name) : ClassName(name) { ClassStore=t; }
   
FCLRealArray::FCLRealArray(void* t, char* o, int s)
  : Operation(o), size(s) { ClassStore=t; }

FCLIntArray::FCLIntArray(void* t, char* o, int s)
  : Operation(o), size(s) { ClassStore=t; }

FreeCheckLink* FreeCheck::next = 0;

void FCLClass::Report()
{ cout << "   " << ClassName << "   " << (unsigned long)ClassStore << "\n"; }

void FCLRealArray::Report()
{
   cout << "   " << Operation << "   " << (unsigned long)ClassStore << 
      "   " << size << "\n";
}

void FCLIntArray::Report()
{
   cout << "   " << Operation << "   " << (unsigned long)ClassStore << 
      "   " << size << "\n";
}

void FreeCheck::Register(void* t, char* name)
{
   FCLClass* f = new FCLClass(t,name);
   if (!f) { cout << "Out of memory in FreeCheck\n"; exit(1); }
#ifdef REG_DEREG
   cout << "Registering   " << name << "   " << (unsigned long)t << "\n";
#endif
}

void FreeCheck::RegisterR(void* t, char* o, int s)
{
   FCLRealArray* f = new FCLRealArray(t,o,s);
   if (!f) { cout << "Out of memory in FreeCheck\n"; exit(1); }
#ifdef REG_DEREG
   cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
}

void FreeCheck::RegisterI(void* t, char* o, int s)
{
   FCLIntArray* f = new FCLIntArray(t,o,s);
   if (!f) { cout << "Out of memory in FreeCheck\n"; exit(1); }
#ifdef REG_DEREG
   cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
}

void FreeCheck::DeRegister(void* t, char* name)
{
   FreeCheckLink* last = 0;
#ifdef REG_DEREG
   cout << "Deregistering " << name << "   " << (unsigned long)t << "\n";
#endif
   for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next)
   {
      if (fcl->ClassStore==t)
      {
	 if (last) last->next = fcl->next; else next = fcl->next;
	 delete fcl; return;
      }
      last = fcl;
   }
   cout << "\nRequest to delete non-existent object of class and location:\n";
   cout << "   " << name << "   " << (unsigned long)t << "\n";
   Exception::PrintTrace(TRUE);
   cout << "\n";
}

void FreeCheck::DeRegisterR(void* t, char* o, int s)
{
   FreeCheckLink* last = 0;
#ifdef REG_DEREG
   cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
   for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next)
   {
      if (fcl->ClassStore==t)
      {
	 if (last) last->next = fcl->next; else next = fcl->next;
	 if (((FCLRealArray*)fcl)->size != s)
	 {
	    cout << "\nArray sizes don't agree:\n";
	    cout << "   " << o << "   " << (unsigned long)t
	       << "   " << s << "\n";
	    Exception::PrintTrace(TRUE);
	    cout << "\n";
	 }
	 delete fcl; return;
      }
      last = fcl;
   }
   cout << "\nRequest to delete non-existent real array:\n";
   cout << "   " << o << "   " << (unsigned long)t << "   " << s << "\n";
   Exception::PrintTrace(TRUE);
   cout << "\n";
}

void FreeCheck::DeRegisterI(void* t, char* o, int s)
{
   FreeCheckLink* last = 0;
#ifdef REG_DEREG
   cout << o << "   " << s << "   " << (unsigned long)t << "\n";
#endif
   for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next)
   {
      if (fcl->ClassStore==t)
      {
	 if (last) last->next = fcl->next; else next = fcl->next;
	 if (((FCLIntArray*)fcl)->size != s)
	 {
	    cout << "\nArray sizes don't agree:\n";
	    cout << "   " << o << "   " << (unsigned long)t
	       << "   " << s << "\n";
	    Exception::PrintTrace(TRUE);
	    cout << "\n";
	 }
	 delete fcl; return;
      }
      last = fcl;
   }
   cout << "\nRequest to delete non-existent int array:\n";
   cout << "   " << o << "   " << (unsigned long)t << "   " << s << "\n";
   Exception::PrintTrace(TRUE);
   cout << "\n";
}

void FreeCheck::Status()
{
   if (next)
   {
      cout << "\nObjects of the following classes remain undeleted:\n";
      for (FreeCheckLink* fcl = next; fcl; fcl = fcl->next) fcl->Report();
      cout << "\n";
   }
   else cout << "\nNo objects remain undeleted\n";
}

#endif


