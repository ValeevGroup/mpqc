//$$ except.h                          Exception handling classes

// A set of classes to simulate exceptions in C++
//
//   Partially copied from Carlos Vidal's article in the C users' journal
//   September 1992, pp 19-28
//
//   Operations defined
//      Try {     }
//      Throw ( exception object )
//      Catch ( exception class ) {      }
//      CatchAll {      }
//      CatchAndThrow
//
//   All catch lists must end with a CatchAll or CatchAndThrow statement
//   but not both.
//
//   When exceptions are finally implemented replace Try, Throw, Catch,
//   CatchAll, CatchAndThrow by try, throw, catch, catch(...), and {}.
//
//   All exception classes must be derived from Exception, have no non-static
//   variables and must include functions
//
//      static long st_type()  { return 2; }
//      long type() const { return 2; }
//
//   where 2 is replaced by a prime number unique to the exception class.
//   See notes for use with levels of exceptions.
//

#ifndef EXCEPTION_LIB
#define EXCEPTION_LIB

#include <setjmp.h>

void Terminate();

/*********** classes for setting up exceptions and reporting ***************/

class Exception;

class Tracer                                    // linked list showing how
{                                               // we got here
   char* entry;
   Tracer* previous;
public:
   Tracer(char*);
   ~Tracer();
   void ReName(char*);
   friend class Exception;
};


class Exception                                 // The base exception class
{
public:
   static Tracer* last;                         // points to Tracer list
   static long st_type() { return 1; }
   virtual long type() const { return 1; }
   static void PrintTrace(Boolean=FALSE);       // for printing trace
   friend class Tracer;
   Exception(int action);
};


inline Tracer::Tracer(char* e)
   : entry(e), previous(Exception::last) { Exception::last = this; }

inline Tracer::~Tracer() { Exception::last = previous; }

inline void Tracer::ReName(char* e) { entry=e; }





/************** the definitions of Try, Throw and Catch *******************/


class JumpItem;
class Janitor;

class JumpBase         // pointer to a linked list of jmp_buf s
{
public:
   static JumpItem *jl;
   static long type;                    // type id. of last exception
   static jmp_buf env;
};

class JumpItem         // an item in a linked list of jmp_buf s
{
public:
   JumpItem *ji;
   jmp_buf env;
   Tracer* trace;                     // to keep check on Tracer items
   Janitor* janitor;                  // list of items for cleanup
   JumpItem() : trace(0), janitor(0), ji(JumpBase::jl)
      { JumpBase::jl = this; }
   ~JumpItem() { JumpBase::jl = ji; }
};

void Throw(const Exception& exc);

void Throw();

#define Try                                             \
      if (!setjmp( JumpBase::jl->env )) {               \
      JumpBase::jl->trace = Exception::last;            \
      JumpItem JI387256156;

#define Catch(EXCEPTION)                                \
   } else if (JumpBase::type % EXCEPTION::st_type() == 0) {


#define CatchAll } else

#define CatchAndThrow  } else Throw();



/******************* cleanup heap following Throw ************************/

class Janitor
{
protected:
   static Boolean do_not_link;                  // set when new is called
   Boolean OnStack;                             // false if created by new
public:
   Janitor* NextJanitor;
   virtual void CleanUp() {}
   Janitor();
#ifdef __GNUC__
   virtual ~Janitor();                          // to stop warning messages
#else
   ~Janitor();
#endif
};




#ifdef DO_FREE_CHECK
// Routines for tracing whether new and delete calls are balanced

class FreeCheck;

class FreeCheckLink
{
protected:
   FreeCheckLink* next;
   void* ClassStore;
   FreeCheckLink();
   virtual void Report()=0;                   // print details of link
   friend class FreeCheck;
};

class FCLClass : public FreeCheckLink         // for registering objects
{
   char* ClassName;
   FCLClass(void* t, char* name);
   void Report();
   friend class FreeCheck;
};

class FCLRealArray : public FreeCheckLink     // for registering real arrays
{
   char* Operation;
   int size;
   FCLRealArray(void* t, char* o, int s);
   void Report();
   friend class FreeCheck;
};

class FCLIntArray : public FreeCheckLink     // for registering int arrays
{
   char* Operation;
   int size;
   FCLIntArray(void* t, char* o, int s);
   void Report();
   friend class FreeCheck;
};


class FreeCheck
{
   static FreeCheckLink* next;
public:
   static void Register(void*, char*);
   static void DeRegister(void*, char*);
   static void RegisterR(void*, char*, int);
   static void DeRegisterR(void*, char*, int);
   static void RegisterI(void*, char*, int);
   static void DeRegisterI(void*, char*, int);
   static void Status();
   friend class FreeCheckLink;
   friend class FCLClass;
   friend class FCLRealArray;
   friend class FCLIntArray;
};

#define FREE_CHECK(Class)                                                  \
public:                                                                    \
   void* operator new(size_t size)                                         \
   {                                                                       \
      void* t = ::operator new(size); FreeCheck::Register(t,#Class);       \
      return t;                                                            \
   }                                                                       \
   void operator delete(void* t)                                           \
   { FreeCheck::DeRegister(t,#Class); ::operator delete(t); }

#define NEW_DELETE(Class)                                                  \
public:                                                                    \
   void* operator new(size_t size)                                         \
   {                                                                       \
      do_not_link=TRUE;                                                    \
      void* t = ::operator new(size); FreeCheck::Register(t,#Class);       \
      return t;                                                            \
   }                                                                       \
   void operator delete(void* t)                                           \
   { FreeCheck::DeRegister(t,#Class); ::operator delete(t); }

#define MONITOR_REAL_NEW(Operation, Size, Pointer)                         \
   FreeCheck::RegisterR(Pointer, Operation, Size);
#define MONITOR_INT_NEW(Operation, Size, Pointer)                          \
   FreeCheck::RegisterI(Pointer, Operation, Size);
#define MONITOR_REAL_DELETE(Operation, Size, Pointer)                      \
   FreeCheck::DeRegisterR(Pointer, Operation, Size);
#define MONITOR_INT_DELETE(Operation, Size, Pointer)                       \
   FreeCheck::DeRegisterI(Pointer, Operation, Size);
#else
#define FREE_CHECK(Class) public:
#define MONITOR_REAL_NEW(Operation, Size, Pointer) {}
#define MONITOR_INT_NEW(Operation, Size, Pointer) {}
#define MONITOR_REAL_DELETE(Operation, Size, Pointer) {}
#define MONITOR_INT_DELETE(Operation, Size, Pointer) {}


#define NEW_DELETE(Class)                                                  \
public:                                                                    \
   void* operator new(size_t size)                                         \
   { do_not_link=TRUE; void* t = ::operator new(size); return t; }         \
   void operator delete(void* t) { ::operator delete(t); }

#endif




#endif




