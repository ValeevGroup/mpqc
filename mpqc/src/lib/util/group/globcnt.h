
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_globcnt_h
#define _util_group_globcnt_h

// A process can create a GlobalCounter using the void CTOR.
// This process can share the string representation of the
// counter with other processes.  They can then use the const
// char * CTOR to create global counters that reference the
// same global counter.

class GlobalCounter {
  private:
    int semid_;
    int controls_release_;
  public:
    GlobalCounter();
    void initialize();
    void initialize(const char *stringrep);
    ~GlobalCounter();

    char *stringrep();

    void wait_for_zero();
    void operator += (int);
    void operator ++();
    void operator --();
    void operator ++(int) { operator++(); }
    void operator --(int) { operator--(); }
    void operator = (int);

    int val();
};

#endif
