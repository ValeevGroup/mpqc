
#include <iostream.h>
#include <util/keyval/keyval.h>
#include <util/misc/bug.h>

double bugtest_global_var_;
double t;

void
y()
{
  t = 1.0/bugtest_global_var_;
}

void
x(const RefDebugger &d)
{
  d->traceback();
  y();
}

int
main(int argc, char **argv)
{
  const char *infile = SRCDIR "/bugtest.in";

  RefKeyVal keyval = new ParsedKeyVal(infile);

  RefDebugger d;

  d = keyval->describedclassvalue("debug");
  if (d.null()) {
      d = new Debugger();
      d->handle_defaults();
      d->set_prefix(999);
      d->set_traceback_on_signal(1);
      d->set_debug_on_signal(1);
      d->set_exit_on_signal(0);
    }
  d->set_exec(argv[0]);

  d->traceback("no particular problem");
  //d->debug("no particular problem");
  x(d);

  StateOutText o("state.dat");
  d.save_state(o);
  o.flush();

  StateInText i("state.dat");
  d.restore_state(i);

  d = 0;
  cout << "bugtest: done" << endl;

  return 0;
}
