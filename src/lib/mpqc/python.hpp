#ifndef PYTHON_HPP
#define PYTHON_HPP

#ifdef HAVE_PYTHON

//#include <Python.h>
#include <stdio.h>
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>

// template<typename T>
// void operator%(const std::string &fmt, const T &t) {
// }

struct Python {
private:

    struct initializer {
        initializer() {
            using namespace boost::python;
            if (!Py_IsInitialized()) {
                Py_InitializeEx(0);
                boost::python::object main = boost::python::import("__main__");
                const char* startup = getenv("PYTHONSTARTUP");
                if (startup) {
                    FILE *fp = fopen(startup, "r");
                    if (fp) PyRun_SimpleFile(fp, startup);
                    fclose(fp);
                }

                // disable catching Ctrl-C
                PyRun_SimpleString("import signal\n"
                                   "signal.signal(signal.SIGINT, signal.SIG_DFL)\n");

                // PyRun_SimpleString("import signal\n"
                //                    "import sys\n"
                //                    "def sigint(signal, frame):  sys.exit(0)\n"
                //                    "signal.signal(signal.SIGINT, sigint)\n");

                PyImport_AddModule("mpqc");
                PyRun_SimpleString("import mpqc");
            }
        }
        // must only call once to avoid crashes with boost::python
        ~initializer() {
            //std::cout << "Py_Finalize" << std::endl;
            Py_Finalize();
        }
    };

public:

    Python() {
        static initializer init;
    }

    struct tuple : boost::python::tuple {
	template<typename T>
	tuple(T t)
	    : boost::python::tuple(boost::python::make_tuple(t)) {}
	template<typename T0, typename T1>
	tuple(T0 t0, T1 t1)
	    : boost::python::tuple(boost::python::make_tuple(t0,t1)) {}
	template<typename T0, typename T1, typename T2, typename T3>
	tuple(T0 t0, T1 t1, T2 t2, T3 t3)
	    : boost::python::tuple(boost::python::make_tuple(t0,t1,t2,t3)) {}
    };

    void error() {
	PyErr_Print();
    };

    void exec(const std::string &cmd) {
	std::cout << cmd << std::endl;
	PyRun_SimpleString(cmd.c_str());
    }
    
    void exec(const boost::python::object &cmd) {
        exec(boost::python::extract<std::string>(cmd));
    }

    void interactive() {
	namespace py = boost::python;
	//update(this->main_module, this->main);
	PyRun_InteractiveLoop(stdin, "<stdin>");
    }

    // void update(boost::python::object ns, boost::python::dict kv) {
    //     namespace py = boost::python;
    //     py::list keys = py::list(kv.iterkeys());
    //     for (size_t i = 0; i < py::len(keys); ++i) {
    //         ns.attr(keys[i]) = kv[keys[i]];
    //     }
    // }

    // 	try {
    // 	    py::object main =
    // 		py::object(py::handle<>(py::borrowed(PyImport_AddModule("__main__"))));
    // 	    py::object ns = main.attr("__dict__");
    // 	    ns["value"] = "fuck";
    // 	    PyRun_SimpleString("import readline\n");
    // 	    // boost::python::exec("import readline\n");
    // 	    PyRun_InteractiveLoop(stdin, "<stdin>");
    // 	}
    // 	catch (...) {
    // 	    PyErr_Print();
    // 	}
    // }

};


#endif // HAVE_PYTHON

#endif /* PYTHON_HPP */
