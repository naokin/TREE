#ifndef __TTNS_THROW_HEADER
#define __TTNS_THROW_HEADER

/* default C++ exceptions */
#include <stdexcept>
#define TTNS_THROW(truth, msg)\
{ if (!(truth)) { throw std::runtime_error(msg); } }

/* redefine for other use-cases */
#ifdef THROW_PYTHON_EXCEPTIONS
#include <Python.h>
#define TTNS_THROW(truth, msg)\
{ if (!(truth)) { PyErr_SetString(PyExc_RuntimeError, msg); return -1; } }
#endif

#ifdef THROW_C_EXCEPTIONS
#include <iostream>
#define TTNS_THROW(truth, msg)\
{ if (!(truth)) { std::cout << msg << std::endl; return -1; } }
#endif

#include <iostream>
#define TTNS_ASSERT(truth, msg) { if (!(truth)) { std::cout << msg << std::endl; } }

#ifdef _ENABLE_TTNS_DEBUG
#define TTNS_DEBUG(msg) { std::cout << "DEBUG $ " << msg << std::endl; }
#else
#define TTNS_DEBUG(msg) { }
#endif

//namespace ttns { using ::operator<<; }

#endif // __TTNS_THROW_HEADER
