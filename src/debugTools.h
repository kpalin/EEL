#ifndef _DEBUGTOOLS_H_
#define _DEBUGTOOLS_H_

#include<Python.h>

#include<stdio.h>
#include<stdarg.h>

inline int getRefCount(PyObject *obj)
{
  return obj->ob_refcnt;
}



#ifndef NDEBUG
void _printDfun(char *file, unsigned int line,char *fmt,...)
{
  va_list argp;
  fprintf(stderr, "%s:%d: ",file,line);
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");

}
   	

#define printDebug(...) _printDfun(__FILE__,__LINE__,__VA_ARGS__)
#else
#define printDebug(...) ;
#endif
   
//
// $Log$
//

#endif
