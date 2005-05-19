#ifndef _ALIGNED_COLS_H_
#define _ALIGNED_COLS_H_


/*
 * $Log$
 * Revision 1.2.2.1  2005/05/19 06:45:07  kpalin
 * Added path for the alignedCols module import
 *
 * Revision 1.2  2005/02/24 11:37:29  kpalin
 * Site annotations.
 *
 * Revision 1.1  2004/07/30 12:11:07  kpalin
 * C api header file for alnColumn objects.
 *
 */

#define PyAln_New_NUM 0
#define PyAln_New_RETURN PyObject*
#define PyAln_New_PROTO (const char *motifName, \
	 int seqX,\
	 int seqY,\
	 int beginX,\
	 int endX,\
	 int beginY,\
	 int endY,\
	 char strand,\
	 double score,\
	 double siteScoreX,\
	 double siteScoreY, \
	 const char *annotX, \
         const char *annotY)

#define PyAln_New_multi_NUM 1
#define PyAln_New_multi_RETURN PyObject*
#define PyAln_New_multi_PROTO (const char *motifName, \
	       PyObject *seqCodes,\
	       PyObject *seqPos,\
	       PyObject *beginEnd,\
	       char strand,\
	       double score,\
	       PyObject  *siteScore, \
	       PyObject *annotation)

#define PyAln_API_pointers 2


#ifdef ALIGNEDCOLS_MODULE

static PyAln_New_RETURN site_new PyAln_New_PROTO;
static PyAln_New_multi_RETURN site_new_multi PyAln_New_multi_PROTO;

#else

static void **PyAln_API;

#define PyAln_New \
(*(PyAln_New_RETURN (*)PyAln_New_PROTO) PyAln_API[PyAln_New_NUM])

#define PyAln_New_Multi \
(*(PyAln_New_multi_RETURN (*)PyAln_New_multi_PROTO) PyAln_API[PyAln_New_multi_NUM])

static int
import_alnCols(void)
{
  PyObject *module = PyImport_ImportModule("eellib.alignedCols");

  if(module!=NULL) {
    PyObject *c_api_object=PyObject_GetAttrString(module,"_C_API");
    if(c_api_object==NULL)
      return -1;
    if(PyCObject_Check(c_api_object))
      PyAln_API=(void**)PyCObject_AsVoidPtr(c_api_object);
    Py_DECREF(c_api_object);
  }
  return 0;

}

#endif

#endif
