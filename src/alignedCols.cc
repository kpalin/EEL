#include <Python.h>

#define ALIGNEDCOLS_MODULE
#include "alignedCols.h"
#include "structmember.h"

/*
 * $Log$
 * Revision 1.2  2004/12/22 11:14:34  kpalin
 * Some fixes for better distributability
 *
 * Revision 1.1  2004/07/30 12:12:21  kpalin
 * Workings for alnColumn objects. Mostly stolen from align.cc
 *
 */



//////////////////////////////////////////////////////////////////////
// Aligned site object


typedef struct {
  PyObject_HEAD

    /* Type-specific fields go here. */

  char *motifName;
  PyObject *seqID; // Tuple of sequence codes
  PyObject *beginEnd; // Tuple of 2-tuples of (begin,end) pairs
  PyObject *siteScore; // Tuple of site scores.
  PyObject *siteSeqPos; // Positions on the site sequence.
  PyObject *annotation; // Additional annotations for the site.
  int seqX;
  int seqY;
  int beginX;
  int endX;
  int beginY;
  int endY;
  char strand;
  double score;
  double siteScoreX;
  double siteScoreY;
} align_siteObject;


static void site_dealloc(align_siteObject* self);


static PyMemberDef site_members[] = {
    {"motif",T_STRING, offsetof(align_siteObject, motifName), 0,
     "Name of the motif."},
    {"seqID",T_OBJECT_EX, offsetof(align_siteObject,seqID),0,
     "Codes of the sequences"},
    {"beginEnd",T_OBJECT_EX, offsetof(align_siteObject,beginEnd),0,
     "2-tuples of (begin,end) positions."},
    {"siteScore",T_OBJECT_EX, offsetof(align_siteObject,siteScore),0,
     "Scores of the sites"},
    {"siteSeqPos",T_OBJECT_EX, offsetof(align_siteObject,siteSeqPos),0,
     "Positions on the site sequence."},
    {"annotation",T_OBJECT, offsetof(align_siteObject,annotation),0,
     "Additional annotation for the sites."},

    {"seqX",T_INT, offsetof(align_siteObject, seqX), 0,
     "Position on site sequence x."},
    {"seqY",T_INT, offsetof(align_siteObject, seqY), 0,
     "Position on site sequence y."},

    {"beginX",T_INT, offsetof(align_siteObject, beginX), 0,
     "Begin position on DNA sequence x."},

    {"beginY",T_INT, offsetof(align_siteObject, beginY), 0,
     "Begin position on DNA sequence y."},

    {"endX",T_INT, offsetof(align_siteObject, endX), 0,
     "End position on DNA sequence x."},

    {"endY",T_INT, offsetof(align_siteObject, endY), 0,
     "End position on DNA sequence y."},
    {"strand",T_CHAR, offsetof(align_siteObject, strand), 0,

     "Strand of the motif."},
    {"score",T_DOUBLE, offsetof(align_siteObject, score), 0,
     "Alignment score this far."},
    {"siteScoreX",T_DOUBLE, offsetof(align_siteObject, siteScoreX), 0,
     "Score of the site on sequence X."},
    {"siteScoreY",T_DOUBLE, offsetof(align_siteObject, siteScoreY), 0,
     "Score of the site on sequence Y."},
    {NULL}  /* Sentinel */
};

static PyTypeObject align_siteType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "alignedCols.alnColumn",             /*tp_name*/
    sizeof(align_siteObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)site_dealloc,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Site pair object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    0,             /* tp_methods */
    site_members,             /* tp_members */
};

static PyMethodDef alnColMethods[] = {
//   {"alignfile",  align_alignfile, METH_VARARGS,
//    "aligns sequences from a gff-file"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};
static PyObject *
site_new(const char *motifName,
	 int seqX,
	 int seqY,
	 int beginX,
	 int endX,
	 int beginY,
	 int endY,
	 char strand,
	 double score,
	 double siteScoreX,
	 double siteScoreY,
	 const char *annotX,
	 const char *annotY)
{
  align_siteObject *ret;
  ret=(align_siteObject*)align_siteType.tp_alloc(&align_siteType,0);

  ret->motifName=new char[strlen(motifName)+1];
  strcpy(ret->motifName,motifName);


  ret->seqID=Py_BuildValue("(ii)",0,1);
  ret->beginEnd=Py_BuildValue("((ii)(ii))",beginX,endX,beginY,endY);
  ret->siteScore=Py_BuildValue("(dd)",siteScoreX,siteScoreY);
  ret->siteSeqPos=Py_BuildValue("(ii)",seqX,seqY);
  ret->annotation=Py_BuildValue("(zz)",annotX,annotY);
  ret->seqX=seqX;
  ret->seqY=seqY;
  ret->beginX=beginX;
  ret->endX=endX;
  ret->beginY=beginY;
  ret->endY=endY;
  ret->strand=strand;
  ret->score=score;
  ret->siteScoreX=siteScoreX;
  ret->siteScoreY=siteScoreY;

  return (PyObject*)ret;
}

static PyObject *
site_new_multi(const char *motifName,
	       PyObject *seqCodes,
	       PyObject *seqPos,
	       PyObject *beginEnd,
	       char strand,
	       double score,
	       PyObject  *siteScore,
	       PyObject *annotation)
{
  align_siteObject *ret;
  ret=(align_siteObject*)align_siteType.tp_alloc(&align_siteType,0);

  ret->motifName=new char[strlen(motifName)+1];
  strcpy(ret->motifName,motifName);


  ret->seqID=seqCodes;
  ret->beginEnd=beginEnd;
  ret->siteScore=siteScore;
  ret->siteSeqPos=seqPos;
  ret->strand=strand;
  ret->score=score;
  ret->annotation=annotation;

  ret->seqX=-1;
  ret->seqY=-1;
  ret->beginX=-1;
  ret->endX=-1;
  ret->beginY=-1;
  ret->endY=-1;
  ret->siteScoreX=-1;
  ret->siteScoreY=-1;

  return (PyObject*)ret;
}








static void site_dealloc(align_siteObject* self)
{
  delete [] self->motifName;

  self->ob_type->tp_free((PyObject*)self);
}

extern "C"
void
initalignedCols(void)
{
  PyObject *m;

  static void *PyAln_API[PyAln_API_pointers];
  PyObject *c_api_object;





  m=Py_InitModule("eellib.alignedCols", alnColMethods);
  
  if(m==NULL)
    return;


  /* Initialize the C API pointer array */
  PyAln_API[PyAln_New_NUM] = (void *)site_new;
  PyAln_API[PyAln_New_multi_NUM] = (void *)site_new_multi;

  /* Create a CObject containing the API pointer array's address */
  c_api_object = PyCObject_FromVoidPtr((void *)PyAln_API, NULL);

  if (c_api_object != NULL)
    PyModule_AddObject(m, "_C_API", c_api_object);



  align_siteType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&align_siteType) < 0)
    return;
  
  
  Py_INCREF(&align_siteType);
  PyModule_AddObject(m, "alnColumn", (PyObject *)&align_siteType);
}


//////////////////////////////////////////////////////////////////////

