#include <Python.h>
#include "structmember.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <list>
#include <deque>
#include <time.h>
#include <stdio.h>
#include <math.h>
using namespace std;


#ifndef SEQ_BUFFER_SIZE
#define SEQ_BUFFER_SIZE 15000000
#endif



unsigned long py_fileLikeTell(PyObject *py_file);
int py_fileLikeSeek(PyObject *py_file, unsigned long pos);


//gets a Python list of lists (Matrix)
vector<vector<double> > parse(PyObject* listoflists)
{
  //do some checking (maybe not necessary)


  if (!listoflists || !PySequence_Check(listoflists)  || 
      PySequence_Size(listoflists)<4 ||
      !PySequence_Check(PyList_GetItem(listoflists, 0)))
    {
      PyErr_SetString(PyExc_ValueError,"Wrong Matrix format");
      vector<vector<double> > ret;
      return ret;
    }
  
  PyObject* line0 = PyList_GetItem(listoflists, 0);
  PyObject* line1 = PyList_GetItem(listoflists, 1);
  PyObject* line2 = PyList_GetItem(listoflists, 2);
  PyObject* line3 = PyList_GetItem(listoflists, 3);

  vector<double> help(PyList_Size(line0));
  vector<vector<double> > ret(4, help);
  
  for(unsigned int i=0; i<help.size(); i++)
    {
      ret[0][i] = PyFloat_AsDouble(PyList_GetItem(line0, i));
      ret[1][i] = PyFloat_AsDouble(PyList_GetItem(line1, i));
      ret[2][i] = PyFloat_AsDouble(PyList_GetItem(line2, i));
      ret[3][i] = PyFloat_AsDouble(PyList_GetItem(line3, i));
      //printf("%g %g %g %g\n",ret[0][i],ret[1][i],ret[2][i],ret[3][i]);
    }
  return ret;
}

void drawvector(vector<vector<double> > ret)
{
  vector<vector<double> >::iterator ret_iter;
  vector<double>::iterator ret_help;
  for(ret_iter=ret.begin(); ret_iter < ret.end();ret_iter++)
    {
      for(ret_help= ret_iter->begin(); ret_help < ret_iter->end(); ret_help++)
	{
	  cout<<*ret_help<<", ";
	}
      cout<<endl;
    }
}  

static PyObject *
matrix_draw(PyObject *self, PyObject *args)
{
  drawvector(parse(PyTuple_GetItem(args, 0)));

  Py_INCREF(Py_None);
  return Py_None;
}

//match a matrix on a sequence
static PyObject *
matrix_match(PyObject *self, PyObject *args)
{
  char* Seq= PyString_AsString (PyTuple_GetItem(args, 1));
  int seq_length=PyString_Size (PyTuple_GetItem(args, 1));
  int i, j;
  int mat_length;
  int nucleotide, compl_nucleotide;

  vector<vector<double> > M= parse(PyTuple_GetItem(args, 0));
  mat_length= M[0].size();

  // histories store the interim results of the matching
  // one for each the 2 complement strangs
  deque<double> history(mat_length,0);
  deque<double> compl_history(mat_length,0);

  // iterators to work with the histories
  // notice that iter is a reverse_iterator
  deque<double>::reverse_iterator iter;
  deque<double>::iterator compl_iter;
  

  cout << "THIS PROBABLY WON'T WORK"<<endl;

  PyObject* ret, *entry;
  if(!(ret=PyList_New(seq_length))) 
    {
      cerr << "could not allocate return list"<<endl;
    }
  
  for (i=0; i<seq_length; i++)
    {
      /*      //the score at position i is the max of the fronts of the histories
      entry=(history.front()<compl_history.front() ?
	     PyFloat_FromDouble(compl_history.front()):
	     PyFloat_FromDouble(history.front()));
      PyList_SET_ITEM(ret, i, entry);
      */
      switch(Seq[i])
	{
	case 'A':
	  //nucleotide stores in which line of the matrix is to search
	  nucleotide=0;        
	  compl_nucleotide=3;
	  break;
	case 'C':
	  nucleotide=1;
	  compl_nucleotide=2;
	  break;
	case 'G':
	  nucleotide=2;
	  compl_nucleotide=1;
	  break;
	case 'T':
	  nucleotide=3;
	  compl_nucleotide=0;
	  break;
	default:
	  cout<<"."<<flush;
	    //cout<<"Wrong letter in Sequence! Reading it like 'N'"<<endl;
	case 'N': 
	case 'X':
	  //histories are used like queues
	  history.pop_front();
	  history.push_back(0);
	  compl_history.pop_front();
	  compl_history.push_back(0);

	  // in case of 'N' a specific value is added to the histories
	  // to make sure, that it becomes no hit
	  iter=history.rbegin();
	  compl_iter=compl_history.begin();
	  for(j=0; j<mat_length; j++)
	    {
	      *iter+=-10;              //what do I add best to the histories?
	      ++iter;  // notice that iter is a reverse_iterator
	      
	      *compl_iter+=-10;        //what do I add best to the histories?
	      ++compl_iter;
	    }
	  
	  continue;
	  break;
	}

      //histories are used like queues
      history.pop_front();
      history.push_back(0);
      compl_history.pop_front();
      compl_history.push_back(0);

      iter=history.rbegin();
      compl_iter=compl_history.begin();
      for(j=0; j<mat_length; j++)
	{
	  *iter+= M[nucleotide][j];
	  ++iter;  // notice that iter is a reverse_iterator

	  *compl_iter+=M[compl_nucleotide][j];
	  ++compl_iter;
	}

      entry=(history.front()<compl_history.front() ?
	     PyFloat_FromDouble(compl_history.front()):
	     PyFloat_FromDouble(history.front()));
      PyList_SET_ITEM(ret, i, entry);

    }
      
  return ret;
}


void addMatch(PyObject *dict,int const pos,char const strand,double const score)
{
  PyDict_SetItem(dict, Py_BuildValue("(ic)",pos,strand),PyFloat_FromDouble(score));

}


void addMatchWithKey(PyObject *dict,PyObject *key,int const pos,char const strand,double const score)
{
  PyObject *subDict;
  if(!PyMapping_HasKey(dict,key)) {
    subDict=PyDict_New();
    PyDict_SetItem(dict,key,subDict);
  } else {
    subDict=PyDict_GetItem(dict,key);
  }
  addMatch(subDict,pos,strand,score);

}



typedef unsigned long int bit32;
bit32 const nucl_A=0,nucl_C=1,nucl_G=2,nucl_T=3;

struct __BGdataCPP {
  vector<unsigned long int> counts;
  vector<unsigned long int> shortCounts;
};

typedef struct {
  PyObject_HEAD

  PyObject *bgSample;
  //unsigned int sampleLen;

  PyObject *strBuf,*py_read,*py_readParam;
  unsigned int bytes_read,buf_p;

  unsigned int order;
  unsigned int qgram;

  bit32 shiftMask;

  int totalCounts;

  bit32 streamHistory;
  unsigned int streamCount;
  struct __BGdataCPP *CP;
} matrix_bgObject;

static void
bg_dealloc(matrix_bgObject* self) 
{
  if(self->bgSample)
    Py_DECREF(self->bgSample);
  if(self->strBuf)
    Py_DECREF(self->strBuf);
  if(self->py_read)
    Py_DECREF(self->py_read);
  if(self->py_readParam)
    Py_DECREF(self->py_readParam);

  delete self->CP;

  self->ob_type->tp_free((PyObject*)self);
}


char bg_getNextChar(matrix_bgObject *self)
{
  char ret=0;

  if(self->buf_p>=self->bytes_read) {
    if(self->strBuf) {
      Py_DECREF(self->strBuf);
    }
    self->strBuf=PyObject_CallObject(self->py_read,self->py_readParam);
    self->bytes_read=PyObject_Size(self->strBuf);
    self->buf_p=0;
  }

  if(self->buf_p<self->bytes_read) {
    ret=PyString_AsString(self->strBuf)[self->buf_p++];
  }
  return ret;
}


char bg_getChar(matrix_bgObject *self,int i)
{
  char ret=0,*ret_str=NULL;
  PyObject *py_ret=NULL;

  py_ret=PySequence_GetItem(self->bgSample,i);
  //py_ret=PySequence_GetSlice(self->bgSample,i,min(i+2,(int)self->sampleLen));
  //py_ret=PyString_FromString("A");

  if(!py_ret) {
    cout<<"Null item from sequence position "<<i<<endl;
    return 0;
  }
  ret_str=PyString_AsString(py_ret);
  ret=ret_str[0];

  Py_DECREF(py_ret);

  return ret;

}


bit32 addNucleotideToGram(bit32 gram,char nucleotide, bit32 shiftMask)
{
  
  switch(nucleotide) {
  case 'A':
    gram=((gram<<2)|nucl_A)&shiftMask;
    break;
  case 'C':
    gram=((gram<<2)|nucl_C)&shiftMask;
    break;
  case 'G':
    gram=((gram<<2)|nucl_G)&shiftMask;
    break;
  case 'T':
    gram=((gram<<2)|nucl_T)&shiftMask;
    break;
  default:
    gram=shiftMask+1;
    break;
  }

  return gram;

}


void bitCodeToStr(char *gram,int len,bit32 code)
{
  for(int j=0;j<len;j++) {
    switch(code&3) {
    case nucl_A:
      gram[len-j-1]='A';
      break;
    case nucl_C:
      gram[len-j-1]='C';
      break;
    case nucl_G:
      gram[len-j-1]='G';
      break;
    case nucl_T:
      gram[len-j-1]='T';
      break;
    }
    code=code>>2;
  }

}



static PyObject*
bg_countGrams(matrix_bgObject *self)
{
  char prevChr=0;
  bit32 bitInd=0;
  
  int gramFeed=self->qgram;

  unsigned long fileStartPos=py_fileLikeTell(self->bgSample);

  for(char nucl=bg_getNextChar(self);nucl!=0;nucl=bg_getNextChar(self)) {
    //nucl=bg_getNextChar(self);
    //if(nucl==0) {
    //  break;
    //}
    //cout<<(int)nucl<<endl;
    switch(nucl) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
      bitInd=addNucleotideToGram(bitInd,nucl,self->shiftMask);
      if(gramFeed>0)
	gramFeed--;

      break;
    case '\n':
      continue;

    case '>':
      if(prevChr==0 || prevChr=='\n') {  //If fasta title, 
	for(char chr=bg_getNextChar(self); chr!='\n' && chr!=0;chr=bg_getNextChar(self)) {    // Skip line
	}
      }
    case 'N': 
    case 'X':
    default:
      gramFeed=self->qgram;
      continue;
    }
    if(gramFeed==0) {
      //cout<<bitInd;
      self->totalCounts++;
      self->CP->counts[bitInd]++;
      self->CP->shortCounts[bitInd>>2]++;
    } 

    prevChr=nucl;

  }

  py_fileLikeSeek(self->bgSample,fileStartPos);

  return (PyObject*)self;
}


extern "C" int
bg_init(matrix_bgObject *self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = {"bgSample","order",NULL};

  if(self->bgSample)
    Py_DECREF(self->bgSample);

  self->order=4;

  if (! PyArg_ParseTupleAndKeywords(args, kwds, "O|i", kwlist,
				    &self->bgSample,&self->order))
    return -1; 


  if(PySequence_Check(self->bgSample)) {   //Given a tupple of grams
    int size=PySequence_Length(self->bgSample);
    double db_qgram=log((double)size)/log(4.0);

    self->qgram=(int)db_qgram;
    self->order=self->qgram-1;
    if(fabs(self->qgram-db_qgram)>0.1) {
      PyErr_SetString(PyExc_ValueError,"Malformed gram count sequence.");
      return -1;
    }

  } else {  // Given file to read the grams
    Py_INCREF(self->bgSample);
    if(self->order>15 || self->order<1) {
      printf("Too high or non positive order (%d). Can only handle up to order 15\n",self->order);
      PyErr_SetString(PyExc_ValueError,"Too high or non positive order. Can only handle up to order 15.");
      return -1;
    }
    self->qgram=self->order+1;

  }
  self->shiftMask=(1<<(self->qgram*2))-1;
  self->CP=new struct __BGdataCPP;
  if(!self->CP) {
    PyErr_NoMemory();
    return -1;
  }


  self->CP->counts.resize(self->shiftMask+1,0);
  self->CP->shortCounts.resize((self->shiftMask>>2)+1,0);




  if(PySequence_Check(self->bgSample)) {   //Given a tupple of grams
    int size=PySequence_Length(self->bgSample);
    for(int i=0;i<size;i++) {
      int value=PyInt_AsLong(PySequence_GetItem(self->bgSample,i));
      self->CP->counts[i]=value;
      self->totalCounts+=value;
      self->CP->shortCounts[i>>2]+=value;
    }
    Py_DECREF(self->bgSample);
    Py_INCREF(Py_None);
    self->bgSample=Py_None;


  } else {


    self->strBuf=NULL;  //Initialize buffering
    self->py_read=PyObject_GetAttrString(self->bgSample,"read");
    if(!self->py_read) {
      PyErr_SetString(PyExc_ValueError,"bgSample must be a file like object.");
      return 0;
    }
    self->py_readParam=Py_BuildValue("(l)",SEQ_BUFFER_SIZE);
    self->buf_p=0;
    self->bytes_read=0;


    bg_countGrams(self);
  }
  return 0;
}

static PyObject*
bg_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  matrix_bgObject *self;

  self = (matrix_bgObject *)type->tp_alloc(type, 0);
  if(self != NULL) {
    self->order=4;
    self->streamHistory=0;
    self->streamCount=0;
    self->shiftMask=0;
    self->totalCounts=0;
    self->qgram=0;
    self->CP=NULL;
    Py_INCREF(Py_None);
    self->bgSample=Py_None;
    cout<<self->order<<endl;
  } else {
    PyErr_NoMemory();
    return NULL;
  }


  if(!bg_init(self,args,kwds)) {
    Py_DECREF(self);
    self=NULL;
  }

  return (PyObject*)self;
}


static PyObject*
bg_stringProb(matrix_bgObject *self, PyObject *args)
{
  PyObject *py_ret;
  double ret;
  unsigned long qgramCnt,contextCnt;
  char *str;
  bit32 gram=0;


  if (!PyArg_ParseTuple(args, "s", &str)) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if(strlen(str)!=self->qgram) {
    PyErr_SetString(PyExc_ValueError,"String must be of length order + 1");
    return 0;
  }

  for(unsigned int i=0;i<self->qgram;i++) {
    gram=addNucleotideToGram(gram,str[i],self->shiftMask);
    if(gram>self->shiftMask) {
      PyErr_SetString(PyExc_ValueError,"String must contain only ACGT");
      return 0;
    }
  }
  qgramCnt=self->CP->counts[gram];
  contextCnt=self->CP->shortCounts[gram>>2];
  ret=((double)qgramCnt)/contextCnt;

  //printf("p=%d/%d=%g\n",qgramCnt,contextCnt,ret);
  
  py_ret=PyFloat_FromDouble(ret);

  return py_ret;
  
}

void bit32toStr(char *bitStr,bit32 val)
{
  for(int i=0;i<32;i++) {
    //cout<<i<<endl;
    if(val&(1<<i)) {
      bitStr[31-i]='1';
    } else {
      bitStr[31-i]='0';
    }
  }
  bitStr[32]=0;
}

double logPnextInStream(matrix_bgObject *self, char nucl)
{   // On success, return value is non positive.
  double ret;
  int qgramCnt,contextCnt;
  static double const ln2=log(2.0);

  self->streamHistory=addNucleotideToGram(self->streamHistory,nucl,self->shiftMask);
  self->streamCount++;


  if(self->streamHistory>self->shiftMask) {
    ret=2.0;
  } else {
    if(self->streamCount>self->order) {
      qgramCnt=self->CP->counts[self->streamHistory];
      contextCnt=self->CP->shortCounts[(self->streamHistory)>>2];;
    } else {



      // Compute the start of the stream cleanly (although slowly)
      qgramCnt=0;
      contextCnt=0;

      bit32 startMask=(1<<(2*(self->qgram - self->streamCount)))-1;
      int startShift=2*self->streamCount;

      /*
      printf("Stream item %d (%c)\n",self->streamCount,nucl);

      char bitStr[33];
      bit32toStr(bitStr,self->shiftMask);
      cout<<"Shift Mask:"<<bitStr<<endl;

      bit32toStr(bitStr,startMask<<startShift);
      cout<<"Start Mask:"<<bitStr<<endl;

      bit32toStr(bitStr,self->streamHistory);
      cout<<"History:   "<<bitStr<<endl;


      char gramStr[self->qgram+1];
      gramStr[self->qgram]=0;
      bitCodeToStr(gramStr,self->qgram,self->streamHistory);
      cout<<"Stream:    "<<gramStr<<endl;
*/

      assert(((startMask<<startShift)&self->streamHistory)==0);




      for(unsigned int i=0;i<= startMask;i++) {
	bit32 gram=(i<<startShift)|self->streamHistory;
	/*
	bitCodeToStr(gramStr,self->qgram,gram);
	bit32toStr(bitStr,gram);
	cout<<gramStr<<" "<<bitStr<<endl;
	*/
	qgramCnt+=self->CP->counts[gram];
	contextCnt+=self->CP->shortCounts[gram>>2];
      }
    }
    ret=( log((double)qgramCnt+0.25)-log((double)contextCnt+1.0) )/ln2;
  }
  return ret;

}

static PyObject*
bg_logPnextInStream(matrix_bgObject *self, PyObject *args)
{
  PyObject *py_ret;
  double ret;
  char *str;


  if (!PyArg_ParseTuple(args, "s", &str)) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  ret=logPnextInStream(self,str[0]);

  if(ret>1.0) {
    PyErr_SetString(PyExc_ValueError,"String must contain only ACGT");
    return 0;
  }



  py_ret=PyFloat_FromDouble(ret);

  return py_ret;
}

static PyObject*
resetStream(matrix_bgObject *self)
{
  self->streamHistory=0;
  self->streamCount=0;

  Py_INCREF(Py_None);
  return Py_None;
}


static PyObject*
giveGramVector(matrix_bgObject *self)
{
  PyObject *gramVec=PyTuple_New(self->shiftMask+1);
  PyObject *itemInt;

  if(!self->CP) {
    PyErr_SetString(PyExc_AttributeError,"Object not initialized.");
    return 0;    
  }

  for(unsigned int i=0;i<=self->shiftMask;i++) {
    itemInt=PyInt_FromLong(self->CP->counts[i]);
    PyTuple_SetItem(gramVec,i,itemInt);
  }

  Py_INCREF(gramVec);
  return gramVec;
}


int startContext(matrix_bgObject *self, char *str)
{

  for(unsigned int i=0;i<self->order;i++) {
    self->streamHistory=addNucleotideToGram(self->streamHistory,str[i],self->shiftMask);
    self->streamCount++;
    if(self->streamHistory>self->shiftMask) {
      return 0;
    }
  }
  return 1;
}



static PyObject*
bg_startContext(matrix_bgObject *self, PyObject *args)
{
  char *str;


  if (!PyArg_ParseTuple(args, "s", &str)) {
    return 0;
  }
  

  if(!startContext(self,str)) {
    PyErr_SetString(PyExc_ValueError,"String must contain only ACGT");
    return 0;
  }

  Py_INCREF(Py_None);
  return Py_None;
  
}



static PyObject*
bg_giveCounts(matrix_bgObject *self)
{
  PyObject *ret=PyDict_New();
    
  char *gram=(char*)malloc(self->qgram+1);


  gram[self->qgram]=0;

  if(!ret) {
    PyErr_NoMemory();
    return 0;
  }
  
  for(unsigned int i=0;i<=self->shiftMask;i++) {
    bitCodeToStr(gram,self->qgram,i);
    PyMapping_SetItemString(ret,gram,PyInt_FromLong(self->CP->counts[i]));
  }

  free(gram);

  return ret;
  
  
}





static PyMethodDef bg_methods[] = {
    {"stringProb", (PyCFunction)bg_stringProb, METH_VARARGS,
     "Return the probability of the last character, given the 'order' first characters."},
    {"logPnextProb", (PyCFunction)bg_logPnextInStream, METH_VARARGS,
     "Return the log2 probability of the next character in stream."},
    {"startContext", (PyCFunction)bg_startContext, METH_VARARGS,
     "'order' first characters of the stream whose probability will be assessed. Only first 'order' characters matter."},
    {"giveGramCounts", (PyCFunction)bg_giveCounts, METH_NOARGS,
     "Returns the gram counts in dictionary."},
    {"resetStream", (PyCFunction)resetStream, METH_NOARGS,
     "Reset markov stream."},
    {"giveGramVector", (PyCFunction)giveGramVector, METH_NOARGS,
     "Gives tupple containing count of grams in order."},

  {NULL}
};

static PyMemberDef bg_members[] = {
    {"totalCounts",T_INT, offsetof(matrix_bgObject,totalCounts), 0,
     "Number of observed q grams."},
    {"order",T_INT, offsetof(matrix_bgObject,order), 0,
     "Order of the model."},
    {"bgSample",T_OBJECT_EX, offsetof(matrix_bgObject, bgSample), 0,
     "Background sample string."},
    {NULL}
};

static PyTypeObject matrix_bgType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "matrix.BackGround",             /*tp_name*/
    sizeof(matrix_bgObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)bg_dealloc,                         /*tp_dealloc*/
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
    "Matrix background object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    bg_methods,             /* tp_methods */
    bg_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)bg_init,      /* tp_init */
    0,                         /* tp_alloc */
    bg_new,                 /* tp_new */
};



static PyObject *
matrix_computeBG(PyObject *self, PyObject *args)
{

  char *seq,prevChr=0;
  int order=3;
  bit32 bitInd=0;
  vector<unsigned long> *counts;

  if (!PyArg_ParseTuple(args, "s|i", &seq,&order)){
    Py_INCREF(Py_None);
    return Py_None;
  }
  cout<<"order: "<<order<<endl;
  int const seqLen=strlen(seq);
  int const qgram=order+1;

  if(order>15) {
    printf("Too high order. Can only handle up to order 15\n");
    Py_INCREF(Py_None);
    return Py_None;
  }
  bit32 const shiftMask=(1<<(qgram*2))-1;


  counts=new vector<unsigned long>;

  counts->resize(shiftMask+1,0);


  for (int i=0; i<seqLen; i++) {
    if(i>0)
      prevChr=seq[i-1];
    switch(seq[i]) { 
    case 'A':
    case 'C':
    case 'G':
    case 'T':
      bitInd=addNucleotideToGram(bitInd,seq[i],shiftMask);
      break;
/*		       
    case 'A':
      bitInd=((bitInd<<2)|nucl_A)&shiftMask;
      break;
    case 'C':
      bitInd=((bitInd<<2)|nucl_C)&shiftMask;
      break;
    case 'G':
      bitInd=((bitInd<<2)|nucl_G)&shiftMask;
      break;
    case 'T':
      bitInd=((bitInd<<2)|nucl_T)&shiftMask;
      break;*/


    case 'N': 
    case 'X':
	  cout<<"."<<flush;
	  i+=qgram;
	  continue;
    case '>':
      if(prevChr==0 || prevChr=='\n') {  //If fasta title, 
	i++;
	while(seq[i]!='\n') {    // Skip line
	  i++;
	}
      }
      continue;
    default:
      continue;
    }
    if(i>=(qgram-1)) {
      //cout<<bitInd;
      (*counts)[bitInd]++;
    } else { 
      //cout<<"Too close"<<endl;
    }
  }
  
#define DEBUG_OUTPUT 1
#ifdef DEBUG_OUTPUT

  cout<<"q: "<<qgram<<endl;
  cout<<"grams: "<<counts->size()<<endl;
  char *gram=(char*)malloc(qgram+1);
  gram[qgram]=0;

  for(unsigned int i=0;i<=shiftMask;i++) {
    bitInd=(bit32)i;
    for(int j=0;j<qgram;j++) {
      switch(bitInd&3) {
      case nucl_A:
	gram[qgram-j-1]='A';
	break;
      case nucl_C:
	gram[qgram-j-1]='C';
	break;
      case nucl_G:
	gram[qgram-j-1]='G';
	break;
      case nucl_T:
	gram[qgram-j-1]='T';
	break;
      }
      bitInd=bitInd>>2;
    }
    cout<<i<<" "<<gram<<" "<<(*counts)[i]<<endl;
  }


#endif

  free(gram);
  delete counts;

  Py_INCREF(Py_None);
  return Py_None;
}


unsigned long py_fileLikeTell(PyObject *py_file)
{  //Execute the python tell function

  unsigned long fileStartPos=0;  // Get start position information
  PyObject  *py_ret;

  py_ret=PyObject_CallMethod(py_file,"tell",NULL);

  if(PyInt_Check(py_ret)){
    //cout<<"It's an int!"<<endl;
    fileStartPos=(unsigned long)PyInt_AsLong(py_ret);

  } else if(PyLong_Check(py_ret)){
    //cout<<"It's a long!"<<endl;
    fileStartPos=PyLong_AsUnsignedLong(py_ret);
  }

  //Py_DECREF(py_callable);
  Py_DECREF(py_ret);

  return fileStartPos;
}


int py_fileLikeSeek(PyObject *py_file, unsigned long pos)
{  //Execute the python seek function
  PyObject *py_input,*py_ret;

  py_input=PyLong_FromUnsignedLong(pos);

  py_ret=PyObject_CallMethod(py_file,"seek","(O)",py_input);


  if(py_ret==NULL) {
    return 0;
  }
  Py_DECREF(py_ret);



  Py_DECREF(py_input);

  return 1;
}


//Returns a map from index to score of possible TFBS.
//The arguments are matrix, sequence and bound.
static PyObject *
matrix_getTFBS(PyObject *self, PyObject *args)
{
  char *Seq=NULL; //= PyString_AsString (PyTuple_GetItem(args, 1));
  //int seq_length;//=PyString_Size (PyTuple_GetItem(args, 1));
  double bound;//=PyFloat_AsDouble(PyTuple_GetItem(args, 2));
  PyObject *py_infile,*py_matrix;
  int i, j;
  int mat_length;
  int nucleotide=-1, compl_nucleotide=-1,loop_status;

  int bytes_read,buf_p;
  int const loop_continue=1,loop_break=2,loop_OK=0;



#ifdef TIME_TFBS
  clock_t before,after;

  // Start timing
  before=clock();
#endif


  if (!PyArg_ParseTuple(args, "OOd" ,&py_matrix,&py_infile,&bound)){
    Py_INCREF(Py_None);
    return Py_None;
  }

  
  //vector<vector<double> > M= parse(PyTuple_GetItem(args, 0));
  vector<vector<double> > M= parse(py_matrix);
  if(M.size()==0) {
    Py_INCREF(Py_None);
    return Py_None;
  }

  mat_length= M[0].size();




  // histories store the interim results of the matching
  // one for each of the 2 complement strands
  deque<double> history(mat_length,0.0);
  deque<double> compl_history(mat_length,0.0);

  // iterators to work with the histories
  deque<double>::reverse_iterator iter;
  deque<double>::iterator compl_iter;
  
  PyObject* ret;
  //float entry;
  //char* strand="+";

  if(!(ret=PyDict_New())) 
    {
      PyErr_NoMemory();
      return 0;
    }
  


  PyObject *py_read=NULL;
  
  
  unsigned long fileStartPos=py_fileLikeTell(py_infile);


  py_read=PyObject_GetAttrString(py_infile,"read");


  PyObject *py_readParam=NULL;
  PyObject *py_strBuf=NULL;

  py_readParam=Py_BuildValue("(l)",SEQ_BUFFER_SIZE);


  //py_strBuf=PyObject_CallObject(py_read,py_readParam);

  //cout<<PyString_AsString(py_strBuf)<<endl;
  //bytes_read=strlen(Seq);
  //cout<<bytes_read;
  //bytes_read=fread(Seq,SEQ_BUFFER_SIZE-1,1,cinfile);
  //buf_p=0;
  loop_status=loop_OK;


  bytes_read=1;
  buf_p=1;
  
   
  for (i=0; bytes_read>0; i++,buf_p++) {
    if(buf_p>=bytes_read) {
      if(py_strBuf)
	Py_DECREF(py_strBuf);
      
      py_strBuf=PyObject_CallObject(py_read,py_readParam);
      Seq=PyString_AsString(py_strBuf);
      bytes_read=PyObject_Size(py_strBuf);
      //bytes_read=fread(Seq,SEQ_BUFFER_SIZE-1,1,cinfile);
      Seq[bytes_read]=0;
      buf_p=0;
      //cout<<"Read "<<bytes_read<<" more bytes"<<endl;

      cout<<"."<<flush;
      if(bytes_read==0) {
	break;
      }
    }

    switch(Seq[buf_p])
      {
      case '\n':
      case ' ':
	i--;
	loop_status=loop_continue;
	break;
      case '>':
	cout<<"Encountered unexpectedly an another sequence!"<<endl;
	PyDict_SetItem(ret,PyString_FromString("NEXT_SEQ"),PyLong_FromUnsignedLong(py_fileLikeTell(py_infile)));
	loop_status=loop_break;
	break;

      case 'A':
	//nucleotide stores in which line of the matrix is to search
	nucleotide=0;        
	compl_nucleotide=3;
	break;
      case 'C':
	nucleotide=1;
	compl_nucleotide=2;
	break;
      case 'G':
	nucleotide=2;
	compl_nucleotide=1;
	break;
      case 'T':
	nucleotide=3;
	compl_nucleotide=0;
	break;
	default:
	  //cout<<"Wrong letter in Sequence! Reading it like 'N'"<<endl;
	case 'N': 
	case 'X':
	  //cout<<"."<<flush;
	  //histories are used like queues
	  history.pop_front();
	  history.push_back(0);
	  compl_history.pop_front();
	  compl_history.push_back(0);

	  // in case of 'N' a specific value is added to the histories
	  // to make sure, that it becomes no hit
	  iter=history.rbegin();
	  compl_iter=compl_history.begin();
	  for(j=0; j<mat_length; j++)
	    {
	      *iter+=-10;              //what do I add best to the histories?
	      ++iter;  
	      
	      *compl_iter+=-10;        //what do I add best to the histories?
	      ++compl_iter;
	    }
	  
	  continue;
	  break;
	}


      


      if(loop_status==loop_break) {
	break;
      } else if(loop_status==loop_continue) {
	loop_status=loop_OK;
	continue;
      }


      // Add the score from nucleotide seq[i] to the lists.
      iter=history.rbegin();
      compl_iter=compl_history.begin();
      for(j=0; j<mat_length; j++)
	{
	  *iter+= M[nucleotide][j];
	  ++iter; 

	  //*compl_iter+=M[compl_nucleotide][mat_length-j-1];
	  *compl_iter+=M[compl_nucleotide][j];
	  ++compl_iter;
	}

      /*
      iter=history.begin();
      cout<<"+"<<Seq[i];
      for(j=0; j<mat_length; j++)
	{
	  printf(" %d %g ",i+2-j-mat_length,*iter);
	  ++iter; 
	}
      cout<<endl;
      iter=compl_history.begin();
      cout<<"-"<<Seq[i];
      for(j=0; j<mat_length; j++)
	{
	  printf(" %d %g ",i+2-j-mat_length,*iter);
	  ++iter; 
	}
      cout<<endl;
      */
      // The front element has the score for match with last nucleotide i.

      int pos=i-mat_length+2;
      if(pos>0 && history.front()>bound) {
	addMatch(ret,pos,'+',history.front());
      //printf("+hit: %d+: %g\n%d-: %g\n\n",pos,history.front(),pos,compl_history.front());
      }
      if(pos>0 && compl_history.front()>bound) {
	//printf("-hit: %d+: %g\n%d-: %g\n\n",pos,history.front(),pos,compl_history.front());
	addMatch(ret,pos,'-',compl_history.front());
      }

      
      //histories are used like queues
      history.pop_front();
      history.push_back(0);
      compl_history.pop_front();
      compl_history.push_back(0);


      /*
      if (history.front()<compl_history.front())
	{
	  entry= compl_history.front();
	  strand="-";
	}
      else
	{
	  entry= history.front();
	  strand="+";
	}

      if(entry>bound)
	{
	  int pos=i-mat_length+1+1;
	  if(pos<=0)
	    continue;
	  PyObject* pair=PyTuple_New(2);
	  PyTuple_SetItem(pair, 0, PyFloat_FromDouble(entry));
	  PyTuple_SetItem(pair, 1, PyString_FromString(strand));
	  //int help=(strand=="+" ?
	  //	    i:
	  //	    i-mat_length+1);
	  PyDict_SetItem(ret, PyInt_FromLong(i-mat_length+1+1), pair);
	}
      */


    }
      



  // Reset file position
  if(!py_fileLikeSeek(py_infile,fileStartPos)) {
    return NULL;
  }

  Py_DECREF(py_readParam);
  Py_DECREF(py_read);

#ifdef TIME_TFBS
  // End timing
  after=clock();

  cout<<"CPU secs: "
      <<((after-before)*1.0/CLOCKS_PER_SEC<<endl;

#endif

  if(py_strBuf)
    Py_DECREF(py_strBuf);

  return ret;
}

//Returns a map from index to score of possible TFBS.
//The arguments are matrix, sequence and bound.
static PyObject *
matrix_getTFBSwithBG(PyObject *self, PyObject *args)
{
  char *Seq=NULL; //= PyString_AsString (PyTuple_GetItem(args, 1));
  //int seq_length;//=PyString_Size (PyTuple_GetItem(args, 1));
  double bound;//=PyFloat_AsDouble(PyTuple_GetItem(args, 2));
  PyObject *py_infile,*py_matrix;
  int i, j;
  int mat_length;
  int nucleotide=-1, compl_nucleotide=-1,loop_status;

  int bytes_read,buf_p;
  int const loop_continue=1,loop_break=2,loop_OK=0;
  matrix_bgObject *bg=NULL;


#ifdef TIME_TFBS
  clock_t before,after;

  // Start timing
  before=clock();
#endif


  if (!PyArg_ParseTuple(args, "OOd|O" ,&py_matrix,&py_infile,&bound,&bg)){
    //Py_INCREF(Py_None);
    //return Py_None;
    return NULL;
  }

  if((PyObject*)bg==Py_None) {
    bg=NULL;
  }

  vector<vector<double> > M= parse(py_matrix);
  if(M.size()==0) {
    PyErr_SetString(PyExc_ValueError,"Malformed matrix");
    return 0;
  }


  mat_length= M[0].size();



  // histories store the interim results of the matching
  // one for each of the 2 complement strands
  deque<double> history(mat_length,0.0);
  deque<double> compl_history(mat_length,0.0);
  deque<double> bg_history(mat_length,0.0);

  // iterators to work with the histories
  deque<double>::reverse_iterator iter;
  deque<double>::iterator compl_iter;
  deque<double>::reverse_iterator bg_iter;
  
  PyObject* ret;

  if(!(ret=PyDict_New())) 
    {
      PyErr_NoMemory();
      return 0;
    }
  


  PyObject *py_read=NULL;
  
  
  unsigned long fileStartPos=py_fileLikeTell(py_infile);


  py_read=PyObject_GetAttrString(py_infile,"read");


  PyObject *py_readParam=NULL;
  PyObject *py_strBuf=NULL;

  py_readParam=Py_BuildValue("(l)",SEQ_BUFFER_SIZE);


  loop_status=loop_OK;


  bytes_read=1;
  buf_p=1;
  
   
  for (i=0; bytes_read>0; i++,buf_p++) {
    if(buf_p>=bytes_read) {
      if(py_strBuf)
	Py_DECREF(py_strBuf);
      
      py_strBuf=PyObject_CallObject(py_read,py_readParam);
      Seq=PyString_AsString(py_strBuf);
      bytes_read=PyObject_Size(py_strBuf);
      //bytes_read=fread(Seq,SEQ_BUFFER_SIZE-1,1,cinfile);
      Seq[bytes_read]=0;
      buf_p=0;
      //cout<<"Read "<<bytes_read<<" more bytes"<<endl;

      cout<<"."<<flush;
      if(bytes_read==0) {
	break;
      }
    }

    switch(toupper(Seq[buf_p]))
      {
      case '\n':
      case ' ':
	i--;
	loop_status=loop_continue;
	break;
      case '>':
	cerr<<"Encountered unexpectedly an another sequence!"<<endl;
	PyDict_SetItem(ret,PyString_FromString("NEXT_SEQ"),PyLong_FromUnsignedLong(py_fileLikeTell(py_infile)));
	loop_status=loop_break;
	break;

      case 'A':
	//nucleotide stores in which line of the matrix is to search
	nucleotide=0;        
	compl_nucleotide=3;
	break;
      case 'C':
	nucleotide=1;
	compl_nucleotide=2;
	break;
      case 'G':
	nucleotide=2;
	compl_nucleotide=1;
	break;
      case 'T':
	nucleotide=3;
	compl_nucleotide=0;
	break;
	default:
	  cerr<<"Wrong letter in Sequence! Reading it like 'N'"<<endl;
	case 'N': 
	case 'X':
	  //cout<<"."<<flush;
	  //histories are used like queues
	  history.pop_front();
	  history.push_back(0);
	  compl_history.pop_front();
	  compl_history.push_back(0);
	  bg_history.pop_front();
	  bg_history.push_back(100000.0);

	  // in case of 'N' a specific value is added to the histories
	  // to make sure, that it becomes no hit
	  iter=history.rbegin();
	  compl_iter=compl_history.begin();
	  bg_iter=bg_history.rbegin();

	  for(j=0; j<mat_length; j++)
	    {
	      *iter+=-10;              //what do I add best to the histories?
	      ++iter;  
	      
	      *compl_iter+=-10;        //what do I add best to the histories?
	      ++compl_iter;
	    }
	  
	  continue;
	  break;
	}


      


      if(loop_status==loop_break) {
	break;
      } else if(loop_status==loop_continue) {
	loop_status=loop_OK;
	continue;
      }


      double bgP=0.0;
      if(bg) {
	bgP=logPnextInStream(bg,Seq[buf_p]);
      }
      // Add the score from nucleotide seq[i] to the lists.
      iter=history.rbegin();
      compl_iter=compl_history.begin();
      bg_iter=bg_history.rbegin();

      for(j=0; j<mat_length; j++)
	{
	  *iter+= M[nucleotide][j];
	  ++iter; 

	  *bg_iter+=bgP;
	  ++bg_iter;

	  //*compl_iter+=M[compl_nucleotide][mat_length-j-1];
	  *compl_iter+=M[compl_nucleotide][j];
	  ++compl_iter;
	}

      // The front element has the score for match with last nucleotide i.

      int pos=i-mat_length+2;

      double WatsonScore=history.front();
      if(bg) {
	WatsonScore-=bg_history.front();

	assert(bg_history.front()<0.0);
	assert(history.front()<0.0);
      }
      //printf("watson: %g\n",WatsonScore);

      if(pos>0 && WatsonScore>bound) {
	addMatch(ret,pos,'+',WatsonScore);
      }

      double CrickScore=compl_history.front();
      if(bg) {
	CrickScore-=bg_history.front();
	assert(compl_history.front()<0.0);
	assert(bg_history.front()<0.0);
      }
      //printf("crick: %g\n",CrickScore);

      if(pos>0 && CrickScore>bound) {
	addMatch(ret,pos,'-',CrickScore);
      }

      
      //histories are used like queues
      history.pop_front();
      history.push_back(0);

      compl_history.pop_front();
      compl_history.push_back(0);

      bg_history.pop_front();
      bg_history.push_back(0);

    }
      



  // Reset file position
  if(!py_fileLikeSeek(py_infile,fileStartPos)) {
    return NULL;
  }

  Py_DECREF(py_readParam);
  Py_DECREF(py_read);

#ifdef TIME_TFBS
  // End timing
  after=clock();

  cout<<"CPU secs: "
      <<((after-before)*1.0/CLOCKS_PER_SEC<<endl;

#endif

  if(py_strBuf)
    Py_DECREF(py_strBuf);

  return ret;
}

//################################################################################

//################################################################################
//################################################################################


struct TFBSscan {
  PyObject *py_matrix;  // Pointer to the matrix itself
  double bound;   // Cutoff
  int mat_length;   //Matrix length
  vector<vector<double> > M;   // Parsed matrix for easy access
  deque<double> history;
  deque<double> compl_history;
  deque<double> bg_history;

//   // iterators to work with the histories
//   deque<double>::reverse_iterator iter;
//   deque<double>::iterator compl_iter;
//   deque<double>::reverse_iterator bg_iter;

};


struct TFBSscan *parseMatricies(int *count,PyObject *mats,PyObject *cutoffs,double cutoff)
{
  *count=PySequence_Length(mats);
  struct TFBSscan *ret=new  struct TFBSscan [*count];
  
  if(!ret) {
    PyErr_SetString(PyExc_MemoryError,"Out of memory.");
    return 0;
  }

  for(int i=0;i<*count;i++) {
    assert(PySequence_Check(mats));
    ret[i].py_matrix=PySequence_GetItem(mats,i);
    //assert(PySequence_Check(ret[i].py_matrix));

    if(cutoffs) {
      ret[i].bound=PyFloat_AsDouble(PyNumber_Float(PySequence_GetItem(cutoffs,i)));
    } else {
      ret[i].bound=cutoff;
    }
    PyObject *matListList=PyObject_GetAttrString(ret[i].py_matrix,"M_weight");
    if(!matListList) {
      return 0;
    }
    ret[i].M=parse(matListList);
    if(ret[i].M.empty()) {
      return 0;
    }
    if(ret[i].M.size()==0) {
      PyErr_SetString(PyExc_ValueError,"Malformed matrix");
      return 0;
    }
    ret[i].mat_length=ret[i].M[0].size();
    ret[i].history.resize(ret[i].mat_length,0.0);
    ret[i].compl_history.resize(ret[i].mat_length,0.0);
    ret[i].bg_history.resize(ret[i].mat_length,0.0);

  }



  return ret;
}


//Returns a map from matrix to index to score of possible TFBS.
//The arguments are matrix, sequence and bound.
static PyObject *
matrix_getAllTFBSwithBG(PyObject *self, PyObject *args)
{
  char *Seq=NULL; 
  PyObject *py_infile,*py_cutoff,*py_matlist;
  double cutoff;
  struct TFBSscan *Mat;
  int matrixCount=-1,seq_i;

  int nucleotide=-1, compl_nucleotide=-1,loop_status;

  int bytes_read,buf_p;
  int const loop_continue=1,loop_break=2,loop_OK=0;
  matrix_bgObject *bg=NULL;


#ifdef TIME_TFBS
  clock_t before,after;

  // Start timing
  before=clock();
#endif


  if (!PyArg_ParseTuple(args, "OOO|O" ,&py_matlist,&py_infile,&py_cutoff,&bg)){
    return NULL;
  }

  // Check background
  if((PyObject*)bg==Py_None) {
    bg=NULL;
  }


  // Check cutoff
  if(PyNumber_Check(py_cutoff)) {
    cutoff=PyFloat_AsDouble(PyNumber_Float(py_cutoff));
    py_cutoff=NULL;
  } else if( PySequence_Check(py_cutoff) && 
	     PySequence_Length(py_cutoff)==PySequence_Length(py_matlist) ) {
    cutoff=0.0;
  } else {
    PyErr_SetString(PyExc_ValueError,"Wrong number of cutoffs/matrices.");
    return 0;
  }

  // Use subroutine to parse matrices
  Mat=parseMatricies(&matrixCount,py_matlist,py_cutoff,cutoff);

  if(!Mat) {
    return 0;
  }

  /*
  printf("matricies: %d cutoffs: %d cutoff: %g \n",
	 PySequence_Length(py_matlist),(py_cutoff?PySequence_Length(py_cutoff):-1),
	 cutoff);

  */
  PyObject* ret;

  if(!(ret=PyDict_New())) 
    {
      PyErr_NoMemory();
      return 0;
    }
  


  PyObject *py_read=NULL;
  
  
  unsigned long fileStartPos=py_fileLikeTell(py_infile);


  py_read=PyObject_GetAttrString(py_infile,"read");


  PyObject *py_readParam=NULL;
  PyObject *py_strBuf=NULL;

  py_readParam=Py_BuildValue("(l)",SEQ_BUFFER_SIZE);


  loop_status=loop_OK;


  bytes_read=1;
  buf_p=1;
  
   
  for (seq_i=0; bytes_read>0; seq_i++,buf_p++) {
    if(buf_p>=bytes_read) {
      //PyObject_Print(PyObject_Repr(py_infile),stdout,0);

      if(py_strBuf)
	Py_DECREF(py_strBuf);
      
      py_strBuf=PyObject_CallObject(py_read,py_readParam);
      Seq=PyString_AsString(py_strBuf);
      bytes_read=PyObject_Size(py_strBuf);

      Seq[bytes_read]=0;
      buf_p=0;
      //cout<<"Read "<<bytes_read<<" more bytes"<<endl;

      cout<<"."<<flush;
      if(bytes_read==0) {
	break;
      }
    }

    switch(toupper(Seq[buf_p]))
      {
      case '\n':
      case ' ':
	seq_i--;
	loop_status=loop_continue;
	break;
      case '>':
	cerr<<"Encountered unexpectedly an another sequence!"<<endl;
	PyDict_SetItem(ret,PyString_FromString("NEXT_SEQ"),PyLong_FromUnsignedLong(py_fileLikeTell(py_infile)));
	loop_status=loop_break;
	break;

      case 'A':
	//nucleotide stores in which line of the matrix is to search
	nucleotide=0;        
	compl_nucleotide=3;
	break;
      case 'C':
	nucleotide=1;
	compl_nucleotide=2;
	break;
      case 'G':
	nucleotide=2;
	compl_nucleotide=1;
	break;
      case 'T':
	nucleotide=3;
	compl_nucleotide=0;
	break;
      default:
	cerr<<"Wrong letter in Sequence! Reading it like 'N'"<<endl;
      case 'N': 
      case 'X':
	//histories are used like queues
	for(int i=0;i<matrixCount;i++) {
	  Mat[i].history.pop_front();
	  Mat[i].history.push_back(-1000.0);
	  Mat[i].compl_history.pop_front();
	  Mat[i].compl_history.push_back(-1000.0);
	  Mat[i].bg_history.pop_front();
	  Mat[i].bg_history.push_back(0.0);

	  // in case of 'N' a specific value is added to the histories
	  // to make sure, that it becomes no hit
	  deque<double>::reverse_iterator iter=Mat[i].history.rbegin();
	  deque<double>::iterator compl_iter=Mat[i].compl_history.begin();
	  deque<double>::reverse_iterator bg_iter=Mat[i].bg_history.rbegin();
	  
	  for(int j=0; j<Mat[i].mat_length; j++)
	    {
	      *iter+=-1000.0;              //what do I add best to the histories?
	      ++iter;  
	      
	      *compl_iter+=-1000.0;        //what do I add best to the histories?
	      ++compl_iter;
	    }
	}
	continue;
	break;
      }


      


    if(loop_status==loop_break) {
      break;
    } else if(loop_status==loop_continue) {
      loop_status=loop_OK;
      continue;
    }
    

    double bgP=0.0;
    if(bg) {
      bgP=logPnextInStream(bg,Seq[buf_p]);
    }
    
    for(int i=0;i<matrixCount;i++) {
      // Add the score from nucleotide seq[seq_i] to the lists.
      deque<double>::reverse_iterator iter=Mat[i].history.rbegin();
      deque<double>::iterator compl_iter=Mat[i].compl_history.begin();
      deque<double>::reverse_iterator bg_iter=Mat[i].bg_history.rbegin();
      
      for(int j=0; j<Mat[i].mat_length; j++)
	{
	  *iter+= Mat[i].M[nucleotide][j];
	  ++iter; 
	  
	  *bg_iter+=bgP;
	  ++bg_iter;
	  
	  //*compl_iter+=M[compl_nucleotide][mat_length-j-1];
	  *compl_iter+=Mat[i].M[compl_nucleotide][j];
	  ++compl_iter;
	}
      
      // The front element has the score for match with last nucleotide i.
      
      int pos=seq_i-Mat[i].mat_length+2;
      
      double WatsonScore=Mat[i].history.front();
      if(bg) {
	WatsonScore-=Mat[i].bg_history.front();
	if(Mat[i].bg_history.front()>0.0) printf("bgP: %g\n",Mat[i].bg_history.front());
	assert(Mat[i].bg_history.front()<0.0);
	assert(Mat[i].history.front()<0.0);
      }
      //printf("watson: %g\n",WatsonScore);
      
      if(pos>0 && WatsonScore>Mat[i].bound) {
	addMatchWithKey(ret,Mat[i].py_matrix,pos,'+',WatsonScore);
      } else {
      }
      
      double CrickScore=Mat[i].compl_history.front();
      if(bg) {
	CrickScore-=Mat[i].bg_history.front();
	assert(Mat[i].compl_history.front()<0.0);
	assert(Mat[i].bg_history.front()<0.0);
      }
      //printf("crick: %g\n",CrickScore);
      
      if(pos>0 && CrickScore>Mat[i].bound) {
	addMatchWithKey(ret,Mat[i].py_matrix,pos,'-',CrickScore);
      }
      
      
      //printf("%g > %g %g\n",Mat[i].bound,WatsonScore,CrickScore);
      

      
      //histories are used like queues
      Mat[i].history.pop_front();
      Mat[i].history.push_back(0.0);
      
      Mat[i].compl_history.pop_front();
      Mat[i].compl_history.push_back(0.0);
      
      Mat[i].bg_history.pop_front();
      Mat[i].bg_history.push_back(0.0);
      
    }
  }


  // Reset file position
  if(!py_fileLikeSeek(py_infile,fileStartPos)) {
    return NULL;
  }

  Py_DECREF(py_readParam);
  Py_DECREF(py_read);

#ifdef TIME_TFBS
  // End timing
  after=clock();
  cout<<"CPU secs: "
      <<((after-before)*1.0/CLOCKS_PER_SEC<<endl;

#endif

  delete [] Mat;

  if(py_strBuf)
    Py_DECREF(py_strBuf);

  return ret;
}



static PyMethodDef matrixMethods[] = {
  {"draw",  matrix_draw, METH_VARARGS,
   "Draws a matrix"},
  {"computeBG",  matrix_computeBG, METH_VARARGS,
   "Kind of computes higher order Background: (seq,order)"},
  {"match", matrix_match, METH_VARARGS,
   "matches a Matrix on a Sequence.\nArguments: Matrix (As List of Lists), Sequence (As String)\nOutput is a sequence of scores. Score[i] is the maximum score (Watson or Crick strand) of BS ending at character i"},
  {"getTFBS", matrix_getTFBS, METH_VARARGS,
   "Returns a map from index to score of possible TFBS.\nThe arguments are matrix, sequence and bound. *depricated* Use withBG instead"},
  {"getTFBSwithBg", matrix_getTFBSwithBG, METH_VARARGS,
   "Returns a map from index to score of possible TFBS.\nThe arguments are matrix, sequence, bound and background object."},
  {"getAllTFBSwithBg", matrix_getAllTFBSwithBG, METH_VARARGS,
   "Returns a map from matrix to index to score of possible TFBS.\nThe arguments are a list of matricies, the sequence, list/float of cutoffs and a background object."},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


//PyMODINIT_FUNC
extern "C"
void
initmatrix(void)
{
  PyObject *m=NULL;

  matrix_bgType.tp_new=PyType_GenericNew;
  if(PyType_Ready(&matrix_bgType)<0)
    return;

  m=Py_InitModule("eellib.matrix", matrixMethods);

  if(m==NULL)
    return;

  Py_INCREF(&matrix_bgType);
  PyModule_AddObject(m,"BackGround",(PyObject*) &matrix_bgType);
}
