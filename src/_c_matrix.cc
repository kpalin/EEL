#include <Python.h>
#include "structmember.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <list>
#include <queue>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <algorithm>

using namespace std;

#include "matrix.h"
#include "debugTools.h"


/*
 * $Log$
 * Revision 1.22  2008/05/21 11:49:08  jazkorho
 * Added comments.
 *
 * Revision 1.21  2008/05/19 08:14:24  jazkorho
 * Rewrote TFBS search code. Scanning with zero-order bg is much faster now.
 *
 * Revision 1.20  2008/02/29 09:10:09  kpalin
 * Report all snps
 *
 * Revision 1.19.2.1  2008/01/21 12:51:02  kpalin
 * Now report all, even weak, SNPs hitting binding sites.
 *
 * Revision 1.19  2007/02/08 09:22:54  kpalin
 * p-value thresholding now giving impossibly high threshold for motifs
 * with too likely highest scoring match.
 *
 * Revision 1.18  2006/12/08 09:49:56  kpalin
 * Fixed a seqfault and added a source file from Pasi rastas which got
 * included in _c_matrix.cc
 *
 * Revision 1.17  2006/11/13 12:38:02  kpalin
 * Added code for p-value threshold computation. The threshold
 * computation part is from Pasi Rastas.
 *
 * Revision 1.16  2006/08/10 11:39:16  kpalin
 * Port to 64bit. Changed custom bit32 type to uint32_t.
 *
 * Revision 1.15  2005/11/25 12:13:49  kpalin
 * Removed direct stdout output from C.
 *
 * Revision 1.14  2005/07/07 09:24:10  kpalin
 * Fixed some compilation problems with Visual C++
 *
 * Revision 1.13  2005/07/05 11:22:01  kpalin
 * A cap for number of SNP:s under observation. Avoids memory problems.
 *
 * Revision 1.12  2005/05/19 07:49:42  kpalin
 * Merged Waterman-Eggert style suboptimal alignments and
 * SNP matching.
 *
 * Revision 1.11.2.1  2005/05/03 08:22:52  kpalin
 * Fixed crashes due to iterator invalidation for erase in deque.
 *
 * Revision 1.11  2005/03/22 13:17:22  kpalin
 * Merged some fixes surfacing from testing the public version.
 *
 * Revision 1.10  2005/03/22 12:29:50  kpalin
 * Flush the markov background context.
 *
 * Revision 1.9  2005/03/03 09:01:45  kpalin
 * Presumably working with proper output.
 *
 * Revision 1.8  2005/03/02 13:36:58  kpalin
 * Just a little bit of testing. Slow but (again ;) working ;) version.
 *
 * Revision 1.7  2005/02/28 08:40:12  kpalin
 * Windows porting.
 *
 * Revision 1.6  2005/02/25 09:27:57  kpalin
 * Little, just little, less inefficient SNP scanner.
 *
 */

unsigned long py_fileLikeTell(PyObject *py_file);
int py_fileLikeSeek(PyObject *py_file, unsigned long pos);


//################################################################################
// Misc. helper methods 

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


char *getAllels(char IUPAC)
{
  char *allels=NULL;
  switch(IUPAC) {
    // Order ACGT
  case 'R':
    allels="AG";
    break;
  case 'Y':
    allels="CT";
    break;
  case 'M':
    allels="AC";
    break;
  case 'K':
    allels="GT";
    break;
  case 'S':
    allels="CG";
    break;
  case 'W':
    allels="AT";
    break;
  }
  //printf("%c=>%s\n",IUPAC,allels);

  return allels;
}



static PyObject *
matrix_draw(PyObject *self, PyObject *args)
{
  drawvector(parse(PyTuple_GetItem(args, 0)));

  Py_INCREF(Py_None);
  return Py_None;
}



void addMatch(PyObject *dict,int const pos,char const strand,double const score,PyObject *snps,double const altScore=0.0)
{
  assert(snps!=NULL);
  assert(dict!=NULL);
  printDebug("refcount(snps of size %d)=%d",PyTuple_Size(snps),getRefCount(snps));
  assert(PyTuple_Check(snps));
  assert(PyTuple_Size(snps)==0 || getRefCount(snps)==1);
  PyObject *hitKey=Py_BuildValue("(icO)",pos,strand,snps);
  assert(getRefCount(hitKey)==1);

  printDebug("refcount(snps of size %d)=%d",PyTuple_Size(snps),getRefCount(snps));
  assert(PyTuple_Check(snps));
  assert(PyTuple_Size(snps)==0 || getRefCount(snps)==1);

  if(PyErr_Occurred()!=NULL) {
#ifndef NDEBUG
    printf("Python error! Probably ran out of memory!");
#endif
    return;
  }

  assert(hitKey!=NULL);
  assert(PyTuple_Check(hitKey));
  assert(PyDict_Check(dict));

  PyObject *scoreT=Py_BuildValue("(dd)",score,altScore);
  assert(getRefCount(scoreT)==1);
  PyDict_SetItem(dict, hitKey,scoreT);
  Py_DECREF(scoreT);
  Py_DECREF(hitKey);
  assert(getRefCount(hitKey)==1);
  assert(getRefCount(scoreT)==1);

}


void addMatchWithKey(PyObject *dict,PyObject *key,int const pos,char const strand,double const score,PyObject *snps,double const altScore)
{
  PyObject *subDict;

  assert(getRefCount(dict)==1);
  if(!PyMapping_HasKey(dict,key)) {
    subDict=PyDict_New();
    printDebug("refcount(subDict)=%d",getRefCount(subDict));
    assert(getRefCount(subDict)==1);
    PyDict_SetItem(dict,key,subDict);
    Py_DECREF(subDict);
    printDebug("refcount(subDict)=%d",getRefCount(subDict));
    assert(getRefCount(subDict)==1);

  } else {
    subDict=PyDict_GetItem(dict,key);
    printDebug("refcount(subDict)=%d",getRefCount(subDict));
    assert(getRefCount(subDict)==1);
  }
  if(PyErr_Occurred()!=NULL) {
    return;
  }
  assert(snps);
  assert(subDict);
  addMatch(subDict,pos,strand,score,snps,altScore);
}


//################################################################################
// BG object definitions and methods 

uint32_t const nucl_A=0,nucl_C=1,nucl_G=2,nucl_T=3;

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
    Py_XDECREF(self->strBuf);
    self->strBuf=PyObject_CallObject(self->py_read,self->py_readParam);
    self->bytes_read=PyObject_Size(self->strBuf);
    self->buf_p=0;
  }

  if(self->buf_p<self->bytes_read) {
    ret=PyString_AsString(self->strBuf)[self->buf_p++];
  }
  return ret;
}


// char bg_getChar(matrix_bgObject *self,int i)
// {
//   char ret=0,*ret_str=NULL;
//   PyObject *py_ret=NULL;

//   py_ret=PySequence_GetItem(self->bgSample,i);
//   //py_ret=PySequence_GetSlice(self->bgSample,i,min(i+2,(int)self->sampleLen));
//   //py_ret=PyString_FromString("A");

//   if(!py_ret) {
//     cout<<"Null item from sequence position "<<i<<endl;
//     return 0;
//   }
//   ret_str=PyString_AsString(py_ret);
//   ret=ret_str[0];

//   Py_DECREF(py_ret);

//   return ret;

// }


uint32_t addNucleotideToGram(uint32_t gram,char nucleotide, uint32_t shiftMask)
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


void bitCodeToStr(char *gram,int len,uint32_t code)
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
  uint32_t bitInd=0;
  
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


int PyGramCount_Check(PyObject *self)
{
  return PyTuple_Check(self)||PyList_Check(self);
}

unsigned int posround(double val)
{ 
  int trunc_val=(int)val;

  if((val-trunc_val)>0.5) {
    trunc_val++;
  }
  return trunc_val;
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


  if(PyGramCount_Check(self->bgSample)) {   //Given a tupple of grams
    int size=PySequence_Length(self->bgSample);
    double db_qgram=log((double)size)/log(4.0);

    self->qgram=(int)posround(db_qgram);
    self->order=self->qgram-1;
    if(fabs(self->qgram-db_qgram)>0.1) {
      char errStr[256];
      sprintf(errStr,"Malformed gram count sequence len=%d dbqgram=%g qgram=%d.",size,db_qgram,self->qgram);
      PyErr_SetString(PyExc_ValueError,errStr);
      return -1;
    }

  } else {  // Given file to read the grams
    Py_INCREF(self->bgSample);
    if(self->order>15 || self->order<1) {
      //printf("Too high or non positive order (%d). Can only handle up to order 15\n",self->order);
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




  if(PyGramCount_Check(self->bgSample)) {   //Given a tupple of grams
    int size=PySequence_Length(self->bgSample);
    for(int i=0;i<size;i++) {
      PyObject *pyValue=PySequence_GetItem(self->bgSample,i);
      assert(pyValue);
      int value=PyInt_AsLong(pyValue);
      Py_DECREF(pyValue);
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
#ifdef DEBUG_OUTPUT
    cout<<self->order<<endl;
#endif
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
  uint32_t gram=0;


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

void uint32_ttoStr(char *bitStr,uint32_t val)
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
    ret=2000.0;
    self->streamCount=0; // Reset the context if got 'N'
  } else {
    if(self->streamCount>self->order) {
      qgramCnt=self->CP->counts[self->streamHistory];
      contextCnt=self->CP->shortCounts[(self->streamHistory)>>2];;
    } else {



      // Compute the start of the stream cleanly (although slowly)
      qgramCnt=0;
      contextCnt=0;

      uint32_t startMask=(1<<(2*(self->qgram - self->streamCount)))-1;
      int startShift=2*self->streamCount;

      /*
      printf("Stream item %d (%c)\n",self->streamCount,nucl);

      char bitStr[33];
      uint32_ttoStr(bitStr,self->shiftMask);
      cout<<"Shift Mask:"<<bitStr<<endl;

      uint32_ttoStr(bitStr,startMask<<startShift);
      cout<<"Start Mask:"<<bitStr<<endl;

      uint32_ttoStr(bitStr,self->streamHistory);
      cout<<"History:   "<<bitStr<<endl;


      char gramStr[self->qgram+1];
      gramStr[self->qgram]=0;
      bitCodeToStr(gramStr,self->qgram,self->streamHistory);
      cout<<"Stream:    "<<gramStr<<endl;
*/

      assert(((startMask<<startShift)&self->streamHistory)==0);




      for(unsigned int i=0;i<= startMask;i++) {
	uint32_t gram=(i<<startShift)|self->streamHistory;
	/*
	bitCodeToStr(gramStr,self->qgram,gram);
	uint32_ttoStr(bitStr,gram);
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

double logBestP(matrix_bgObject *self)
{   // Compute the best log P value possible with this background object.

  double ret=0.0;
  int qgramCnt,contextCnt;
  static double const ln2=log(2.0);
  char *nucls="ACGT";
  uint32_t gram;

  for(uint32_t context=0;context<(self->shiftMask>>2);context++) {
    contextCnt=self->CP->shortCounts[context];
    for(int i=0;i<4;i++) {
      gram=addNucleotideToGram(context,nucls[i],self->shiftMask);
      qgramCnt=self->CP->counts[gram];

      ret=min(ret,(log((double)qgramCnt+0.25)-log((double)contextCnt+1.0))/ln2);
      
    }
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
    PyObject *pyCnt=PyInt_FromLong(self->CP->counts[i]);
    PyMapping_SetItemString(ret,gram,pyCnt);
    assert(getRefCount(pyCnt)==2); 
    Py_DECREF(pyCnt);
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
  uint32_t bitInd=0;
  vector<unsigned long> *counts;

  if (!PyArg_ParseTuple(args, "s|i", &seq,&order)){
    Py_INCREF(Py_None);
    return Py_None;
  }
  //cout<<"order: "<<order<<endl;
  int const seqLen=strlen(seq);
  int const qgram=order+1;

  if(order>15) {
    printf("Too high order. Can only handle up to order 15\n");
    Py_INCREF(Py_None);
    return Py_None;
  }
  uint32_t const shiftMask=(1<<(qgram*2))-1;


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
#ifdef DEBUG_OUTPUT
      cout<<"."<<flush;
#endif
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
  
#ifdef DEBUG_OUTPUT

  cout<<"q: "<<qgram<<endl;
  cout<<"grams: "<<counts->size()<<endl;
  char *gram=(char*)malloc(qgram+1);
  gram[qgram]=0;

  for(unsigned int i=0;i<=shiftMask;i++) {
    bitInd=(uint32_t)i;
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

  free(gram);

#endif

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


//################################################################################
// c++ p-value code generously donated by Pasi Rastas under GPL

#ifdef  _GLIBCXX_DEBUG
#include <ext/hash_map>
using __gnu_debug::hash_map;
#else
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#endif
typedef hash_map< int, double > myHashMap;

//this one uses table
int tresholdFromP(const intMatrix &mat, const double &p,const doubleArray bgDist)
{
    int numA = mat.size();
    int n = mat[0].size();

    int maxT = 0;
    int minV = INT_MAX;

    for (int i = 0; i < n; ++i) {
	int max = mat[0][i];
	int min = max;
	for (int j = 1; j < numA; ++j) {
	    int v = mat[j][i];
	    if (max < v)
		max = v;
	    else if (min > v)
		min = v;
	}
	maxT += max;
	if (minV > min)
	    minV = min;
    }
    int R = maxT - n * minV;

    //cout << "Size " << R << "\n" ;

    // use indexes, map, own implementation...
    doubleArray table0(R+1, 0.0);
    doubleArray table1(R+1, 0.0);

    for (int j = 0; j < numA; ++j)
	table0[mat[j][0] - minV] += bgDist[j]; // change this to use own background model

    for (int i = 1; i < n; ++i) {
	for (int j = 0; j < numA; ++j) {
	    int s = mat[j][i] - minV;
	    for (int r = s; r <= R; ++r)
		table1[r] += bgDist[j] * table0[r - s]; // change this to use own background model
	}
	for (int r = 0; r <= R; ++r) {
	    table0[r] = table1[r];
	    table1[r] = 0.0;
	}
    }

    double sum = 0.0;
    int prevNonZero=R;

    for (int r = R; r >= 0; --r) {
      sum += table0[r];
      //cout << "sum = " << sum << "\n" ;
      if (sum > p) {
	//cout << "tol = " << r << "\n" ;
	return max(r+1,prevNonZero-1) + n * minV;
      }
      if(table0[r]>0.0) {
	prevNonZero=r;
      }

    }
#ifdef DEBUG_OUTPUT
    cerr << "Error: No treshold found!";
#endif
    return INT_MAX;
    //cout << "sum = " << sum << "\n";
}

// same as above, uses hashmap instead of arrays...
// somewhat slower...
int tresholdFromP2(const intMatrix &mat, const double &p, const doubleArray bgDist)
{
  int numA = mat.size();  // Size of the alphabet
  int n = mat[0].size();  // Length of the matrix

    int maxT = 0;
    int minV = INT_MAX;

    for (int i = 0; i < n; ++i) {
	int max = mat[0][i];
	int min = max;
	for (int j = 1; j < numA; ++j) {
	    int v = mat[j][i];
	    if (max < v)
		max = v;
	    else if (min > v)
		min = v;
	}
	maxT += max;
	if (minV > min)
	    minV = min;
    }

    myHashMap table0;
    myHashMap table1;

    for (int j = 0; j < numA; ++j)
	table0[ mat[j][0] ] += bgDist[j]; // change this to use own background model

    for (int i = 1; i < n; ++i) {
	//cout << "Size = " << table0.size() << "\n";
	for (int j = 0; j < numA; ++j) {
	    int s = mat[j][i];
	    for (myHashMap::iterator it = table0.begin(); 
		 it != table0.end(); ++it) {
		table1[ it->first + s ] += bgDist[j] * (it->second); // change this to use own background model
	    }
	}
	table0 = table1;
	table1.clear();
    }

    //cout << "maxT = " << maxT << " minV = " << minV << "\n";
    double sum = 0.0;
    int prevNonZero=2*(maxT+1);  // Basically infinite

    for (int r = maxT; r >= n * minV; --r) {
      sum += table0[r];

      //cout << "sum = " << sum << "\n" ;
      if (sum > p) {
	//cout << "tol = " << r << "\n" ;
	return ((r + 1)+prevNonZero)/2;
      }

      if(table0[r]>0.0) {
	prevNonZero=r;
      }

    }
#ifdef DEBUG_OUTPUT
    cerr << "Error: No treshold found!";
#endif
    return INT_MAX;
}

// same as above, uses hashmaps with customisations
// somewhat slower...
int tresholdFromP3(const intMatrix &mat, const double &p, const doubleArray bgDist)
{
  int numA = mat.size();  // Size of the alphabet
  int n = mat[0].size();  // Length of the matrix

    int maxT = 0;
    int minV = INT_MAX;

    for (int i = 0; i < n; ++i) {
	int max = mat[0][i];
	int min = max;
	for (int j = 1; j < numA; ++j) {
	    int v = mat[j][i];
	    if (max < v)
		max = v;
	    else if (min > v)
		min = v;
	}
	maxT += max;
	if (minV > min)
	    minV = min;
    }

    myHashMap table0;
    myHashMap table1;

    for (int j = 0; j < numA; ++j)
	table0[ mat[j][0] ] += bgDist[j]; // change this to use own background model

    for (int i = 1; i < n; ++i) {
	//cout << "Size = " << table0.size() << "\n";
	for (int j = 0; j < numA; ++j) {
	    int s = mat[j][i];
	    for (myHashMap::iterator it = table0.begin(); 
		 it != table0.end(); ++it) {
		table1[ it->first + s ] += bgDist[j] * (it->second); // change this to use own background model
	    }
	}
	table0 = table1;
	table1.clear();
    }

    //cout << "maxT = " << maxT << " minV = " << minV << "\n";
    double sum = 0.0;
    map<int,double> sortedTable(table0.begin(),table0.end());
    int prevNonZero=maxT;

    for (map<int,double>::const_reverse_iterator r=sortedTable.rbegin(); r->first >= n * minV; ++r) {
      sum += r->second;

      //cout << "sum = " << sum << " r = "<<r->first<<","<<r->second<<"\n" ;
      if (sum > p) {
	//cout << "tol = " << r << "\n" ;
	//--r;
	//cout<<" r = "<<r->first<<","<<r->second<<" ehkä = "<<(r->first+1)<<"\n";
	return (r->first+1+prevNonZero)/2;
      }
      if(r->second>0.0) {
	prevNonZero=r->first;
      }

    }
#ifdef DEBUG_OUTPUT
    cerr << "Error: No treshold found!";
#endif
    return INT_MAX;
}

// Wrappers from Kimmo Palin
intMatrix *pyMatrix2IntMatrix(const PyObject  *py_matrix,double *multiplier)
{ // Parse the python array of arrays to type of vector of vectors that is good for Pasis code.
  //Also round the floating point values to integers in range ROUNDING_RANGE
  intMatrix *mat=new intMatrix();

  //const int matLen=PySequence_Fast_GET_SIZE(PySequence_Fast_GET_ITEM(py_matrix,0));
  double minV=DBL_MAX,maxV=-DBL_MAX;


  for(int i=0;
      i<PySequence_Fast_GET_SIZE(py_matrix);i++) {
    PyObject *pyNuc=PySequence_Fast_GET_ITEM(py_matrix,i);
    for(int j=0;
	j<PySequence_Fast_GET_SIZE(pyNuc);j++) {
      double thisVal=PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pyNuc,j));
      minV=min(minV,thisVal);
      maxV=max(maxV,thisVal);
    }
  }
  (*multiplier)=ROUNDING_RANGE/(maxV-minV);

#ifndef NDEBUG
  double err,SSE=0.0,SAE=0.0,maxAbsErr=0.0;
#endif
  for(int i=0;
      i<PySequence_Fast_GET_SIZE(py_matrix);i++) {
    intArray nucl;
    PyObject *pyNuc=PySequence_Fast_GET_ITEM(py_matrix,i);
    for(int j=0;
	j<PySequence_Fast_GET_SIZE(pyNuc);j++) {
      double val=PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pyNuc,j));
      int intVal=(int)round(val*(*multiplier));
      nucl.push_back(intVal);
      //cerr<<val<<"~"<<val*(*multiplier)<<"~"<<(int)round(PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pyNuc,j))*(*multiplier))<<'\t';
      //nucl.push_back((int)round(PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pyNuc,j))*(*multiplier)));
#ifndef NDEBUG
      err=((float)intVal)/(*multiplier)-val;
      SSE+=err*err;
      SAE+=(err>0?err:-err);
      maxAbsErr=max(maxAbsErr,(err>0?err:-err));
#endif      
    }
#ifndef NDEBUG
    cerr<<endl;
#endif
    mat->push_back(nucl);
  }

#ifndef NDEBUG
  int nm=PySequence_Fast_GET_SIZE(py_matrix)*PySequence_Fast_GET_SIZE(PySequence_Fast_GET_ITEM(py_matrix,0));
  cerr<<"Multiplier:"<<*multiplier<<" SSE:"<<SSE<<" SAE:"<<SAE<<" RMSE:"<<sqrt(SSE/nm)<<" MAE:"<<SAE/nm<<" MaxAbsErr:"<<maxAbsErr<<endl;
#endif
  return mat;

}

// Compute the score threshold for given matrix, p-value and 0th order background
static PyObject *
matrix_thresholdFromP(PyObject *self, PyObject *args)
{
  PyObject *bgdist,*py_matrix;
  double pval,multiplier;
  int cutoff;
  double cutoffScore;
  doubleArray bgdistA;

#ifdef TIME_TFBS
  clock_t before,after;

  // Start timing
  before=clock();
#endif


  if (!PyArg_ParseTuple(args, "OdO" ,&py_matrix,&pval,&bgdist)) {
    return NULL;
  }

  if(!PySequence_Check(py_matrix)) {
    PyErr_SetString(PyExc_ValueError,"No matrix list given.");
    return 0;
  }
  if(!PySequence_Check(bgdist)) {
    PyErr_SetString(PyExc_ValueError,"No background distribution given.");
    return 0;
  }
  if(PySequence_Length(bgdist)!=4) {
    PyErr_SetString(PyExc_ValueError,"Wrong size of background distribution.");
    return 0;
  }


  // Cast background distribution fit for Pasis code
  double bgSum=0.0;
  for(int i=0;i<PySequence_Length(bgdist);i++) {
    PyObject *myFloat=PySequence_Fast_GET_ITEM(bgdist,i);

    myFloat=PyNumber_Float(myFloat);
    if(!myFloat) {
      PyErr_SetString(PyExc_ValueError,"Invalid background distribution value.");
      return 0;
    }

    double myDouble=PyFloat_AsDouble(myFloat);
    Py_DECREF(myFloat);

    bgdistA.push_back(myDouble);
    bgSum+=myDouble;
  }
  // Normalize to distribution
  for(int i=0;i<PySequence_Length(bgdist);i++) {
      
    bgdistA[i]/=bgSum;
  }

  intMatrix *myIntMat=pyMatrix2IntMatrix(py_matrix,&multiplier);

  cutoff=tresholdFromP2(*myIntMat,pval,bgdistA);
  delete myIntMat;
  cutoffScore=cutoff/multiplier;
  
#ifdef TIME_TFBS
  // End timing
  after=clock();
  cout<<"CPU secs: "<<((after-before)*1.0/CLOCKS_PER_SEC)<<endl;
		       
#endif
	 
  return PyFloat_FromDouble(cutoffScore);

}



//################################################################################
// TFBS search

// Code-character transformation tables
static const int code_to_char[] = {'A', 'C', 'G', 'T', 'N', 'R', 'Y', 'M', 'K', 'S', 'W'};
static const int snptable[][2] = {{0,2},  // R -> AG
    {1,3},  // Y -> CT
    {0,1},  // M -> AC
    {2,3},  // K -> GT
    {1,2},  // S -> GC
    {0,3}}; // W -> AT


// TFBS search with zero-order background
static PyObject * matrix_getAllTFBSzeroOrderBG(PyObject *self, PyObject *args){

    // Check parameters
    PyObject *py_sequence,*py_cutoff,*py_matlist, *py_bgA, *py_bgC, *py_bgG, *py_bgT;

    if (!PyArg_ParseTuple(args, "OOOOOOO" ,&py_matlist,&py_sequence,&py_cutoff, &py_bgA, &py_bgC, &py_bgG, &py_bgT)){
        return NULL;
    }

    // Get background
    doubleArray bg(4,0);
    bg[0] = PyFloat_AsDouble(PyNumber_Float(py_bgA));
    bg[1] = PyFloat_AsDouble(PyNumber_Float(py_bgC));
    bg[2] = PyFloat_AsDouble(PyNumber_Float(py_bgG));
    bg[3] = PyFloat_AsDouble(PyNumber_Float(py_bgT));

    // We need to know how many matrices there are
    int matrixcount = PySequence_Length(py_matlist);

    // Check cutoff
    double cutoffs[matrixcount];

    if(PyNumber_Check(py_cutoff)) { // Absolute cutoff, same for each matrix
        double cutoff=PyFloat_AsDouble(PyNumber_Float(py_cutoff));
        for (int i = 0; i < matrixcount; ++i){
            cutoffs[i] = cutoff;
        }
    }
    else if( PySequence_Check(py_cutoff) && PySequence_Length(py_cutoff)==PySequence_Length(py_matlist) )
    {
        for (int i = 0; i < matrixcount; ++i){
            cutoffs[i] = PyFloat_AsDouble(PyNumber_Float(PySequence_GetItem(py_cutoff,i)));
        }
    }
    else {
        PyErr_SetString(PyExc_ValueError,"Wrong number of cutoffs/matrices.");
        return NULL;
    }

    PyObject* ret;

    if(!(ret=PyDict_New())) {
        PyErr_NoMemory();
        return NULL;
    }

    // -----------------
    // Read the sequence
    // This supposes that we can read the whole sequence to memory at once

    charArray sequence;
    vector<int> start_pos;
    vector<int> end_pos;
    vector<int> snp_pos;

    unsigned long fileStartPos=py_fileLikeTell(py_sequence);

    // Get the length of the sequence
    int seq_length = PySequence_Length(py_sequence);

    // Length as a Python object
    PyObject *py_readParam = Py_BuildValue("(l)",seq_length);

    // We want to call "read" on py_sequence
    PyObject *py_read = PyObject_GetAttrString(py_sequence,"read");

    // Calling "read" and accessing the resulting string in C++
    PyObject *py_seq_as_string = PyObject_CallObject(py_read,py_readParam);
    char * raw_sequence = PyString_AsString(py_seq_as_string);

    int position = 0;
    int clean = 0;
    try {
        sequence.reserve(seq_length);
    }
    catch (...) {
        PyErr_SetString(PyExc_MemoryError,"Ran out of memory trying to load sequence. This probably means that the sequence was too large for the program to handle.");
        Py_DECREF(py_readParam);
        Py_DECREF(py_read);
        Py_DECREF(py_seq_as_string);
        Py_DECREF(ret);
        return NULL;
    }

    for (int i = 0; i < seq_length; ++i){
        char nuclChr = toupper(raw_sequence[i]);
        int code = -1;

        switch (nuclChr){
            case 'A': code = 0; break;
            case 'C': code = 1; break;
            case 'G': code = 2; break;
            case 'T': code = 3; break;
            case 'R': code = 5; break;
            case 'Y': code = 6; break;
            case 'M': code = 7; break;
            case 'K': code = 8; break;
            case 'S': code = 9; break;
            case 'W': code = 10; break;
            case '\n':
                case ' ': code = -1; break;
                case '>': // Beginning of an another sequence. Shouldn't happen.
                    PyErr_SetString(PyExc_ValueError,"Encountered unexpectedly an another sequence!");
                    {
                        PyObject *pyNextSeq=PyString_FromString("NEXT_SEQ");
                        PyObject *pyFilePos=PyLong_FromUnsignedLong((unsigned long)sequence.size()+1);
                        PyDict_SetItem(ret,pyNextSeq,pyFilePos);
                        Py_DECREF(pyNextSeq);
                        Py_DECREF(pyFilePos);
                    }
                    return ret;
            case 'N':
            case 'X':
                default:  code = 4; break; // Wrong characters are read as N
        }

        if (code != -1){
            sequence.push_back(code);

            if (code < 4 && !clean){
                clean = 1;
                start_pos.push_back(position);
            }
            else if (code >= 4 && clean ){
                clean = 0;
                end_pos.push_back(position);
            }
            
            if (code >= 5){
                snp_pos.push_back(position);
            }
            ++position;
        }

    }
    end_pos.push_back(position);

    Py_DECREF(py_readParam);
    Py_DECREF(py_read);
    Py_DECREF(py_seq_as_string);

    // Reset file position
    if(!py_fileLikeSeek(py_sequence,fileStartPos)) {
        return NULL;
    }
    
    // The actual scanning part!
    // (with trivial algorithm)

    // Process the matrices for scanning
    for (int matrix = 0; matrix < matrixcount; ++matrix){
        PyObject * py_matrix = PySequence_GetItem(py_matlist,matrix);
        PyObject * matListList=PyObject_GetAttrString(py_matrix,"freq");
        if(!matListList) {
            PyErr_SetString(PyExc_ValueError,"Malformed matrix");
            return NULL;
        }
        doubleMatrix freq_matrix = parse(matListList);

        assert(freq_matrix.size() == 4);
        int mat_len = freq_matrix[0].size();
        for (int i = 1; i < 4; ++i){
            assert(freq_matrix[i].size() == mat_len);
        }

        // Matrix as given, for scanning Watson strand
        doubleMatrix pssm;

        for (int i = 0; i < 4;++i){
            doubleArray row;
            for (int j = 0; j < mat_len; ++j){
                row.push_back( (log(freq_matrix[i][j]) - log(bg[i])) / log(2) );
            }
            pssm.push_back(row);
        }

        naiveAlgorithm(sequence, start_pos, end_pos, pssm, cutoffs[matrix], ret, py_matrix, '+');

        // Check positions overlapping SNPs . 
        // Assuming there are relatively few SNPs, this shouldn't take too long
        getHitsWithSNPs(sequence, snp_pos, pssm, bg, cutoffs[matrix], ret, py_matrix, '+');

        // Reverse complement, for Crick strand
        doubleMatrix pssm_inverse;

        for (int i = 0; i < 4;++i){
            doubleArray row;
            for (int j = 0; j < mat_len; ++j){
                row.push_back( (log(freq_matrix[3 - i][mat_len - j - 1]) - log(bg[i])) / log(2) );
            }
            pssm_inverse.push_back(row);
        }

        naiveAlgorithm(sequence, start_pos, end_pos, pssm_inverse, cutoffs[matrix], ret, py_matrix, '-');

        // SNP positions for reverse complement matrix
        getHitsWithSNPs(sequence, snp_pos, pssm_inverse, bg, cutoffs[matrix], ret, py_matrix, '-');
    }

//     // The actual scanning part!
//     // (with MMACF)
// 
//     vector<doubleMatrix> matrices;
//     doubleArray cutoffs_parameter;
//     vector<PyObject*> py_matrices;
//     charArray strands;
// 
//     // Process the matrices for scanning
//     for (int matrix = 0; matrix < matrixcount; ++matrix){
//         PyObject * py_matrix = PySequence_GetItem(py_matlist,matrix);
//         PyObject * matListList=PyObject_GetAttrString(py_matrix,"freq");
//         if(!matListList) {
//             PyErr_SetString(PyExc_ValueError,"Malformed matrix");
//             return NULL;
//         }
//         doubleMatrix freq_matrix = parse(matListList);
// 
//         assert(freq_matrix.size() == 4);
//         int mat_len = freq_matrix[0].size();
//         for (int i = 1; i < 4; ++i){
//             assert(freq_matrix[i].size() == mat_len);
//         }
// 
//         // Matrix as given, for scanning Watson strand
//         doubleMatrix pssm;
// 
//         for (int i = 0; i < 4;++i){
//             doubleArray row;
//             for (int j = 0; j < mat_len; ++j){
//                 row.push_back( (log(freq_matrix[i][j]) - log(bg[i])) / log(2) );
//             }
//             pssm.push_back(row);
//         }
// 
//         matrices.push_back(pssm);
//         cutoffs_parameter.push_back(cutoffs[matrix]);
//         py_matrices.push_back(py_matrix);
//         strands.push_back('+');
// 
//         // Check positions overlapping SNPs first. 
//         // Assuming there are relatively few SNPs, this shouldn't take too long
//         getHitsWithSNPs(sequence, snp_pos, pssm, bg, cutoffs[matrix], ret, py_matrix, '+');
// 
//         // Reverse complement, for Crick strand
//         doubleMatrix pssm_inverse;
// 
//         for (int i = 0; i < 4;++i){
//             doubleArray row;
//             for (int j = 0; j < mat_len; ++j){
//                 row.push_back( (log(freq_matrix[3 - i][mat_len - j - 1]) - log(bg[i])) / log(2) );
//             }
//             pssm_inverse.push_back(row);
//         }
// 
//         matrices.push_back(pssm_inverse);
//         cutoffs_parameter.push_back(cutoffs[matrix]);
//         py_matrices.push_back(py_matrix);
//         strands.push_back('-');
// 
//         // SNP positions for reverse complement matrix
//         getHitsWithSNPs(sequence, snp_pos, pssm_inverse, bg, cutoffs[matrix], ret, py_matrix, '-');
//     }
// 
//     // Scan "clean" sections of sequence with Aho-Corasick filter algorithm
//     multipleMatrixAhoCorasickLookaheadFiltration(sequence, start_pos, end_pos, matrices, bg, cutoffs_parameter, ret, py_matrices, strands);

    return ret;
}


// The naive algorithm
void naiveAlgorithm(const charArray &s, const intArray &start_pos, const intArray &end_pos, const doubleMatrix &p,  const double tol, PyObject *ret_dict, PyObject * py_matrix, const char strand){
    
    //const int n = s.size();
    const int m = p[0].size();

    for (int slice = 0; slice < (int) start_pos.size(); ++slice){

        int start = start_pos[slice];
        int end = end_pos[slice];

        if (end - start + 1 >= m){
            

            for (int i = start; i < end - m + 1; ++i )
            {

                double score = 0;
                for (int j = 0; j < m; ++j)
                {
                    score += p[ s[i + j] ][j];
                }
                if (score >= tol)
                {
                    PyObject *snps=PyTuple_New(0);
                    addMatchWithKey(ret_dict, py_matrix, i + 1, strand, score, snps, score);

                }
            }
        }
    }
}




doubleArray expectedDifferences(const doubleMatrix &mat, const doubleArray &bg)
{
    int numA = mat.size();
    int m = mat[0].size();
    doubleArray ret(m);
    
    for (int i = 0; i < m; ++i)
    {
        double max = DBL_MIN;
        for (int j = 0; j < numA; ++j)
        {
            if (max < mat[j][i])
                max = mat[j][i];
        }
        
        ret[i] = max;
        
        for (int j = 0; j < numA; ++j)
        {
            ret[i] -= bg[j] * mat[j][i];
        }
    }
    
    return ret;
}


// Actual scanning subroutine for zero-order background
// Uses a Aho-Corasick-based filter to find potential matches
// Scans only the parts of the sequence marked as "clean", ie.
// no SNPs or symbol N. 
void multipleMatrixAhoCorasickLookaheadFiltration(const charArray &s, const intArray &start_pos, const intArray &end_pos, const vector<doubleMatrix> &matrices, const doubleArray &bg, const doubleArray &tol, PyObject *ret_dict, vector<PyObject*> py_matrices, const charArray &strands)
{

    const int numA = 4;
    const int q = 8; 

    intArray m(matrices.size(), 0);
    int min_length = INT_MAX;

    for (int i = 0; i < (int) matrices.size(); ++i)
    {
        m[i] = matrices[i][0].size();
        min_length = min(min_length, m[i]);
    }
        
    // Calculate entropies for all matrices
    vector<doubleArray> goodnesses;
    goodnesses.reserve(matrices.size());
    
    for (int i = 0; i < (int)matrices.size(); ++i)
    {
        goodnesses.push_back(expectedDifferences(matrices[i], bg));
    }
    
    // Find best window positions for all matrices
    intArray window_positions;
    window_positions.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            window_positions.push_back(0);
        }
        else
        {
            double current_goodness = 0;
            for (int i = 0; i < q; ++i)
            {
                current_goodness += goodnesses[k][i];
            }
            
            double max_goodness = current_goodness;
            int window_pos = 0;
            
            for (int i = 0; i < m[k] - q; ++i)
            {
                current_goodness -= goodnesses[k][i];
                current_goodness += goodnesses[k][i+q];
                if (current_goodness > max_goodness)
                {
                    max_goodness = current_goodness;
                    window_pos = i+1;
                }
            }
            window_positions.push_back(window_pos);
        }
    }
        
    // Calculate lookahead scores for all matrices
    doubleMatrix T;
    T.reserve(matrices.size());
    
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        doubleArray C(m[k],0);
        for (int j = m[k] - 1; j > 0; --j) {
            double max = DBL_MIN;
            for (int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            C[j - 1] = C[j] + max;
        }
        T.push_back(C);
    }
    
    // Pre-window scores
    doubleArray P;
    P.reserve(matrices.size());
    for (int k = 0; k < (int)matrices.size(); ++k)
    {
        double B = 0;
        for (int j = 0; j < window_positions[k]; ++j)
        {
            double max = DBL_MIN;
            for (int i = 0; i < numA; ++i) {
                if (max < matrices[k][i][j])
                    max = matrices[k][i][j];
            }
            B += max;
        }
        P.push_back(B);
    }
    
    // Arrange matrix indeces not in window by entropy, for use in scanning
    intMatrix orders;
    orders.reserve(matrices.size());
    doubleMatrix L;
    L.reserve(matrices.size());
    
    for (int k = 0; k < (int) matrices.size(); ++k)
    {
        if (q >= m[k])
        {
            intArray temp_int;
            orders.push_back(temp_int);
            doubleArray temp_dbl;
            L.push_back(temp_dbl);
        }
        else
        {
            intArray order(m[k]-q, 0);
            for (int i = 0; i < window_positions[k]; ++i)
            {
                order[i] = i;
            }
            for (int i = window_positions[k]+q; i < m[k]; ++i)
            {
                order[i-q] = i;
            }
            
            compareRows comp;
            comp.goodness = &(goodnesses[k]);
            
            sort(order.begin(), order.end(), comp);
            
            // Scanning lookahead scores
            orders.push_back(order);
            
            doubleArray K(m[k]-q, 0); 
            for (int j = m[k]-q-1; j > 0; --j)
            {
                double max = DBL_MIN;
                for (int i = 0; i < numA; ++i)
                {
                    if (max < matrices[k][i][order[j]])
                        max = matrices[k][i][order[j]];
                }
                K[j - 1] = K[j] + max;
            }
            L.push_back(K);
        }
    }
    
    // Actual automaton construction begins
    
    // We first construct a temporary transition table that uses integers to indicate transitions
    // Pointers cannot be used because vector resizes would mess them up
    vector<ConstructionACStateMulti> tempACMachine;
    queue<ConstructionQueueElementMulti> stateQ;
    
    // Construct states
    
    ConstructionACStateMulti newState;
    for (int i = 0; i < numA; ++i)
    {
        newState.transition[i] = -1;
    }
    
    tempACMachine.push_back(newState);
    
    ConstructionQueueElementMulti qElement;
    qElement.prev = 0;
    qElement.i = 0;
    {
        OutputListElementMulti temp;
        temp.score = 0;
        temp.full = false;
        for (int i = 0; i < (int)matrices.size(); ++i)
        {
            temp.matrix = i;
            qElement.scores.push_back(temp);
        }
    }
    stateQ.push(qElement);
    
    while (!stateQ.empty())
    {
        qElement = stateQ.front();
        stateQ.pop();
        
        if (qElement.i == q)
        {
            tempACMachine[qElement.prev].output.insert(tempACMachine[qElement.prev].output.end(), qElement.scores.begin(), qElement.scores.end());
        }
        else if (qElement.i < q)
        {
            for (int j = 0; j < numA; ++j)
            {
                list<OutputListElementMulti> expand;
                list<OutputListElementMulti> finished;
                for (list<OutputListElementMulti>::iterator x = qElement.scores.begin(); x != qElement.scores.end(); ++x)
                {
                    if (x->score + matrices[x->matrix][j][qElement.i+window_positions[x->matrix]] + T[x->matrix][qElement.i+window_positions[x->matrix]] + P[x->matrix] >= tol[x->matrix])
                    {
                        if (m[x->matrix] == qElement.i + 1)
                        {
                            finished.push_back(*x);
                            finished.back().score += matrices[x->matrix][j][qElement.i];
                            finished.back().full = true;
                        }
                        else 
                        {
                            expand.push_back(*x);
                            expand.back().score += matrices[x->matrix][j][qElement.i + window_positions[x->matrix]];
                        }
                    }
                }
                
                
                if (!expand.empty() || !finished.empty())
                {
                    
                    ConstructionQueueElementMulti newElement;
                    newState.output.clear();
                    newState.output.insert(newState.output.end(), finished.begin(), finished.end());
                    tempACMachine.push_back(newState);
                    tempACMachine[qElement.prev].transition[j] = tempACMachine.size() - 1;
                    if (!expand.empty())
                    {
                        newElement.prev = tempACMachine.size() - 1;
                        newElement.scores.insert(newElement.scores.end(), expand.begin(), expand.end());
                        newElement.i = qElement.i + 1;
                        stateQ.push(newElement);
                    }
                }
            }
        }
    }
    
    for (int i = 0; i < numA; ++i)
    {
        if (tempACMachine[0].transition[i] == -1)
            tempACMachine[0].transition[i] = 0;
    }
    
    // Construct fail function and final output function
    
    int * fail;
    fail = new int [tempACMachine.size()]; // temp array for fail function
    fail[0] = 0; // probably unnecessary
    
    queue<int> failQ;
    for (int i = 0; i < numA; ++i)
    {
        if (tempACMachine[0].transition[i] != 0)
        {
            failQ.push(tempACMachine[0].transition[i]);
            fail[tempACMachine[0].transition[i]] = 0;
        }
    }
    while(!failQ.empty())
    {
        int current = failQ.front();
        failQ.pop();
        for (int i = 0; i < numA; ++i)
        {
            if (tempACMachine[current].transition[i] != -1)
            {
                int next = tempACMachine[current].transition[i];
                failQ.push(next);
                int state = fail[current];
                while (tempACMachine[state].transition[i] == -1)
                    state = fail[state];
                fail[next] = tempACMachine[state].transition[i];
                tempACMachine[next].output.insert(tempACMachine[next].output.end(), tempACMachine[fail[next]].output.begin(), tempACMachine[fail[next]].output.end());
                
            }
        }
    }
    
    // The final AC automaton that will be used in scanning
    // Storing transitions as pointers to next state saves a significant amount of time during scanning, as there's no need to do pointer arithmetics
    FinalACStateMulti * FinalACMachine;
    FinalACMachine = new FinalACStateMulti[tempACMachine.size()];
    
    // Convert the temporary AC automaton to final one using pointers and precompute all transitions
    
    for (int j = 0; j < (int) tempACMachine.size(); ++j)
    {
        for (int i = 0; i < numA; ++i)
        {
            int state = j;
            while (tempACMachine[state].transition[i] == -1)
            {
                state = fail[state];
            }
            FinalACMachine[j].transition[i] = &(FinalACMachine[tempACMachine[state].transition[i]]);
        }
        FinalACMachine[j].output = tempACMachine[j].output;
    }
    
    delete[] fail;
    int number_of_states = (int) tempACMachine.size();
    tempACMachine.clear();
    
    // AC automaton is now ready
    
    
    
    if (number_of_states <= 1)
    {
        delete[] FinalACMachine;
        return;
    }
    
    // Scan the given sequence with the AC automaton
    
    FinalACStateMulti * state;
    
    charArray::const_iterator x;
    
    for (int slice = 0; slice < (int) start_pos.size(); ++slice){

        int start = start_pos[slice];
        int end = end_pos[slice];

        if (end - start + 1 >= min_length){
            
            x = s.begin() + start;
            state  = FinalACMachine;

            for (int i = start; i < end; ++i )
            {
                state = (*state).transition[*x];
                ++x;
                if (!(state->output.empty()))
                {
                    for (list<OutputListElementMulti>::iterator y = state->output.begin(); y != state->output.end(); ++y)
                    {
                        if (y->full == true)
                        {
                            PyObject *snps=PyTuple_New(0);
                            addMatchWithKey(ret_dict,py_matrices[y->matrix],(i-m[y->matrix] + 2),strands[y->matrix],y->score,snps,y->score);

                            continue;
                        }
                        if (i - q - window_positions[y->matrix] + 1 >= start && i + m[y->matrix] - q - window_positions[y->matrix] < end)
                        {
                            double score = y->score;
                            int k = y->matrix;
                            int limit = m[k] - q;
                            int ii = i - q - window_positions[k] + 1;
                            double tolerance = tol[k];
                            intArray::iterator z = orders[k].begin();   
                            for (int j = 0; j < limit  ;++j)
                            {
                                score += matrices[k][s[ii+(*z)]][*z];
                                if (score + L[k][j] < tolerance)
                                    break;
                                ++z;
                            }
                            if (score >= tolerance)
                            {
                                PyObject *snps=PyTuple_New(0);
                                addMatchWithKey(ret_dict,py_matrices[k],(i-q-window_positions[k]+2),strands[k],score,snps,score);
                            }
                        }
                    }
                }
            }
        }
    }
    
    delete[] FinalACMachine;
}


// Finds matches for a single matrix near given SNP positions
// A bit inefficient, but does the job.
void getHitsWithSNPs(const charArray &sequence, const intArray &snp_pos, const doubleMatrix &pssm, const doubleArray &bg, const double cutoff, PyObject *ret_dict, PyObject *py_matrix, const char strand)
{
    if (snp_pos.empty())
        return;

    int mat_len = pssm[0].size();

    vector<int>::const_iterator snp_first = snp_pos.begin();
    vector<int>::const_iterator snp_last = snp_pos.begin();
    int snps_in_window = 1;

    int position = max(0, (*snp_last) - mat_len + 1);
    int end = min((int)sequence.size() - mat_len + 1, (*snp_first)+1);

    double score;
    double minscore;
    double maxscore;

    int code;

    // Iterates over all sequence positions overlapping SNPs

    while (snp_last != snp_pos.end()){

        while (position < end){
            if (position + mat_len - 1 == *(snp_first+1)){ // A new SNP just entered the scanning window
                snps_in_window += 1;
                snp_first += 1;
                end = min((int)sequence.size() - mat_len + 1, (*snp_first)+1);
            }

            score = 0;
            
            if (snps_in_window <= MAX_SNP_COUNT){

                // Gets the maximum score at given position
        
                for (int i = 0; i < mat_len; ++i){
                    code = sequence[position + i];
                    if (code < 4)
                        score += pssm[code][i];
                    else if (code == 4){
                        score = DBL_MIN;  // Fail.
                        break;
                    }
                    else {
                        score += max(pssm[snptable[code-5][0]][i], pssm[snptable[code-5][1]][i]);
                    }
                }
        
                // There is a hit iff the maximum score at the position exceeds the cutoff bound
        
                if (score >= cutoff){
        
                    // This messy bit finds out which alleles may occur at which positions
                    // and constructs the appropriate entry to the return dictionary
                    // Yes, it's messy. I'm sorry.
        
                    maxscore = score;
                    minscore = score;
        
                    charArray allele(snps_in_window, 0);
        
                    // Iterate over possible allele combinations to see which ones score enough
        
                    for (int snpcode = 0; snpcode < (1 << snps_in_window); ++snpcode){
                        int current_snpcode = snpcode;
                        score = 0;
        
                        for (int i = 0; i < mat_len; ++i){
                            code = sequence[position + i];
                            if (code < 4)
                                score += pssm[code][i];
                            else {
                                score += pssm[snptable[code-5][current_snpcode & 1]][i];
        
                                current_snpcode = current_snpcode >> 1;
                            }
                        }
        
                        if (score >= cutoff){
                            minscore = min(score, minscore);
                            current_snpcode = snpcode;
        
        
        
                            for (int i = 0; i < snps_in_window; ++i){
                                code = sequence[*(snp_last+i)];
                                if(allele[i] == 0){
                                    allele[i] = code_to_char[snptable[code-5][current_snpcode & 1]];
                                }
                                else if (allele[i] != code_to_char[snptable[code-5][current_snpcode & 1]]){
                                    allele[i] = 'N'; // If both alleles are possible, mark it as N
                                }
                                current_snpcode = current_snpcode >> 1;
                            }
                        }
                    }
        
                    // If both alleles are possible and they score almost the same,
                    // the SNP is omited from data
                    // LARGE_AFFY_DELTA is declared as preprocessor constant
        
                    int real_size = 0;
        
                    for (int i = 0; i < snps_in_window; ++i){
                        code = sequence[*(snp_last+i)];
                        double scorediff = pssm[snptable[code-5][0]][*(snp_last+i)-position] - pssm[snptable[code-5][1]][*(snp_last+i)-position];
                        if (allele[i] != 'N' || fabs(scorediff) > LARGE_AFFY_DELTA) ++real_size;
                    }
        
                    // This part finally constructs the Python tuple for SNP data
        
                    PyObject *snps = PyTuple_New( real_size );
                    int tuple_pos = 0;
        
                    for (int i = 0; i < snps_in_window; ++i){
                        code = sequence[*(snp_last+i)];
                        double scorediff = pssm[snptable[code-5][0]][*(snp_last+i)-position] - pssm[snptable[code-5][1]][*(snp_last+i)-position];
                        if (allele[i] != 'N' || fabs(scorediff) > LARGE_AFFY_DELTA){
                            PyObject *py_snp = Py_BuildValue("(ccid)", code_to_char[code], allele[i],  *(snp_last+i) - position + 1, scorediff);
        
                            PyTuple_SetItem(snps, tuple_pos, py_snp);
                            ++tuple_pos;
                        }
                    }
        
                    addMatchWithKey(ret_dict, py_matrix, position + 1, strand, maxscore, snps, minscore);
                }
            }

            position += 1;
            if (position == (*snp_last) + 1){ // Last SNP in window just fell out of the window
                snps_in_window -= 1;
                snp_last += 1;
            }

        }

        // No more SNPs in current window, jump to next SNP

        snp_first += 1;
        snp_last = snp_first;
        snps_in_window = 1;

        position = max(0, (*snp_last) - mat_len + 1);
        end = min((int)sequence.size() - mat_len + 1, (*snp_first)+1);
    }
}

// TFBS search with higher-order markov background
static PyObject * matrix_getAllTFBSMarkovBG(PyObject *self, PyObject *args){

    // Check parameters
    PyObject *py_sequence,*py_cutoff,*py_matlist;
    matrix_bgObject *bg; 

    if (!PyArg_ParseTuple(args, "OOOO" ,&py_matlist,&py_sequence,&py_cutoff, &bg)){
        return NULL;
    }

    // We need to know how many matrices there are
    int matrixcount = PySequence_Length(py_matlist);

    // Check cutoff
    double cutoffs[matrixcount];

    if(PyNumber_Check(py_cutoff)) { // Absolute cutoff, same for each matrix
        double cutoff=PyFloat_AsDouble(PyNumber_Float(py_cutoff));
        for (int i = 0; i < matrixcount; ++i){
            cutoffs[i] = cutoff;
        }
    }
    else if( PySequence_Check(py_cutoff) && PySequence_Length(py_cutoff)==PySequence_Length(py_matlist) )
    {
        for (int i = 0; i < matrixcount; ++i){
            cutoffs[i] = PyFloat_AsDouble(PyNumber_Float(PySequence_GetItem(py_cutoff,i)));
        }
    }
    else {
        PyErr_SetString(PyExc_ValueError,"Wrong number of cutoffs/matrices.");
        return NULL;
    }

    PyObject* ret;

    if(!(ret=PyDict_New())) {
        PyErr_NoMemory();
        return NULL;
    }

    // -----------------
    // Read the sequence
    // This supposes that we can read the whole sequence to memory at once
    // Also, we precalculate the probabilities given by the background
    // so we only need to do bg calculations near SNP positions when
    // scanning. This, of course, requires quite much memory.

    charArray sequence;
    doubleArray bgProps;


    unsigned long fileStartPos=py_fileLikeTell(py_sequence);

    // Get the length of the sequence
    int seq_length = PySequence_Length(py_sequence);

    // Length as a Python object
    PyObject *py_readParam = Py_BuildValue("(l)",seq_length);

    // We want to call "read" on py_sequence
    PyObject *py_read = PyObject_GetAttrString(py_sequence,"read");

    // Calling "read" and accessing the resulting string in C++
    PyObject *py_seq_as_string = PyObject_CallObject(py_read,py_readParam);
    char * raw_sequence = PyString_AsString(py_seq_as_string);
    
    try {
        sequence.reserve(seq_length);
        bgProps.reserve(seq_length);
    }
    catch (...) {
        PyErr_SetString(PyExc_MemoryError,"Ran out of memory trying to load sequence. This probably means that the sequence was too large for the program to handle.");
        Py_DECREF(py_readParam);
        Py_DECREF(py_read);
        Py_DECREF(py_seq_as_string);
        Py_DECREF(ret);
        return NULL;
    }
    logPnextInStream(bg, 'N');

    for (int i = 0; i < seq_length; ++i){
        char nuclChr = toupper(raw_sequence[i]);
        int code = -1;

        switch (nuclChr){
            case 'A': code = 0; break;
            case 'C': code = 1; break;
            case 'G': code = 2; break;
            case 'T': code = 3; break;
            case 'R': code = 5; break;
            case 'Y': code = 6; break;
            case 'M': code = 7; break;
            case 'K': code = 8; break;
            case 'S': code = 9; break;
            case 'W': code = 10; break;
            case '\n':
                case ' ': code = -1; break;
                case '>': // Beginning of an another sequence. Shouldn't happen.
                    PyErr_SetString(PyExc_ValueError,"Encountered unexpectedly an another sequence!");
                    {
                        PyObject *pyNextSeq=PyString_FromString("NEXT_SEQ");
                        PyObject *pyFilePos=PyLong_FromUnsignedLong((unsigned long)sequence.size()+1);
                        PyDict_SetItem(ret,pyNextSeq,pyFilePos);
                        Py_DECREF(pyNextSeq);
                        Py_DECREF(pyFilePos);
                    }
                    return ret;
            case 'N':
                case 'X': code = 4; break;
                default:  code = 4; break; // Wrong characters are read as N
        }

        if (code != -1){
            sequence.push_back(code);
            if (code < 4){
                bgProps.push_back(logPnextInStream(bg, nuclChr));
            }
            else {
                bgProps.push_back(logPnextInStream(bg, 'N'));
            }
        }
    }

    Py_DECREF(py_readParam);
    Py_DECREF(py_read);
    Py_DECREF(py_seq_as_string);

    // Reset file position
    if(!py_fileLikeSeek(py_sequence,fileStartPos)) {
        return NULL;
    }

    // The actual scanning part!

    for (int matrix = 0; matrix < matrixcount; ++matrix){
        PyObject * py_matrix = PySequence_GetItem(py_matlist,matrix);
        PyObject * matListList=PyObject_GetAttrString(py_matrix,"M_weight");
        if(!matListList) {
            PyErr_SetString(PyExc_ValueError,"Malformed matrix");
            return NULL;
        }

        // Matrix as given, for scanning Watson strand
        // As BG is higher order, background is not included in weight matrix
        doubleMatrix pssm = parse(matListList);

        assert(pssm.size() == 4);
        int mat_len = pssm[0].size();
        for (int i = 1; i < 4; ++i){
            assert(pssm[i].size() == mat_len);
        }

        // Find all hits for the matrix
        getHitsWithMarkovBG(sequence, pssm, bg, bgProps, cutoffs[matrix], ret, py_matrix, '+');

        // Reverse complement, for Crick strand
        doubleMatrix pssm_inverse;

        for (int i = 0; i < 4;++i){
            doubleArray row;
            for (int j = 0; j < mat_len; ++j){
                row.push_back( pssm[3 - i][mat_len - j - 1]);
            }
            pssm_inverse.push_back(row);
        }

        // Find all hits for the matrix
        getHitsWithMarkovBG(sequence, pssm_inverse, bg, bgProps, cutoffs[matrix], ret, py_matrix, '-');
    }

    return ret;
}

// Actual scanning part for markov bg
// Finds hits for single TFBS matrix
void getHitsWithMarkovBG(const charArray &sequence, const doubleMatrix &pssm, matrix_bgObject * bg, const doubleArray &bgProps, const double cutoff, PyObject *ret_dict, PyObject *py_matrix, const char strand){

    const int n = sequence.size();
    const int m = pssm[0].size();

    int snps_in_window = 0;
    double score = 0;
    char allele [MAX_SNP_COUNT];
    double score_diffs [MAX_SNP_COUNT];
    int positions [MAX_SNP_COUNT];


    for (int i = 0; i < m - 1; ++i){
        if (sequence[i] > 4)
            ++snps_in_window;
    }

    for (int i = 0; i < n - m + 1; ++i){
        if (sequence[i + m - 1] > 4)
            ++snps_in_window;

        if (snps_in_window == 0){ // No SNPs
            score = 0;

            // Calculate the score
            for (int j = 0; j < m; ++j){
                if (sequence[i+j] == 4){
                    score = DBL_MIN;
                    break;
                }
                score += pssm[sequence[i+j]][j];
                score -= bgProps[i+j];
            }
            if (score > cutoff){
                PyObject *snps=PyTuple_New(0);
                addMatchWithKey(ret_dict,py_matrix, i+1, strand, score, snps, score);
            }
        }


        else if (snps_in_window <= MAX_SNP_COUNT){ // SNPs 
            double maxscore = DBL_MIN;
            double minscore = DBL_MAX;
            char N_in_window = 0;

    // Resets the alleles after previous SNP
            for (int j = 0; j < snps_in_window; ++j){
                allele[j] = 0;
            }
            int code;

    // Iterate over possible SNP combinations
    // Gets quite messy again.
            for (int snpcode = 0; snpcode < (1 << snps_in_window); ++snpcode){
                int current_snpcode = snpcode;
                score = 0;

        // The BG object is used to give the background score
        // at this position for the given allele combination.
        // As we have to do this for each combination at each
        // position, this takes quite a lot unnecessary time
        // if there are lots of SNPs. Improving this should
        // be a priority.

                logPnextInStream(bg, 'N'); // Resets the stream

        // Position just before the actual window
        // Take alleles into account in calculation of BG prop.
                for (int j = min((int)bg->order, i); j > 0; --j){
                    code = sequence[i-j];
                    if (code <= 4)
                        logPnextInStream(bg, code_to_char[code]);
                    else {
                        logPnextInStream(bg, code_to_char[snptable[code-5][current_snpcode & 1]]);
                        current_snpcode = current_snpcode >> 1;
                    }
                }

        // Calculate the score
                for (int j = 0; j < m; ++j){
                    code = sequence[i + j];
                    if (code < 4){
                        score += pssm[code][j];
                        score -= logPnextInStream(bg, code_to_char[code]);
                    }
                    else if (code > 4) {
                        score += pssm[snptable[code-5][current_snpcode & 1]][j];
                        score -= logPnextInStream(bg, code_to_char[snptable[code-5][current_snpcode & 1]]);
                        current_snpcode = current_snpcode >> 1;
                    }
                    else {
                        N_in_window = 1;
                        break;
                    }
                }
                if (N_in_window){
                    maxscore = DBL_MIN;
                    break;
                }

        // Do bookkeeping stuff if given combination scores well enough
                if (score > cutoff){
                    maxscore = max(score, maxscore);
                    minscore = min(score, minscore);
                    current_snpcode = snpcode;

                    int allele_n = 0;

                    for (int j = -min((int)bg->order, i); j < m; ++j){
                        code = sequence[i + j];
                        if (code > 4){
                            if(allele[allele_n] == 0){
                                allele[allele_n] = code_to_char[snptable[code-5][current_snpcode & 1]];
                            }
                            else if (allele[allele_n] != code_to_char[snptable[code-5][current_snpcode & 1]]){
                                allele[allele_n] = 'N'; // If both alleles are possible, mark it as N
                            }
                            current_snpcode = current_snpcode >> 1;
                            ++allele_n;
                        }
                    }
                }
            }

            if (maxscore > cutoff){ // We actually got a hit - store it to the dictionary

        // If both alleles are possible and they score almost the same,
        // the SNP is omited from data
        // LARGE_AFFY_DELTA is declared as preprocessor constant

                int real_size = 0;
                int k = 0;

                for (int j = -min((int)bg->order, i); j < m; ++j){
                    code = sequence[i+j];
                    if (code > 4){
                        if (j < 0){ // SNP in background only
                            score_diffs[k] = 0;
                        }
                        else{
                            score_diffs[k] = pssm[snptable[code-5][0]][j] - pssm[snptable[code-5][1]][j];
                        }
                        positions[k] = j;
                        if (allele[k] != 'N' || fabs(score_diffs[k]) > LARGE_AFFY_DELTA) ++real_size;
                        ++k;
                    }
                }

        // This part finally constructs the Python tuple for SNP data

                PyObject *snps = PyTuple_New( real_size );
                int tuple_pos = 0;

                for (k = 0; k < snps_in_window; ++k){
                    if (allele[k] != 'N' || fabs(score_diffs[k]) > LARGE_AFFY_DELTA){
                        code = sequence[i+positions[k]];
                        PyObject *py_snp = Py_BuildValue("(ccid)", code_to_char[code], allele[k],  positions[k] + 1, score_diffs[k]);

                        PyTuple_SetItem(snps, tuple_pos, py_snp);
                        ++tuple_pos;
                    }
                }

                addMatchWithKey(ret_dict, py_matrix, i + 1, strand, maxscore, snps, minscore);
            }

        }

        if ((i - (int)(bg->order)) >= 0){
            if (sequence[i - (int)(bg->order)] > 4){
                --snps_in_window;
            }
        }
    }
}

static PyMethodDef matrixMethods[] = {
    {"draw",  matrix_draw, METH_VARARGS,
    "Draws a matrix"},
    {"computeBG",  matrix_computeBG, METH_VARARGS,
    "Kind of computes higher order Background: (seq,order)"},
    {"getAllTFBSzeroOrderBG", matrix_getAllTFBSzeroOrderBG, METH_VARARGS,
    "Returns a map from matrix to index to score of possible TFBS given zero-order background.\nThe arguments are a list of matricies, the sequence, list/float of cutoffs and A/C/G/T frequences as floats."},
    {"getAllTFBSMarkovBG", matrix_getAllTFBSMarkovBG, METH_VARARGS,
    "Returns a map from matrix to index to score of possible TFBS given a markov background\nThe arguments are a list of matricies, the sequence, list/float of cutoffs and a background object."},
    {"thresholdFromP", matrix_thresholdFromP, METH_VARARGS,
    "Returns a score threshold for given matrix corresponding to a p-value\nArguments: a [[matrix]], p-value, [background distribution]."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

//PyMODINIT_FUNC
extern "C"
void
init_c_matrix(void)
{
  PyObject *m=NULL;

  matrix_bgType.tp_new=PyType_GenericNew;
  if(PyType_Ready(&matrix_bgType)<0)
    return;

  m=Py_InitModule("eellib._c_matrix", matrixMethods);

#ifndef NDEBUG
  cerr<<"Loaded matrix.cc"<<endl;
#endif
  if(m==NULL)
    return;

  Py_INCREF(&matrix_bgType);
  PyModule_AddObject(m,"BackGround",(PyObject*) &matrix_bgType);
}
