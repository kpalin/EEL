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
#include <float.h>
#include <algorithm>
using namespace std;

#include "matrix.h"

#ifndef SEQ_BUFFER_SIZE
#define SEQ_BUFFER_SIZE 15000000
#endif

/*
 * $Log$
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



PyObject *SNPdat::buildPySNP(int refPos)
{
  return Py_BuildValue("(ccid)",this->ambig,this->allele,refPos-this->pos,this->scoreDif);
}

PyObject *TFBShit::buildPySNPs()
{
  int size=this->sigGenotype.size();
  int realSize=0;
  PyObject *ret=PyTuple_New(size);

  for(unsigned int i=0;i<this->sigGenotype.size();i++) {
    if(this->sigGenotype[i].allele!='N' || this->sigGenotype[i].scoreDif>LARGE_AFFY_DELTA) {
      PyObject *obj=this->sigGenotype[i].buildPySNP(this->mat->length());
#ifndef NDEBUG
      assert(PyTuple_Check(obj));
      char *al=PyString_AsString(PyTuple_GetItem(obj,1));
      assert(*al=='A'||*al=='C'||*al=='G'||*al=='T'||*al=='N');
#endif
      PyTuple_SetItem(ret,realSize++,obj);
    }
  }
  if(size!=realSize) {
    int rval=_PyTuple_Resize(&ret,realSize);
    assert(rval==0);
  }
  return ret;
}

int SNPdat::diffAllele(SNPdat &other)
{
  return this->pos==other.pos && this->ambig==other.ambig && this->allele!=other.allele;
}
int SNPdat::operator==(SNPdat &other)
{
  return this->pos==other.pos && this->allele==other.allele && this->ambig==other.ambig;
}

bool operator<(const class TFBSscan &t1, const class TFBSscan &t2){
  return t1.length()<t2.length() || (t1.length()==t2.length() && t1.py_matrix<t2.py_matrix);
}

bool lessThan(TFBSscan* t1,TFBSscan* t2){
  return *t1<*t2;
}


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
  PyObject *hitKey=Py_BuildValue("(icO)",pos,strand,snps);

  if(PyErr_Occurred()!=NULL) {
#ifndef NDEBUG
    printf("Python error! Probably ran out of memory!");
#endif
    return;
  }

  assert(hitKey!=NULL);
  assert(PyTuple_Check(hitKey));
  assert(PyDict_Check(dict));
  PyDict_SetItem(dict, hitKey,Py_BuildValue("(dd)",score,altScore));

}


void addMatchWithKey(PyObject *dict,PyObject *key,int const pos,char const strand,double const score,PyObject *snps,double const altScore)
{
  PyObject *subDict;
  if(!PyMapping_HasKey(dict,key)) {
    subDict=PyDict_New();
    PyDict_SetItem(dict,key,subDict);
  } else {
    subDict=PyDict_GetItem(dict,key);
  }
  if(PyErr_Occurred()!=NULL) {
    return;
  }
  assert(snps);
  assert(subDict);
  addMatch(subDict,pos,strand,score,snps,altScore);
}



bit32 const nucl_A=0,nucl_C=1,nucl_G=2,nucl_T=3;

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
    ret=2000.0;
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

double logBestP(matrix_bgObject *self)
{   // Compute the best log P value possible with this background object.

  double ret=0.0;
  int qgramCnt,contextCnt;
  static double const ln2=log(2.0);
  char *nucls="ACGT";
  bit32 gram;

  for(bit32 context=0;context<(self->shiftMask>>2);context++) {
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

//################################################################################
//################################################################################


void TFBSscan::halfHistories()
{
  // Remove the SNPs that we have past and left behind

//   cout<<"halving"<<endl;
#ifndef NDEBUG
  int before=this->history.size();
#endif

  deque<deque<double> >::iterator Iter=this->history.begin();
  deque<deque<double> >::iterator compl_Iter=this->compl_history.begin();
  // Remove every second value
  while(Iter!=this->history.end()) {
    compl_Iter=this->compl_history.erase(compl_Iter);
    compl_Iter++;
    Iter=this->history.erase(Iter);
    Iter++;
  }
  assert((before>>1)==this->history.size());
  assert(this->history.size()>0);
  assert(this->compl_history.size()>0);
}

void printIntBits(int val)
{
   char str[33];
   bit32toStr(str,(bit32)val);
   printf("%s %d\n",str,val);
}

vector<int> snpIndexFromCode(int snpCode,int snpCount)
{
  vector<int> ret;
  int  ind;


  for(int i=0;i<snpCount  ;i+=1) {
    ind=snpCode & 1;
    ret.push_back((i<<1)|ind);
    snpCode=snpCode>>1;
  }
  return ret;
}

vector<double> TFBSscan::WatsonScore() 
{
  vector<double> ret;
  
  for(unsigned int i=0;i<this->history.size();i++) {
    ret.push_back(this->history[i].front());
  }

  return ret;
}


vector<double> TFBSscan::CrickScore() 
{
  vector<double> ret;
  for(unsigned int i=0;i<this->compl_history.size();i++) {
    ret.push_back(this->compl_history[i].front());
  }

  return ret;
}

TFBSscan::TFBSscan(PyObject *mat,double cutoff)
{
  this->py_matrix=mat;
  this->bound=cutoff;

  PyObject *matListList=PyObject_GetAttrString(this->py_matrix,"M_weight");
  if(!matListList) {
    PyErr_SetString(PyExc_ValueError,"Malformed matrix");
    return;
  }
  this->M=parse(matListList);
  if(this->M.empty()||this->M.size()==0) {
    PyErr_SetString(PyExc_ValueError,"Malformed matrix");
    return;
  }
  this->mat_length=this->M[0].size();
  this->history.push_front( deque<double>(this->mat_length,0.0) );
  this->compl_history.push_front(deque<double>(this->mat_length,0.0));
}
  
void TFBSscan::nextChar(char chr)
{

  char *allels=getAllels(chr);


  if(allels==NULL) {
    this->nextACGT(chr);
  } else {
    for(int allel_p=0;allel_p<2;allel_p++) {
      //this->SNPs.push_back(SNPdat(chr,allels[allel_p],0));
      int newStart;
      if(allel_p==0) {
	newStart=this->history.size();

	this->history.insert(this->history.end(),
				this->history.begin(),this->history.end());
	this->compl_history.insert(this->compl_history.end(),
			this->compl_history.begin(),this->compl_history.end());
	this->nextACGT(allels[allel_p],0,newStart);
      } else if(allel_p==1) {
	this->nextACGT(allels[allel_p],newStart);
      } else {
	cout<<"Nobody expects the Spanish Inquisition!"<<endl;
      }
    }
//      printf("added snp buffers: %d\n",this->history.size());
  }


}


void TFBSscan::nextACGT(char chr,int fromCode,int toCode)
{
  //printf("chr: %c fromCode: %d history.size(): %d\n",chr,fromCode,this->history.size());
  if(toCode<=fromCode) {
    toCode=this->history.size();
  }
#ifndef NDEBUG
  if(toCode>8) {
    printf("toCode: %d\n",toCode);
  }
#endif
  for(unsigned int i=fromCode;i<toCode;i++) {
    this->nextACGTsingle(chr,i);
    //cout<<i<<endl;
  }

}

void TFBSscan::nextACGTsingle(char chr,int snpCode)
{
  int nucleotide=-1,compl_nucleotide=-1;

  assert(snpCode>=0 && snpCode<this->history.size());
  //histories are used like queues
  this->history[snpCode].pop_front();
  this->history[snpCode].push_back(0.0);
  
  this->compl_history[snpCode].pop_front();
  this->compl_history[snpCode].push_back(0.0);
  

  switch(chr) {
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
    //cerr<<"Wrong letter in Sequence! Reading it like 'N'"<<endl;
  case 'N': 
  case 'X':
    //histories are used like queues
    this->history[snpCode].pop_front();
    this->history[snpCode].push_back(-1000.0);
    this->compl_history[snpCode].pop_front();
    this->compl_history[snpCode].push_back(-1000.0);

    // in case of 'N' a specific value is added to the histories
    // to make sure, that it becomes no hit
    deque<double>::reverse_iterator iter=this->history[snpCode].rbegin();
    deque<double>::iterator compl_iter=this->compl_history[snpCode].begin();
      
    for(int j=0; j<this->mat_length; j++) {
      *iter+=-1000.0;              //what do I add best to the histories?
      ++iter;  
      
      *compl_iter+=-1000.0;        //what do I add best to the histories?
      ++compl_iter;
    }
    nucleotide=-1;
    compl_nucleotide=-1;
    return;
    break;
  }
  // If we got a proper nucleotide
  // Add the score from nucleotide seq[seq_i] to the lists.
  deque<double>::reverse_iterator iter=this->history[snpCode].rbegin();
  deque<double>::iterator compl_iter=this->compl_history[snpCode].begin();

  for(int j=0; j<this->mat_length; j++) {
    //printf("pos:%d score:%g\n",j,this->M[nucleotide][j]);      
    *iter+= this->M[nucleotide][j];
    ++iter; 
	  
    //*compl_iter+=M[compl_nucleotide][mat_length-j-1];
    *compl_iter+=this->M[compl_nucleotide][j];
    ++compl_iter;
  }

 //  printf("\n%c %d\n",chr,snpCode);
//   for(unsigned int j=0;j<this->history[snpCode].size();j++) {
//     printf("%g,",this->history[snpCode][j]);
//   }
//   printf("\n");
  

}

TFBShit::TFBShit(TFBSscan* mat,unsigned int seqPos,char strand)
{
  this->mat=mat;
  this->pos=seqPos-mat->length()+1; // seqPos is the sequence number of the last nucleotide of the site.
  this->strand=strand;
  this->score=-DBL_MAX;
  this->minScore=DBL_MAX;
}

void TFBShit::addHit(double Score,vector<class SNPdat> &genotype)
{

  // Be nice and give the maximum score if both alleles are accepted.
  this->score=max(this->score,Score);
  this->minScore=min(this->score,Score);


//   printf("%g: ",Score);
//   for(int i=0;i<genotype.size();i++) {
//     printf("%d%c%c ",genotype[i].pos,genotype[i].ambig,genotype[i].allele);
//   }
//   cout<<endl;


  if(this->sigGenotype.size()==0 && genotype.size()>0) {
    this->sigGenotype=genotype;
  } else {
    assert(this->sigGenotype.size()==genotype.size());
    for(unsigned int i=0;i<this->sigGenotype.size();i++) {
      assert(this->sigGenotype[i].pos==genotype[i].pos);
      assert(this->sigGenotype[i].ambig==genotype[i].ambig);
      if(this->sigGenotype[i].allele!=genotype[i].allele) {
// 	printf("N:ing. %d%c%c %d%c%c\n",this->sigGenotype[i].pos,this->sigGenotype[i].ambig,this->sigGenotype[i].allele,genotype[i].pos,genotype[i].ambig,genotype[i].allele);
	this->sigGenotype[i].allele='N';
      }
    }

  
  }

#ifndef NDEBUG
  if(this->sigGenotype.size()>0) {
    int kokoe=this->sigGenotype.size();
    for(unsigned int i=0;i<this->sigGenotype.size();i++) {
      SNPdat apu=this->sigGenotype[i];
      if(!(this->sigGenotype[i].allele=='A' ||
	   this->sigGenotype[i].allele=='C' ||
	   this->sigGenotype[i].allele=='G' ||
	   this->sigGenotype[i].allele=='T' ||
	   this->sigGenotype[i].allele=='N'))
	printf("i=%d size:%d pos,ambig,al %d %d %d\n",i,this->sigGenotype.size(),this->sigGenotype[i].pos,this->sigGenotype[i].ambig,this->sigGenotype[i].allele);
    }
  }

#endif
}


char *SNPdat::alleles()
{
  return getAllels(this->ambig);
}

double TFBSscan::setSNPscoreDif(SNPdat &snp,int crick)
{

  char *allels=snp.alleles();
  int pos=this->length()-snp.pos-1;

  assert(allels);
  assert(snp.pos>=0);


  if(pos>=0) {  // SNPs in the background context does not affect this.
    assert(pos>=0);
    assert(pos<this->length());

    snp.scoreDif=this->matItem(pos,allels[0],crick)-this->matItem(pos,allels[1],crick);
  } else {
    snp.scoreDif=0.0;
  }

  return snp.scoreDif;


}

double TFBSscan::matItem(int i,char nucl,int crick)
{
  int code=TFBSscan::ACGTtoCode(nucl);

  if(crick) {
    code=3-code;
    i=this->length()-1-i;
  }

  assert(i>=0);
  assert(i<this->length());
  assert(code<4 && code>=0);

  return this->M[code][i];
}

int TFBSscan::ACGTtoCode(char nucl)
{
  int nucleotide;

  switch(nucl) {
  case 'A':
    //nucleotide stores in which line of the matrix is to search
    nucleotide=0;        
    break;
  case 'C':
    nucleotide=1;
    break;
  case 'G':
    nucleotide=2;
    break;
  case 'T':
    nucleotide=3;
    break;
  default:
    nucleotide=-1;
  }


  // Complement is 3-nucleotide
  return nucleotide;
}
 
vector<SNPdat> TFBShelper::getSNPs(int snpCode,int matInd,int crick)
{

  // Get SNPs for the matrix matInd.  This might have less SNP matches 
  // than the recorded background region.

  vector<SNPdat> ret;


  if(this->SNPcount()>0) {

    // Correct the snpCode for the larger number of SNPs
    if(matInd>=0 && this->matricies[matInd]->allelCount()<this->allelCount()) {
      snpCode*=this->allelCount()/this->matricies[matInd]->allelCount();
    }

    vector<int> ind=snpIndexFromCode(snpCode,this->SNPcount());
//     if(this->SNPcount()>2) {
//       for(int i=0;i<this->SNPs.size();i++){
// 	printf("%d%c%c\n",this->SNPs[i].pos,this->SNPs[i].ambig,this->SNPs[i].allele);
//       }
//     }
    for(unsigned int i=0;i<ind.size();i++) {
      assert(ind[i]<(int)this->allelCount());
      if(matInd<0 || this->SNPs[ind[i]].pos<(int)(this->matricies[matInd]->length()+this->bgOrder())) {
	// Don't record the SNPs that are not on the matrices region.

	// Make a copy
	SNPdat curSnp=this->SNPs[ind[i]];
	this->matricies[matInd]->setSNPscoreDif(curSnp,crick);
	ret.push_back(curSnp);
      }
    }
  }
  
  return ret;
  
}




vector<TFBShit*> TFBShelper::getMatches()
{

  vector<TFBShit*> ret;
#ifndef NDEBUG
  if(this->SNPcount()>3) {
    printf("SNPs: %d\n",this->SNPcount());
  }
#endif

  for(unsigned int matInd=0;matInd<this->matrixCount();matInd++) {
    vector<double> WatsonScores=this->matricies[matInd]->WatsonScore();
    vector<double> CrickScores=this->matricies[matInd]->CrickScore();
    double scoreBound=this->matricies[matInd]->bound;

    class TFBShit *watson=NULL,*crick=NULL;

    for(unsigned int snpCode=0;snpCode<this->matricies[matInd]->allelCount();snpCode++) {
      assert(WatsonScores.size()==CrickScores.size());
      assert(WatsonScores.size()==this->matricies[matInd]->allelCount());
      double bgProb=this->getBGprob(matInd,snpCode);
      double WatsonScore=WatsonScores[snpCode]-bgProb;
      double CrickScore=CrickScores[snpCode]-bgProb;
      if(WatsonScore>scoreBound) {
	if(!watson) {
	  watson=new TFBShit(this->matricies[matInd],this->seqPos(),'+');
	}
// 	printf("bgprob:%g, watsonscore:%g %d\n",bgProb,WatsonScore,snpCode);
	vector<SNPdat> mySNPs=this->getSNPs(snpCode,matInd);
	watson->addHit(WatsonScore,mySNPs);
      }
      if(CrickScore>scoreBound) {
	if(!crick) {
	  crick=new TFBShit(this->matricies[matInd],this->seqPos(),'-');
	}
	vector<SNPdat> mySNPs=this->getSNPs(snpCode,matInd,1);
	crick->addHit(CrickScore,mySNPs);
      }
    }
    if(watson) {
      ret.push_back(watson);
    }
    if(crick) {
      ret.push_back(crick);
    }
  }
  return ret;
}

double TFBShelper::getBGprob(int matInd,int snpCode)
{

  if(!this->haveBG) {
    return 0.0;
  }
 
  TFBSscan *mat=this->matricies[matInd];
  unsigned int matAlleles=mat->allelCount();

  int pos=mat->length()-1;

  double bgP=-DBL_MAX;
  int alleleFactor=this->allelCount()/matAlleles;


#ifndef NDEBUG
  if(this->allelCount()!=this->bg.size()) {
    printf("allels=%d  bg.size()=%d bgProb.size()=%d\n",this->allelCount(),this->bg.size(),this->probBuffer.size());
  }
  assert(this->allelCount()==this->bg.size());
  double oldBgP=this->probBuffer[snpCode*alleleFactor][pos];


  // Return the most likely background probability i.e. Be conservative.
  for(int sCode=snpCode*alleleFactor;sCode<(snpCode+1)*alleleFactor;sCode++) {
    //printIntBits(sCode);
    bgP=max(bgP,this->probBuffer[sCode][pos]);
    // Since we report the SNPs in the markov context, we should get equal values here. 
    // i.e. this should be innecessary loop.
    if(fabs(oldBgP-bgP)>=1e-5) {
      printf("ERO: %g\n",fabs(oldBgP-bgP));
      vector<int> ind=snpIndexFromCode(snpCode,this->SNPcount());
      for(unsigned int i=0;i<ind.size();i++) {
	printf("%d%c%c\n",this->SNPs[i].pos,this->SNPs[i].ambig,this->SNPs[i].allele);
      }
      assert(fabs(oldBgP-bgP)<1e-5);

    }
  }
#else
  bgP=this->probBuffer[snpCode*alleleFactor][pos];

#endif

  return bgP;
}


vector<double> TFBShelper::getBGprobs(int matInd)
{

  vector<double> P;

  for(unsigned int i=0;i<this->probBuffer.size();i++) {
    P.push_back(this->getBGprob(matInd,i));
  }
  //printf("Giving %g for mat length %d(%d)\n",P,mat->length(),pos);

  return P;
}

TFBShelper::TFBShelper(matrix_bgObject *bgIn,vector<TFBSscan*> &mat) : matricies(mat)
{ 
  this->seqCount=0;

  this->maxLen=0;

  for(unsigned int i=0;i<this->matricies.size();i++) {
#ifndef NDEBUG
    cout<<this->matricies[i]->length()<<endl;
#endif
    this->maxLen=max(this->maxLen,this->matricies[i]->length());
  }


  this->haveBG=(bgIn!=NULL);
  if(this->haveBG) {
#ifndef NDEBUG
    printf("Best log2(p)=%g\n",logBestP(bgIn));
#endif
    this->bg.push_back(*bgIn);  //Copy !!!
    this->maxLen+=this->bgOrder();
    this->probBuffer.push_back(deque<double>(this->maxLen,0.0));
  }

}


void TFBShelper::removeScannerHistories()
{


  int snpId=0;
  unsigned int lookBackDist=0;

  // Longest matrix first
  for(int i=this->matricies.size()-1;i>=0 && snpId<this->SNPs.size();i--) {
    // Loop min(number of matricies, number of SNPs) times.


    lookBackDist=this->matricies[i]->length()+this->bgOrder();
    

    // Furthest SNP first
    while(snpId<this->SNPs.size() && this->SNPs[snpId].pos>(int)lookBackDist) {
      snpId++;
    }
    if(snpId<this->SNPs.size() && this->SNPs[snpId].pos==(int)lookBackDist) {
      this->matricies[i]->halfHistories();
    }
  }




}


void TFBShelper::nextChar(char chr)
{
  this->seqCount++;

  for(unsigned int i=0;i<this->SNPs.size();i++) {
    this->SNPs[i].pos++;
  }
  //static int erot=0;


  // Remove the SNPs that we have past and left behind
  if(this->SNPs.size()>0) {

    this->removeScannerHistories();

    if( this->SNPs[0].pos>this->maxLen) {



      // Remove the passed SNP alleles
      this->SNPs.pop_front();
      this->SNPs.pop_front();
      
      if(this->haveBG) {
	// Remove record keeping for the alleles for the background
	deque<deque<double> >::iterator probIter=this->probBuffer.begin();
	deque<matrix_bgObject>::iterator bgIter=this->bg.begin();
	// Remove every second value
	while(probIter!=this->probBuffer.end()) {
	  probIter=this->probBuffer.erase(probIter);
	  probIter++;
	  bgIter=this->bg.erase(bgIter);
	  bgIter++;
	}

      }
      
    }
  }


  
  
  char *allels=getAllels(chr);
  if(allels) {
    for(int allel_p=0;allel_p<2;allel_p++) {
      this->SNPs.push_back(SNPdat(chr,allels[allel_p],0));
    }
  }
  for(unsigned int i=0;i<this->matricies.size();i++) {
    this->matricies[i]->nextChar(chr);
  }


  if(this->haveBG) {
    if(allels==NULL) {
      this->nextACGT(chr);
    } else {
      int newStart=this->probBuffer.size();
      this->probBuffer.insert(this->probBuffer.end(),
			      this->probBuffer.begin(),this->probBuffer.end());
      this->bg.insert(this->bg.end(),
		      this->bg.begin(),this->bg.end());
      for(int allel_p=0;allel_p<2;allel_p++) {
	if(allel_p==0) {
	  this->nextACGT(allels[allel_p],0,newStart);
	} else if(allel_p==1) {
	  this->nextACGT(allels[allel_p],newStart);
	} else {
	  cout<<"Now something completely different!"<<endl;
	}
      }
    }
//     printf("added BG snps: %d buffers: %d\n",this->SNPcount(),this->probBuffer.size());
  }
  
//   if(chr!='A' && chr!='C' && chr!='G' && chr!='T' ) {
//     printf("snps: %d buffers: %d bg=%d\n",this->SNPcount(),this->probBuffer.size(),this->haveBG);
//   }

}

void TFBShelper::nextACGT(char chr,unsigned int startFrom,unsigned int upTo)
{

  if(upTo<=startFrom) {
    upTo=this->probBuffer.size();
  }

  for(unsigned int i=startFrom;i<upTo;i++) {  // ACG or T with previous snps

    double bgP=logPnextInStream(&this->bg[i],chr);

    // Iterate the buffer
    this->probBuffer[i].pop_back();
    this->probBuffer[i].push_front(0.0);

  // Add the probabilities for this character
    for(int j=0;j<this->maxLen;j++) {
      this->probBuffer[i][j]+=bgP;
    }
  }
}

vector<TFBSscan*> parseMatricies(int *count,PyObject *mats,PyObject *cutoffs,double cutoff)
{
  double bound;
  *count=PySequence_Length(mats);
  vector<TFBSscan*> ret;
  


  assert(PySequence_Check(mats));
  for(int i=0;i<*count;i++) {
    //assert(PySequence_Check(ret[i].py_matrix));

    if(cutoffs) {
      bound=PyFloat_AsDouble(PyNumber_Float(PySequence_GetItem(cutoffs,i)));
    } else {
      bound=cutoff;
    }
    ret.push_back(new TFBSscan(PySequence_GetItem(mats,i),bound));


    if(PyErr_Occurred()!=NULL) {
      return vector<TFBSscan*>();
    }
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
  int matrixCount=-1,seq_i;

  int loop_status;

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
  vector<TFBSscan*> Mat=parseMatricies(&matrixCount,py_matlist,py_cutoff,cutoff);

  if(Mat.size()==0 || PyErr_Occurred()!=NULL) {
    return 0;
  }
  sort(Mat.begin(),Mat.end(),lessThan);
  TFBShelper scanner(bg,Mat);
  /*
  printf("matricies: %d cutoffs: %d cutoff: %g \n",
	 PySequence_Length(py_matlist),(py_cutoff?PySequence_Length(py_cutoff):-1),
	 cutoff);

  */
  PyObject* ret;

  if(!(ret=PyDict_New())) {
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

      // Reading from a file like object.
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
    } // End of read.
    
    char nuclChr;
    switch(nuclChr=toupper(Seq[buf_p]))
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
      default:
	scanner.nextChar(nuclChr);
	break;
      }


      


    if(loop_status==loop_break) {
      break;
    } else if(loop_status==loop_continue) {
      loop_status=loop_OK;
      continue;
    }
    

    vector<TFBShit*> hits=scanner.getMatches();
//     if(hits.size()>0) {
//       printf("hits: %d pos=%d\n",hits.size(),scanner.seqPos());
//     }
    for(unsigned int i=0;i<hits.size();i++) {
      PyObject *snps=hits[i]->buildPySNPs();
      assert(PyTuple_Check(snps));
      assert(ret);
      addMatchWithKey(ret,hits[i]->mat->py_matrix,hits[i]->pos,hits[i]->strand,hits[i]->score,snps,hits[i]->minScore);

#ifndef NDEBUG
      if((hits[i]->score-hits[i]->minScore)>1.0) {
	char *str=PyString_AsString(PyObject_GetAttrString(hits[i]->mat->py_matrix,"name"));

	printf("pos=%d %s score_delta=%g\n",hits[i]->pos,str,hits[i]->score-hits[i]->minScore);
      }
#endif
      delete hits[i];
      hits[i]=NULL;
      if(PyErr_Occurred()!=NULL) {
	return NULL;
      }
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

  if(py_strBuf)
    Py_DECREF(py_strBuf);

  return ret;
}



static PyMethodDef matrixMethods[] = {
  {"draw",  matrix_draw, METH_VARARGS,
   "Draws a matrix"},
  {"computeBG",  matrix_computeBG, METH_VARARGS,
   "Kind of computes higher order Background: (seq,order)"},
  {"getAllTFBSwithBg", matrix_getAllTFBSwithBG, METH_VARARGS,
   "Returns a map from matrix to index to score of possible TFBS.\nThe arguments are a list of matricies, the sequence, list/float of cutoffs and a background object."},
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
