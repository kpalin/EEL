#include <Python.h>
#include "structmember.h"

#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#define PI 3.1415926

#define ALIGN_OUTPUT
#ifdef ALIGN_OUTPUT
#include <iostream>
#endif

#include <sys/times.h>
#include <unistd.h>

#define CHECKING_DECREF(X) if((PyObject*)(X)==Py_None) printf("none decref line %d",__LINE__); Py_DECREF(X);

#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <limits>



/*
 *
 * $Log$
 * Revision 1.2  2004/06/22 12:29:05  kpalin
 * Presumably working. No return values back to python
 * but the alignment looks like working as planned.
 *
 * Revision 1.1  2004/06/22 08:41:18  kpalin
 * A full version that compiles. Doesn't probably work
 * and is totally untested and slow but it's just a first trial.
 *
 *
 */



/* The binding sites have to be less than MAX_BP_DIST apart */
#ifndef MAX_BP_DIST
#define MAX_BP_DIST 1000
#endif

//#undef NDEBUG
#include <assert.h>




using namespace std;


#include "multiAlign.h"



void ABset::addChar(int seq,int chr,char strand)
{
  if(strand=='+') {
    chr=chr*2;
  } else if(strand=='-') {
    chr=chr*2+1;
  } else {
    abort();
  }
  B[chr].push_back(seq);
  if(B[chr].size()>=(uint)limit) {
    A.push_back(chr);
  }
}

vector<int> ABset::nextB()
{

  return B[A[pointer++]];
}

int ABset::moreAlphas()
{
  if(pointer<(int)A.size())
    return 1;
  return 0;
}

ABset Inputs::getAB(PointerVec &p,int limit)
{
  ABset out=ABset(factors(),limit);

  for(int i=0;i<sequences();i++) {
    out.addChar(i,seq[i][p[i]].ID,seq[i][p[i]].strand);
  }

  return out;
}


PointerVec::PointerVec(vector<int> &p, vector<int> &dims,vector<int> *dimFactors)
{
  dimlen=dims;
  limData=NULL;
  m=dims.size();
  bound.resize(m,0);
  this->dimFactors=dimFactors;


  setValue(p);
}


posind PointerVec::difference(int i) {
  assert(limData!=NULL);
  return limData->getSite(dimlen[i],i).pos-
    limData->getSite(*this,i).epos;
}



const PointerVec& PointerVec::operator--(int dummy)
{
  uint i=0;

  do {
    if(bound[i]==0) {
      if( p[i]>0 &&
	  ( limData==NULL ||  // Bound decreasing to box within limitBP
	    difference(i)<limitBP) ) {
	p[i]--;
	matrix_p-=(*dimFactors)[i];
	break;
      } else {
	p[i]=dimlen[i]-1;
	matrix_p+=p[i]*(*dimFactors)[i];
	i++;
      }
    } else {
      i++;
    }
  } while(i<m);

  if(i>=m) {
    ok=0;
  }

  assert(matrix_p==dataPoint());

  return *this;
}

const PointerVec& PointerVec::operator++(int dummy)
{
  uint i=0;

  do {
    if(bound[i]==0) {
      p[i]++;
      matrix_p+=(*dimFactors)[i];
      if(p[i]>=dimlen[i]) {
	p[i]=0;
	matrix_p-=dimlen[i]*(*dimFactors)[i];
	i++;
      } else {
	break;
      }
    } else {
      i++;
    }
  } while(i<m);

  if(i>=m) {
    ok=0;
  }

  assert(matrix_p==dataPoint());
  return *this;
}




Matrix::Matrix(vector<int> &dims)
{
  matSize=1;

  dim=dims.size();
  dimLen=dims;

  dimFactors.resize(dim);


  for(int i=0;i<dim;i++) {
    matSize*=dimLen[i];
    if(i==0) {
      dimFactors[i]=1;
    } else {
      dimFactors[i]=dimLen[i-1]*dimFactors[i-1];
    }
  }
  try {
    data.resize(matSize);
  }
  catch (std::bad_alloc const&) {
    cout << "Memory allocation fail!" << endl;
    PyErr_SetString(PyExc_MemoryError,"Out of memory!");
  }

}

PointerVec Matrix::getOrigin()
{
  vector<int> zeros=vector<int> (dim,0);

  return PointerVec(zeros,dimLen,&dimFactors);
}


int PointerVec::dataPoint()
{ // Compute the mapping from multi to single dimensional array

  int pos=0;

  for(uint i=0;i<m;i++) {
    pos+=p[i]*(*dimFactors)[i];
    assert(p[i]<dimlen[i]);
  }
  return pos;
}


matrixentry Matrix::getValue(PointerVec &p)
{
  int pos=p.getMatrixCoord();

  return data[pos];
}

void Matrix::setValue(PointerVec &p,matrixentry &val)
{
  int pos=p.getMatrixCoord();

  data[pos]=val;

}



Inputs::Inputs(PyObject *inpSeq)
{
  // inpSeq is pySequence of pySequences as represented in the GFF file

  PyObject *inpIter=PyObject_GetIter(inpSeq);
  PyObject *inSite;

  if(inpIter==NULL || PyErr_Occurred()!=NULL) {
    return;
  }
  assert(inpIter!=NULL);

  while((inSite=PyIter_Next(inpIter))!=NULL) {
    if(addSite(inSite)==0 || PyErr_Occurred()) {
      return;
    }
  }


}

int Inputs::addSite(PyObject *site)
{
  // site should be a sequence with fields as in GFF file

  int n=PySequence_Size(site);
  int tf_id,seq_id;


  if(n<0) {
    PyErr_SetString(PyExc_TypeError,"A non sequence site data");
    return 0;
  } else if(n<7) {
    PyErr_SetString(PyExc_ValueError,"Not enough fields in the site data");
    return 0;
  }
 
  // Sequence name and id
  string seqName=string(PyString_AsString(PySequence_GetItem(site,0)));
  if(PyErr_Occurred()) return 0;
  
  map<string,uint>::iterator seq_iter=SEQ_to_id.find(seqName);

  if(seq_iter==SEQ_to_id.end()) {
    seq_id=SEQ_to_id.size();
    SEQ_to_id[seqName]=seq_id;
    seq.resize(seq_id+1);
  } else {
    seq_id=seq_iter->second;
  }

  // Transcription factor name name and id
  string TFname=string(PyString_AsString(PySequence_GetItem(site,2)));
  if(PyErr_Occurred()) return 0;
  
  map<string,uint>::iterator id_iter=TF_to_id.find(TFname);

  if(id_iter==TF_to_id.end()) {
    tf_id=TF_to_id.size();
    TF_to_id[TFname]=tf_id;
  } else {
    tf_id=id_iter->second;
  }


  // Create a struct for the new site
  id_triple new_site;

  new_site.ID=tf_id;
  new_site.pos=(posind)PyInt_AsLong(PyNumber_Int(PySequence_GetItem(site,3)));
  if(PyErr_Occurred()) return 0;
  new_site.epos=(posind)PyInt_AsLong(PyNumber_Int(PySequence_GetItem(site,4)));
  if(PyErr_Occurred()) return 0;
  new_site.weight=(double)PyFloat_AsDouble(PyNumber_Float(PySequence_GetItem(site,5)));
  if(PyErr_Occurred()) return 0;

  char *strand_p=PyString_AsString(PySequence_GetItem(site,6));
  if(PyErr_Occurred()) return 0;
  new_site.strand=(char)strand_p[0];


  // Finally add the new site to the input object
  seq[seq_id].push_back(new_site);

  return 1;
}


vector<id_triple> Inputs::getSites(PointerVec &p)
{
  vector<id_triple> out=vector<id_triple>();

  for(int i=0;i<sequences();i++) {
    out.push_back(seq[i][p[i]]);
  }

  return out;
}

vector<int> Inputs::sequenceLens() 
{
  vector<int> dims=vector<int>();
  for(int i=0;i<sequences();i++) {
    dims.push_back(seq[i].size());
  }
  return dims;
}



matrixentry::matrixentry(store v,PointerVec &p)
{
  backTrack=p;
  value=v;
}

store matrixentry::getValue()
{
  return value;
}

PointerVec& matrixentry::getBacktrace()
{
  return backTrack;
}


PointerVec PointerVec::project(vector<int> &boundDims)
{
  PointerVec ret=PointerVec();

  if(boundDims.size()>=m) {
    return ret;
  }

  ret.dimlen=dimlen;
  ret.dimFactors=dimFactors;
  ret.bound.resize(m,0);
  ret.p=p;

  for(uint i=0;i<m;i++) {
    ret.p[i]=ret.p[i]-1;
  }

  for(uint i=0;i<boundDims.size();i++) {
    ret.bound[boundDims[i]]=1;
    ret.p[boundDims[i]]=p[boundDims[i]];    
  }


  matrix_p=dataPoint();
  ok=1;
  return ret;
}


PointerVec PointerVec::projectAndLimit(vector<int> &boundDims,Inputs *indata,
				       int maxbp)
{
  PointerVec ret=project(boundDims);
  ret.limData=indata;
  ret.dimlen=p;
  ret.limitBP=maxbp;

  return ret;
}
  


void PointerVec::setValue(vector<int> &np)
{
  p=np;
  matrix_p=dataPoint();
  ok=1;
}

void PointerVec::clearProject()
{
  for(uint i=0;i<m;i++) {
    bound[i]=0;
  }
}


//returns the square of a mod 2 PI double
inline double squaremod(double val)
{
  double f=fabs(val);
  f-=2*PI*trunc(f/(2*PI));

#ifndef NDEBUG
  double apuf=fabs(val);
  apuf-=2*PI*round(apuf/(2*PI));
  if(apuf>PI) apuf-= 2*PI;

  if(f>(PI+1e9) || f<(-PI-1e9)) {
    printf("f=%g!=%g=PI\n",f,PI);
  }
#endif
  assert(fabs(f-apuf)<1e9);
  assert(f<(PI+1e9));
  assert(f>(-PI-1e9));
  return f * f;
}





typedef struct {
  unsigned long int x;
  unsigned long int y;
} matCoord;

bool operator<(const pair<matCoord,store> &t1,store& t2) {
  return t1.second<t2;
}

bool operator<(const matCoord& t1,const matCoord& t2){
  return t1.x < t2.x || (t1.x==t2.x && t1.y < t2.y);
}

typedef struct {
  store value;
  matCoord bgin;
  matCoord end;
} MS_res;  // Memory Save Result




bool operator<(const MS_res& t1,const MS_res& t2){
  return t1.value < t2.value || (t1.value==t2.value && t1.bgin < t2.bgin) ||
   (t1.value==t2.value && t1.bgin.x==t2.bgin.x  && t1.bgin.y==t2.bgin.y && t1.end < t2.end) ;
}



bool operator>(const MS_res& t1,const MS_res& t2){
  return !(t1<t2);  // Not really, but I'm too lazy to do it correctly for now
}


//////////////////////////////////////////////////////////////////////
// Aligned site object

typedef struct {
  PyObject_HEAD

    /* Type-specific fields go here. */

  char *motifName;
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

static PyMemberDef column_members[] = {
    {"motif",T_STRING, offsetof(align_siteObject, motifName), 0,
     "Name of the motif."},

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


static void column_dealloc(align_siteObject* self)
{
  delete [] self->motifName;

  self->ob_type->tp_free((PyObject*)self);
}

static PyTypeObject malign_colType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "align.SitePair",             /*tp_name*/
    sizeof(align_siteObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)column_dealloc,                         /*tp_dealloc*/
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
    column_members,             /* tp_members */
};

static PyObject *
column_new(const char *motifName,
	 int seqX,
	 int seqY,
	 int beginX,
	 int endX,
	 int beginY,
	 int endY,
	 char strand,
	 double score,
	 double siteScoreX,
	 double siteScoreY)
{
  align_siteObject *ret;
  ret=(align_siteObject*)malign_colType.tp_alloc(&malign_colType,0);

  ret->motifName=new char[strlen(motifName)+1];
  strcpy(ret->motifName,motifName);

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

//////////////////////////////////////////////////////////////////////




static 
void
malignment_dealloc(malign_AlignmentObject* self)
{
    
    delete self->CP;
    

    self->ob_type->tp_free((PyObject*)self);
}



extern "C" int
malignment_init(malign_AlignmentObject *self, PyObject *args, PyObject *kwds)
{


  double lambda, xi, mu, nuc_per_rotation,nu;
  string firstSeqName,secondSeqName,sequence;
  int result_ask;
  PyObject *data;
  static char *kwlist[] = {"data","results_ask","lambda","xi","mu","nu","nuc_per_rotation",NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oiddddd", kwlist, 
				  &data,&result_ask,
			&lambda, &xi, &mu, &nu, &nuc_per_rotation)){
    return -1;
  }

  //cout<<"NYT AJETAAN INIT_IÄ"<<endl;
  self->CP=new struct __CPSTUF;
  if(!self->CP) {
    PyErr_SetString(PyExc_MemoryError,"Out of memory!");
    return -1;
  }
  assert(self->CP!=NULL);

  if(PySequence_Check(data)==0) {
    PyErr_SetString(PyExc_ValueError,"First parameter must be sequence!");
    return -1;
  }

  malign_alignCommon(self,data,result_ask,lambda,xi,mu,nu,nuc_per_rotation);

  if(PyErr_Occurred()) {
    Py_XDECREF(self);
    return -1;
  }

  return 0;
}

static 
PyObject *
malignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    malign_AlignmentObject *self;

    self = (malign_AlignmentObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
      self->secs_to_align = 0.0;
      self->nu=0.0;
      self->names = (PyTupleObject*)PyTuple_New(0);
      if (self->names == NULL)
	{
	  CHECKING_DECREF(self);
	  return NULL;
	}
        
      self->bestAlignments=PyList_New(0);
      if (self->bestAlignments == NULL)
	{
	  CHECKING_DECREF(self);
	  return NULL;
	}

      self->CP=new struct __CPSTUF;
      cout<<"ASETETAAN CP UUTEEN MALIGN OBJEKTIIN!"<<endl<<flush;

      if (self->CP==NULL) {
	return NULL;
      }
    }

    
    self->item_count=0;
    assert(self->CP!=NULL);
    cout<<"CP ON EPÄNULL"<<endl<<flush;
    return (PyObject *)self;
}









static PyMethodDef malignment_methods[] = {
//     {"nextBest", (PyCFunction)malignment_nextBest, METH_NOARGS,
//      "Return the next best alignment from the matrix."
//     },
    {NULL}  /* Sentinel */
};

static PyMemberDef malignment_members[] = {
    {"bestAlignments",T_OBJECT_EX, offsetof(malign_AlignmentObject, bestAlignments), 0,
     "List of best alignments"},
    {"secs_to_align",T_DOUBLE, offsetof(malign_AlignmentObject, secs_to_align), 0,
     "CPU time for Alignment"},
    {"names",T_OBJECT_EX, offsetof(malign_AlignmentObject, names), 0,
     "Names of the sequences"},
    {"Lambda",T_DOUBLE, offsetof(malign_AlignmentObject,lambda), 0,
     "Parameter for bonus"},
    {"Xi",T_DOUBLE, offsetof(malign_AlignmentObject,xi), 0,
     "Parameter for distance penalty"},
    {"Nu",T_DOUBLE, offsetof(malign_AlignmentObject,nu), 0,
     "Parameter for distance difference penalty"},
    {"Mu",T_DOUBLE, offsetof(malign_AlignmentObject,mu), 0,
     "Parameter for rotation penalty"},
    {"nuc_per_rotation",T_DOUBLE, offsetof(malign_AlignmentObject,nuc_per_rotation), 0,
     "Parameter for nucleotides per 360 deg rotation of DNA"},
    {"item_count",T_INT, offsetof(malign_AlignmentObject,item_count), 0,
     "Number of filled cells."},
    {"askedResults",T_INT, offsetof(malign_AlignmentObject,askedresults), 0,
     "Number of filled cells."},
    {NULL}  /* Sentinel */
};


static PyTypeObject malign_AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "multiAlign.MAlignment",             /*tp_name*/
    sizeof(malign_AlignmentObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)malignment_dealloc,                         /*tp_dealloc*/
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
    "Multiple alignment object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    malignment_methods,             /* tp_methods */
    malignment_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)malignment_init,      /* tp_init */
    0,                         /* tp_alloc */
    malignment_new,                 /* tp_new */
};






//////////////////////////////////////////////////////////////////////





inline double anglepenalty(double const d,double const D,double const nucl_per_rotation)
{
  double const theta=(d-D)*2.0*PI/nucl_per_rotation;
  return squaremod(theta)/(d+D);
}


inline double penalty(malign_AlignmentObject *dat,posind d1,posind d2)
{
  double val=dat->mu*(d1+d2)/2.0;

  if( (d1+d2)>0.0) {
    val=val+ dat->nu*(d1-d2)*(d1-d2)/(d1+d2)+
      dat->xi*anglepenalty(d1,d2,dat->nuc_per_rotation);
  }
  return val;
}





void outputMemory(double bytes)
{
  if(bytes>(1024*1024)) {
    cout<<(bytes/(1024*1024.0))<<" megabytes";
  } else {
    cout<<(bytes/(1024))<<" kilobytes";
  }
  
}




static PyObject *
malignObject(malign_AlignmentObject *self)
{ 
  store Score;
  // The main multiple alignment algorithm.

  int Empty=0,nonEmpty=0;

  for(PointerVec entry=self->CP->dynmat->getOrigin();entry.isOK();entry++) {
#ifndef NDEBUG
    //entry.output();
#endif
    ABset AB=self->CP->indata->getAB(entry);

    vector<id_triple> sites=self->CP->indata->getSites(entry);

    if(AB.moreAlphas()) {
      nonEmpty++;
//       for(int i=0;i<self->CP->indata.sequences();i++) {
// 	cout<<sites[i].ID<<sites[i].strand<<",";
//       }
//       cout<<" non empty A";
//       cout<<endl;
    } else {
      Empty++;
    }

    store maxScore=0;
    PointerVec maxSource=PointerVec();
    while(AB.moreAlphas()) {
      vector<int> Bak=AB.nextB();

      store base=0;
      int n=Bak.size();
      for(int i=0;i<n;i++) {
	base+=sites[Bak[i]].weight;
      }
      base*=self->lambda*n*(n-1)/2;


      for(PointerVec k=entry.projectAndLimit(Bak,self->CP->indata,MAX_BP_DIST);
	  k.isOK(); k--) {  // Exponential loop

	Score=base+(*self->CP->dynmat)[k].getValue();
	for(uint i=0;i<Bak.size();i++) {
	  for(uint j=i+1;j<Bak.size();j++) {
	    Score=Score - penalty(self,k.difference(i),k.difference(j));
	  }
	}



	if(Score>maxScore) {
	  maxScore=Score;
	  maxSource=k;
	}

	if(PyErr_Occurred()) return NULL;
      }


    }
    (*self->CP->dynmat)[entry]=matrixentry(maxScore,maxSource);

  }

#ifdef DEBUG_OUTPUT
  cout<<"Matrix dimensions: ";
  for(int i=0;i<self->CP->dynmat->dims();i++) {
    cout<<" "<<self->CP->dynmat->dims(i);
  }
  cout<<endl
      <<"Matrix size: "<<self->CP->dynmat->maxSize()<<endl
      <<"Empty:       "<<Empty<<endl 
      <<"Non empty:   "<<nonEmpty<<endl 
      <<"share of non empty: "<<nonEmpty/(1.0*Empty+nonEmpty)<<endl
      <<"Using matrix of "<<self->CP->dynmat->size()<<" items."<<endl;
#endif
  return (PyObject*)self;
}



////////////////////////////////////////////////////////////

static PyObject *
malign_alignCommon( malign_AlignmentObject *self,PyObject *data,int result_ask,
		   double lambda,double xi,double mu,double nu,
		   double nuc_per_rotation)
{
  PyObject *ret_obj;
  self->CP->indata=new Inputs(data);

  if(PyErr_Occurred()) return NULL;


  if(self->CP->indata->sequences()<3) {
    PyErr_SetString(PyExc_EOFError,"Too few sequences in input");
    return NULL;
  }




  vector<int> dimensions=self->CP->indata->sequenceLens();
  self->CP->dynmat=new Matrix(dimensions);

  if(PyErr_Occurred()) return NULL;
  

  self->lambda=lambda;
  self->xi=xi;
  self->mu=mu;
  self->nu=nu;
  self->nuc_per_rotation=nuc_per_rotation;
  self->askedresults=result_ask;
  self->secs_to_align=0;


  tms before,after;
  long ticks_per_sec=sysconf(_SC_CLK_TCK);
  long ticks_to_align;
  // Start timing
  times(&before);


  ret_obj=(PyObject*)malignObject(self);


  // End timing
  times(&after);

  ticks_to_align=((after.tms_utime-before.tms_utime)+
		 (after.tms_stime-before.tms_stime));
  self->secs_to_align+=((double)ticks_to_align)/ticks_per_sec;
  


  return ret_obj;

}





static PyObject *
malign_aligndata(PyObject *self, PyObject *args)
{
  double lambda, xi, mu, nuc_per_rotation,nu;
  string firstSeqName,secondSeqName,sequence;
  int result_ask;
  PyObject *data;

  PyErr_SetString(PyExc_ValueError,"First parameter must be sequence!");
  return NULL;


  // Trying to avoid this.
  if (!PyArg_ParseTuple(args, "Oiddddd", &data,&result_ask,
			&lambda, &xi, &mu, &nu, &nuc_per_rotation)){
    return Py_BuildValue("s", "");
  }

  if(PySequence_Check(data)==0) {
    return Py_BuildValue("s","");
  }
  

  return malign_alignCommon((malign_AlignmentObject*)self,data,result_ask,lambda,xi,mu,nu,nuc_per_rotation);
}




//////////////////////////////////////////////////////////////////////





static PyMethodDef malignMethods[] = {
  {"aligndata",  malign_aligndata, METH_VARARGS,
   "aligns computed sequences\nArguments: data,result_ask,lambda,xi,mu,nu,nuc_per_rotation"},
//   {"alignfile",  align_alignfile, METH_VARARGS,
//    "aligns sequences from a gff-file"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


//PyMODINIT_FUNC
extern "C"
void
initmultiAlign(void)
{
    PyObject* m=NULL;

    malign_AlignmentType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&malign_AlignmentType) < 0)
      return;

    malign_colType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&malign_colType) < 0)
      return;
    
    m=Py_InitModule("multiAlign", malignMethods);

    if(m==NULL)
      return;

    Py_INCREF(&malign_AlignmentType);
    PyModule_AddObject(m, "MultiAlignment", (PyObject *)&malign_AlignmentType);

    Py_INCREF(&malign_colType);
    PyModule_AddObject(m, "SiteColumn", (PyObject *)&malign_colType);

#ifndef NDEBUG
    cout<<"MultiAlignLoaded"<<endl;
#endif
}

