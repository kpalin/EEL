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


#include "alignedCols.h"

// GDB command for breaking after importing dynamic module:
//   br _PyImport_LoadDynamicModule

/*
 *
 * $Log$
 * Revision 1.4  2004/07/23 11:53:31  kpalin
 * Compiles, but I don't think it works. The projection stuff seems more
 * difficult than I expected.
 *
 * Revision 1.3  2004/07/12 12:49:25  kpalin
 * Now a working connection *FROM* python2.2. TODO: Sparse matrix and
 * connection *TO* python.
 *
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


siteCode fromSiteAndStrand(int chr,char strand)
{
  if(chr<0) {
    cerr<<"Invalid site code "<<chr<<endl;
    abort();
  }
  if(strand=='+') {
    chr=chr*2;
  } else if(strand=='-') {
    chr=chr*2+1;
  } else {
    cerr<<"Invalid strand "<<chr<<endl;
    abort();
  }

  return (siteCode)chr;
}

char toStrand(siteCode code)
{
  if(code%2==0) {
    return '+';
  } else {
    return '-';
  }
}

int toSite(siteCode code)
{
  return code/2;
}


void ABset::addChar(int seq,int chr,char strand)
{
  chr=fromSiteAndStrand(chr,strand);

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


  //p.output();
  for(int i=0;i<sequences();i++) {
    out.addChar(i,seq[i][p[i]].ID,seq[i][p[i]].strand);
  }
  //cout<<(out.moreAlphas()?"more alphas":"no more alphas")<<endl;
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


posind PointerVec::difference(int i) const {
  assert(limData!=NULL);
  return abs(limData->getSite(dimlen[i],i).pos-
    limData->getSite(*this,i).epos);
}



const PointerVec& PointerVec::operator--(int dummy)
{
  uint i=0;
  posind ero=-1;

  do {
    ero=difference(i);
    if( p[i]>0 &&
	( limData==NULL || bound[i]==1 || // Bound decreasing to box within limitBP
	  difference(i)<limitBP) ) {
      p[i]--;
      matrix_p-=(*dimFactors)[i];
      // 	cout<<"("<<i<<":"<<ero<<")"<<flush;
      break;
    } else{
      int delta_p=max((dimlen[i]-1)-p[i]-1,0);
      p[i]+=delta_p;
      matrix_p+=delta_p*(*dimFactors)[i];
      i++;
    } 
    assert(i==m|| p[i]>=0 );
    assert(i==m|| p[i]<dimlen[i]);
  } while(i<m);

  //this->output();
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
    if(bound[i]==0) {  // If the dimension is free
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

PointerVec Matrix::argMax()
{
  PointerVec p=this->getOrigin();
  PointerVec maxstore,*backer;
  store maxval=0.0;

  while(p.isOK()) {
    if(this->getValue(p).getValue()>maxval) {
      for(backer=this->getValue(p).getBacktraceP();
	  backer->isOK() && this->getValue(*backer).getValue()>=0.0;
	  backer=this->getValue(*backer).getBacktraceP()) {
	/* pass*/
      }
      if(!backer->isOK()) { // Only store, if we do not backtrack to negative.
	maxstore=p;
	maxval=this->getValueP(p)->getValue();
      }
    }
    p++;
  }
#ifndef NDEBUG
  if(!maxstore.isOK()) {
    cout<<"Busted matrix, all zeros"<<endl;
  }
#endif
    
  return maxstore;
}


int PointerVec::dataPoint() const
{ // Compute the mapping from multi to single dimensional array

  int pos=0;

  for(uint i=0;i<m;i++) {
    pos+=p[i]*(*dimFactors)[i];
    if(p[i]>=dimlen[i]) {
      cout<<"p[i] = "<<p[i]<<endl<<dimlen[i]<<endl;
      assert(p[i]<dimlen[i]);
    }
  }
  assert(pos>=0);
  return pos;
}


matrixentry const &Matrix::getValue(PointerVec const &p) const
{
  int pos=p.getMatrixCoord();

  return data[pos];
}

matrixentry *Matrix::getValueP(PointerVec const &p)
{
  int pos=p.getMatrixCoord();

  return &data[pos];
}

void Matrix::setValue(PointerVec &p,matrixentry &val)
{
  int pos=p.getMatrixCoord();

  data[pos]=val;

}




matrixentry &Matrix::operator[](PointerVec &p)
{ 
  return data[p.getMatrixCoord()]; 
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
    Py_DECREF(inSite);
  }
  Py_DECREF(inpIter);
  
  map<string,uint>::iterator seqName=sequenceNames();
  for(int i=0;i<sequences();i++) {
    cout<<seqName->first<<",";
    seqName++;
    sort(seq[i].begin(),seq[i].end());
  }
  cout<<endl;



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
  string &seqName=*(new string(PyString_AsString(PySequence_GetItem(site,0))));
  if(PyErr_Occurred()) return 0;
  
  map<string,uint>::iterator seq_iter=SEQ_to_id.find(seqName);

  if(seq_iter==SEQ_to_id.end()) {
    seqNames.push_back(seqName);
    seq_id=SEQ_to_id.size();
    SEQ_to_id[seqName]=seq_id;
    seq.resize(seq_id+1);
  } else {
    seq_id=seq_iter->second;
  }

  // Transcription factor name name and id
  string &TFname=*(new string(PyString_AsString(PySequence_GetItem(site,2))));
  if(PyErr_Occurred()) return 0;
  
  map<string,uint>::iterator id_iter=TF_to_id.find(TFname);

  if(id_iter==TF_to_id.end()) {
    factorNames.push_back(TFname);
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


void matrixentry::setInitData(store v,int site,char strand,PointerVec *p)
{
  assert(&p);

  if(p && p->isOK()) {
    backTrack=new PointerVec(*p);
  } else {
    backTrack=NULL;
  }
  value=v;
  siteEtStrand=fromSiteAndStrand(site,strand);
}


matrixentry::matrixentry(store v,int site,char strand,PointerVec &p)
{
  this->setInitData(v,site,strand,&p);
}

matrixentry::matrixentry(store v,int site,char strand,PointerVec *p)
{
  this->setInitData(v,site,strand,p);
//   if(p && p->isOK()) {
//     backTrack=new PointerVec(*p);
//   } else {
//     backTrack=NULL;
//   }
//   value=v;
//   siteEtStrand=site;
}

store matrixentry::getValue() const
{
  return value;
}

int matrixentry::getSite() const
{
  return toSite(siteEtStrand);
}

char matrixentry::getStrand() const
{
  assert(siteEtStrand>=0);

  return toStrand(siteEtStrand);
}

void PointerVec::output()
{
  if(this->isOK()) {
   // cout<<"valid("<<m<<"): ";
    for(uint i=0;i<m;i++) {
      cout<<p[i]<<(bound[i]?"b":"f")<<dimlen[i]<<",";
    } } else { cout<<"Invalid!"; }
  //  cout<<endl; 
  cout<<flush;
}



/* Projections:

   A dimension is "free" if we wish to explore the cube in that
   dimension. 
   The "bound" dimensions are irrelevant and are not explored
   because they have no effect on the scoring.
*/
PointerVec PointerVec::project(vector<int> &boundDims)
{
  PointerVec ret=PointerVec();

  if(p[0]==3 && p[1]==5 && p[2]==0) {
    cerr<<"JOTAIN JÄNNÄÄ"<<endl;
  }


  if(boundDims.size()>m) {  
    return ret;
  }

  ret.dimlen=dimlen;
  for(uint i=0;i<ret.dimlen.size();i++) {
    ret.dimlen[i]=min(this->dimlen[i],p[i]+1);
  }

//   for(uint i=0;i<dimlen.size();i++) {
//     cout<<dimlen[i]<<"="<<ret.dimlen[i]<<endl;
//   }


  ret.dimFactors=dimFactors;
  ret.bound.resize(m,1);
  //ret.p.resize(p.size());
  ret.p=p;
  ret.m=m;

//   for(uint i=0;i<m;i++) {
//     ret.p[i]=max(0,p[i]);
//   }
// #ifndef NDEBUG
//   ret.ok=1;ret.output();
// #endif

  for(uint i=0;i<boundDims.size();i++) {
    ret.bound[boundDims[i]]=0;
    ret.p[boundDims[i]]=max(0,ret.p[boundDims[i]]-1);
  }


  ret.matrix_p=ret.dataPoint();
  ret.ok=1;
  return ret;
}


PointerVec PointerVec::projectAndLimit(vector<int> &boundDims,Inputs *indata,
				       int maxbp)
{
  PointerVec ret=this->project(boundDims);
  if(!ret.isOK()) {
    cout<<"IHMEJAKUMMA"<<endl;
    assert(ret.isOK());
  }
  ret.limData=indata;
  ret.limitBP=maxbp;

  for(uint i=0;i<boundDims.size();i++) {
    int bd=boundDims[i];
    while(ret.p[bd]>0 && 
	  indata->getSite(ret.p[bd],bd).epos >= indata->getSite(this->p[bd],bd).pos) {
      ret.p[bd]--;
      ret.matrix_p-=(*ret.dimFactors)[bd];
    }
  }

  //cout<<"not proj";this->output();
  //cout<<"projected";ret.output();
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
inline double squaremodpi(double val)
{
  double f=fabs(val);
  f-=2*PI*round(f/(2*PI));

#ifndef NDEBUG
  if(f>(PI+1e-6) || f<(-PI-1e-6)) {
    printf("f=%g!=%g=PI\n",f,PI);
  }
#endif

  if(fabs(val-round(val))>1.0) {
    cout<<"round error: abs("<<fabs(val)<<"-"<<round(val)<<")=~"<<
      abs(fabs(val)-round(val))<<"<1.0"<<endl;
  }
  assert(f<(PI+1e-9));
  assert(f>(-PI-1e-9));
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

  self->bestAlignments=PyList_New(0);
  if (self->bestAlignments == NULL || PyErr_Occurred())
    {
      CHECKING_DECREF(self);
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
      self->names = PyTuple_New(0);
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
      {"nextBest", (PyCFunction)malignment_nextBest, METH_NOARGS,
       "Return the next best alignment from the matrix."
      },
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
    {"memSaveUsed",T_INT, offsetof(malign_AlignmentObject,memSaveUsed), 0,
     "Placeholder for compatibility with 2d alignment."},
    {NULL}  /* Sentinel */
};


static PyTypeObject malign_AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "multiAlign.MultiAlignment",             /*tp_name*/
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



static PyObject*
malignment_nextBest(malign_AlignmentObject *self)
{
  /*
    return value:
    goodAlign= [ ( (seq1..seqk), Score, Motif, 
                 ( (start1,end1)..(startk,endk )) , Strand )..]
   */

  PointerVec *p,start=self->CP->dynmat->argMax();

  PyObject *goodAlign=PyList_New(0);
  //  PyObject *colTuple;
  PyObject *seqs;
  PyObject *coords;
  PyObject *seqPos,*siteScore;
  int seqC=0;

#ifndef NDEBUG
  p=&start;  
  cout<<"eka:"<<flush;
  start.output();
#endif

  matrixentry *matVal=NULL;


  for(p=&start;p && p->isOK();p=matVal->getBacktraceP()) {


    matVal=self->CP->dynmat->getValueP(*p);

    char strand=matVal->getStrand();
    siteCode motifID=matVal->getSite();

    // Tuple for one column of the alignment
    //    colTuple=PyTuple_New(5);


    // Tuple for sequence codes
    seqs=PyTuple_New(self->CP->indata->sequences());
    // Tuple for coordinates
    coords=PyTuple_New(self->CP->indata->sequences());

    seqPos=PyTuple_New(self->CP->indata->sequences());
    siteScore=PyTuple_New(self->CP->indata->sequences());




    seqC=0;
    vector<id_triple>  sites=self->CP->indata->getSites(*p);

    for(uint i=0;i<sites.size();i++) {
      if((siteCode)sites[i].ID==motifID && sites[i].strand==strand) {
	// Sequence
	PyTuple_SetItem(seqs,seqC,PyInt_FromLong(i));


	// Coordinates on the sequence
	PyObject *coord;
	coord=Py_BuildValue("(ii)",(int)sites[i].pos,(int)sites[i].epos);
	PyTuple_SetItem(coords,seqC,coord);

	// Position on the site sequence
	PyTuple_SetItem(seqPos,seqC,PyInt_FromLong((*p)[i]));

	PyTuple_SetItem(siteScore,seqC,PyFloat_FromDouble(sites[i].weight));

	seqC++; 
      }
    }

    // Resizing might change the location of the tuples
    _PyTuple_Resize(&seqs,seqC);
    _PyTuple_Resize(&coords,seqC);
    _PyTuple_Resize(&seqPos,seqC);
    _PyTuple_Resize(&siteScore,seqC);

//     PyTuple_SetItem(colTuple,0,seqs);
//     // The alignment score up till here
//     PyTuple_SetItem(colTuple,1,PyFloat_FromDouble(self->CP->dynmat->getValue(*p).getValue()));
//     PyTuple_SetItem(colTuple,2,PyInt_FromLong(motifID));
//     PyTuple_SetItem(colTuple,3,coords);
//     PyTuple_SetItem(colTuple,4,PyString_FromStringAndSize(&strand,1));

    
    PyObject *alnRow=PyAln_New_Multi(self->CP->indata->factor(motifID),
				    seqs,seqPos,coords,strand,
				    self->CP->dynmat->getValue(*p).getValue(),
				    siteScore);

    PyList_Append(goodAlign,alnRow);


    matVal->negate();

#ifndef NDEBUG
    cout<<"Kierretään :"<<flush;
    p->output();
    cout<<"="<<self->CP->dynmat->getValue(*p).getValue()<<endl;
#endif

  }

  return goodAlign;
}



inline double anglepenalty(double const d,double const D,double const nucl_per_rotation)
{
  double const theta=(d-D)*2.0*PI/nucl_per_rotation;
  double const totlen=d+D;

  if(((PI*PI)/totlen)<FLT_EPSILON) { // Check that we can compute it.
    return 0.0;
  } else {
    return squaremodpi(theta)/totlen;
  }
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

  int Empty=0,nonEmpty=0,unNecessary=0;

  for(PointerVec entry=self->CP->dynmat->getOrigin();entry.isOK();entry++) {
#ifndef NDEBUG
    //entry.output();
#endif
    ABset AB=self->CP->indata->getAB(entry);

    vector<id_triple> sites=self->CP->indata->getSites(entry);

#ifndef NDEBUG
    if(AB.moreAlphas()) {
      nonEmpty++;
    } else {
      Empty++;
    }
#endif

    store maxScore=SCORE_FAIL;
    siteCode maxSite=SITE_FAIL;
    char maxStrand=STRAND_FAIL;

    PointerVec *maxSourceP=NULL;
    while(AB.moreAlphas()) {
      maxScore=-DBL_MAX;
      vector<int> Bak=AB.nextB();

      store base=0;
      int n=Bak.size();
      for(int i=0;i<n;i++) {
	base+=sites[Bak[i]].weight;
      }
      base*=self->lambda*n*(n-1)/2;


      if(base>maxScore) {
	maxScore=base;  // On our own. The first site of the alignment.
	maxSite=sites[Bak[0]].ID;
	maxStrand=sites[Bak[0]].strand;
	delete maxSourceP;
	maxSourceP=NULL;
      }

#ifndef NDEBUG
      cout<<"."<<flush;
#endif
      for(PointerVec k=entry.projectAndLimit(Bak,self->CP->indata,MAX_BP_DIST);
	  k.isOK(); k--) {  // Exponential loop

	assert(k.getMatrixCoord()>=0);

#ifndef NDEBUG
	cout<<"*"<<flush;
#endif

	// Continuing alignment.
#ifndef NDEBUG
	if((*self->CP->dynmat)[k].getValue()==0.0) {
	  unNecessary++;
	}
#endif
	Score=base+(*self->CP->dynmat)[k].getValue();
	for(uint i=0;i<Bak.size();i++) {
	  for(uint j=i+1;j<Bak.size();j++) {
	    Score=Score - penalty(self,k.difference(i),k.difference(j));
	  }
	}


	if(Score>maxScore) {
	  maxScore=Score;
	  delete maxSourceP;
	  maxSite=sites[Bak[0]].ID;
	  maxStrand=sites[Bak[0]].strand;
	  maxSourceP=new PointerVec(k);
	  //k.output();
	}

	if(PyErr_Occurred()) return NULL;
      }


    }

    if(maxScore>0.0) {
#ifndef NDEBUG
      entry.output();
      cout<<"="<<maxScore<<"<-";
      if(maxSourceP) {
	maxSourceP->output();
      }else { cout<<"NULL";  }
      cout<<endl;
#endif
      (*self->CP->dynmat)[entry]=matrixentry(maxScore,maxSite,maxStrand,maxSourceP);
    } else {
      (*self->CP->dynmat)[entry]=matrixentry();
    }
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
      <<"Un necessary lookups: "<<unNecessary<<endl
      <<"Using matrix of "<<self->CP->dynmat->size()<<" items."<<endl;
#endif
  return (PyObject*)self;
}



////////////////////////////////////////////////////////////

int setSeqNames(malign_AlignmentObject *self,Inputs &data)
{
  int n=data.sequences();
  map<string,uint>::iterator iter=data.sequenceNames();

  self->names=PyTuple_New(n);

  if(PyErr_Occurred()) return 0;

  for(int i=0;i<n;i++,iter++) {
    PyTuple_SET_ITEM(self->names,i,
		       PyString_FromString(iter->first.c_str()));
    if(PyErr_Occurred()) return 0;
  }
  return 1;

}

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
  self->memSaveUsed=0;
  self->secs_to_align=0;

  if(!setSeqNames(self,*self->CP->indata)){
    return NULL;
  }

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

    
    m=Py_InitModule("multiAlign", malignMethods);

    if(m==NULL)
      return;

    if(import_alnCols()<0)
      return;

    Py_INCREF(&malign_AlignmentType);
    PyModule_AddObject(m, "MultiAlignment", (PyObject *)&malign_AlignmentType);


#ifndef NDEBUG
    cout<<"MultiAlignLoaded"<<endl;
#endif
}

