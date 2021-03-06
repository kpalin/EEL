#include <Python.h>
#include <stdlib.h>
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

#include "alignPenalties.h"



#include "alignedCols.h"

// GDB command for breaking after importing dynamic module:
//   br _PyImport_LoadDynamicModule

/*
 *
 * $Log$
 * Revision 1.29  2006/10/19 07:57:30  kpalin
 * Fixed a seqfault.
 *
 * Revision 1.28  2006/08/31 10:09:34  kpalin
 * Minimum remembered pairwise alignments and speed improvements for multi align.
 *
 * Revision 1.27  2006/08/29 06:08:40  kpalin
 * Speed improvements.
 *
 * Revision 1.26  2006/08/25 12:18:54  kpalin
 * Caching distances to pointervec limiter.
 *
 * Revision 1.25  2006/05/12 06:50:54  kpalin
 * Close to perfection.
 *
 * Revision 1.24  2006/05/11 12:44:43  kpalin
 * Works. But leaks a bit of memory in matrixentry.
 *
 * Revision 1.23  2006/05/04 07:58:38  kpalin
 * Speed improvements.
 *
 * Revision 1.22  2006/05/03 11:00:13  kpalin
 * Merged. Merged a fix with the main devel. branch.
 *
 * Revision 1.21  2006/05/03 10:11:17  kpalin
 * Fixed a bug resulting wrong alignments depending on order of the input
 * sequences.
 *
 * Revision 1.20  2006/01/05 10:33:45  kpalin
 * Fast and working version.
 *
 * Revision 1.19  2006/01/05 09:21:54  kpalin
 * Fast and working version. Still needs cleanup.
 *
 * Revision 1.18  2006/01/04 12:28:30  kpalin
 * Works with getLimited.
 *
 * Revision 1.16  2006/01/04 11:16:04  kpalin
 * Maybe working but debugging version.
 *
 * Revision 1.15  2006/01/03 14:03:39  kpalin
 * Again improved. Should be faster.
 *
 * Revision 1.14  2006/01/02 12:44:34  kpalin
 * Cleaning up coordinate caching.
 *
 * Revision 1.13  2005/12/30 12:43:18  kpalin
 * Significantly (should be factor of n) faster.
 *
 * Revision 1.12  2005/11/25 12:11:02  kpalin
 * Fixed few bugs (input etc.)
 *
 * Revision 1.11  2005/11/24 13:28:34  kpalin
 * Seamingly working version of multi D alignment.
 *
 * Revision 1.10  2005/10/03 10:18:41  kpalin
 * Presumably working multiple alignment version. Not yet usable though.
 *
 * Revision 1.9  2005/03/22 13:24:02  kpalin
 * Totally different approach. Doing multi-D matrix.
 *
 * Revision 1.8  2005/02/09 11:08:09  kpalin
 * Not really working version. As all others before this too.
 *
 * Revision 1.7  2004/12/22 11:14:34  kpalin
 * Some fixes for better distributability
 *
 * Revision 1.6  2004/12/14 14:07:22  kpalin
 * *** empty log message ***
 *
 * Revision 1.5  2004/07/30 12:09:18  kpalin
 * Slow but presumably working version.
 *
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


// Return zero if one of the sequences do not have this factor.
const int PointerVec::allHasFactor() const {
  return this->myMat->allHasFactor(this->getMotif());
}



PointerVec::PointerVec(Matrix *mat,Inputs const *indata)
{
  this->ok=0;
  this->limiterPvec=NULL;
  limData=indata;
  m=mat->dims();
  this->myMat=mat;

  this->matrix_p.resize(this->m,0);

  this->resetMotifCode();
  if(!this->allHasFactor()) {
    (*this)++;
  }
  // Fail if there are no common motifs.
  if(this->allHasFactor()) {
    ok=1;
  }


}


/* set matrix coordinate for site just before the limiter 
   that was previously used with this pointing to motif of type tfID */
void PointerVec::setPrevMatrixCoord(motifCode const tfID,seqCode const i) {
  assert(this->limiterPvec);
  assert(tfID==this->getMotif());

  int newCoord=this->myMat->CoordCacheGet(tfID,i);
  if(this->matrix_p[i]<newCoord) {
    this->myMat->CoordCacheReset(i);
    this->matrixIndexSet(i,0);
    //this->matrix_p[i]=0;
  } else {
    this->matrixIndexSet(i,newCoord);
    //this->matrix_p[i]=newCoord;
  }
  
  motifCode limiterTF=this->limiterPvec->getMotif();
  while(this->matrix_p[i]<(this->myMat->countTFinSeq(i,tfID)) && 
	this->difference(*this->limiterPvec,i,tfID,limiterTF)>=0) {
    //this->matrixIndexInc(i);
    this->matrix_p[i]++;  // BE VERY VERY CAREFULL!!! 
    //this->matrix_p[i]++;
  }

  if( this->matrix_p[i] > 0) {  
    this->matrixIndexDec(i);
    //this->matrix_p[i]--;
  } else {
    this->matrixIndexSet(i,this->matrixIndex(i));  // BE VERY VERY CAREFULL!!! 
  }


  this->myMat->CoordCacheSet(tfID,i,this->matrix_p[i]);

}



bool PointerVec::updateRestCoords()
{ /* For fixed matrix_p[0], update matrix_p[i] i>0 such that
     matrix_p[i] is pointing to a site just before limiterPvec[i].  
     Return 0 if such site does not exist or it is further than limitBP. 
     Return 1 otherwise. */
  motifCode const tfid=this->getMotif();
  //motifCode const limiterTFid=this->limiterPvec->getMotif();

  
  for(uint i=1;i<this->m;i++) {
    this->setPrevMatrixCoord(tfid,i);
    //this->diff2limiter[i]=this->difference(*this->limiterPvec,i,tfid,limiterTFid);
    if(this->difference2limiter(i)>=this->limitBP || this->difference2limiter(i)<0) {
#ifdef DEBUG_OUTPUT
      printf("afterBreak %d\n",this->diff2limiter[i]);
#endif

      // this is out of bounds
      return 0;
    }
    assert(this->checkAtBorder(i));
  

  }

  return 1;
}

// TODO:: THIS MUST BE MADE MUCH MUCH MORE EFFICIENT
// Return true if pointer is valid with respect to limiterPvec

// Return:
//    True  if this is within limits
//    Set this->ok to False, if this->difference(*this->limiterPvec,0)>=this->limitBP
bool PointerVec::decFirst()
{
#ifdef DEBUG_OUTPUT
  cout<<"dF()"<<endl;
#endif


  // BE VERY CAREFULL
  do {
    //this->matrixIndexDec(0);
    this->matrix_p[0]--;
    if(this->matrix_p[0]<0) {
      break;
    }
    this->resetMotifCode();
  } while(!this->allHasFactor());

  this->matrixIndexSet(0,this->matrixIndex(0));
  // CAN EASE UP NOW

  if(this->matrix_p[0]<0) {
    this->ok=0;
    return 0;
  }

  //this->curMotifCode=this->getMotifLaborous();
  // Here we dont need to know the new TF id because we only use dimension 0.
  
  //this->diff2limiter[0]=this->difference(*this->limiterPvec,0);
  
  if(this->diff2limiter[0]<0) {
    //this->ok=0;
#ifdef DEBUG_OUTPUT
    printf("negDist %d\n",this->diff2limiter[0]);
#endif
    // this is ahead of limiterPvec
    return 0;
  }
  if(this->diff2limiter[0]>=this->limitBP) {
    // this is out of bounds on sequence 0.
#ifdef DEBUG_OUTPUT
    printf("HighDist: %d\n",this->diff2limiter[0]);
#endif
    this->ok=0;
    return 0;
  } else {
    return this->updateRestCoords();
  }
  assert(0);
} 



bool PointerVec::checkWithinLimits() const 
{
  bool ret=1;
  
  if(!this->limiterPvec) {
    ret=0;
  } else {
    //motifCode tfID=this->getMotif();

    for(seqCode i=0;(unsigned int)i<(unsigned int)this->m;i++) {
      //int dist=this->difference(*this->limiterPvec,i,tfID);
      int dist=this->difference2limiter(i);
      ret&=(dist>=0);
      ret&=(dist<this->limitBP);
    }
  }
  return ret;
}

bool PointerVec::checkAtBorder()
{
  bool ret=1;

  for(unsigned int i=1;i<this->m;i++) {
    ret&=this->checkAtBorder(i);
  };

  return ret;
}

bool PointerVec::checkAtBorder(seqCode i)
{
  bool ret=1;
  motifCode tfID=this->getMotif();

  if(this->matrix_p[i]<(this->myMat->countTFinSeq(i,tfID)-1)){
    ret&=(this->getSite(i).epos<this->limiterPvec->getSite(i).pos);
    this->matrixIndexInc(i);
    //this->matrix_p[i]++;
    ret&=(this->getSite(i).epos>=this->limiterPvec->getSite(i).pos);
    this->matrixIndexDec(i);
    //this->matrix_p[i]--;

    if(!ret) {
      this->output();
    }
  }

  return ret;
}



bool PointerVec::checkLT() const
{
  if(limiterPvec && !(*this<=*limiterPvec)) {
    cout<<"###"<<endl;
     
    this->output();
    limiterPvec->output();
    for(uint i=0; i<this->m;i++) {
      printf("\n%d: %g <=? %g",i,this->getSite(i).epos,limiterPvec->getSite(i).pos);
    }
    cout<<"###"<<endl;
    return 0;
  }
  return 1;
}

// procedure nextLookBack():
// move this to the 'next' position in the hypercube whose side length is limitBP 
// whose right hand corner is at limiterPvec
//   

void PointerVec::nextLookBack() 
{
  seqCode seq;
  //int tmpDelta;

  assert(this->checkWithinLimits());

#ifdef DEBUG_OUTPUT
  cout<<"nLB()"<<endl;
#endif
  motifCode tfid=this->getMotif();

  for(seq=this->m-1;seq>0;seq--) {
    //this->matrix_p[seq]--;
    this->matrixIndexDec(seq);
    //this->diff2limiter[seq]=this->difference(*this->limiterPvec,seq,tfid);


    if(this->matrix_p[seq]>=0 && this->diff2limiter[seq]<this->limitBP) {
      assert(this->diff2limiter[seq]>0);
      return;
    } else {
      // Mentiin liian kauas rajoittajasta tai alun etupuolelle

      this->matrixIndexSet(seq,this->myMat->CoordCacheGet(tfid,seq));
      //this->matrix_p[seq]=this->myMat->CoordCacheGet(tfid,seq);
      //this->diff2limiter[seq]=this->difference(*this->limiterPvec,seq,tfid);
      assert(this->ok && this->checkAtBorder(seq));

    }
  }

   if(seq==0) {
     while(!this->decFirst() && this->isOK())
       assert(this->difference2limiter(seq)>=0);
   }
   


#ifdef DEBUG_OUTPUT
   if(ok) {
     motifCode tfid=this->getMotif();

     assert(this->matrix_p[0]<this->myMat->dims(0));
     printf("dist[0]= %d\n",this->difference2limiter(0));
     for(uint i=1;i<this->m;i++) {
       cout<<this->matrix_p[i]<<"<<"<<this->myMat->countTFinSeq(i,tfid)<<" "<<this->limData->factor(tfid)<<endl;
       assert(this->matrix_p[i]<this->myMat->countTFinSeq(i,tfid));
     }
   }
#endif
  assert(!this->ok || this->checkWithinLimits());


}

// BasicPointerVec::BasicPointerVec(class PointerVec &p)
// {
//   this->matrix_p=p.matrix_p;
//   this->ok=p.ok;
// }


bool PointerVec::operator<=(const PointerVec &other) const
{
  bool lt=1;

  for(uint i=0;lt && i<this->m;i++) {
    //    lt=lt&&(this->matrix_p[i]<=other.matrix_p[i]);
    lt=lt&&(this->getSite(i).epos<=other.getSite(i).pos);
  }
  return lt;
}



const PointerVec& PointerVec::operator++(int dummy)
{
  seqCode seq;

  for(seq=this->m-1;seq>0;seq--) {
    int sites=this->myMat->countTFinSeq(seq,this->getMotif());
    if(this->matrix_p[seq]<(sites-1)) {
      //this->matrix_p[seq]++;
      this->matrixIndexInc(seq);
      break;
    } else {
      this->matrixIndexSet(seq,0);
      //this->matrix_p[seq]=0;
    }
  }
  if(seq==0) {
    do {
      this->matrixIndexInc(0);
      //this->matrix_p[0]++;
      if(this->matrix_p[0]<this->myMat->dims(0)) {
	this->resetMotifCode();
      } else {
	break;
      } 
    } while(!this->allHasFactor());
    
    if(this->matrix_p[0]>=this->myMat->dims(0)) {
      this->matrixIndexSet(0,-1);
      //this->matrix_p[0]=-1; // THIS INVALIDATES THIS PointerVec!!!!
      this->ok=0;
    }
  }


#ifdef DEBUG_OUTPUT
  else if(ok) { // Need condition to avoid messing up with loop in getOrigin()
    motifCode tfid=this->getMotif();
    assert(this->matrix_p[0]<this->myMat->dims(0));
    for(uint i=1;i<this->m;i++) {
      cout<<this->matrix_p[i]<<"<"<<this->myMat->countTFinSeq(i,tfid)<<" "<<this->limData->factor(tfid)<<endl;
      assert(this->matrix_p[i]<this->myMat->countTFinSeq(i,tfid));
    }
  }
#endif

  return *this;
}

// Recursively free the allocated memory
void Matrix::freeData(vector<matrixentry> *d)
{
  delete d;
  d=NULL;
}

// Recursively free the allocated memory
void Matrix::freeData(vector<void*> *d,int dimI)
{
  if(dimI==(this->dim-1)) {
    freeData((vector<matrixentry>*)d);
  } else if(dimI<(this->dim-1)) {
    for(unsigned int i=0;i<d->size();i++) {
      freeData((vector<void*> *)d->at(i),dimI+1);
      d->at(i)=NULL;
    }
    if(dimI>0) {
      delete d;
      d=NULL;
    }
  } else {
    cerr<<"It's just a flesh wound!"<<endl;
    abort();
  }

}

Matrix::~Matrix()
{
  freeData(&this->data,0);
}

// Return either *vector<void*>  or *vector<matrixentry> for the last level
void * Matrix::allocateData(seqCode level,motifCode mot)
{
  unsigned int const TFcount=this->countTFinSeq(level,mot);


  if(level<(this->dim-1)) {
    vector<void*> *curDat=new vector<void*>;

    curDat->reserve(TFcount);


    for(uint i=0;i<TFcount;i++) {
      curDat->push_back(this->allocateData(level+1,mot));
    }
    this->matSize+=sizeof(void*);
    return (void*)curDat;
  } else if(level==(this->dim-1)) {
    vector<matrixentry> *curDat=new vector<matrixentry>;
    this->cells+=TFcount;

    curDat->resize(TFcount);

    this->matSize+=sizeof(matrixentry)*curDat->size();
    return (void*)curDat;
  }

  return (void*) NULL; // SHOULDN'T COME HERE!!
}


void Matrix::initAllHaveFactor()
{
  for(motifCode tfID=0;tfID<this->indata->factors();tfID++) {
    int sites=1;
    for(seqCode i=1;(unsigned int)i<(unsigned int)this->dim;i++) {
      sites*=this->countTFinSeq(i,tfID);
    }
    this->allHaveFactor.push_back((sites>0));
  }
    
}

Matrix::Matrix(Inputs *indata)
{
  this->indata=indata;

  this->dim=indata->sequences();
  this->dimLen=indata->sequenceLens();


  this->matSize=0;
  


  try {
    // Allocate the helpers for the sparse matrix

    // TODO: Get rid of TFs that do not occur in some of the sequences.
    this->tfIndex.resize(indata->factors());
    for(motifCode tf=0;tf<indata->factors();tf++) {
      this->tfIndex[tf].resize(indata->sequences());
    }
    for(seqCode seq=0;seq<indata->sequences();seq++) {
#ifndef NDEBUG
      cerr<<seq<<":"<<this->dimLen[seq]<<endl;
#endif
      for(posCode pos=0;pos<this->dimLen[seq];pos++) {
#ifndef NDEBUG
	id_triple const &debugSite=indata->getSite(pos,seq);
	motifCode debugTF=debugSite.ID;
	this->tfIndex[debugTF][seq].push_back(pos);
#else
	this->tfIndex[indata->getSite(pos,seq).ID][seq].push_back(pos);
#endif
      }
    }

    
#ifdef DEBUG_OUTPUT
    for(int i=0;i<indata->factors();i++) {
      int motifCount=countTFinSeq(0,i);
      cout<<i<<" "<<indata->factor(i)<<" "<<motifCount<<flush;
      for(int j=1;j<indata->sequences();j++) {
	//assert(motifCount==countTFinSeq(j,i));
	motifCount=countTFinSeq(j,i);
	cout<<" "<<motifCount<<flush;
      }
      cout<<endl;
    }
#endif

    // Allocate the sparse matrix itself

    // First dimension straight away.
    this->data.reserve(this->dimLen[0]);
    this->cells=0; //this->dimLen[0];

    for(posCode basePos=0;basePos<this->dimLen[0];basePos++) {
      motifCode baseMot=indata->getSite(basePos,0).ID;
      this->data.push_back(this->allocateData(1,baseMot));
    }
  }
  catch (std::bad_alloc const&) {
    PyErr_SetString(PyExc_MemoryError,"Out of memory!");
  }

  this->CoordCacheInit();
  this->initAllHaveFactor();
  
}

PointerVec Matrix::getOrigin()
{
  PointerVec orig=PointerVec(this,this->indata);

  return orig;
}

BasicPointerVec Matrix::argMax()
{
  BasicPointerVec maxstore;
  const BasicPointerVec *backer;
  store maxval=0.0;

  for(PointerVec p=this->getOrigin();p.isOK();p++) {
    assert( p.getSite(0).ID==p.getSite(1).ID );
    /*assert( (!p.isOK()) 
	    || (this->getValue(p).getValue()>=0.0)
	    || ( (backer=this->getValue(p).getBacktraceP())==NULL) 
	    || (this->getValue(*backer).getValue()<0.0) ); */ 
    if(this->getValue(p).getValue()>maxval) {
      for(backer=this->getValue(p).getBacktraceP();
	  backer->isOK() && this->getValue(*backer).getValue()>=0.0;
	  backer=this->getValue(*backer).getBacktraceP()) {
	/* pass*/
      }
      if(!backer->isOK()) { // Store only if we do not backtrack to negative.
	maxstore=p;
	maxval=this->getValueP(p)->getValue();
      }
    }
  }
#ifdef DEBUG_OUTPUT
  if(!maxstore.isOK()) {
    cout<<"Busted matrix, all zeros"<<endl;
  }
#endif
    
  return maxstore;
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
    inSite=NULL;
  }
  Py_DECREF(inpIter);
  inpIter=NULL;
  
  map<string,motifCode>::iterator seqName=sequenceNames();
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
  motifCode tf_id;
  seqCode seq_id;


  if(n<0) {
    PyErr_SetString(PyExc_TypeError,"A non sequence site data");
    return 0;
  } else if(n<7) {
    PyErr_SetString(PyExc_ValueError,"Not enough fields in the site data");
    return 0;
  }

  
 
  // Sequence name and id
  PyObject *tmp_obj;
  if(!(tmp_obj=PySequence_GetItem(site,0))) {
    PyErr_SetString(PyExc_ValueError,"No sequence name!!");
    return 0;
  }
  string seqName=string(PyString_AsString(tmp_obj));
  Py_DECREF(tmp_obj);
  tmp_obj=NULL;
  if(PyErr_Occurred()) return 0;
  
  map<string,motifCode>::iterator seq_iter=SEQ_to_id.find(seqName);

  if(seq_iter==SEQ_to_id.end()) {
    seqNames.push_back(seqName);
    seq_id=SEQ_to_id.size();
    SEQ_to_id[seqName]=seq_id;
    seq.resize(seq_id+1);
  } else {
    seq_id=seq_iter->second;
  }

  // Transcription factor name name and id
  if(!(tmp_obj=PySequence_GetItem(site,2))) {
    PyErr_SetString(PyExc_ValueError,"No TF name!!");
    return 0;
  }
  string TFname=string(PyString_AsString(tmp_obj));
  Py_DECREF(tmp_obj);
  tmp_obj=NULL;
  if(PyErr_Occurred()) return 0;


  
  map<string,motifCode>::iterator id_iter=TF_to_id.find(TFname);

  if(id_iter==TF_to_id.end()) {
    factorNames.push_back(TFname);
    tf_id=TF_to_id.size();
    TF_to_id[TFname]=tf_id;
  } else {
    tf_id=id_iter->second;
  }


  // Create a struct for the new site
  id_triple new_site;


  if( (tmp_obj=PyNumber_Int(PySequence_GetItem(site,3))) ) {
    new_site.pos=(posind)PyInt_AsLong(tmp_obj);
    Py_DECREF(tmp_obj);
    tmp_obj=NULL;
  }
  if(PyErr_Occurred()) {
    PyObject *errStr=PyObject_Str(PySequence_GetItem(site,4));
    PyErr_Format(PyExc_TypeError,"Invalid format for start position '%s' (%s:%d)",PyString_AsString(errStr),__FILE__,__LINE__);
    Py_DECREF(errStr);
    errStr=NULL;
    return 0;
  }

  if( (tmp_obj=PyNumber_Int(PySequence_GetItem(site,4))) ) {
    new_site.epos=(posind)PyInt_AsLong(tmp_obj);
    Py_DECREF(tmp_obj);
    tmp_obj=NULL;
  }

  if(PyErr_Occurred()) {
    PyObject *errStr=PyObject_Str(PySequence_GetItem(site,4));
    PyErr_Format(PyExc_TypeError,"Invalid format for end position '%s' (%s:%d)",PyString_AsString(errStr),__FILE__,__LINE__);
    Py_DECREF(errStr);
    errStr=NULL;
    return 0;
  }

  if( (tmp_obj=PyNumber_Float(PySequence_GetItem(site,5))) ) {
    new_site.weight=(double)PyFloat_AsDouble(tmp_obj);
    Py_DECREF(tmp_obj);
    tmp_obj=NULL;
  }


  if(PyErr_Occurred()) {
    PyObject *errStr=PyObject_Str(PySequence_GetItem(site,5));
    PyErr_Format(PyExc_TypeError,"Invalid format for weight '%s' (%s:%d)",PyString_AsString(errStr),__FILE__,__LINE__);
    Py_DECREF(errStr);
    errStr=NULL;
    return 0;
  }

  PyObject *annotSlice=PySequence_GetSlice(site,8,99999);
  PyObject *sep=PyString_FromString("\t");
  PyObject *jStr=_PyString_Join(sep,annotSlice);
  new_site.annot=string(PyString_AsString(jStr));

  Py_DECREF(jStr);
  jStr=NULL;
  Py_DECREF(sep);
  sep=NULL;
  Py_DECREF(annotSlice);
  annotSlice=NULL;


			

  if(PyErr_Occurred()) {
    //cout<<seqName<<" "<<TFname<<" "<<new_site.pos<<endl;
    return 0;
  }

  if(!(tmp_obj=PySequence_GetItem(site,6))) {
    PyErr_SetString(PyExc_ValueError,"No strand!!");
    return 0;
  }
  char *strand_p=PyString_AsString(tmp_obj);
  Py_DECREF(tmp_obj);
  tmp_obj=NULL;
  if(PyErr_Occurred()) return 0;
  new_site.strand=(char)strand_p[0];

  new_site.ID=2*tf_id+(new_site.strand=='+'?0:1);


  // Finally add the new site to the input object
  seq[seq_id].push_back(new_site);

  return 1;
}


vector<id_triple> Inputs::getSites(PointerVec const &p) const
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





// matrixentry::matrixentry(store v,BasicPointerVec &p):value(v),backTrack2(p)
// {
// }

// matrixentry::matrixentry(store v,BasicPointerVec *p):value(v),backTrack2(*p)
// {
// }




void PointerVec::output() const 
{ 
  if(this->isOK()) {
   // cout<<"valid("<<m<<"): ";
    cout<<"(";
    for(uint i=0;i<m;i++) {
      cout<<this->operator[](i);
      if(i<(m-1)) {
	 cout<<",";
      } else {
	cout<<"["<<this->limData->factor(this->getMotif())<<this->getSite(0).strand<<"])";
      }
    } 
  } else { cout<<"Invalid)"; }
  //  cout<<endl; 
  cout<<flush;
}


void PointerVec::setValue(vector<int> &np)
{
  // TODO: Something for this. Remove or fix.
  //p=np;
  abort();

  ok=0;
}

// Assist functions.
//store PointerVec::getValue() const { return myMat->getValueP(*this)->getValue(); }
vector<id_triple>  PointerVec::getSites() const { return limData->getSites(*this); }



//motifCode const PointerVec::getMotifLaborous() const {  return this->limData->getSite(this->matrix_p[0],0).ID; }
//motifCode const PointerVec::getMotif() const {  return this->curMotifCode; }
//motifCode const PointerVec::getMotif() const {  return this->limData->getSite(this->matrix_p[0],0).ID; }




void Matrix::CoordCacheInit() {
  this->coordUpdateCache.clear();
  this->coordUpdateCache.resize(this->indata->sequences());
  for(seqCode seq=0;(unsigned int)seq<this->coordUpdateCache.size();seq++) {
    this->coordUpdateCache[seq].resize(this->indata->factors(),0);
  }
}

PointerVec PointerVec::getLimited(int limitbp) const
{
  PointerVec ret=PointerVec(*this);


  ret.limitBP=limitbp;

  ret.limiterPvec=this;
  ret.diff2limiter.resize(this->m,0);

  for(seqCode i=1;(unsigned int)i<(unsigned int)this->m;i++) {
    if(ret.myMat->CoordCacheGet(ret.getMotif(),i)>ret.matrix_p[i]) {
#ifdef DEBUG_OUTPUT
      cout<<"Reseting "<<i<<" at "<<ret.myMat->CoordCacheGet(ret.getMotif(),i)<<endl;
      ret.output();
      cout<<endl;
#endif
      ret.myMat->CoordCacheReset(i);
      break;
      
    }
  }

  while(ret.isOK() && !ret.decFirst());

  if(!ret.isOK()) {
    goto finally;
  }

  

  for(seqCode i=1;(unsigned int)i<(unsigned int)this->m;i++) {
    while(ret.isOK() && ret.difference2limiter(i)<0) {
      //ret.matrix_p[i]--;
      ret.matrixIndexDec(i);
      if(ret.matrix_p[i]<0 || ret.difference2limiter(i)>=ret.limitBP) {
	ret.ok=0;
	goto finally;
      }
    }
    ret.myMat->CoordCacheSet(ret.getMotif(),i,ret.matrix_p[i]);
  }

//   if(ok) {
//     if(!ret.checkAtBorder()) {
//       cout<<"Not at border:";
//       ret.output();
//     }
//     if(!ret.checkWithinLimits()) {
//       cout<<"Not within limits:";
//       ret.output();
//       ret.limiterPvec->output();
//       for(uint i=0;i<this->m;i++) {
// 	cout<<ret.difference(*ret.limiterPvec,i)<<",";
//       }
//       cout<<endl;
//     }
//   }
 finally:
  assert(!ret.ok || ret.checkWithinLimits());
  assert(!ret.ok || ret.checkAtBorder());


  return ret;
}

// PointerVec PointerVec::getLimited(int limitbp)
// {
//   PointerVec ret=PointerVec(*this);
//   ret.setLimit(limitbp);

//   assert(!ret.ok || ret.checkWithinLimits());


//   return ret;
// }




// Return true, if the sites are the same in all dimensions.
int PointerVec::allSame() 
{
  motifCode ID=this->getSite(0).ID;

  for(uint i=1;i<m;i++) {
    if(this->getSite(i).ID!=ID) {
      return 0;
    }
  }
  return 1;
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
    self->CP=NULL;

    if(self->names) {
      Py_DECREF(self->names);
      self->names=NULL;
    }
    if(self->bestAlignments) {
      Py_DECREF(self->bestAlignments);
      self->bestAlignments=NULL;
    }

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

  //cout<<"NYT AJETAAN INIT_I�"<<endl;
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
      //CHECKING_DECREF(self);
      return -1;
    }
  self->alpha=0;
  self->beta=0;
  self->Rsquared=-1.0;
  self->RMSE=0.0;

  malign_alignCommon(self,data,result_ask,lambda,xi,mu,nu,nuc_per_rotation);

  if(PyErr_Occurred()) {
    //Py_XDECREF(self);
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
    cout<<"CP ON EP�NULL"<<endl<<flush;
    return (PyObject *)self;
}










static PyMethodDef malignment_methods[] = {
      {"nextBest", (PyCFunction)malignment_nextBest, METH_NOARGS,
       "Return the next best alignment from the matrix."},
      {"suboptimal", (PyCFunction)malignment_nextBest, METH_NOARGS,
       "Return the next best alignment from the matrix. TODO:: Fix this to Waterman-Eggert like"
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
    {"alpha",T_DOUBLE, offsetof(malign_AlignmentObject,alpha), 0,
     "Intercept of the expectation model."},
    {"beta",T_DOUBLE, offsetof(malign_AlignmentObject,beta), 0,
     "Slope of the expectation model."},
    {"Rsquared",T_DOUBLE, offsetof(malign_AlignmentObject,Rsquared), 0,
     "Coefficient of determination of the expectation model"},
    {"RMSE",T_DOUBLE, offsetof(malign_AlignmentObject,RMSE), 0,
     "Root mean squared error of the expectation model."},
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

  const BasicPointerVec *p;
  BasicPointerVec start=self->CP->dynmat->argMax();

  PyObject *goodAlign=PyList_New(0);
  //  PyObject *colTuple;
  PyObject *seqs;
  PyObject *coords;
  PyObject *seqPos,*siteScore,*annotations;
  int seqC=0;

#ifdef DEBUG_OUTPUT
  cout<<"eka:"<<flush;
  //start.output();
#endif

  matrixentry *matVal=NULL;


  for(p=&start;p && p->isOK();p=matVal->getBacktraceP()) {

    matVal=self->CP->dynmat->getValueP(*p);

    
    PointerVec pp=PointerVec(*p,self->CP->dynmat,self->CP->indata);
    char strand=pp.getSite(0).strand;
    siteCode motifID=pp.getSite(0).ID;

    // Tuple for one column of the alignment
    //    colTuple=PyTuple_New(5);


    // Tuple for sequence codes
    seqs=PyTuple_New(self->CP->indata->sequences());
    // Tuple for coordinates
    coords=PyTuple_New(self->CP->indata->sequences());

    seqPos=PyTuple_New(self->CP->indata->sequences());
    siteScore=PyTuple_New(self->CP->indata->sequences());
    annotations=PyTuple_New(self->CP->indata->sequences());




    seqC=0;
    vector<id_triple>  sites=self->CP->indata->getSites(pp);

    for(uint i=0;i<sites.size();i++) {
      if((siteCode)sites[i].ID==motifID && sites[i].strand==strand) {
	// Sequence
	PyTuple_SetItem(seqs,seqC,PyInt_FromLong(i));


	// Coordinates on the sequence
	PyObject *coord;
	coord=Py_BuildValue("(ii)",(int)sites[i].pos,(int)sites[i].epos);
	PyTuple_SetItem(coords,seqC,coord);

	// Position on the site sequence
	PyTuple_SetItem(seqPos,seqC,PyInt_FromLong(pp[i]));

	PyTuple_SetItem(siteScore,seqC,PyFloat_FromDouble(sites[i].weight));

	// Annotation
	PyTuple_SetItem(annotations,seqC,PyString_FromString(sites[i].annot.c_str()));
	seqC++; 
      }
    }

    // Resizing might change the location of the tuples
    _PyTuple_Resize(&seqs,seqC);
    _PyTuple_Resize(&coords,seqC);
    _PyTuple_Resize(&seqPos,seqC);
    _PyTuple_Resize(&siteScore,seqC);
    _PyTuple_Resize(&annotations,seqC);

//     PyTuple_SetItem(colTuple,0,seqs);
//     // The alignment score up till here
//     PyTuple_SetItem(colTuple,1,PyFloat_FromDouble(self->CP->dynmat->getValue(*p).getValue()));
//     PyTuple_SetItem(colTuple,2,PyInt_FromLong(motifID));
//     PyTuple_SetItem(colTuple,3,coords);
//     PyTuple_SetItem(colTuple,4,PyString_FromStringAndSize(&strand,1));

    
    PyObject *alnRow=PyAln_New_Multi(self->CP->indata->factor(motifID),
				    seqs,seqPos,coords,strand,
				    self->CP->dynmat->getValue(pp).getValue(),
				    siteScore,annotations);

    PyList_Append(goodAlign,alnRow);
    Py_DECREF(alnRow);
    alnRow=NULL;


    matVal->negate();

#ifdef DEBUG_OUTPUT
    cout<<"Kierret��n :"<<flush;
    pp.output();
    cout<<"="<<self->CP->dynmat->getValue(pp).getValue()<<endl;
#endif

  }

  return goodAlign;
}



// inline double anglepenalty(double const d,double const D,double const nucl_per_rotation)
// {
//   double const theta=(d-D)*2.0*PI/nucl_per_rotation;
//   double const totlen=d+D;

//   if(((PI*PI)/totlen)<FLT_EPSILON) { // Check that we can compute it.
//     return 0.0;
//   } else {
//     return squaremodpi(theta)/totlen;
//   }
// }



inline const double penalty(const malign_AlignmentObject* const dat,const posind d1,const posind d2)
{
  double val=dat->mu*(d1+d2)/2.0;

  assert(d1>=0);
#ifndef NDEBUG
  if(d2<0) {printf("d2=%g\n",d2);}
#endif
  assert(d2>=0);
  if( (d1+d2)>0.0) {
    val=val+ dat->nu*(d1-d2)*(d1-d2)/(d1+d2)+
      dat->xi*anglepenalty(d1,d2,dat->nuc_per_rotation);
  }
  return val;
}

// Compute the sum of pairs penalty for adjacent columns 
// left and rigtht with n rows (sequences)
inline const double multiPenaltyLR(const malign_AlignmentObject* const self,const PointerVec &left,const PointerVec &right,const int n) 
{
  double Score=0.0;

  const motifCode leftTF=left.getMotif();
  const motifCode rightTF=right.getMotif();

  for(seqCode i=0;i<n;i++) {
    for(seqCode j=i+1;j<n;j++) {
      Score+=penalty(self,left.difference(right,i,leftTF,rightTF),left.difference(right,j,leftTF,rightTF));
    }
  }
  return Score;
}

// Compute the sum of pairs penalty for adjacent columns 
// left and left->limiterPvec with n rows (sequences)
inline const double multiPenalty(const malign_AlignmentObject* const self,const PointerVec &left,const PointerVec &right,const int n) 
{
  double Score=0.0;

  for(seqCode i=0;i<n;i++) {
    for(seqCode j=i+1;j<n;j++) {
      Score+=penalty(self,
		     left.difference2limiter(i),
		     left.difference2limiter(j));
    }
  }
  return Score;
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
malignObject(malign_AlignmentObject * const self)
{ 
  store Score=0.0,absBest=0.0;
  store maxScore=SCORE_FAIL;
  // The main multiple alignment algorithm.

  int Empty=0,nonEmpty=0;
#ifndef NDEBUG
  int unNecessary=0;
#endif


  // Iterate over the whole multi D matrix with the pointer entry.
  PointerVec entry=self->CP->dynmat->getOrigin();

  int n=self->CP->dynmat->dims();

  for(;entry.isOK();entry++) {
#ifdef DEBUG_OUTPUT
    entry.output();
    assert(entry.getSite(0).ID==entry.getSite(1).ID);
    assert(entry.getSite(0).strand==entry.getSite(1).strand);
    
    int a=entry[0];
    int b=entry[1];
    if(a==b) {
      cout<<"."<<flush; 
    } else {
      cout<<"x"<<flush;
    }
#endif

    //vector<id_triple> sites=entry.getSites();


    //siteCode maxSite=SITE_FAIL;
    //char maxStrand=STRAND_FAIL;

    BasicPointerVec maxSource;
    maxScore=-DBL_MAX;
    //const motifCode entryTF=entry.getMotif().

    store base=0;

    // Compute the fixed bonus part of the penalty function
    for(int i=0;i<n;i++) {
      const id_triple &site=entry.getSite(i);
      base+=site.weight;
    }
    base*=self->lambda*n*(n-1)/2;



    maxScore=base;  // On our own. The first site of the alignment.


    // Look for the previous alignments
    //PointerVec k=PointerVec(entry);


    for(PointerVec k=entry.getLimited(MAX_BP_DIST);
	k.isOK(); k.nextLookBack()) {  // Exponential loop
	 
      assert(k.checkWithinLimits());
#ifdef DEBUG_OUTPUT
      k.output();cout<<"s ";
    cout<<"dists: (";
    for(int i=0;i<n;i++) {
      if(i>0) cout<<",";
      cout<<k.difference(entry,i);
    }
    cout<<")"<<endl;
#endif
	 
      // Continuing alignment.
#ifndef NDEBUG
      if((*self->CP->dynmat)[k].getValue()==0.0) {
	unNecessary++;
      }
#endif
      
      // Add the score from the lookup table
      //Score=base+k.getValue();

      // Compute the negative penalty
      //for(int i=0;i<n;i++) {
      //	for(int j=i+1;j<n;j++) {
      //	  Score=Score - penalty(self,k.difference(entry,i),k.difference(entry,j));
      //	}
      //}
      Score=base+k.getValue()-multiPenalty(self,k,entry,n);
      assert(k.getValue()>=0.0);
      //printf("k.ok=%d  maxSource.ok=%d maxScore=%g Score=%g k.value=%g\n",k.isOK(),maxSource.isOK(),maxScore,Score,k.getValue());
      assert(maxScore>=base);
      if(Score>maxScore) {
	maxScore=Score;
	maxSource=BasicPointerVec(k);
	//if(!k.isOK() || !maxSource.isOK()) {
	//printf("k.ok=%d  maxSource.ok=%d maxScore=%g\n",k.isOK(),maxSource.isOK(),maxScore);
	  //}
	assert(k.isOK()==maxSource.isOK());
      }

      if(PyErr_Occurred()) return NULL;
    }
    if(maxScore>absBest)
      absBest=maxScore;
    if(maxScore<0.0)
      printf("BUSTED\n");
    assert(maxScore>0.0);
    if(maxScore>0.0) {
#ifdef DEBUG_OUTPUT
      //entry.output();
      cout<<"="<<maxScore<<"<-";
      //maxSource.output();
      cout<<endl;
#endif
      //printf("maxscore OK %d %g  > %g > %g\n",maxSource.isOK(),Score,maxScore,base);
      (*self->CP->dynmat)[entry]=matrixentry(maxScore,maxSource);
      //printf("value=%g  %g\n\n",(*self->CP->dynmat)[entry].getValue(),entry.getValue());
      nonEmpty++;
    } else {
      (*self->CP->dynmat)[entry]=matrixentry();
      Empty++;
    }
    maxSource=BasicPointerVec();
  }
#ifdef DEBUG_OUTPUT
  cout<<"Matrix dimensions: ";
  for(int i=0;i<self->CP->dynmat->dims();i++) {
    cout<<" "<<self->CP->dynmat->dims(i);
  }
  cout<<endl
      <<"Matrix size: "<<self->CP->dynmat->maxSize()<<endl
      <<"Used Cells:  "<<self->CP->dynmat->usedCells()<<endl
      <<"Max Value:   "<<absBest<<endl
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
  map<string,motifCode>::iterator iter=data.sequenceNames();

  self->names=PyTuple_New(n);

  if(PyErr_Occurred()) return 0;
  for(int i=0;i<n;i++,iter++) {
    PyTuple_SET_ITEM(self->names,(int)iter->second,
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


  if(self->CP->indata->sequences()<2) {
    PyErr_SetString(PyExc_EOFError,"Too few sequences in input");
    return NULL;
  }




  vector<int> dimensions=self->CP->indata->sequenceLens();
  self->CP->dynmat=new Matrix(self->CP->indata);

  if(PyErr_Occurred()) return NULL;
  

  self->lambda=lambda;
  self->xi=xi;
  self->mu=mu;
  self->nu=nu;
  self->nuc_per_rotation=nuc_per_rotation;
  self->askedresults=result_ask;
  self->memSaveUsed=0;
  self->secs_to_align=0;

#ifdef DEBUG_OUTPUT
  cout
    <<"Lambda:"<<lambda<<endl
    <<"Xi    :"<<xi<<endl
    <<"Mu    :"<<mu<<endl
    <<"Nu    :"<<nu<<endl
    <<"npr   :"<<nuc_per_rotation<<endl
    <<"result:"<<result_ask<<endl;
#endif

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

  // SUDDEN DEATH!!! THIS IS JUNK FUNCTION OR METHOD
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

    
    m=Py_InitModule("eellib.multiAlign", malignMethods);

    if(m==NULL)
      return;

    if(import_alnCols()<0)
      return;

    Py_INCREF(&malign_AlignmentType);
    PyModule_AddObject(m, "MultiAlignment", (PyObject *)&malign_AlignmentType);


#ifdef DEBUG_OUTPUT
    cout<<"MultiAlignLoaded"<<endl;
#endif
}

