#ifndef _MULTIALIGN_H__

#define _MULTIALIGN_H__

/*
 *
 *$Log$
 *Revision 1.17  2006-12-08 09:49:56  kpalin
 *Fixed a seqfault and added a source file from Pasi rastas which got
 *included in _c_matrix.cc
 *
 *Revision 1.16  2006/08/31 10:09:34  kpalin
 *Minimum remembered pairwise alignments and speed improvements for multi align.
 *
 *Revision 1.15  2006/08/29 06:08:28  kpalin
 *Speed improvements.
 *
 *Revision 1.14  2006/08/25 12:18:31  kpalin
 *Caching distances to the pointervec limiter.
 *
 *Revision 1.13  2006/05/12 06:51:05  kpalin
 *Close to perfection.
 *
 *Revision 1.12  2006/05/11 12:45:05  kpalin
 *Works. But leaks a bit of memory in matrixentry.
 *
 *Revision 1.11  2006/05/04 07:58:53  kpalin
 *Speed improvements.
 *
 *Revision 1.10  2006/05/03 11:00:17  kpalin
 *Merged. Merged a fix with the main devel. branch.
 *
 *Revision 1.9  2006/05/03 10:11:22  kpalin
 *Fixed a bug resulting wrong alignments depending on order of the input
 *sequences.
 *
 *Revision 1.8  2006/01/05 09:22:09  kpalin
 *Fast and working version.
 *
 *Revision 1.7  2006/01/04 12:29:39  kpalin
 *Almost compatible with multiALign.cc:1.18
 *
 *Revision 1.6  2006/01/02 12:44:16  kpalin
 *Cleaning up limit coordinate caching.
 *
 *Revision 1.5  2005/11/24 13:29:45  kpalin
 *Possibly irrelevant cast to int added.
 *
 *Revision 1.4  2005/10/03 10:18:41  kpalin
 *Presumably working multiple alignment version. Not yet usable though.
 *
 *Revision 1.3  2005/03/22 13:24:02  kpalin
 *Totally different approach. Doing multi-D matrix.
 *
 *Revision 1.2  2004/07/30 12:06:55  kpalin
 *First basically working version. The problem is the immense
 *time and space complexity.
 *
 *Revision 1.1  2004/07/12 12:48:18  kpalin
 *Class definitions and some other standard header stuff
 *
 *
 */

/*
#ifndef NDEBUG
#ifndef DEBUG_OUTPUT
#define DEBUG_OUTPUT
#endif
#endif
*/
#include <algorithm>
#include <functional>
#include <math.h>
#include <float.h>
#include <limits.h>

typedef unsigned int uint;

typedef double posind;

typedef double store;
typedef signed short int siteCode;



#define SCORE_FAIL -DBL_MAX
#define SITE_FAIL -1
#define STRAND_FAIL 'f'

typedef int motifCode;
typedef int seqCode;
typedef int posCode;



struct id_triple
{
  motifCode ID;
  posind pos;
  posind epos;
  //uint pos;
  double weight;
  char strand;
  string annot;
  bool operator<(const id_triple &other) const { 
    return (int)this->epos<(int)other.epos || 
      ((int)this->epos==(int)other.epos && this->weight<other.weight) ||
      ((int)this->epos==(int)other.epos && this->weight==other.weight && this->strand<other.strand);
  }

  //Sanna's addition
  bool operator==(const id_triple &other) const { 
    return ((int)this->ID==(int)other.ID && (int)this->pos==(int)other.pos && (int)this->epos==(int)other.epos && this->weight==other.weight && this->strand==other.strand);
  }
};

class BasicPointerVec {
 protected:
 public:
  vector<int> matrix_p;
  bool ok;
  BasicPointerVec(): ok(0) { };
    BasicPointerVec(class PointerVec &p);
    int matrixIndex(seqCode const i) const {return this->matrix_p[i];}
    int isOK() const { return this && this->ok; }
};


class Matrix;

class PointerVec: public BasicPointerVec {
  //vector<int> p;
  //vector<int> matrix_p;
  //vector<int> dimlen;
  //static 
  uint m;

  //static 
  class Inputs const * limData;
  //static 
  int limitBP;
  int dataPoint() const ;

  // Distances from this to the limiter. Used as a cache.
  vector<int> diff2limiter;

  // Highest point of the limited lookback.
  const class PointerVec *limiterPvec;
  //static const class PointerVec *limiterPvec;
  // Pointer to the matrix this pointer is associated with:
  //static 
  class Matrix * myMat;

  

  motifCode curMotifCode;
  bool decFirst();
  bool updateRestCoords();
  const int allHasFactor() const;

  void setPrevMatrixCoord(motifCode const tfID,seqCode const dimension); 


  // Set the matrix index for sequence seq to value val
  int matrixIndexSet(seqCode seq,int val) {
    this->matrix_p[seq]=val;
    if(val>=0 && this->limiterPvec && this->isOK()) {
      this->diff2limiter[seq]=this->difference(*this->limiterPvec,seq);
    }
    return this->matrix_p[seq];
  }

  // Decrease matrix index for sequence seq by 1
  int matrixIndexDec(seqCode seq) {
    return this->matrixIndexSet(seq,this->matrix_p[seq]-1);
  }
  // Increase matrix index for sequence seq by 1
  int matrixIndexInc(seqCode seq) {
    return this->matrixIndexSet(seq,this->matrix_p[seq]+1);
  }
  
  // Update the internal motif code cache.
  void resetMotifCode() {
    this->curMotifCode=this->getMotifLaborous();
  }

public:
  PointerVec() {ok=0; limiterPvec=NULL;limData=NULL; }
  //  PointerVec(const PointerVec &other);
  //PointerVec(vector<int> &,vector<int> &,vector<int>*);
  PointerVec(Matrix* mat,Inputs const *inp);

  PointerVec(const BasicPointerVec &p,class Matrix * mat,class Inputs  *inp);
  void setValue(vector<int> &np);


  PointerVec getLimited(int limitbp) const;
  int allSame();
  //vector<int> &getValue() { return p; }

  

  // Jump to the "next" entry in the matrix that has the same site on each dimension
  const PointerVec&  operator++(int);
  bool operator<=(const PointerVec &other) const; 
  bool checkLT() const ;
  bool checkAtBorder();
  bool checkAtBorder(seqCode i);
  bool checkWithinLimits() const ;

  int difference2limiter(seqCode const i) const {
    // Return the distance FROM the end of this to the beginning of limiterPvec
    assert(this->limiterPvec!=NULL);
    assert(this->diff2limiter[i]==this->difference(*this->limiterPvec,i));
    return this->diff2limiter[i];
  }

  // Return the distance FROM the end of this TO the beginning of other on sequence i when this is pointer to motifs of type thisTFid and other is pointer to motifs of type otherTFid
  int difference(PointerVec const &other,seqCode const i,motifCode const thisTFid, motifCode const otherTFid) const {
    assert(limData!=NULL);
    assert(thisTFid==this->getMotif());
    int otherRight=(int)(other.getSite(i,otherTFid).pos-this->getSite(i,thisTFid).epos);

    //assert((otherRight-1)==this->difference(other,i,thisTFid));
    return otherRight-1; //max(otherRight,otherLeft);
  } ;


  // Return the distance FROM the end of this TO the beginning of other on sequence i when this is pointer to motifs of type thisTFid
  int difference(PointerVec const &other,seqCode const i,motifCode const thisTFid) const {
#ifndef NDEBUG
    assert(limData!=NULL);
    assert(thisTFid==this->getMotif());

    id_triple otherSite=other.getSite(i);
    id_triple thisSite=this->getSite(i,thisTFid);

    int ero=otherSite.pos;
    ero-=thisSite.epos;
    int otherRight=(int)(other.getSite(i).pos-this->getSite(i,thisTFid).epos);
    //int otherLeft=(int)(this->getSite(i).pos-other.getSite(i).epos);
    //return otherRight-1; //max(otherRight,otherLeft);
    //assert((otherRight-1)==this->difference(other,i,thisTFid,other.getMotif()));
#endif
    return this->difference(other,i,thisTFid,other.getMotif());
  } ;

  // Return the distance FROM the end of this TO the beginning of other on sequence i
  int difference(PointerVec const &other,seqCode const i) const  {
    return this->difference(other,i,this->getMotif());
  }
  int operator[](seqCode const i) const;

  void nextLookBack();

  void output() const;
  


  // Assist functions:
  store getValue()  const; // inline

  motifCode const getMotif() const { 
#ifndef NDEBUG
    const motifCode thisMotif=this->getMotifLaborous();
    if(this->curMotifCode!=thisMotif) {
      fprintf(stderr,"VÄÄRÄ MOTIFCODE!!\n");
      //abort();
    }
#endif
    return this->curMotifCode;
  }
  motifCode const getMotifLaborous() const; // inline

  vector<id_triple>  getSites() const;
  id_triple const &getSite(seqCode const i) const;
  id_triple const &getSite(seqCode const  i,motifCode const tfID) const;

};



class Inputs {
  vector< vector<id_triple> > seq;
  map<string,motifCode> TF_to_id;
  map<string,seqCode> SEQ_to_id;
  vector<string> seqNames;
  vector<string> factorNames;


public:
  Inputs() { return; }
  Inputs(PyObject* data);
  int addSite(PyObject* item);
  int sequences() const  { return SEQ_to_id.size(); }
  int factors() const { return TF_to_id.size()*2; }
  char const *sequence(seqCode i) const {return seqNames[(int)i].c_str(); }
  char const *factor(motifCode i) const {return factorNames[((int)i)/2].c_str(); }
  vector<int> sequenceLens();
  int sequenceLens(seqCode seqC) const {return seq[seqC].size(); }


  map<string,seqCode>::iterator sequenceNames() { return SEQ_to_id.begin(); }
  map<string,motifCode>::iterator transfacNames() { return TF_to_id.begin(); }

  vector<id_triple>  getSites(PointerVec const &here) const;

  id_triple const &getSite(PointerVec const &p,seqCode i) const { 
#if 1
    posCode my_p=p[i];
    vector<id_triple> my_seq=seq.at(i);
    return my_seq.at(my_p); 
#else
    return seq[i][p[i]];
#endif
  }

  id_triple const &getSite(posCode const pos,seqCode const i) const { 
#if 1
    vector<id_triple> const &my_seq=seq.at(i);
    return my_seq.at(pos);
#else
    return seq[i][pos]; 
#endif
  }
};





class matrixentry
{
  store value;

  //following values are for the backtracing
  BasicPointerVec backTrack2;
public:
  matrixentry() { value=-1.0;  }
  matrixentry(store v,BasicPointerVec &p):value(v),backTrack2(p) { };
  matrixentry(store v,BasicPointerVec *p):value(v),backTrack2(*p) { } ;
  store getValue()  const {
    return this->value;
  }
  const BasicPointerVec *getBacktraceP() const {
    return &backTrack2; }

  void negate() { value*=-1.0; }

};





/* class matrixentry */
/* { */
/*   store value; */

/*   //following values are for the backtracing */
/*   BasicPointerVec backTrack; */
  
/*   //void setInitData(store,BasicPointerVec*); */
/*  public: */
/*   matrixentry(): value(0.0) { }  */
/*     matrixentry(store v ,BasicPointerVec &p): value(0.0),backTrack(p) { assert(p.isOK()); printf("value=%g\n",v); }  */
/*       //matrixentry(store,BasicPointerVec*); */
/*       store getValue()  const  { return value;} */
/*       BasicPointerVec *getBacktraceP() {return &backTrack; } */
      
/*       void negate() { value*=-1.0; } */
      
/* }; */


class Matrix {

  int dim,matSize;
  vector<int> dimLen;

  vector<int> dimFactors;

  vector<void*> data;

  Inputs const *  indata;

  int cells;

  vector<bool> allHaveFactor;

  // p=tfIndex[tfID][seqID][posCode] <=>  indata->seq[seqID][p] is the posCode:th 
  //                                      occurrence of motif tfID in sequence seqID
  vector< vector< vector<int> > > tfIndex;

  void* allocateData(seqCode level,motifCode mot);

  vector< vector<int> > coordUpdateCache;
  void initAllHaveFactor();

  void freeData(vector<void*> *d,int dimI);
  void freeData(vector<matrixentry> *d);
public:
  Matrix() {return; }
  Matrix(Inputs* indata);

  ~Matrix(); 
  PointerVec getOrigin();
  BasicPointerVec argMax();

  void CoordCacheInit();
  int CoordCacheGet(motifCode const tfID,seqCode const i) const{
    return this->coordUpdateCache[i][tfID];
    //return 0;
  }


  const bool allHasFactor(motifCode const tfID) const {
    return this->allHaveFactor[tfID];
  }
  void CoordCacheSet(motifCode const tfID,seqCode const i,int value) {
    //assert(value>=this->coordUpdateCache[tfID][i]);
    this->coordUpdateCache[i][tfID]=value;
  }

  void CoordCacheReset(seqCode const i) {
    this->coordUpdateCache[i].clear();
    this->coordUpdateCache[i].resize(this->indata->factors(),0);
  }

  int usedCells() {return this->cells;}

  matrixentry &getValue(BasicPointerVec const &p) const { return *this->getValueP(p);};
  matrixentry *getValueP(BasicPointerVec const &) const;
  void setValue(BasicPointerVec &p,matrixentry &val) { *(this->getValueP(p))=val; };

  matrixentry &operator[](BasicPointerVec &p) { return this->getValue(p); };

  int size() { return data.size(); }
  int maxSize() { int s=1;for(int i=0;i<dim;i++) s*=dimLen[i]; return s; }
  int dims() const {return dim; }
  int dims(int i) const { return dimLen[i]; }
  const int countTFinSeq(seqCode const seq,motifCode const tfID) const { 
    return this->tfIndex[tfID][seq].size(); }


  const int by_seq_tf_pos(int const seqID, motifCode const tfID,posCode const pos) const { 
#ifndef NDEBUG
    vector< vector<int> > seqTF=this->tfIndex.at(tfID);
    vector<int>  TFid=seqTF.at(seqID);
    int tfCode = TFid.at(pos);
    return tfCode;
#else
    return this->tfIndex[tfID][seqID][pos]; 
#endif
  }

};



struct __CPSTUF {
  __CPSTUF() { 
    indata=NULL;
    dynmat=NULL;
    return; }

  ~__CPSTUF() {
    if(indata) {
      delete indata;
    }
    if(dynmat) {
      delete dynmat;
    }
  }
  Inputs *indata;

  Matrix *dynmat;



  friend class PointerVec;
  //priority_queue<MS_res> bestAligns;

  //map<store,pair<matCoord,matCoord> > bestAlignsTmp;
}  ;

typedef struct {
  PyObject_HEAD

    /* Type-specific fields go here. */
  PyObject *names;  //PyTupleObject 

  PyObject *bestAlignments;

  double secs_to_align;
  double lambda;
  double xi; 
  double mu;
  double nu;
  double nuc_per_rotation;
  int askedresults;
  int memSaveUsed;
  // Expected value model
  // ln(Evalue)=alpha+beta*score
  double alpha;
  double beta;
  double Rsquared;
  double RMSE;

  int item_count;

  struct __CPSTUF *CP;

} malign_AlignmentObject;


static PyObject *
malign_alignCommon(  malign_AlignmentObject *self,PyObject *data,int result_ask,
		   double lambda,double xi,double mu,double nu,
		   double nuc_per_rotation);

static PyObject*
malignment_nextBest(malign_AlignmentObject *self);

inline BasicPointerVec::BasicPointerVec(class PointerVec &p):matrix_p(p.matrix_p), ok(p.ok) 
{ 
}

inline store PointerVec::getValue()  const { return this->myMat->getValueP(*this)->getValue(); } ;
inline motifCode const PointerVec::getMotifLaborous() const {  
#if 1
  posCode p=this->matrix_p.at(0);
  id_triple const &trip= this->limData->getSite(p,0);
    
  return trip.ID; 

#else
  return this->limData->getSite(this->matrix_p[0],0).ID; 
#endif
};


inline matrixentry *Matrix::getValueP(BasicPointerVec const &p) const
{
  vector<void*> const *d=&this->data;
  //void *d=&this->data;

  for(int dimI=0;dimI<(this->dim-1);dimI++) {
    int mat_ind=p.matrixIndex(dimI);
    //d=((vector<void*>*)d)->at(mat_ind);
    //d=(vector<void*>*)d[mat_ind];
    d=(vector<void*>*)d->at(mat_ind);
  }

  matrixentry *ret=&((vector<matrixentry>*)d)->at(p.matrixIndex(this->dim-1));

  return ret;
}

inline int PointerVec::operator[](seqCode const i) const {
  int out;
  if(i==0) {
    out=this->matrix_p[i];
  } else{
    out=this->myMat->by_seq_tf_pos(i,this->getMotif(),this->matrix_p[i]); 
  }
  return out;
} 


inline id_triple const &PointerVec::getSite(seqCode const i) const 
{ 
  assert(this->isOK());
  if(i==0) {
    return limData->getSite(this->matrix_p[i],i);
  } else {
    return this->getSite(i,this->getMotif());
  }
}

inline id_triple const &PointerVec::getSite(seqCode const i,motifCode const tfID) const 
{ 
  // !!! DOES NOT CHECK FOR CORRECT tfID 

  assert(tfID==this->getMotif());

  assert(this->isOK());
  if(i==0) {
    return limData->getSite(this->matrix_p[i],i);
  } else {
    return limData->getSite(myMat->by_seq_tf_pos(i,tfID,this->matrix_p[i]),i); 
  }
}

inline PointerVec::PointerVec(const BasicPointerVec &p,class Matrix * mat,class Inputs  *inp):BasicPointerVec(p),limData(inp),myMat(mat)
{
  this->m=mat->dims();
  this->resetMotifCode();
}


#endif
