#ifndef _MULTIALIGN_H__

#define _MULTIALIGN_H__

/*
 *
 *$Log$
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

#ifndef NDEBUG
#ifndef DEBUG_OUTPUT
#define DEBUG_OUTPUT
#endif
#endif

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
    return (int)this->epos<(int)other.epos || ((int)this->epos==(int)other.epos && this->weight<other.weight);
  }
};

class PointerVec {
  //vector<int> p;
  vector<int> matrix_p;
  //vector<int> dimlen;
  uint m;
  int ok;

  class Inputs *limData;
  int limitBP;
  int dataPoint() const ;

  // Highest point of the limited lookback.
  class PointerVec *limiterPvec;
  // Pointer to the matrix this pointer is associated with:
  class Matrix* myMat;

  

  bool decFirst();
  bool updateRestCoords();
  int allHasFactor();

  int getPrevMatrixCoord(motifCode const tfID,seqCode const dimension); 

public:
  PointerVec() {ok=0; limiterPvec=NULL;limData=NULL; }
  //  PointerVec(const PointerVec &other);
  //PointerVec(vector<int> &,vector<int> &,vector<int>*);
  PointerVec(Matrix* mat,Inputs* inp);
  void setValue(vector<int> &np);

  void setLimit(PointerVec &p,int limitbp);
  int isOK() const { return this && ok; }
  int allSame();
  //vector<int> &getValue() { return p; }


  // Jump to the "next" entry in the matrix that has the same site on each dimension
  const PointerVec&  operator++(int);
  bool operator<=(const PointerVec &other) const; 
  bool checkLT();
  bool checkAtBorder();
  bool checkAtBorder(seqCode i);


  int difference(PointerVec const &other,seqCode const i,motifCode const tfID) const ;
  int difference(PointerVec const &other,seqCode const i) const  {
    return this->difference(other,i,this->getMotif());
  }
  int operator[](seqCode const i) const;
  int matrixIndex(seqCode const i) const {return this->matrix_p[i];}

  void nextLookBack();

  void output();
  


  // Assist functions:
  store getValue()  const;

  vector<id_triple>  getSites();
  id_triple getSite(seqCode const i) const;
  id_triple getSite(seqCode const  i,motifCode const tfID) const;
  motifCode const getMotif() const;

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
  int sequences() { return SEQ_to_id.size(); }
  int factors() { return TF_to_id.size()*2; }
  char const *sequence(seqCode i) {return seqNames[(int)i].c_str(); }
  char const *factor(motifCode i) {return factorNames[((int)i)/2].c_str(); }
  vector<int> sequenceLens();
  int sequenceLens(seqCode seqC) {return seq[seqC].size(); }


  map<string,seqCode>::iterator sequenceNames() { return SEQ_to_id.begin(); }
  map<string,motifCode>::iterator transfacNames() { return TF_to_id.begin(); }

  vector<id_triple>  getSites(PointerVec &here);

  id_triple getSite(PointerVec const &p,seqCode i) const { return seq[i][p[i]]; }
  id_triple getSite(posCode pos,seqCode i) const { return seq[i].at(pos); }
};








class matrixentry
{
  store value;

  //following values are for the backtracing
  PointerVec *backTrack;

  void setInitData(store,PointerVec*);
public:
  matrixentry() { value=-1.0;backTrack=(PointerVec*)2;}
  matrixentry(store,PointerVec&);
  matrixentry(store,PointerVec*);
  store getValue()  const;
  PointerVec *getBacktraceP() const{return backTrack; }

  void negate() { value*=-1.0; }

};


class Matrix {

  int dim,matSize;
  vector<int> dimLen;

  vector<int> dimFactors;

  vector<void*> data;

  Inputs *indata;

  int cells;

  // p=tfIndex[seqID][tfID][posCode] <=>  indata->seq[seqID][p] is the posCode:th 
  //                                      occurrence of motif tfID in sequence seqID
  vector< vector< vector<int> > > tfIndex;

  void* allocateData(seqCode level,motifCode mot);

  vector< vector<int> > coordUpdateCache;
  void CoordCacheInit();


public:
  Matrix() {return; }
  Matrix(Inputs* indata);

  PointerVec getOrigin();
  PointerVec argMax();

  int CoordCacheGet(motifCode const tfID,seqCode const i) const{
    return this->coordUpdateCache[tfID][i];
  }

  void CoordCacheSet(motifCode const tfID,seqCode const i,int value) {
    //assert(value>=this->coordUpdateCache[tfID][i]);
    this->coordUpdateCache[tfID][i]=value;
  }

  int usedCells() {return this->cells;}

  matrixentry &getValue(PointerVec const &) const;
  matrixentry *getValueP(PointerVec const &) const;
  void setValue(PointerVec &,matrixentry &);

  matrixentry &operator[](PointerVec &p);

  int size() { return data.size(); }
  int maxSize() { int s=1;for(int i=0;i<dim;i++) s*=dimLen[i]; return s; }
  int dims() const {return dim; }
  int dims(int i) const { return dimLen[i]; }
  int countTFinSeq(seqCode seq,motifCode tfID) const { 
    return this->tfIndex[seq][tfID].size(); }


  int by_seq_tf_pos(int seqID, motifCode tfID,posCode pos) const { 
#ifndef NDEBUG
    vector< vector<int> > seqTF=this->tfIndex.at(seqID);
    vector<int>  TFid=seqTF.at(tfID);
    int tfCode = TFid.at(pos);
    return tfCode;
#else
    return this->tfIndex[seqID][tfID][pos]; 
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

  int item_count;

  struct __CPSTUF *CP;

} malign_AlignmentObject;


static PyObject *
malign_alignCommon(  malign_AlignmentObject *self,PyObject *data,int result_ask,
		   double lambda,double xi,double mu,double nu,
		   double nuc_per_rotation);

static PyObject*
malignment_nextBest(malign_AlignmentObject *self);

#endif
