#ifndef _MULTIALIGN_H__

#define _MULTIALIGN_H__

/*
 *
 *$Log$
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

struct id_triple
{
  uint ID;
  posind pos;
  posind epos;
  //uint pos;
  double weight;
  char strand;
  bool operator<(const id_triple &other) const { return this->pos<other.pos; }
};

class PointerVec {
  vector<int> p;
  vector<int> dimlen;
  vector<char> bound;  //Bit vector
  uint m;
  int ok;

  vector<int> *dimFactors;
  int matrix_p;

  class Inputs *limData;
  int limitBP;
  int dataPoint() const ;

public:
  PointerVec() {ok=0; limData=NULL; }
  PointerVec(vector<int> &,vector<int> &,vector<int>*);
  void setValue(vector<int> &np);

  int isOK() const { return this && ok; }
  vector<int> &getValue() { return p; }

  const PointerVec&  operator++(int);
  const PointerVec & operator--(int);


  PointerVec project(vector<int> &);
  PointerVec projectAndLimit(vector<int> &,Inputs *,int maxbp);
  posind difference(int i) const ;
  void clearProject();
  int operator[](int i) const {return p[i]; }

  vector<char> &getBound() { return bound; }
  int getMatrixCoord() const { assert(matrix_p==dataPoint()); return matrix_p; }

  void output();
  
};

class ABset {
  vector<int> A; // Set of characters over limit
  vector<vector<int> > B; // Set of sequences for given character
  int limit;
  int pointer;

public:
  ABset(int dim,int limit=2) { B.resize(dim*2); this->limit=limit; pointer=0; }
  void addChar(int seq,int chr,char strand);
  int moreAlphas();
  vector<int> nextB();
};









class Inputs {
  vector< vector<id_triple> > seq;
  map<string,uint> TF_to_id;
  map<string,uint> SEQ_to_id;
  vector<string> seqNames;
  vector<string> factorNames;
public:
  Inputs() { return; }
  Inputs(PyObject* data);
  int addSite(PyObject* item);
  int sequences() { return SEQ_to_id.size(); }
  int factors() { return TF_to_id.size(); }
  char const *sequence(int i) {return seqNames[i].c_str(); }
  char const *factor(int i) {return factorNames[i].c_str(); }
  vector<int> sequenceLens();

  map<string,uint>::iterator sequenceNames() { return SEQ_to_id.begin(); }
  map<string,uint>::iterator transfacNames() { return TF_to_id.begin(); }

  vector<id_triple>  getSites(PointerVec &here);

  id_triple getSite(PointerVec const &p,int i) const { return seq[i][p[i]]; }
  id_triple getSite(int pos,int i) const { return seq[i][pos]; }
  ABset getAB(PointerVec &p,int limit=2);
};








class matrixentry
{
  store value;
  // Site and strand code.
  siteCode siteEtStrand;


  //following values are for the backtracing
  PointerVec *backTrack;

  void setInitData(store,int,char,PointerVec*);
public:
  matrixentry() { value=0.0;siteEtStrand=SITE_FAIL;backTrack=NULL;}
  matrixentry(store,int,char,PointerVec&);
  matrixentry(store,int,char,PointerVec*);
  store getValue()  const;
  int getSite() const;
  char getStrand() const ;
  PointerVec *getBacktraceP() const{return backTrack; }

  void negate() { value*=-1.0; }

};


class Matrix {

  int dim,matSize;
  vector<int> dimLen;

  vector<int> dimFactors;

  vector<matrixentry> data;


public:
  Matrix() {return; }
  Matrix(vector<int> &);

  PointerVec getOrigin();
  PointerVec argMax();

  matrixentry const &getValue(PointerVec const &) const ;
  matrixentry *getValueP(PointerVec const &);
  void setValue(PointerVec &,matrixentry &);

  matrixentry &operator[](PointerVec &p);

  int size() { return data.size(); }
  int maxSize() { int s=1;for(int i=0;i<dim;i++) s*=dimLen[i]; return s; }
  int dims() {return dim; }
  int dims(int i) { return dimLen[i]; }
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
