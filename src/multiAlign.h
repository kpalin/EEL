#ifndef _MULTIALIGN_H__

#define _MULTIALIGN_H__

/*
 *
 *$Log$
 *
 */

#ifndef NDEBUG
#ifndef DEBUG_OUTPUT
#define DEBUG_OUTPUT
#endif
#endif

typedef unsigned int uint;

typedef double posind;

typedef double store;

struct id_triple
{
  uint ID;
  posind pos;
  posind epos;
  //uint pos;
  double weight;
  char strand;
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
  int dataPoint();

public:
  PointerVec() {ok=0; limData=NULL; }
  PointerVec(vector<int> &,vector<int> &,vector<int>*);
  void setValue(vector<int> &np);

  int isOK() { return ok; }
  vector<int> &getValue() { return p; }
  const PointerVec&  operator++(int);
  const PointerVec & operator--(int);
  PointerVec project(vector<int> &);
  PointerVec projectAndLimit(vector<int> &,Inputs *,int maxbp);
  posind difference(int i);
  void clearProject();
  int operator[](int i) {return p[i]; }

  int getMatrixCoord() { assert(matrix_p==dataPoint()); return matrix_p; }

  void output() {
    for(uint i=0;i<m;i++) {
      cout<<p[i]<<",";
    } cout<<endl; }
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

public:
  Inputs() { return; }
  Inputs(PyObject* data);
  int addSite(PyObject* item);
  int sequences() { return SEQ_to_id.size(); }
  int factors() { return TF_to_id.size(); }
  vector<int> sequenceLens();

  vector<id_triple>  getSites(PointerVec &here);

  id_triple getSite(PointerVec &p,int i) { return seq[i][p[i]]; }
  id_triple getSite(int pos,int i) { return seq[i][pos]; }
  ABset getAB(PointerVec &p,int limit=2);
};

class matrixentry
{
  store value;

  //following values are for the backtracing
  PointerVec backTrack;
public:
  matrixentry() { return ;}
  matrixentry(store,PointerVec&);
  store getValue();
  PointerVec &getBacktrace();
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
  matrixentry getValue(PointerVec &);
  void setValue(PointerVec &,matrixentry &);
  matrixentry &operator[](PointerVec &p) { return data[p.getMatrixCoord()]; }
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
  PyTupleObject *names;

  PyObject *bestAlignments;

  double secs_to_align;
  double lambda;
  double xi; 
  double mu;
  double nu;
  double nuc_per_rotation;
  int askedresults;

  int item_count;

  struct __CPSTUF *CP;

} malign_AlignmentObject;


static PyObject *
malign_alignCommon(  malign_AlignmentObject *self,PyObject *data,int result_ask,
		   double lambda,double xi,double mu,double nu,
		   double nuc_per_rotation);


#endif
