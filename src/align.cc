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

#include <time.h>
//#include <unistd.h>



#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <limits>



#ifdef HAVE_GZSTREAM
#include <gzstream.h>
#endif


#include "alignedCols.h"


/*
 *
 *  $Log$
 *  Revision 1.18  2005/01/25 09:29:29  kpalin
 *  Removed usage of long double to make it compile on MacOSX
 *
 *  Revision 1.17  2005/01/13 13:16:50  kpalin
 *  Moved the requesting of sequences to be aligned to Python side
 *  of stuff. Much better.
 *
 *  Revision 1.16  2005/01/12 13:35:21  kpalin
 *  Allowing #-comments in GFF files.
 *
 *  Revision 1.15  2005/01/07 11:58:46  kpalin
 *  Windows port. Should work with Visual C++ 7.1
 *
 *  Revision 1.14  2004/12/22 11:14:34  kpalin
 *  Some fixes for better distributability
 *
 *  Revision 1.13  2004/12/14 14:07:22  kpalin
 *  *** empty log message ***
 *
 *  Revision 1.12  2004/08/26 10:39:47  kpalin
 *  Improved documentation and proper return value from aligndata if
 *  incorrect number of arguments. i.e. if ParseTuple fails.
 *
 *  Revision 1.11  2004/07/30 12:07:50  kpalin
 *  Output changed to use alnColumn objects.
 *
 *  Revision 1.10  2004/02/26 11:27:36  kpalin
 *  Changed the output from nextBest() function to be new SitePair objects
 *  instead a large tuple of numbers. It's easier to change in the future
 *  and easier to use on python side.
 *
 *  Revision 1.9  2004/02/13 09:40:25  kpalin
 *  Corrected bugs
 *
 *  Revision 1.8  2004/02/11 09:39:41  kpalin
 *  Enabled memory saving features.
 *  Breaks 'more' command but saves a lot of memory.
 *
 *  First find the good matches and find the alignment
 *  only when asked. Increases the running time but drops
 *  the memory consuption. SAVE_MEM_LIMIT is threshold
 *  for memory saving feature.
 *
 *  Revision 1.7  2004/02/04 10:29:46  kpalin
 *  Corrected memory usage reporting for large inputs.
 *
 *  Revision 1.6  2004/01/14 10:06:31  kpalin
 *  Checkings for decref:ing Py_None
 *
 *  Revision 1.4  2003/12/30 11:21:36  kpalin
 *  Added a conditional dependency on gzstream and zlib to allow processing of
 *  gziped files.
 *
 *  Revision 1.3  2003/12/29 12:43:32  kpalin
 *  Ilmeisesti jotain uutta. En tiedä mitä.
 *
 *
 */


#define CHECKING_DECREF(X) if((PyObject*)(X)==Py_None) printf("none decref line %d",__LINE__); Py_DECREF(X);


/* The binding sites have to be less than MAX_BP_DIST apart */
#ifndef MAX_BP_DIST
#define MAX_BP_DIST 1000
#endif

#ifndef OUTPUTFREQ
#define OUTPUTFREQ 30000
#endif

#undef SMALLX   /* Define SMALLX if you want smaller sequence to be x. 
		 This (as of 22.10.2003) speeds up the mehtod a little*/
#define SMALLX 1

// #define SAVE_MEM 1

#ifndef SAVE_MEM_LIMIT
#define SAVE_MEM_LIMIT 1024*1024*1024  // Memory limit after which we use (slower) memory frugal algorithm.
#endif


#ifdef HAS_NO_TRUNC

double trunc(double x)
{
  return (x>0?1:-1)*floor(fabs(x));
}

#endif


#ifdef SAVE_MEM        // Under memory constraints, define SAVE_MEM. Only a tiny speed loss.
typedef double store;
#else
typedef double store;
#endif



#undef _ORIG_ALIGN_
//#undef NDEBUG
#include <assert.h>


using namespace std;

typedef unsigned int uint;

class triple
{
public:
  string TF;
  uint pos;
  uint epos;
  double weight;
  char strand;
  string annot;
//   ~triple();  
};

// triple::~triple()
// {
//   if(this->annot) {
//     free(this->annot);
//     this->annot=NULL;
//   }
// }
struct id_triple
{
  string TF;
  uint ID;
  double pos;
  double epos;
  //uint pos;
  double weight;
  char strand;
  string annot;
  //  ~id_triple();
};

// id_triple::~id_triple()
// {
//   if(this->annot) {
//     free(this->annot);
//     this->annot=NULL;
//   }
// }
// Needed for sort
bool operator<(const triple& t1, const triple& t2){
  return t1.epos<t2.epos || (t1.epos==t2.epos && t1.weight<t2.weight);
}

typedef map<string, vector<triple> > matchlisttype;

struct matrixentry
{
  store value;

  //following values are for the backtracing
  int x;
  int y;
};













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

inline double anglepenalty(double const d,double const D,double const nucl_per_rotation)
{
  double const theta=(d-D)*2.0*PI/nucl_per_rotation;
  return squaremod(theta)/(d+D);
}






//splits a string at special key word
vector<string> str_split(const string& str, const string& delimiters = " ")
{
    vector<string> tokens;
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
}

// Dirty trick. Requires references.
int parseLine(string const &line,triple &tri,string &sequence)
{
  string TF,junk;
  int pos,epos;
  double weight;
  char strand;
  istringstream help(line);
  if(line[0]=='#') {
    return 2;
  }
  //if (! (help >> sequence))
  if( ! getline(help,sequence,'\t'))
    {
#ifdef ALIGN_OUTPUT
      cerr<<"wrong format in sequence part:"<<endl
	  <<line<<endl;
#endif
      return 0;
    }

  //  if (! (help >> TF) || ! (help >> TF))
  if (! getline(help,TF,'\t') || ! getline(help,TF,'\t'))
    {
#ifdef ALIGN_OUTPUT
      cerr<<"wrong format in feature part:"<<endl
	  <<line<<endl;
#endif
      return 0;
      }
  if (! (help >> pos))
    {
#ifdef ALIGN_OUTPUT
      cerr<<"wrong format in start position part:"<<endl
	  <<line<<endl;
#endif
      return 0;
    }
  if (! (help >> epos))
    {
#ifdef ALIGN_OUTPUT
      cerr<<"wrong format in end position part:"<<endl
	  <<line<<endl;
#endif
      return 0;
      }
  if ( ! (help>>weight))
    {
#ifdef ALIGN_OUTPUT
      cerr<<"wrong format in score part:"<<endl
	  <<line<<endl;
#endif
      return 0;
    }

  if ( ! (help>>strand)) 
    {
#ifdef ALIGN_OUTPUT
      cerr<<"wrong format in strand part:"<<endl
	  <<line<<endl;
#endif
      return 0;
    }
  tri.TF = TF;
  tri.pos=pos;
  tri.weight=weight;
  tri.epos=epos;
  tri.strand=strand;
  

  // Take spacer
  
  help>>junk;

  // Get annotation
  string ends;
  getline(help,ends,'\n');
  if(ends.length()>2) {
    tri.annot=ends;
  } else {
    tri.annot=string();
  }
  return 1;
}




// this function parses a file in gff format
matchlisttype* parseStream(istream *in)
{
  matchlisttype* matchlist=new matchlisttype;
  if (!*in) 
    {
      #ifdef ALIGN_OUTPUT
      cerr << "Error while parsing GFF data"<< endl;
      #endif
      return matchlist;      
    }


  string line,sequence;
  triple tri;
  while(getline(*in, line))
  {
    if(parseLine(line,tri,sequence)==1) {
      (*matchlist)[sequence].push_back(tri);
    }
  }

  //cout<<"The whole stream parsed successfully!"<<endl;

  // the lists must be sorted
  matchlisttype::iterator iter;
  for(iter=matchlist->begin(); iter!= matchlist->end(); iter++){
    sort(iter->second.begin(), iter->second.end());
  }
  return matchlist;
}



#ifdef ALIGN_OUTPUT
static PyObject *
align_draw(PyObject *self, PyObject *args)
{ 
  char* file;
  if (!PyArg_ParseTuple(args, "s", &file)){
    Py_INCREF(Py_None);
    return Py_None;
  }

  ifstream inData(file);

  //cout<<"filename: "<<file<<endl;

  matchlisttype* ml=parseStream(&inData);
  matchlisttype::iterator iter;
  for(iter=ml->begin(); iter!= ml->end(); iter++)
    {
      cout<<iter->first<<":"<<endl;
      for(uint i=0; i<iter->second.size(); i++)
	{
	  cout<< (iter->second)[i].TF <<", "<<(iter->second)[i].pos
	      <<", "<<(iter->second)[i].weight<<endl;
	}
    }
  
  delete ml;
      
  Py_INCREF(Py_None);
  return Py_None;

}
#endif



inline int indexBeforeOrAtBp(vector<int> const ind,vector<id_triple> const seq, double const pos)
{
  // Much worse version of indexing. Assymptotically n*log n whereas n suffices
  // new version run for 1500x1370 matrix took 6.0 CPU seconds
  // for the old version, the run took 132 CPU seconds

  int a=0,b=ind.size()-1,m, ret=-2;


#ifdef EXTRADEBUG
  cout<<"Target:"<<pos<<endl;
#endif
  while((b-a)>1) {
    m=(a+b)/2;
    if(seq[ind[m]].epos>=pos) {
      b=m;
#ifdef EXTRADEBUG
      cout<<"Going b="<<m<<":"<<seq[ind[m]].epos<<" is larger"<<endl;
#endif
    } else {
#ifdef EXTRADEBUG
      cout<<"Going a="<<m<<":"<<seq[ind[m]].epos<<" is smaller"<<endl;
#endif
      a=m;
    }
  } 

  if(b<0) // Empty list
    ret=-1;
  else if(seq[ind[b]].epos<pos)
    ret=b;
  else if (seq[ind[a]].epos<pos)
    ret=a;
  else if(a==0 && seq[ind[a]].epos>=pos) // no match before
    ret=-1;


#ifdef EXTRADEBUG
  cout<<"Returning: "<<ret;
  if(ret>=0)
    cout<<" ind= "<<ind[ret]<<" epos="<<seq[ind[ret]].epos;
  else if(ret==-2){
    cout<<" ="<<a<<" ~ "<<ind[a]<<" ~ "<<seq[ind[a]].epos;
    cout<<" b="<<b<<" ~ "<<ind[b]<<" ~ "<<seq[ind[b]].epos;
  }
  cout<<endl;
#endif

  assert(ret==-1 || seq[ind[ret]].epos<pos);
      
  return ret;
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


struct __CPSTUF {
  vector<id_triple> seq_x;
  vector<id_triple> seq_y;
  vector<string> ID_to_TF;
  vector<vector<int> > index;
  vector<vector<matrixentry> > matrix;

  vector<int> prevYindex;
  vector<double> prevPos;

  priority_queue<MS_res> bestAligns;

  //map<store,pair<matCoord,matCoord> > bestAlignsTmp;
}  ;

typedef struct {
  PyObject_HEAD

    /* Type-specific fields go here. */
  PyStringObject *x_name,*y_name;
  PyObject *names;
  PyStringObject *output;
  PyObject *bestAlignments;

  double secs_to_align;
  double lambda;
  double xi;
  double mu;
  double nu;
  double nuc_per_rotation;
  int askedresults;
  int mem_usage;
  int item_count;
  unsigned int expectedMemUsage;
  double fill_factor;
  int memSaveUsed;

  struct __CPSTUF *CP;

} align_AlignmentObject;




static 
void
alignment_dealloc(align_AlignmentObject* self)
{
    CHECKING_DECREF(self->x_name);
    CHECKING_DECREF(self->y_name);
    CHECKING_DECREF(self->output);
    CHECKING_DECREF(self->bestAlignments);
    
    delete self->CP;
    

    self->ob_type->tp_free((PyObject*)self);
}



extern "C" int
alignment_init(align_AlignmentObject *self, PyObject *args, PyObject *kwds)
{

    double secs=0.0;
    PyObject  *x_name=NULL,*y_name=NULL,*output=NULL;

    static char *kwlist[] = {"x_name","y_name","secs_to_align","output", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|SSdS", kwlist,
				      &x_name,&y_name,&secs,&output))
        return -1; 

    self->secs_to_align = secs;

    if (x_name) {
        CHECKING_DECREF(self->x_name);
        Py_INCREF(x_name);
        self->x_name = (PyStringObject*)x_name;
    }

    if (y_name) {
        CHECKING_DECREF(self->y_name);
        Py_INCREF(y_name);
        self->y_name = (PyStringObject*)y_name;
    }


    Py_INCREF(x_name);
    Py_INCREF(y_name);
    self->names=Py_BuildValue("(SS)",x_name,y_name);

    if (output) {
        CHECKING_DECREF(self->output);
        Py_INCREF(output);
        self->output = (PyStringObject*)output;
    }

    

    self->nu=1.0;
    self->mu=1.0;
    self->lambda=0.5;
    self->xi=1.0;
    self->expectedMemUsage=0;
    return 0;
}

static 
PyObject *
alignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    align_AlignmentObject *self;

    self = (align_AlignmentObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
      self->secs_to_align = 0.0;
      self->nu=0.0;
      self->x_name = (PyStringObject*)PyString_FromString("");
      if (self->x_name == NULL)
	{
	  CHECKING_DECREF(self);
	  return NULL;
	}
      self->y_name = (PyStringObject*)PyString_FromString("");
      if (self->y_name == NULL)
	{
	  CHECKING_DECREF(self);
	  return NULL;
	}
      self->output = (PyStringObject*)PyString_FromString("");
      if (self->output == NULL)
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
      if (self->CP==NULL) {
	return NULL;
      }
      self->memSaveUsed=0;

    }

    self->item_count=0;
    alignment_init(self,args,kwds);

    return (PyObject *)self;
}



double findMax(int &pos_x,int &pos_y,align_AlignmentObject *self)
{
  // Return maximum element that does not have a link to negative matrix position
  store max=-1.0;
  pos_x=-1;
  pos_y=-1;
 
  for (uint x=0; x< self->CP->matrix.size(); x++) {
    for(uint y=0; y<self->CP->matrix[x].size(); y++) {

      if(self->CP->matrix[x][y].value>0) {
	int tmp_x=self->CP->matrix[x][y].x;
	int tmp_y=self->CP->matrix[x][y].y;
	if(tmp_x>=0 && tmp_y>=0 && self->CP->matrix[tmp_x][tmp_y].value<0.0) {
	  self->CP->matrix[x][y].value*=-1.0;
	} else if (self->CP->matrix[x][y].value > max) {
	  max= self->CP->matrix[x][y].value;
	  pos_x=x;
	  pos_y=y;
	}
      }
    }
  }

  return max;
}







static PyObject*
alignment_nextBest(align_AlignmentObject *self);


static PyMethodDef alignment_methods[] = {
    {"nextBest", (PyCFunction)alignment_nextBest, METH_NOARGS,
     "Return the next best alignment from the matrix."
    },
    {NULL}  /* Sentinel */
};

static PyMemberDef alignment_members[] = {
    {"bestAlignments",T_OBJECT_EX, offsetof(align_AlignmentObject, bestAlignments), 0,
     "List of best alignments"},
    {"secs_to_align",T_DOUBLE, offsetof(align_AlignmentObject, secs_to_align), 0,
     "CPU time for Alignment"},
    
    {"names",T_OBJECT_EX, offsetof(align_AlignmentObject, names), 0,
     "Names of the sequences, x=0, y=1"},
    
    {"x_name",T_OBJECT_EX, offsetof(align_AlignmentObject, x_name), 0,
     "Name of the X string"},
    {"y_name",T_OBJECT_EX, offsetof(align_AlignmentObject, y_name), 0,
     "Name of the Y string"},
    {"output",T_OBJECT_EX, offsetof(align_AlignmentObject, output), 0,
     "The string output from align. Going away"},
    {"Lambda",T_DOUBLE, offsetof(align_AlignmentObject,lambda), 0,
     "Parameter for bonus"},
    {"Xi",T_DOUBLE, offsetof(align_AlignmentObject,xi), 0,
     "Parameter for distance penalty"},
    {"Nu",T_DOUBLE, offsetof(align_AlignmentObject,nu), 0,
     "Parameter for distance difference penalty"},
    {"Mu",T_DOUBLE, offsetof(align_AlignmentObject,mu), 0,
     "Parameter for rotation penalty"},
    {"nuc_per_rotation",T_DOUBLE, offsetof(align_AlignmentObject,nuc_per_rotation), 0,
     "Parameter for nucleotides per 360 deg rotation of DNA"},
    {"fill_factor",T_DOUBLE, offsetof(align_AlignmentObject,fill_factor), 0,
     "Factor of filled cells out of all possible"},
    {"mem_usage",T_INT, offsetof(align_AlignmentObject,mem_usage), 0,
     "Mem usage of alignment matrix in bytes."},
    {"memSaveUsed",T_INT, offsetof(align_AlignmentObject,memSaveUsed), 0,
     "True if alignment used low memory feature."},
    {"item_count",T_INT, offsetof(align_AlignmentObject,item_count), 0,
     "Number of filled cells."},
    {"askedResults",T_INT, offsetof(align_AlignmentObject,askedresults), 0,
     "Number of filled cells."},
    {NULL}  /* Sentinel */
};


static PyTypeObject align_AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "align.Alignment",             /*tp_name*/
    sizeof(align_AlignmentObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)alignment_dealloc,                         /*tp_dealloc*/
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
    "Alignment object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    alignment_methods,             /* tp_methods */
    alignment_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)alignment_init,      /* tp_init */
    0,                         /* tp_alloc */
    alignment_new,                 /* tp_new */
};






//////////////////////////////////////////////////////////////////////


inline int indexYBeforeOrAtBp(int xx,struct __CPSTUF *CP, double const pos)
{
  // MUCH better version of indexBeforeOrAtBp
  // run for 1500x1370 matrix took 6.0 CPU seconds
  // for the old version the run took 132 CPU seconds
  int yy;
  vector<int> const &index=CP->index[CP->seq_x[xx].ID];


  assert(xx>=0 && xx<(int)CP->prevYindex.size());


  //printf("\n\nAsking for xx=%d  pos=%f (prevPos=%f)\n",xx,pos,CP->prevPos[xx]);
  if(CP->prevYindex[xx]<0 || pos < CP->prevPos[xx]) {  // Not going forward. Start from beginning.
    //printf("STARTING FROM BEGINNING\n");
    yy=0;
  } else {  // Going forward. Start where we left off.
    yy=CP->prevYindex[xx];
  }

  assert(yy>=0);

  //printf("Searching.. yy=%d < %d\n",yy,index.size());



  while(yy < (int)index.size() && CP->seq_y[index.at(yy)].epos <= pos) {
    //printf("Skipping in a loop.. yy=%d (epos=%f)\n",yy,CP->seq_y[index.at(yy)].epos);
    yy++;
  }


  yy--; // We went too far. Let's go back.


  
  /*
  printf("xx=%d yy=%d < len(index)=%d\n",xx,yy, index.size());
  printf("epos=%f\n",
	 ((yy>=0 && (yy)<(int)index.size())?CP->seq_y[index.at(yy)].epos:-1.0));
  printf("pos=%f\n",pos);
  printf("next epos=%f\n",((yy+1)<(int)index.size()?CP->seq_y[index.at(yy+1)].epos:-1.0));
  */


  CP->prevYindex[xx]=yy;
  CP->prevPos[xx]=pos;



  assert(yy==-1 || CP->seq_y[index.at(yy)].epos<pos);
  assert(yy==-1 || yy==((int)index.size()-1)  ||  CP->seq_y[index.at(yy+1)].epos > pos  );
#ifndef NDEBUG
  double prev=-1.0;
  double next=1e9;

  if(yy>0) {
    prev=CP->seq_y[index[yy-1]].epos;
  }
  if(yy>=0 && yy<((int)index.size()-1)) {
    next=CP->seq_y[index[yy+1]].epos;
  }
  if(yy>=0) {
    if(!(prev<CP->seq_y[index[yy]].epos && CP->seq_y[index[yy]].epos<pos && CP->seq_y[index[yy]].epos<next))
      printf("prev=%g epos=%g pos=%g next=%g\n",prev,CP->seq_y[index[yy]].epos,pos,next);
    assert(prev<CP->seq_y[index[yy]].epos && CP->seq_y[index[yy]].epos<pos && CP->seq_y[index[yy]].epos<next);
  }
#endif

  return yy;
}



void outputMemory(double bytes)
{
  if(bytes>(1024*1024)) {
    cout<<(bytes/(1024*1024.0))<<" megabytes";
  } else {
    cout<<(bytes/(1024))<<" kilobytes";
  }
  
}

void memReport(align_AlignmentObject *self)
{
  self->fill_factor=(self->item_count*1.0/(self->CP->seq_x.size() *1.0* self->CP->seq_y.size()));


  double normSize=((self->CP->seq_x.size() * 1.0)*self->CP->seq_y.size() * sizeof(matrixentry));
		   
  cout <<"Sequence length: x: "<<self->CP->seq_x.size()<<", y: "<<self->CP->seq_y.size()<<endl
       <<"Expected memory usage ";
  outputMemory(self->expectedMemUsage);
  
  cout<<endl<<"used memory ~ ";
  outputMemory(self->mem_usage);
  cout<<endl<<"normal matrix size would be ~ ";
     //<< (unsigned long)(((self->CP->seq_x.size() * 1.0)*self->CP->seq_y.size() * sizeof(matrixentry))/1024.0) 

  outputMemory(normSize);
  cout<<endl<<"So using only "<<(100.0*(self->mem_usage/normSize))<<" percent of the matrix"<<endl
      <<"Filling only "<<(self->fill_factor*100.0) <<" percent of the cells"<<endl;
  
      
}

int indexYafterOrAtRealY(int sx, int x,struct __CPSTUF *CP,int real_y)
{
  int s=0,e=CP->matrix[x].size()-1,m=(s+e)>>1;
  int real_m,ret;

  if(real_y==0) { //Quick escape
    return 0;
  }
  while((e-s)>1){
    m=(s+e)>>1;
    real_m=CP->index[CP->seq_x[sx+x].ID][m];
    if(real_m<real_y){
      s=m;
    } else if(real_m>real_y) {
      e=m;
    } else {
      break;
    }
  };
  
  ret=CP->matrix[x].size(); 
  if(e>=s && CP->index[CP->seq_x[sx+x].ID][e]>=real_y) {
    ret=e;
  } 
  if(e>=s && CP->index[CP->seq_x[sx+x].ID][m]>=real_y) {
    ret=m;
  } 
  if(e>=s && CP->index[CP->seq_x[sx+x].ID][s]>=real_y) {
    ret=s;
  } 

  return ret;

}


static PyObject *
alignObject(align_AlignmentObject *self,unsigned long int const sx, unsigned long int const sy,unsigned long int const ex, unsigned long int const ey)
{ 
  /* Compute the alignment matrix in rectangle (sx,sy) -> (ex,ey)  ((ey,ex) exclusive) */

  assert(sx<ex);
  assert(ex<=self->CP->seq_x.size());
  assert(sy<ey);
  assert(ey<=self->CP->seq_y.size());

  // allocating the matrix
  for (uint i=sx; i<ex; i++){
    //cout<<"row:"<<i<<" len:"<<self->CP->index[self->CP->seq_x[i].ID].size()<<endl;
    self->CP->matrix.push_back(vector<matrixentry> (self->CP->index[self->CP->seq_x[i].ID].size()));
    self->item_count+=self->CP->index[self->CP->seq_x[i].ID].size();
    self->mem_usage += self->CP->index[self->CP->seq_x[i].ID].size()*sizeof(matrixentry);	
  }

  memReport(self);


//   for(int i=0;i<self->CP->matrix.size();i++) {
//     cout<<"matrix["<<i<<"].size()="<<self->CP->matrix[i].size()<<endl;
//   }
      

  // filling the matrix using dynamic programming
  matrixentry entry;
  int x, y;

  int cells_filled=0;
  double help;

  for (x=0; x<(int)(ex-sx); x++) {

    cells_filled+=1+self->CP->matrix[x].size();
    if(cells_filled>OUTPUTFREQ) {
      cout<< "filling line "<<x+1<<" of "<<(ex-sx) << " of length "<<self->CP->matrix[x].size()
#ifndef NDEBUG
	  <<" ("<<self->CP->ID_to_TF[self->CP->seq_x[sx+x].ID/2].c_str()<<")"
#endif
	  <<endl;
      cells_filled-=OUTPUTFREQ;
    }

    //for(y=0; y<(int)self->CP->matrix[x].size(); y++){
    for(y=indexYafterOrAtRealY(sx,x,self->CP,sy); y<(int)self->CP->matrix[x].size(); y++){
      int const real_y= self->CP->index[self->CP->seq_x[sx+x].ID][y];
//        if((unsigned int)real_y<sy || (unsigned int)real_y>=ey) {
// 	 cout<<"skipping "<<y<<endl;
//  	continue;
//        }
      if((unsigned int)real_y>=ey) {
      	break;
      }
      entry.value= (store)self->lambda*(self->CP->seq_x[sx+x].weight + self->CP->seq_y[real_y].weight);
      if(entry.value<0.0) entry.value=0.0;
      entry.x=-1;
      entry.y=-1;
      //if(x>250) cout<<x<<","<<y<<endl;
      int xx=x-1;
      while(xx>=0 && self->CP->seq_x[sx+x].pos<(self->CP->seq_x[sx+xx].epos+0.5)) {
	xx--;
	
      }

      //printf("x,y=%d,%d (val=%g=%g*(%g+%g))\n",x,real_y,entry.value,self->lambda,self->CP->seq_x[x].weight,self->CP->seq_y[real_y].weight);

//       cout<<"Before:"<<x<<"x"<<y<<":"<<entry.value<<endl;

      int real_yy;
      double delta_x,delta_y;
      for(; 
	  xx>=0 && 
	    (delta_x=(self->CP->seq_x[sx+x].pos-self->CP->seq_x[sx+xx].epos-1.0))<MAX_BP_DIST; 
	  xx--){
#ifdef HEAVY_DEBUG
	assert(indexBeforeOrAtBp(self->CP->index[self->CP->seq_x[sx+xx].ID],self->CP->seq_y,self->CP->seq_y[real_y].pos-0.5) ==
	       indexYBeforeOrAtBp(sx+xx,self->CP,self->CP->seq_y[real_y].pos-0.5));
#endif
	//for (int yy=indexBeforeOrAtBp(self->CP->index[self->CP->seq_x[xx].ID],self->CP->seq_y,self->CP->seq_y[real_y].pos-0.5);
	for (int yy=indexYBeforeOrAtBp(sx+xx,self->CP,self->CP->seq_y[real_y].pos-0.5);
	     yy>=0 && 
	       (delta_y=(self->CP->seq_y[real_y].pos - self->CP->seq_y[(real_yy=self->CP->index[self->CP->seq_x[sx+xx].ID][yy])].epos-1.0))<MAX_BP_DIST 
	       && real_yy>=(int)sy; 
	     yy--){

	  assert(delta_x>-0.5);
	  assert(delta_y>-0.5);
	  assert(xx<x);
	  //if(delta_y<0.5)
	  //  printf("delta_x %g, delta_y %g\n",delta_x,delta_y);
	  assert((int)real_yy<(int)real_y);

	  // The actual recursion formula
	  /* TEX SYNTAX */
	  /*
	    \STATE $ D_{i,j}=\max_{0\le k <i, 0\le l<j} \left\{0,%
	    D_{k,l}+\lambda(w_i+v_j)-\mu \frac{\delta p+\delta q}{2}-
	    \nu \frac{(\delta p-\delta q)^2}{\delta p + \delta q}-
	    \xi \phi^2/(\delta p+\delta q)
	    \right\}$
	    
	  */


	  /* 
	  // This formula gives non finite results due to division by zero.
	  help= self->CP->matrix[xx][yy].value +
	    self->lambda *(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight) -
	    self->mu *(delta_x + delta_y)/2.0-
	    self->nu * (delta_x-delta_y)*(delta_x-delta_y)/(delta_x+delta_y)-
	    self->xi *anglepenalty(delta_x,delta_y,self->nuc_per_rotation);
	  */

	  // Be carefull for correctness and speed!!!
	  //if(delta_x==0.0 && delta_y==0.0) {  // Security check for division by zero
// 	  if(real_yy<sy || real_yy>=ey) 
// 	    cout<<xx<<"=xx:yy="<<yy<<":[xx].size()="<<self->CP->matrix[xx].size()<<endl;
	  if((delta_x+delta_y)==0.0) {
	    help=(double)self->CP->matrix[xx][yy].value +
	      self->lambda *(self->CP->seq_x[sx+x].weight + self->CP->seq_y[real_y].weight) -
	      self->mu *(delta_x + delta_y)/2.0;
	  } else {
	    help=(double)self->CP->matrix[xx][yy].value +
	      self->lambda *(self->CP->seq_x[sx+x].weight + self->CP->seq_y[real_y].weight) -
	      self->mu *(delta_x + delta_y)/2.0 -
	      self->nu * (delta_x-delta_y)*(delta_x-delta_y)/(delta_x+delta_y)-
	      self->xi *anglepenalty(delta_x,delta_y,self->nuc_per_rotation);
	  }
	  
	  assert(finite(help));

	    //xi *squaremod((delta_x - delta_y) * 2.0 * PI/nuc_per_rotation) / 
	    //(delta_x + delta_y);
	  

	  //if(help>entry.value && xx!=x && yy!=y){
	  //if(help>entry.value && xx<x && real_yy<real_y){
	  if(help>entry.value){
	    entry.value= (store)help;
	    //if((sx+xx)>0 && yy>0) {  // Make sure that the local alignment ending at (0,0) is properly terminated
	      entry.x=sx+xx;
	      entry.y=yy;
	      //}
	  }
	}
      }
      
//       cout<<"After:"<<x<<"x"<<y<<":"<<entry.value<<endl;
      self->CP->matrix[x][y]= entry;
    }
  }
  




  


  CHECKING_DECREF(self->output);
  self->output=(PyStringObject*)Py_BuildValue("s","No output. Use nextBest()");



  return (PyObject*)self;
}


////////////////////////////////////////////////////////////

unsigned long maxSitesInLimitedDNA(vector<id_triple> &seq,double const limit)
{
  // Returns the maximum number of binding sites on seq
  // that fall within limit of each other.

  unsigned int s=0,e=0;
  unsigned long maxSites=0,curSites=0;

  while(e<seq.size()) {
    if( (seq[e].pos-seq[s].epos)<=limit) {
      e++;
      curSites++;
    } else {
      maxSites=(maxSites>curSites?maxSites:curSites);
      s++;
      curSites--;
    } 
  }

  return maxSites;
}


void clearMatrix(align_AlignmentObject *self)
{
  // Need to clear matrix. It's temporary if we use memory saving.
  for(vector<vector<matrixentry> >::iterator iter=self->CP->matrix.begin();
      iter!=self->CP->matrix.end();iter++) {
    self->item_count-=(*iter).size();
    self->mem_usage-=(*iter).size()*sizeof(matrixentry);
    (*iter).clear();
      
  }
  self->CP->matrix.clear();
}

static PyObject *
alignMemorySaveObject(align_AlignmentObject *self)
{ 
 
  unsigned long maxRowLen=0;
  double storingLimit=0.0;
  map<matCoord,pair<matCoord,store> > maxses;
  priority_queue<MS_res,vector<MS_res>,greater<MS_res> > localBestAligns;



  self->memSaveUsed=1;


  for(vector<vector<int> >::iterator iter=self->CP->index.begin();
      iter!=self->CP->index.end();iter++) {
    maxRowLen=(maxRowLen>(*iter).size()?maxRowLen:(*iter).size());
  }

  // allocating the matrix (only constant number of rows maximum)
  // constant = all sites on all posible bases.

  // unsigned long const slizeWidth=MAX_BP_DIST*self->CP->index.size()/2;

  // Better have it more strict:
  unsigned long const slizeWidth=maxSitesInLimitedDNA(self->CP->seq_x,MAX_BP_DIST+0.5)+1;
  
  for (uint i=0; i<slizeWidth; i++){
    self->CP->matrix.push_back(vector<matrixentry> (maxRowLen));
    self->item_count+=maxRowLen;
    self->mem_usage += maxRowLen*sizeof(matrixentry);	
  }
  
#ifndef NDEBUG
  cout<<"Matrix size: "<<slizeWidth<<"x"<<maxRowLen<<endl;
#endif

  memReport(self);

  // filling the matrix using dynamic programming
  matrixentry entry;
  int x, y;

  int cells_filled=0;
  double help;

  for (x=0; x<(int)self->CP->seq_x.size(); x++) {

    cells_filled+=1+self->CP->matrix[x%slizeWidth].size();
    if(cells_filled>OUTPUTFREQ) {
      cout<< "filling line "<<x+1<<" of "<<self->CP->seq_x.size() << " of length "<<self->CP->index[self->CP->seq_x[x].ID].size()<<endl;
      cells_filled-=OUTPUTFREQ;
    }


    for(y=0; y<(int)self->CP->index[self->CP->seq_x[x].ID].size(); y++){
      int const real_y= self->CP->index[self->CP->seq_x[x].ID][y];

      assert((unsigned int)real_y<self->CP->seq_y.size());

      entry.value= (store)self->lambda*(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight);
      if(entry.value<0.0) entry.value=0.0;
      entry.x=x;
      entry.y=y;
      //if(x>250) cout<<x<<","<<y<<endl;
      int xx=x-1;
      while(xx>=0 && self->CP->seq_x[x].pos<(self->CP->seq_x[xx].epos+0.5)) {
	xx--;
	
      }

      //printf("x,y=%d,%d (val=%g=%g*(%g+%g))\n",x,real_y,entry.value,self->lambda,self->CP->seq_x[x].weight,self->CP->seq_y[real_y].weight);

      int real_yy;
      double delta_x,delta_y;
      for(; 
	  xx>=0 && 
	    (delta_x=(self->CP->seq_x[x].pos-self->CP->seq_x[xx].epos-1.0))<MAX_BP_DIST; 
	  xx--){
#ifdef HEAVY_DEBUG
	assert(indexBeforeOrAtBp(self->CP->index[self->CP->seq_x[xx].ID],self->CP->seq_y,self->CP->seq_y[real_y].pos-0.5) ==
	       indexYBeforeOrAtBp(xx,self->CP,self->CP->seq_y[real_y].pos-0.5));
#endif

	for (int yy=indexYBeforeOrAtBp(xx,self->CP,self->CP->seq_y[real_y].pos-0.5);
	     yy>=0 && 
	       (delta_y=(self->CP->seq_y[real_y].pos - self->CP->seq_y[(real_yy=self->CP->index[self->CP->seq_x[xx].ID][yy])].epos-1.0))<MAX_BP_DIST; 
	     yy--){

	  assert(delta_x>-0.5);
	  assert(delta_y>-0.5);
	  assert(xx<x);
	  //if(delta_y<0.5)
	  //  printf("delta_x %g, delta_y %g\n",delta_x,delta_y);
	  assert((int)real_yy<(int)real_y);

	  // The actual recursion formula
	  /* TEX SYNTAX */
	  /*
	    \STATE $ D_{i,j}=\max_{0\le k <i, 0\le l<j} \left\{0,%
	    D_{k,l}+\lambda(w_i+v_j)-\mu \frac{\delta p+\delta q}{2}-
	    \nu \frac{(\delta p-\delta q)^2}{\delta p + \delta q}-
	    \xi \phi^2/(\delta p+\delta q)
	    \right\}$
	    
	  */


	  /* 
	  // This formula gives non finite results due to division by zero.
	  help= self->CP->matrix[xx][yy].value +
	    self->lambda *(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight) -
	    self->mu *(delta_x + delta_y)/2.0-
	    self->nu * (delta_x-delta_y)*(delta_x-delta_y)/(delta_x+delta_y)-
	    self->xi *anglepenalty(delta_x,delta_y,self->nuc_per_rotation);
	  */

	  // Be carefull for correctness and speed!!!
	  //if(delta_x==0.0 && delta_y==0.0) {  // Security check for division by zero
	  if((delta_x+delta_y)==0.0) {
	    help=(double)self->CP->matrix[xx%slizeWidth][yy].value +
	      self->lambda *(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight) -
	      self->mu *(delta_x + delta_y)/2.0;
	  } else {
	    help=(double)self->CP->matrix[xx%slizeWidth][yy].value +
	      self->lambda *(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight) -
	      self->mu *(delta_x + delta_y)/2.0 -
	      self->nu * (delta_x-delta_y)*(delta_x-delta_y)/(delta_x+delta_y)-
	      self->xi *anglepenalty(delta_x,delta_y,self->nuc_per_rotation);
	  }
	  
	  assert(finite(help));

	    //xi *squaremod((delta_x - delta_y) * 2.0 * PI/nuc_per_rotation) / 
	    //(delta_x + delta_y);
	  

	  if(help>entry.value){
	    entry.value= (store)help;
	    entry.x=self->CP->matrix[xx%slizeWidth][yy].x;
	    entry.y=self->CP->matrix[xx%slizeWidth][yy].y;

	    // Now the entry.(x,y) contains the _starting_ coordinate of
	    // the best local alignment ending here.
	  }
	}
      }
      
      self->CP->matrix[x%slizeWidth][y] = entry;
      // The best local alignment ending to (x,y) starts at (entry.x,entry.y).
      // add the current maximum to MAXSES if getting better alignment
      // ending at that position

      if(entry.value>storingLimit) {
	matCoord sCoord={entry.x,entry.y};
	map<matCoord,pair<matCoord,store> >::iterator prevGood=maxses.find(sCoord);
 	if(prevGood==maxses.end() ||
 	   prevGood->second < entry.value ){  // If new or improved
 	  matCoord eCoord={x,y};
	  if(prevGood==maxses.end()) { // If new
	    MS_res balg={entry.value,sCoord,eCoord};
	    localBestAligns.push(balg);
	  }
 	  maxses[sCoord]=pair<matCoord,store> (eCoord,entry.value);
// 	  storingLimit=(storingLimit>entry.value?storingLimit:entry.value);
 	}
	if(localBestAligns.size()>(unsigned int)self->askedresults) {
	  MS_res maybeWorst=localBestAligns.top();
	  map<matCoord,pair<matCoord,store> >::iterator maybeBetter=maxses.find(maybeWorst.bgin);
	  while(fabs(maybeBetter->second.second-maybeWorst.value) > numeric_limits<store>::epsilon()) {
	    // priority queue and maxses do not match.
	    maybeWorst.value=maybeBetter->second.second;
	    maybeWorst.end=maybeBetter->second.first;
	    localBestAligns.pop();
	    localBestAligns.push(maybeWorst);


	    maybeWorst=localBestAligns.top();
	    maybeBetter=maxses.find(maybeWorst.bgin);
	  }
	  // Take off the worst.
	  maxses.erase(localBestAligns.top().bgin);
	  localBestAligns.pop();
	  storingLimit=(double)localBestAligns.top().value;
	}

      }


      // END TO MAXSES
    }
  }
  

  //cout<<"best: "<<storingLimit <<" maxs: "<<maxses.size() <<endl;

  // Flip the maximums in other format over to the self.
//   map<matCoord,pair<matCoord,store> >::iterator iter=maxses.begin();
//   map<matCoord,pair<matCoord,store> >::iterator prev;
//   while(iter!=maxses.end()) {
//     //
//    self->CP->bestAlignsTmp[iter->second.second]=pair<matCoord,matCoord> (iter->first,iter->second.first);
//     prev=iter;
//     iter++;
//     maxses.erase(prev);
//   }


   assert(maxses.size()==localBestAligns.size());
   for(map<matCoord,pair<matCoord,store> >::iterator iter=maxses.begin();
       iter!=maxses.end();iter++) {
     MS_res balg={iter->second.second,iter->first,iter->second.first};
     self->CP->bestAligns.push(balg);
   }

  clearMatrix(self);
//   for(vector<vector<matrixentry> >::iterator iter=self->CP->matrix.begin();
//       iter!=self->CP->matrix.end();iter++) {
//     (*iter).clear();
//   }
//   self->CP->matrix.clear();

  CHECKING_DECREF(self->output);
  self->output=(PyStringObject*)Py_BuildValue("s","No output. Use nextBest()");



  return (PyObject*)self;
}



void estimateMemoryConsumption(align_AlignmentObject *self)
{
  // Compute expected memory usage.

  //self->expectedMemUsage=self->CP->seq_x.size()*self->CP->seq_y.size()*sizeof(matrixentry)/(self->CP->ID_to_TF.size());

  self->expectedMemUsage=0;
  for(uint i=0;i<self->CP->seq_x.size();i++) {
    self->expectedMemUsage+=self->CP->index[self->CP->seq_x[i].ID].size()*sizeof(matrixentry);
  }
}

////////////////////////////////////////////////////////////

void align_WrongSeqErr(matchlisttype *matchlist,char *message)
{
  
  PyObject *value;
  PyObject *seqList=PyList_New(matchlist->size());
  PyObject *messageStr;
  matchlisttype::iterator iter;
  int i;


  messageStr=PyString_FromString(message);
  
  for(iter=matchlist->begin(),i=0; iter!=matchlist->end(); iter++,i++) {
    PyList_SetItem(seqList,i,PyString_FromString(iter->first.c_str()));
  }

  value=Py_BuildValue("(OO)",messageStr,seqList);
  PyErr_SetObject(PyExc_AttributeError,value);
}


int getSequences(matchlisttype *matchlist, string &firstSeqName,string &secondSeqName,char *firstSeq,char *secondSeq)
{
  matchlisttype::iterator iter;


  if(!secondSeq && matchlist->size()!=2){
    align_WrongSeqErr(matchlist,"Wrong number of sequences in the data!");
    return 0;
  }else if(firstSeq) {
    // Get first user given sequence name
    firstSeqName=string(firstSeq);
    if(matchlist->find(firstSeqName)==matchlist->end()) {
      PyErr_SetString(PyExc_AttributeError,firstSeq);
      return 0;
    }
    // Get second user given sequence name
    if(secondSeq) {
      secondSeqName=string(secondSeq);
      if(matchlist->find(secondSeqName)==matchlist->end()) {
	PyErr_SetString(PyExc_AttributeError,secondSeq);
	return 0;
      }
    } else {  //  if(!secondSeq)

      // Set the second sequence to be something else than the first sequence
      for(iter=matchlist->begin();
	  iter!=matchlist->end() &&
	    iter->first.compare(firstSeqName)==0;iter++) {
	// pass
      }
      if(iter==matchlist->end()) {
	PyErr_SetString(PyExc_AssertionError,"Shoud not be able to get here!");
	return 0;
      }
      secondSeqName=string(iter->first);
    }
  } else {
    iter=matchlist->begin();
    firstSeqName=string(iter->first);
    iter++;
    secondSeqName=string(iter->first);
  }


  return 1;
}


static PyObject *
align_alignCommon(PyObject *self, PyObject *args,istream *data)
{
  PyObject *ret_obj;
  char* stub;
  double lambda, xi, mu, nuc_per_rotation,nu;
  int result_ask;
  char *firstSeq=NULL,*secondSeq=NULL;
  string firstSeqName,secondSeqName;

  if (!PyArg_ParseTuple(args, "siddddd|zz", &stub,&result_ask,
			&lambda, &xi, &mu, &nu,&nuc_per_rotation,&firstSeq,&secondSeq)){
    PyErr_SetString(PyExc_AttributeError, "Couldn't parse arguments!");
    return NULL;
  }


  matchlisttype* matchlist=parseStream(data);

  if(!getSequences(matchlist,firstSeqName,secondSeqName,firstSeq,secondSeq)) {
    delete matchlist;
    return (PyObject*)NULL;
  }
  
  vector<triple> seq_x,seq_y,tmp_seq;
  seq_x=(*matchlist)[firstSeqName];
  seq_y=(*matchlist)[secondSeqName];
  

  // Check that we get short rows.
  string tmp_str;
  if(seq_x.size()<seq_y.size()) {
#ifndef SMALLX
    tmp_str=firstSeqName;
    firstSeqName=secondSeqName;
    secondSeqName=tmp_str;

    tmp_seq=seq_y;
    seq_y=seq_x;
    seq_x=tmp_seq;
#endif
  } else {
#ifdef SMALLX
    tmp_str=firstSeqName;
    firstSeqName=secondSeqName;
    secondSeqName=tmp_str;

    tmp_seq=seq_y;
    seq_y=seq_x;
    seq_x=tmp_seq;
#endif
  }      

  //vector<id_triple> id_seq_x, id_seq_y;
  map<string, uint> TF_to_ID;
  uint id_counter=0;
  for(uint i=0; i< seq_x.size(); i++){
    if (TF_to_ID.find(seq_x[i].TF) == TF_to_ID.end()){
      TF_to_ID[seq_x[i].TF]=id_counter;
      id_counter++;
    }
  }
  for(uint i=0; i< seq_y.size(); i++){
    if (TF_to_ID.find(seq_y[i].TF) == TF_to_ID.end()){
      TF_to_ID[seq_y[i].TF]=id_counter;
      id_counter++;
    }
  }





  PyObject *retargs=Py_BuildValue("(ss)",firstSeqName.c_str(),secondSeqName.c_str());


  align_AlignmentObject *ret_self=(align_AlignmentObject *)alignment_new(&align_AlignmentType,retargs,NULL);
  CHECKING_DECREF(retargs);


  ret_self->CP->ID_to_TF.resize(TF_to_ID.size());
  //  vector<string> ID_to_TF(TF_to_ID.size());
  for(map<string, uint>::iterator i=TF_to_ID.begin();i!= TF_to_ID.end(); i++){
    ret_self->CP->ID_to_TF[i->second]=i->first;
  }

  id_triple idtr;
  for(uint i=0; i< seq_x.size(); i++){
    idtr.ID=TF_to_ID[seq_x[i].TF];
    if(seq_x[i].strand=='+') {
      idtr.ID=idtr.ID*2;
    } else {
      idtr.ID=idtr.ID*2+1;
    }


    idtr.pos=(double)seq_x[i].pos;
    idtr.epos=(double)seq_x[i].epos;
    idtr.weight=seq_x[i].weight;
    idtr.strand=seq_x[i].strand;
    idtr.annot=seq_x[i].annot;

    //id_seq_x.push_back(idtr);
    ret_self->CP->seq_x.push_back(idtr);
    ret_self->CP->prevYindex.push_back(0);
    ret_self->CP->prevPos.push_back(-1.0);

  }

  
  for(uint i=0; i< seq_y.size(); i++){
    idtr.ID=TF_to_ID[seq_y[i].TF];
    if(seq_y[i].strand=='+') {
      idtr.ID=idtr.ID*2;
    } else {
      idtr.ID=idtr.ID*2+1;
    }

    idtr.pos=(double)seq_y[i].pos;
    idtr.epos=(double)seq_y[i].epos;
    idtr.weight=seq_y[i].weight;
    idtr.strand=seq_y[i].strand;
    idtr.annot=seq_y[i].annot;

    //id_seq_y.push_back(idtr);
    ret_self->CP->seq_y.push_back(idtr);
  }


  delete matchlist;



  ret_self->lambda=lambda;
  ret_self->xi=xi;
  ret_self->mu=mu;
  ret_self->nu=nu;
  ret_self->nuc_per_rotation=nuc_per_rotation;
  ret_self->askedresults=result_ask;
  ret_self->secs_to_align=0;

  clock_t before,after;
  clock_t clocks_to_align;
  // Start timing
  before=clock();

  ret_self->mem_usage=0;

  // building the index map
  ret_self->CP->index.resize(ret_self->CP->ID_to_TF.size()*2);

  for(uint id=0; id<ret_self->CP->index.size(); id++){
    for(uint y=0; y<ret_self->CP->seq_y.size(); y++) {
      if (id == ret_self->CP->seq_y[y].ID){
	ret_self->CP->index[id].push_back(y);
	ret_self->mem_usage += sizeof(int);
      }
    }
  }

  estimateMemoryConsumption(ret_self);

#ifdef SAVE_MEM
  if(ret_self->expectedMemUsage>SAVE_MEM_LIMIT) {
    cout<<"Saving memory!"<<endl;
    ret_obj=(PyObject*)alignMemorySaveObject(ret_self);
  } else {
    ret_obj=(PyObject*)alignObject(ret_self,0,0,ret_self->CP->seq_x.size(),ret_self->CP->seq_y.size());
  }
#else
  ret_obj=(PyObject*)alignObject(ret_self,0,0,ret_self->CP->seq_x.size(),ret_self->CP->seq_y.size());
#endif

  // End timing
  after=clock();

  clocks_to_align=after-before;
  ret_self->secs_to_align+=((double)clocks_to_align)/CLOCKS_PER_SEC;
  


  return ret_obj;

}







static PyObject *
align_alignfile(PyObject *self, PyObject *args)
{
  char* file;
  PyObject *ret;

  ret=PySequence_GetItem(args,0);
  if (!ret) {
    PyErr_SetString(PyExc_AttributeError,"Invalid file argument");
    return  (PyObject*)NULL;
  } else {
    file=PyString_AsString(ret);
  }

  istream *inData;
  ifstream clearData;

#ifdef HAVE_GZSTREAM   // If we have zlib
  
  igzstream gzData(file);
  inData=&gzData;

  if(!gzData.good()) {
    clearData.open(file);
    inData=&clearData;
  }
#else // If we do not have zlib
  clearData.open(file);
  inData=&clearData;
#endif
  ret=align_alignCommon(self,args,inData);


#ifdef HAVE_GZSTREAM
  gzData.close();
#endif


  clearData.close();

  return ret;
}


static PyObject *
align_aligndata(PyObject *self, PyObject *args)
{
  char* data;

  if (!(data=PyString_AsString(PySequence_GetItem(args, 0)))){
    PyErr_SetString(PyExc_AttributeError,"Invalid data argument!");
    return 0;
  }
  
  istringstream inData(data);

  return align_alignCommon(self,args,&inData);
}




//////////////////////////////////////////////////////////////////////





static PyMethodDef alignMethods[] = {
  #ifdef ALIGN_OUTPUT
  {"draw",  align_draw, METH_VARARGS,
   "Draws the file"},
  #endif
  {"aligndata",  align_aligndata, METH_VARARGS,
   "aligns computed sequences. Input: data,#results,lambda,xi,mu,nu,nuc_per_rot"},
  {"alignfile",  align_alignfile, METH_VARARGS,
   "aligns sequences from a gff-file. Input: filename,,#results,lambda,xi,mu,nu,nuc_per_rot"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


//PyMODINIT_FUNC
extern "C"
void
initalign(void)
{
    PyObject* m=NULL;

    align_AlignmentType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&align_AlignmentType) < 0)
      return;

    m=Py_InitModule("eellib.align", alignMethods);

    if(m==NULL)
      return;


    if(import_alnCols()<0)
      return;

    Py_INCREF(&align_AlignmentType);
    PyModule_AddObject(m, "Alignment", (PyObject *)&align_AlignmentType);

}


static PyObject*
alignment_nextBest(align_AlignmentObject *self)
{

  // backtracing
  int pos_x=0, pos_y=0;
  int sx=0,sy=0,real_sy=0,real_ey=0;

  PyObject *ret=PyList_New(0);


  if(self->memSaveUsed) {
    
    MS_res bounds=self->CP->bestAligns.top();
    self->CP->bestAligns.pop();
    
    sx=bounds.bgin.x;
    sy=bounds.bgin.y;
    uint ex=bounds.end.x;
    uint ey=bounds.end.y;
    real_sy= self->CP->index[self->CP->seq_x[sx].ID][sy];
    real_ey= self->CP->index[self->CP->seq_x[ex].ID][ey];

    //cout<<"sx"<<sx<<" sy"<<sy<<" ex"<<ex<<" ey"<<ey<<":len="<<self->CP->index[self->CP->seq_x[sx].ID].size()<<endl;

    assert(self->CP->seq_x[sx].ID==self->CP->seq_y[real_sy].ID);
    assert(self->CP->seq_x[ex].ID==self->CP->seq_y[real_ey].ID);


    clearMatrix(self);

    clock_t before,after;
    clock_t clocks_to_align;
    // Start timing
    before=clock();

    alignObject(self,sx,real_sy,ex+1,real_ey+1);

    // End timing
    after=clock();
    
    clocks_to_align=after-before;
    self->secs_to_align+=((double)clocks_to_align)/CLOCKS_PER_SEC;


    pos_x=ex;
    pos_y=ey;

  } else {
    if(findMax(pos_x,pos_y,self)<0.0) {  // No more alignments
      CHECKING_DECREF(ret);
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

    

  
  
  // Format the output.
  //string ret;
  assert(pos_x>=0);
  assert(pos_y>=0);
  assert((pos_x-sx)<(int)self->CP->matrix.size() );
  assert( pos_y<(int)self->CP->matrix[pos_x-sx].size());
  if (pos_x>=0 && pos_y>=0 && (pos_x-sx)<(int)self->CP->matrix.size() && pos_y<(int)self->CP->matrix[pos_x-sx].size()){
    do {
      uint const real_y= self->CP->index[self->CP->seq_x[pos_x].ID][pos_y]; 
     
      // don't look at this position again
      if (self->CP->matrix[pos_x-sx][pos_y].value >0.0) {
	self->CP->matrix[pos_x-sx][pos_y].value *= -1.0;
      }
      

#ifndef NDEBUG
      if(self->CP->seq_x[pos_x].strand!=self->CP->seq_y[real_y].strand)
	printf("strand conflict %d%c %d%c\n",self->CP->seq_x[pos_x].ID,self->CP->seq_x[pos_x].strand,self->CP->seq_y[real_y].ID,self->CP->seq_y[real_y].strand);
#endif


      PyObject *ret_item=PyAln_New(self->CP->ID_to_TF[self->CP->seq_x[pos_x].ID/2].c_str(), // motif name
				  pos_x, // on site sequence X
				  real_y, // on site sequence Y
				  (int)self->CP->seq_x[pos_x].pos, 
				  (int)self->CP->seq_x[pos_x].epos,
				  (int)self->CP->seq_y[real_y].pos, // start DNA pos on Y
				  (int)self->CP->seq_y[real_y].epos, //end DNA pos on Y
				  (char)self->CP->seq_y[real_y].strand, //strand
				  (double)fabs((double)self->CP->matrix[pos_x-sx][pos_y].value), // Total align score this far
				  (double)self->CP->seq_x[pos_x].weight,
				   (double)self->CP->seq_y[real_y].weight,
				   self->CP->seq_x[pos_x].annot.c_str(),
				   self->CP->seq_y[real_y].annot.c_str());

      PyList_Append(ret,ret_item);
      CHECKING_DECREF(ret_item);
		  

      int new_pos_x= self->CP->matrix[pos_x-sx][pos_y].x;
      pos_y= self->CP->matrix[pos_x-sx][pos_y].y;
      pos_x= new_pos_x;
    } while (pos_x>=0 && pos_y>=0); 
  }
    
  return ret;

}
