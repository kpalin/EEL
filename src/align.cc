/*
 * Local variables:
 *  compile-command: "../make -k"
 * End:
 */



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


#include <vector>
#include <map>
#include <algorithm>


/* The binding sites have to be less than MAX_BP_DIST apart */
#ifndef MAX_BP_DIST
#define MAX_BP_DIST 1000
#endif

#ifndef OUTPUTFREQ
#define OUTPUTFREQ 3000
#endif

#undef SMALLX   /* Define SMALLX if you want smaller sequence to be x. 
		 This (as of 22.10.2003) speeds up the mehtod a little*/
#define SMALLX 1

// #define SAVE_MEM 1

#ifdef SAVE_MEM        // Under memory constraints, define SAVE_MEM. Only a tiny speed loss.
typedef float store;
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
};

struct id_triple
{
  string TF;
  uint ID;
  double pos;
  double epos;
  //uint pos;
  double weight;
  char strand;
};

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
int parseLine(string const line,triple &tri,string &sequence)
{
  string TF;
  int pos,epos;
  double weight;
  char strand;
  istringstream help(line);
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
    if(parseLine(line,tri,sequence)) {
      (*matchlist)[sequence].push_back(tri);
    }
  }



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




struct __CPSTUF {
  vector<id_triple> seq_x;
  vector<id_triple> seq_y;
  vector<string> ID_to_TF;
  vector<vector<int> > index;
  vector<vector<matrixentry> > matrix;

  vector<int> prevYindex;
  vector<double> prevPos;
}  ;

typedef struct {
  PyObject_HEAD

    /* Type-specific fields go here. */
  PyStringObject *x_name,*y_name;
  PyStringObject *output;
  PyObject *bestAlignments;

  double secs_to_align;
  double lambda;
  double xi;
  double mu;
  double nu;
  double nuc_per_rotation;
  uint num_of_align;
  int mem_usage;
  int item_count;
  double fill_factor;

  struct __CPSTUF *CP;

} align_AlignmentObject;






static 
void
alignment_dealloc(align_AlignmentObject* self)
{
    Py_XDECREF(self->x_name);
    Py_XDECREF(self->y_name);
    Py_XDECREF(self->output);
    Py_XDECREF(self->bestAlignments);
    
    delete self->CP;
    

    self->ob_type->tp_free((PyObject*)self);
}



extern "C" int
alignment_init(align_AlignmentObject *self, PyObject *args, PyObject *kwds)
{

    double secs;
    PyObject  *x_name=NULL,*y_name=NULL,*output=NULL;

    static char *kwlist[] = {"x_name","y_name","secs_to_align","output", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|SSdS", kwlist,
				      &x_name,&y_name,&secs,&output))
        return -1; 

    self->secs_to_align = secs;

    if (x_name) {
        Py_XDECREF(self->x_name);
        Py_INCREF(x_name);
        self->x_name = (PyStringObject*)x_name;
    }

    if (y_name) {
        Py_XDECREF(self->y_name);
        Py_INCREF(y_name);
        self->y_name = (PyStringObject*)y_name;
    }


    if (output) {
        Py_XDECREF(self->output);
        Py_INCREF(output);
        self->output = (PyStringObject*)output;
    }

    

    self->num_of_align=10;
    self->nu=1.0;
    self->mu=1.0;
    self->lambda=0.5;
    self->xi=1.0;
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
	  Py_DECREF(self);
	  return NULL;
	}
      self->y_name = (PyStringObject*)PyString_FromString("");
      if (self->y_name == NULL)
	{
	  Py_DECREF(self);
	  return NULL;
	}
      self->output = (PyStringObject*)PyString_FromString("");
      if (self->output == NULL)
	{
	  Py_DECREF(self);
	  return NULL;
	}
        
      self->bestAlignments=PyList_New(0);
      if (self->bestAlignments == NULL)
	{
	  Py_DECREF(self);
	  return NULL;
	}

      self->CP=new struct __CPSTUF;
      if (self->CP==NULL) {
	return NULL;
      }

    }

    self->item_count=0;
    alignment_init(self,args,kwds);

    return (PyObject *)self;
}



double findMax(int &pos_x,int &pos_y,align_AlignmentObject *self)
{
  store max=-1.0;
  pos_x=-1;
  pos_y=-1;
 
  for (uint x=0; x< self->CP->matrix.size(); x++) {
    for(uint y=0; y<self->CP->matrix[x].size(); y++) {
      if (self->CP->matrix[x][y].value > max) {
	max= self->CP->matrix[x][y].value;
	pos_x=x;
	pos_y=y;
      }
    }
  }

  return max;
}

static PyObject*
alignment_nextBest(align_AlignmentObject *self)
{

  // backtracing
  int pos_x=0, pos_y=0;
  PyObject *ret=PyList_New(0);

  do {  // Until find distinct new (sub)optimal alignment
    int c=0;

    if(findMax(pos_x,pos_y,self)<0.0) {  // No more alignments
      Py_DECREF(ret);
      Py_INCREF(Py_None);
      return Py_None;
    }
  
    // Do not return matches ending to same position.
    if (pos_x>=0 && pos_y>=0 && pos_x<(int)self->CP->matrix.size() && pos_y<(int)self->CP->matrix[pos_x].size()){
      int tmp_x=pos_x,tmp_y=pos_y,new_tmp_x;
      do {
	new_tmp_x=self->CP->matrix[tmp_x][tmp_y].x;
	tmp_y=self->CP->matrix[tmp_x][tmp_y].y;
	tmp_x=new_tmp_x;
      } while(tmp_x>=0 && tmp_y>=0 && self->CP->matrix[tmp_x][tmp_y].value>0.0);
      
      if(tmp_x>=0 && tmp_y>=0 && self->CP->matrix[tmp_x][tmp_y].value<0.0) {
	tmp_x=pos_x;
	tmp_y=pos_y;
	while(self->CP->matrix[tmp_x][tmp_y].value>0.0) {
	  c++;
	  self->CP->matrix[tmp_x][tmp_y].value *= -1.0;
	  new_tmp_x=self->CP->matrix[tmp_x][tmp_y].x;
	  tmp_y=self->CP->matrix[tmp_x][tmp_y].y;
	  tmp_x=new_tmp_x;
	} 

      }
    
    }
#ifdef EXTRADEBUG
    cout<<"Going "<<c<<" steps from ("<<pos_x<<","<<pos_y<<")"<<endl;
#endif

  } while(pos_x>=0 && pos_y>=0 && self->CP->matrix[pos_x][pos_y].value<0.0);

  // End overlap matching check.
    

  
  
  // Format the output.
  //string ret;
  if (pos_x>=0 && pos_y>=0 && pos_x<(int)self->CP->matrix.size() && pos_y<(int)self->CP->matrix[pos_x].size()){
    //while (matrix[pos_x][pos_y].x || matrix[pos_x][pos_y].y)
    do {
      uint const real_y= self->CP->index[self->CP->seq_x[pos_x].ID][pos_y];
      
      // don't look at this position again
      if (self->CP->matrix[pos_x][pos_y].value >0.0) {
	self->CP->matrix[pos_x][pos_y].value *= -1.0;
      }
      

#ifndef NDEBUG
      if(self->CP->seq_x[pos_x].strand!=self->CP->seq_y[real_y].strand)
	printf("strand conflict %d%c %d%c\n",self->CP->seq_x[pos_x].ID,self->CP->seq_x[pos_x].strand,self->CP->seq_y[real_y].ID,self->CP->seq_y[real_y].strand);
#endif
      PyObject *ret_item=Py_BuildValue("(iids(ii)(ii)c)",
				       pos_x, real_y, (double)abs((double)self->CP->matrix[pos_x][pos_y].value), 
				       self->CP->ID_to_TF[self->CP->seq_x[pos_x].ID/2].c_str(),
				       (int)self->CP->seq_x[pos_x].pos, (int)self->CP->seq_x[pos_x].epos,
				       (int)self->CP->seq_y[real_y].pos, (int)self->CP->seq_y[real_y].epos, (char)self->CP->seq_y[real_y].strand);
      PyList_Append(ret,ret_item);
      Py_XDECREF(ret_item);
		  
// 	  snprintf(result,255,"D[%d][%d]=%.2f %s (%d,%d) <=> (%d,%d)\n",
// 		   pos_x, real_y, abs(self->CP->matrix[pos_x][pos_y].value), 
// 		   self->CP->ID_to_TF[self->CP->seq_x[pos_x].ID].c_str(),
// 		   (uint)self->CP->seq_x[pos_x].pos, (uint)self->CP->seq_x[pos_x].epos,
// 		   (uint)self->CP->seq_y[real_y].pos, (uint)self->CP->seq_y[real_y].epos);
      
      //ret=result + ret;

      int new_pos_x= self->CP->matrix[pos_x][pos_y].x;
      pos_y= self->CP->matrix[pos_x][pos_y].y;
      pos_x= new_pos_x;
    } while (pos_x>=0 && pos_y>=0); 
  }
    
  return ret;

}



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
    {"item_count",T_INT, offsetof(align_AlignmentObject,item_count), 0,
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

static PyObject *
alignObject(align_AlignmentObject *self)
{ 
 
  tms before,after;
  long ticks_per_sec=sysconf(_SC_CLK_TCK);
  long ticks_to_align;
  uint mem_usage=0;
  

  // Start timing
  times(&before);


  // building the index map
  self->CP->index.resize(self->CP->ID_to_TF.size()*2);

  for(uint id=0; id<self->CP->index.size(); id++){
    for(uint y=0; y<self->CP->seq_y.size(); y++) {
      if (id == self->CP->seq_y[y].ID){
	self->CP->index[id].push_back(y);
	mem_usage += sizeof(int);
      }
    }
  }

  // allocating the matrix
  for (uint i=0; i<self->CP->seq_x.size(); i++){
    self->CP->matrix.push_back(vector<matrixentry> (self->CP->index[self->CP->seq_x[i].ID].size()));
    self->item_count+=self->CP->index[self->CP->seq_x[i].ID].size();
    mem_usage += self->CP->index[self->CP->seq_x[i].ID].size()*sizeof(matrixentry);	
  }
  self->mem_usage=mem_usage;
  self->fill_factor=(self->item_count*1.0/(self->CP->seq_x.size() * self->CP->seq_y.size()));


  cout <<"Sequence length: x: "<<self->CP->seq_x.size()<<", y: "<<self->CP->seq_y.size()<<endl
       <<"used memory ~ "<< (mem_usage>>10)  <<" kilobytes"<<endl
       <<"normal matrix size would be ~ "
       << ((self->CP->seq_x.size() * self->CP->seq_y.size() * sizeof(matrixentry))>>10) <<" kilobytes"<<endl
    <<"So using only "<<(100.0*mem_usage)/(self->CP->seq_x.size() * self->CP->seq_y.size() * sizeof(matrixentry))<<" percent of the matrix"<<endl
       <<"Filling only "<<(self->fill_factor*100.0) <<" percent of the cells"<<endl;
  
      

  // filling the matrix using dynamic programming
  matrixentry entry;
  int x, y;

  int cells_filled=0;
  double help;

  for (x=0; x<(int)self->CP->seq_x.size(); x++) {

    cells_filled+=1+self->CP->matrix[x].size();
    if(cells_filled>OUTPUTFREQ) {
      cout<< "filling line "<<x+1<<" of "<<self->CP->seq_x.size() << " of length "<<self->CP->matrix[x].size()<<endl;
      cells_filled-=OUTPUTFREQ;
    }


    for(y=0; y<(int)self->CP->matrix[x].size(); y++){
      int const real_y= self->CP->index[self->CP->seq_x[x].ID][y];

      entry.value= (store)self->lambda*(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight);
      if(entry.value<0.0) entry.value=0.0;
      entry.x=-1;
      entry.y=-1;
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
	//for (int yy=indexBeforeOrAtBp(self->CP->index[self->CP->seq_x[xx].ID],self->CP->seq_y,self->CP->seq_y[real_y].pos-0.5);
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
	    help=(double)self->CP->matrix[xx][yy].value +
	      self->lambda *(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight) -
	      self->mu *(delta_x + delta_y)/2.0;
	  } else {
	    help=(double)self->CP->matrix[xx][yy].value +
	      self->lambda *(self->CP->seq_x[x].weight + self->CP->seq_y[real_y].weight) -
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
	    entry.x=xx;
	    entry.y=yy;
	  }
	}
      }
      
      self->CP->matrix[x][y]= entry;
    }
  }
  



  /*
  // Output preparation
  string alignment;
  char result[255];
  snprintf(result,255,"### lambda=%f mu=%f nu=%f xi=%f Nucleotides per rotation=%f\n",self->lambda,self->mu,self->nu,self->xi,self->nuc_per_rotation);
  alignment.append(string(result));

  snprintf(result,255,"### D[%s][%s]\n",PyString_AsString((PyObject*)self->x_name),
	   PyString_AsString((PyObject*)self->y_name));
  alignment.append(string(result));


  

  // backtracing
  int pos_x=0, pos_y=0;
  for(uint i=1; i<=self->num_of_align; i++){
    double max=0.0;
    for (uint x=0; x< self->CP->matrix.size(); x++) {
      for(uint y=0; y<self->CP->matrix[x].size(); y++) {
	if (self->CP->matrix[x][y].value > max) {
	  max= self->CP->matrix[x][y].value;
	  pos_x=x;
	  pos_y=y;
	}
      }
    }

    string ret;
    if (pos_x>=0 && pos_y>=0 && pos_x<(int)self->CP->matrix.size() && pos_y<(int)self->CP->matrix[pos_x].size()){
      //while (matrix[pos_x][pos_y].x || matrix[pos_x][pos_y].y)
	do {
	  uint const real_y= self->CP->index[self->CP->seq_x[pos_x].ID][pos_y];
	  
	  // don't look at this position again
	  if (self->CP->matrix[pos_x][pos_y].value >0) {
	    self->CP->matrix[pos_x][pos_y].value *= -1;
	  }
	  
	  snprintf(result,255,"D[%d][%d]=%.2f %s (%d,%d) <=> (%d,%d)\n",
		   pos_x, real_y, abs(self->CP->matrix[pos_x][pos_y].value), 
		   self->CP->ID_to_TF[self->CP->seq_x[pos_x].ID/2].c_str(),
		   (uint)self->CP->seq_x[pos_x].pos, (uint)self->CP->seq_x[pos_x].epos,
		   (uint)self->CP->seq_y[real_y].pos, (uint)self->CP->seq_y[real_y].epos);

	  ret=result + ret;
	  int new_pos_x= self->CP->matrix[pos_x][pos_y].x;
	  pos_y= self->CP->matrix[pos_x][pos_y].y;
	  pos_x= new_pos_x;
	} while (pos_x>=0 && pos_y>=0); 
    }
    snprintf(result,255,"###  Alignment No %d  ###\n",i);
    ret= result + ret;
    
    alignment.append(ret);
  }


  */

  // End timing
  times(&after);

  ticks_to_align=((after.tms_utime-before.tms_utime)+
		 (after.tms_stime-before.tms_stime));
  self->secs_to_align=((double)ticks_to_align)/ticks_per_sec;
  
  //snprintf(result,255,"### Alignment took %.1f CPU seconds.\n",self->secs_to_align);
  //alignment.append(string(result));

  


  Py_DECREF(self->output);
  self->output=(PyStringObject*)Py_BuildValue("s","No output. Use nextBest()");



  return (PyObject*)self;
}



static PyObject *
align_alignCommon(PyObject *self, PyObject *args,istream *data)
{
  char* stub;
  int num_of_align;
  double lambda, xi, mu, nuc_per_rotation,nu;
  string firstSeqName,secondSeqName;

  if (!PyArg_ParseTuple(args, "siddddd", &stub, &num_of_align, 
			&lambda, &xi, &mu, &nu,&nuc_per_rotation)){
    return  Py_BuildValue("s", "");
  }


  matchlisttype* matchlist=parseStream(data);
  matchlisttype::iterator iter;

  if(matchlist->size()<2)
    {
      #ifdef ALIGN_OUTPUT
      cerr<<"Error: too few sequences in data"<<endl;
      #endif
      return Py_BuildValue("s", "");
    }

 if(matchlist->size()==2){
    iter=matchlist->begin();
    firstSeqName=string(iter->first);
    iter++;
    secondSeqName=string(iter->first);
  }

  if(matchlist->size()>2)
    {
      cout<<"Please select two sequences  you want to be aligned."<<endl;
      vector<string> seqNames;
      for(iter=matchlist->begin(); iter!= matchlist->end(); iter++){
        seqNames.push_back(iter->first);
      }
      for(uint i=0; i<seqNames.size(); i++){
        cout<<"("<<i<<") "<<seqNames[i]<<endl;
      }
      int name1, name2;
      cin>>name1;
      if(!cin.good()) {
        cin.clear();
        cin.ignore(INT_MAX,'\n');
        name1=-1;
      }

      while(name1<0 || name1 >= (int)seqNames.size()){
        cout<<"wrong number!"<<endl;
        cin>>name1;
        if(!cin.good()) {
          cin.clear();
          cin.ignore(INT_MAX,'\n');
          name1=-1;
        }
      }
      firstSeqName=seqNames[name1];
      cout<<"Using sequence "<<firstSeqName<<" as first sequence"<<endl
          <<"select second sequence"<<endl;
      cin>>name2;
      if(!cin.good()) {
        cin.clear();
        cin.ignore(INT_MAX,'\n');
        name2=-1;
      }
      while(name2<0 || name2>=(int)seqNames.size() || name1==name2){
        cout<<"wrong number or number equals first one!"<<endl;
        cin>>name2;
        if(!cin.good()) {
          cin.clear();
          cin.ignore(INT_MAX,'\n');
          name2=-1;
        }
      }
      secondSeqName=seqNames[name2];
      cout<<"Using sequence "<<secondSeqName<<" as second sequence"<<endl;
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
  Py_DECREF(retargs);


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
    //id_seq_y.push_back(idtr);
    ret_self->CP->seq_y.push_back(idtr);
  }






  ret_self->lambda=lambda;
  ret_self->xi=xi;
  ret_self->mu=mu;
  ret_self->nu=nu;
  ret_self->nuc_per_rotation=nuc_per_rotation;
  ret_self->num_of_align=num_of_align;


  return (PyObject*)alignObject(ret_self);
  //id_seq_x, id_seq_y, ID_to_TF, 
  //	       num_of_align, lambda, xi, mu, nuc_per_rotation,firstSeqName,secondSeqName);
}




static PyObject *
align_alignfile(PyObject *self, PyObject *args)
{
  char* file;
  int num_of_align;
  double lambda, xi, mu, nuc_per_rotation,nu;
  PyObject *ret;
  string firstSeqName,secondSeqName;

  if (!PyArg_ParseTuple(args, "siddddd", &file, &num_of_align, 
			&lambda, &xi, &mu, &nu,&nuc_per_rotation)){
    return  Py_BuildValue("s", "");
  }

  ifstream inData(file);

  ret=align_alignCommon(self,args,&inData);

  inData.close();

  return ret;
}


static PyObject *
align_aligndata(PyObject *self, PyObject *args)
{
  char* data;
  int num_of_align;
  double lambda, xi, mu, nuc_per_rotation,nu;
  string firstSeqName,secondSeqName,sequence;


  if (!PyArg_ParseTuple(args, "siddddd", &data, &num_of_align, 
			&lambda, &xi, &mu, &nu,&nuc_per_rotation)){
    return Py_BuildValue("s", "");
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
   "aligns computed sequences"},
  {"alignfile",  align_alignfile, METH_VARARGS,
   "aligns sequences from a gff-file"},
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
    
    m=Py_InitModule("align", alignMethods);

    
    if(m==NULL)
      return;

    Py_INCREF(&align_AlignmentType);
    PyModule_AddObject(m, "Alignment", (PyObject *)&align_AlignmentType);
}

