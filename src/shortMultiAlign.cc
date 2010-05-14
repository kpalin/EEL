#include <Python.h>
#include <stdlib.h>
#include "structmember.h"

#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <iostream>

#include <sys/times.h>
#include <unistd.h>

#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <queue>

#include <limits>

#include <assert.h>

#include <functional>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "alignedCols.h"

#include "multiAlign.cc"

#define PI 3.1415926
#ifndef MAX_BP_DIST
#define MAX_BP_DIST 1000
#endif

typedef int seqCode;
typedef int posCode;
typedef int motifCode;
typedef signed short int siteCode;

typedef double posind;
typedef unsigned int uint;
typedef unsigned long int ulint;

using namespace std;



class runningInd {
  int ind;

  public:
  runningInd() {
    ind=-1;
  }

  int giveInd() {
    ind++;
    return ind;
  }
};

//Class copied from Kimmo's code with most of the methods
class m_Inputs {
  vector<vector<id_triple> > seq;
  map<string,motifCode> TF_to_id;
  map<string,seqCode> SEQ_to_id;
  vector<string> seqNames;
  vector<string> factorNames;


public:
  m_Inputs() { return; }
  m_Inputs(PyObject* data);
  //Sanna's constructor, for aligning partial sequences
  m_Inputs(vector<id_triple> first, vector<id_triple> second);

  int addSite(PyObject* item);
  int sequences() const  { return SEQ_to_id.size(); }
  int factors() const { return TF_to_id.size()*2; }
  char const *sequence(seqCode i) const {return seqNames[(int)i].c_str(); }
  char const *factor(motifCode i) const {return factorNames[((int)i)/2].c_str(); }
  vector<int> sequenceLens();
  int sequenceLens(seqCode seqC) const {return seq[seqC].size(); }


  map<string,seqCode>::iterator sequenceNames() { return SEQ_to_id.begin(); }
  map<string,motifCode>::iterator transfacNames() { return TF_to_id.begin(); }

  id_triple const &getSite(posCode const pos,seqCode const i) const {
    vector<id_triple> const &my_seq=seq.at(i);
    return my_seq.at(pos);
  }

  //Sanna's methods
  vector<id_triple> getSites(vector<posCode> places) {
    vector<id_triple> triples;
    for(int i=0;i<seq.size();i++) {
      triples.push_back(seq.at(i).at(places[i]));
    }
    return triples;
  }

  vector<id_triple> getSeq(seqCode seqID) {
    return seq.at(seqID);
  }

  seqCode seqIdFromName(string name) {
    return SEQ_to_id.find(name)->second;
  }

  seqCode motifIdFromName(string name) {
    return (TF_to_id.find(name)->second*2);
  }

  posCode getPosCode(posind pos,seqCode i) {
    posCode j;
    for(j=0;seq.at(i).at(j).pos!=pos;j++);
    return j;
  }

  posCode getPosCode(posind pos,seqCode i,string motif) {
    posCode j;
    for(j=0;seq.at(i).at(j).pos!=pos;j++);
    while(factor(seq.at(i).at(j).ID)!=motif && seq.at(i).at(j).pos==pos) {
      j++;
      for(;seq.at(i).at(j).pos!=pos;j++);
    }
    if(factor(seq.at(i).at(j).ID)==motif && seq.at(i).at(j).pos==pos) {
      return j;
    } else {
      cout<<"Something wrong, searched position not found. pos="<<pos<<", i="<<i<<", factor name "<<motif<<"\n";
      return -1;
    }
  }

  vector<id_triple> const getPrevSites(posCode const pos,seqCode const i) const {
    id_triple curSite=getSite(pos,i);
    vector<id_triple> sites;
    posCode k=pos-1;
    id_triple site;
    bool ok=1;
    while(ok && k>=0) {
      site=getSite(k,i);
      if(site.epos<curSite.pos) {
        if(curSite.pos-site.epos<MAX_BP_DIST) {
          sites.push_back(site);
        } else {
          ok=0;
        }
      }
      k--;
    }
    return sites;
  }

  uint const getIndgap(posCode const pos,seqCode const i) const {
    id_triple curSite=getSite(pos,i);
    posCode k=pos-1;
    id_triple site;
    bool ok=1;
    uint indgap=0;
    while(ok && k>=0) {
      site=getSite(k,i);
      if(site.epos>=curSite.pos) {
        indgap++;
      }
      if(curSite.pos-site.epos>=MAX_BP_DIST) {
        ok=0;
      }
      k--;
    }
    return indgap;
  }

  m_Inputs *getPart(vector<vector<posCode> > inds) {
    m_Inputs *retm_Inputs=new m_Inputs();

    retm_Inputs->TF_to_id=this->TF_to_id;
    retm_Inputs->factorNames=this->factorNames;

    int count=0;
    for(int i=0;i<sequences();i++) {
      retm_Inputs->seq.resize(count+1);
      if(inds.at(i).at(1)>0) {
        retm_Inputs->SEQ_to_id[this->seqNames.at(i)]=count;
        retm_Inputs->seqNames.push_back(this->seqNames.at(i));
        for(posCode j=inds.at(i).at(0);j<=inds.at(i).at(1);j++) {
          retm_Inputs->seq.at(count).push_back(getSite(j,i));
        }
        count++;
      }
    }

    return retm_Inputs;
  }


  m_Inputs *getPart(seqCode seq1, seqCode seq2) {
    m_Inputs *retInputs=new m_Inputs();

    retInputs->TF_to_id=this->TF_to_id;
    retInputs->factorNames=this->factorNames;

    retInputs->seq.resize(2);
    retInputs->SEQ_to_id[this->seqNames.at(seq1)]=seq1;
    retInputs->SEQ_to_id[this->seqNames.at(seq2)]=seq2;
    retInputs->seqNames.push_back(this->seqNames.at(seq1));
    retInputs->seqNames.push_back(this->seqNames.at(seq2));
    for(posCode j=0;j<sequenceLens(seq1);j++) {
      retInputs->seq.at(0).push_back(getSite(j,seq1));
    }
    for(posCode j=0;j<sequenceLens(seq2);j++) {
      retInputs->seq.at(1).push_back(getSite(j,seq2));
    }

    return retInputs;
  }


  int findSeq(string seqName) {
    map<string,seqCode>::iterator iter=SEQ_to_id.find(seqName);
    if(iter!=SEQ_to_id.end()) {
      return iter->second;
    } else {
      return -1;
    }
  }
};

m_Inputs::m_Inputs(vector<id_triple> first, vector<id_triple> second) {
  this->seq.push_back(first);
  this->seq.push_back(second);
  this->SEQ_to_id["0"]=0;
  this->SEQ_to_id["1"]=1;
}

m_Inputs::m_Inputs(PyObject *inpSeq) {
  // inpSeq is pySequence of pySequences as represented in the GFF file

  PyObject *inpIter=PyObject_GetIter(inpSeq);
  PyObject *inSite;

  if(inpIter==NULL || PyErr_Occurred()!=NULL) {
    return;
  }
  assert(inpIter!=NULL);
  int iter=0;

  while((inSite=PyIter_Next(inpIter))!=NULL) {
    iter++;
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


int m_Inputs::addSite(PyObject *site) {
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


vector<int> m_Inputs::sequenceLens() {
  vector<int> dims=vector<int>();
  for(int i=0;i<sequences();i++) {
    dims.push_back(seq[i].size());
  }
  return dims;
}



class m_Matrix {
  int dims;
  vector<int> dimsizes;
  ulint jumps[11];
  ulint size;
  float *scores;
  uint *backtracks;

public:
  m_Matrix() { return; }

  m_Matrix(int dims, vector<int> dimsizes) {
    this->dims=dims;
    this->dimsizes=dimsizes;
    size=1;
    for(int i=0;i<dims;i++) {
      size*=(ulint)dimsizes.at(i);
    }
    ulint jump=1;
    for(int i=0;i<dims;i++) {
      jumps[i]=jump;
      jump*=dimsizes.at(i);
    }
    jumps[dims]=jump;

    scores=new float[size];
    backtracks=new uint[size*dims];
  }

  ~m_Matrix() {
    if(scores) {
      delete [] scores;
      scores=NULL;
    }
    if(backtracks) {
      delete [] backtracks;
      backtracks=NULL;
    }
  }

  void inScore(float score, int *place) {
    ulint ind=0;
    for(int i=0;i<dims;i++) {
      ind+=place[i]*jumps[i];
    }
    if(ind<size) {
      scores[ind]=score;
    } else {
      cout<<"inScore - Out of range\n";
      cout<<"Size "<<size<<"\n";
      cout<<"Index "<<ind<<"\n";
      cout<<"Place";
      for(int i=0;i<dims;i++) {
        cout<<" "<<place[i];
      }
      cout<<"\n";
      cout<<"Jumps";
      for(int i=0;i<=dims;i++) {
        cout<<" "<<jumps[i];
      }
      cout<<"\n";
    }
  }

  float outScore(int *place) {
    ulint ind=0;
    for(int i=0;i<dims;i++) {
      ind+=place[i]*jumps[i];
    }
    if(ind<size) {
      return scores[ind];
    } else {
      cout<<"outScore - Out of range\n";
      cout<<"Size "<<size<<"\n";
      cout<<"Index "<<ind<<"\n";
      cout<<"Place";
      for(int i=0;i<dims;i++) {
        cout<<" "<<place[i];
      }
      cout<<"\n";
      return NULL;
    }
  }

  void inBT(uint *BT, int *place) {
    ulint ind=0;
    for(int i=0;i<dims;i++) {
      ind+=place[i]*jumps[i];
    }
    if(ind<size*dims) {
      for(ulint i=0;i<(ulint)dims;i++) {
        if(ind+i*jumps[dims]>=size*dims) {
          cout<<"inBT - Out of range\n";
          cout<<""<<ind<<"+"<<i<<"*"<<jumps[dims]<<"="<<ind+i*jumps[dims]<<"\n";
        } else {
          backtracks[ind+i*jumps[dims]]=BT[i];
        }
      }
    } else {
      cout<<"inBT - Out of range\n";
      cout<<"Size "<<size<<"\n";
      cout<<"Index "<<ind<<"\n";
      cout<<"Place";
      for(int i=0;i<dims;i++) {
        cout<<" "<<place[i];
      }
      cout<<"\n";
      cout<<"Jumps";
      for(int i=0;i<=dims;i++) {
        cout<<" "<<jumps[i];
      }
      cout<<"\n";
    }
  }

  void outBT(int *place, uint *BT) {
    ulint ind=0;
    for(int i=0;i<dims;i++) {
      ind+=place[i]*jumps[i];
    }
    if(ind<size*dims) {
      for(ulint i=0;i<(ulint)dims;i++) {
        BT[i]=backtracks[ind+i*jumps[dims]];
      }
    } else {
      cout<<"outBT - Out of range\n";
      cout<<"Size "<<size<<"\n";
      cout<<"Index "<<ind<<"\n";
      cout<<"Place";
      for(int i=0;i<dims;i++) {
        cout<<" "<<place[i];
      }
      cout<<"\n";
    }
  }
};

class SimpleMultiAlign {
  double lambda;
  double xi;
  double mu;
  double nu;
  double nuc_per_rot;

  int gaplimit;

  int dims;

  m_Inputs *data;

  int curind[10];

  m_Matrix *dynm_Matrix;

  float bestScore;
  uint bestInd[10];

  int validsquares;

  list<vector<posCode> > fulls;

  vector<vector<posCode> > oldpaths;
  posind limit[10];

public:
  SimpleMultiAlign() { return; }
  SimpleMultiAlign(PyObject* inpSeq, int gaplimit, double lambda, double mu, double nu, double xi, double nuc_per_rot) {

    this->gaplimit=gaplimit;
    this->lambda=lambda;
    this->mu=mu;
    this->nu=nu;
    this->xi=xi;
    this->nuc_per_rot=nuc_per_rot;

    this->validsquares=0;

    this->data=new m_Inputs(inpSeq);

    this->dims=data->sequences();

    //Going through indexes is started from 0
    for(int i=0;i<dims;i++) {
      curind[i]=0;
      bestInd[i]=0;
    }

    //Empty the list of previous cells, just in case
    fulls.clear();

    //Formating the matrix
    this->dynm_Matrix=new m_Matrix(dims,data->sequenceLens());

    bestScore=-DBL_MAX;

    //Multialign
    fillm_Matrix(0);

    vector<int> seqLens=data->sequenceLens();
    ulint size=1;
    for(uint i=0;i<seqLens.size();i++) {
      size*=seqLens.at(i);
    }
    cout<<"Filled "<<this->validsquares<<" squares of "<<size<<"\n";
  }

  SimpleMultiAlign(m_Inputs *curdata, int gaplimit, double lambda, double mu, double nu, double xi, double nuc_per_rot) {

    this->gaplimit=gaplimit;
    this->lambda=lambda;
    this->mu=mu;
    this->nu=nu;
    this->xi=xi;
    this->nuc_per_rot=nuc_per_rot;

    this->validsquares=0;

    this->data=curdata;

    this->dims=data->sequences();

    //Going through indexes is started from 0
    for(int i=0;i<dims;i++) {
      curind[i]=0;
      bestInd[i]=0;
    }

    //Formating the matrix
    this->dynm_Matrix=new m_Matrix(dims,data->sequenceLens());

    bestScore=-DBL_MAX;

    //Multialign
    fillm_Matrix(0);

    vector<int> seqLens=data->sequenceLens();
    ulint size=1;
    for(uint i=0;i<seqLens.size();i++) {
      size*=seqLens.at(i);
    }
    cout<<"Filled "<<this->validsquares<<" squares of "<<size<<"\n";
  }

  ~SimpleMultiAlign() {
    if(data) {
      delete data;
      data=NULL;
    }
    if(dynm_Matrix) {
      delete dynm_Matrix;
      dynm_Matrix=NULL;
    }
  }

  map<string,seqCode>::iterator sequenceNames() {
    return data->sequenceNames();
  }

  int giveDims() {
    return dims;
  }

  void fillScore() {
    //Fetching info of current cell
    vector<id_triple> curSites;
    curSites.reserve(dims);
    map<motifCode,int> numberOfMotifs;
    numberOfMotifs.clear();

    for(int k=0;k<dims;k++) {
      curSites.push_back(data->getSite(curind[k],k));
      if(numberOfMotifs[curSites.at(k).ID]) {
        numberOfMotifs[curSites.at(k).ID]++;
      } else {
        numberOfMotifs[curSites.at(k).ID]=1;
      }
    }
    float score;
    uint BT[10];
    for(int k=0;k<dims;k++) {
      BT[k]=0;
    }

    //Finding the motif that is represented most in the cell
    map<motifCode,int>::iterator nomIter;
    for(nomIter=numberOfMotifs.begin();nomIter!=numberOfMotifs.end();nomIter++) {
      if(nomIter->second>=dims-gaplimit) {
        break;
      }
    }

    //Cheking if the cell is good enough
    if(nomIter==numberOfMotifs.end()) {
      //If the same motif is not represented in enough seqs in this place, the score is -infinity
      score=-DBL_MAX;
    } else {
      //Saving the code of the most represented motif
      motifCode curID=nomIter->first;

      int backind[10];

      //Searching for earlier motifs in the sequences
      vector<vector<id_triple> > prevSites;
      prevSites.clear();
      vector<int> seqsToHandle;
      seqsToHandle.clear();
      int indgap[10];
      for(int k=0;k<dims;k++) {
        //If the sequence has the motif we are looking at in this spot and it's not the first one, the previous motifs are checked
        if(curSites.at(k).ID==curID && curind[k]>0) {
          vector<id_triple> sites=data->getPrevSites(curind[k],k);
          //If earlier motifs are found, they will be saved and the seq marked as to be handeld
          if(sites.size()>0) {
            seqsToHandle.push_back(k);
            prevSites.push_back(sites);
            backind[k]=1;
            indgap[k]=data->getIndgap(curind[k],k);
          } else {
          //If there are no previous motifs, the seq will not be handeld more at this instance
            backind[k]=0;
            indgap[k]=0;
          }
        } else {
          //If the motif in the seq is not the one we are looking at or the current motif is the first of the seq, the seq will not be handeld more at this instance
          backind[k]=0;
          indgap[k]=0;
        }
      }

      float curmaxscore=0;
      int maxbackind[10];
      for(int k=0;k<dims;k++) {
        maxbackind[k]=0;
      }
      //Previous cells will only be looked at if there are enough seqs to be handeld
      //No holes allowed, a list is used to remember the previous cells
      if((int)seqsToHandle.size()==dims) {
        //First, go through the cells in the list
        list<vector<posCode> >::iterator fullit=fulls.begin();
        if(fullit!=fulls.end()) {
          //At the beginning of the list, there might be some cells that are too far away, those will be dropped out
          int toofar=0;
          
          vector<id_triple> cell=data->getSites(*fullit);
          if((curSites.at(dims-1).pos-cell.at(dims-1).epos)>=MAX_BP_DIST) {
            toofar=1;
          } else if(cell.at(dims-1).epos>curSites.at(dims-1).pos) {
            cout<<"ERROR! Prev. site in the list is overlapping/after the current site, shouldn't be in the list.\n";
            cout<<"Start of current site curSites.at(dims-1).pos="<<curSites.at(dims-1).pos<<", end of prev. site in the list cell.at(dims-1).epos="<<cell.at(dims-1).epos<<"\n";
          }
          while(fullit!=fulls.end() && toofar) {
            fulls.pop_front();
            fullit=fulls.begin();
            if(fullit!=fulls.end()) {
              cell=data->getSites(*fullit);
              toofar=0;
              if((curSites.at(dims-1).pos-cell.at(dims-1).epos)>=MAX_BP_DIST) {
                toofar=1;
              } else if(cell.at(dims-1).epos>curSites.at(dims-1).pos) {
                cout<<"ERROR! Prev. site in the list is overlapping/after the current site, shouldn't be in the list.\n";
                cout<<"Start of current site curSites.at(dims-1).pos="<<curSites.at(dims-1).pos<<", end of prev. site in the list cell.at(dims-1).epos="<<cell.at(dims-1).epos<<"\n";
              }
            }
          }
        }

        //The rest will be handeld
        float curscore=0;
        while(fullit!=fulls.end()) {

          vector<id_triple> cell=data->getSites(*fullit);
          for(int i=0; i<dims; i++) {
            if(cell.at(i).epos>curSites.at(i).pos) {
              cout<<"ERROR! Prev. site in the list is overlapping/after the current site, shouldn't be in the list.\n";
              cout<<"Start of current site curSites.at("<<i<<").pos="<<curSites.at(i).pos<<", end of prev. site in the list cell.at("<<i<<").epos="<<cell.at(i).epos<<"\n";
            }
          }

          int place[10];
          for(int x=0;x<dims;x++) {
            place[x]=(*fullit)[x];
          }
          curscore=dynm_Matrix->outScore(place);
          vector<id_triple> curTriples=data->getSites(*fullit);
          float penalty=0;
          float deltaphi;

          //Calculating the part of the score depending on current and previous cell
          for(uint k=0;k<dims;k++) {
            for(uint l=0;l<dims;l++) {
              if(k<l) {
                int x;
                x=(int)(data->getSite(curind[k],k).pos-curTriples.at(k).epos-1);
                int y;
                y=(int)(data->getSite(curind[l],l).pos-curTriples.at(l).epos-1);
                if(x-y>=0) {
                  deltaphi=(x-y)*2*PI/nuc_per_rot;
                  deltaphi+=PI;
                  deltaphi=fmod(deltaphi,2*PI);
                  deltaphi-=PI;
                } else {
                  deltaphi=(x-y)*2*PI/nuc_per_rot;
                  deltaphi-=PI;
                  deltaphi=fmod(deltaphi,2*PI);
                  deltaphi+=PI;
                }
                float curpenalty=-(mu*(x+y)/2+nu*(x-y)*(x-y)/(x+y)+xi*deltaphi*deltaphi/(x+y));
                penalty+=curpenalty;
              }
            }
          }

          if(penalty>0) {
            cout<<"Positive penalty, something went wrong\n";
          }

          curscore+=penalty;

          if(curscore>curmaxscore) {
            curmaxscore=curscore;
            for(int k=0;k<dims;k++) {
              maxbackind[k]=curind[k]-(*fullit)[k];
            }
          }

          fullit++;
        }

        //Finding and handlign possible new cells
        int i=1;
        if(!fulls.empty()) {
          //In the last dimension, the cells have been handeld up to lastCell
          //So lastCell is found in the list of the previous cells of the last dimension
          //The list begins from the closest eligible previous cell
          list<vector<posCode> >::reverse_iterator revit=fulls.rbegin();
          while(i<=(int)prevSites.at(dims-1).size() && (prevSites.at(dims-1).at(i-1).epos!=data->getSites(*revit).at(dims-1).epos || prevSites.at(dims-1).at(i-1).ID!=data->getSites(*revit).at(dims-1).ID)) {
            i++;
          }
          if(i>(int)prevSites.at(dims-1).size()) {
            cout<<"ERROR! i="<<i<<", prevSites.at(dims-1).size()="<<prevSites.at(dims-1).size()<<"\n";
          }
          //One before is the fisrt unhandeld cell
          i--;
        } else {
          //If the list is empty, all the previous cells will be gone through
          i=prevSites.at(dims-1).size();
        }

        while(i>0) {
          backind[dims-1]=i;

          int ok=1;
          while(ok) {
            //The part of the score that depends on the previous cell

            //The scor of the previous cell is fetched
            int place[10];
            for(int k=0;k<dims;k++) {
              place[k]=curind[k]-backind[k]-indgap[k];
            }
            curscore=dynm_Matrix->outScore(place);

            float penalty=0;
            float deltaphi;

            //If the score of the previous cell is larger than 0, the calculation will continue
            if(curscore>0) {
              //Adding the cell to the list
              vector<id_triple> triple;
              vector<posCode> tabplace;
              for(int k=0;k<dims;k++) {
                triple.push_back(prevSites.at(k).at(backind[k]-1));
                tabplace.push_back(place[k]);
              }

              fulls.push_back(tabplace);
              //Calculating the part of the score that depends on the previous cell
              for(uint k=0;k<seqsToHandle.size();k++) {
                for(uint l=0;l<seqsToHandle.size();l++) {
                  if(k<l) {
                    int x;
                    x=(int)(data->getSite(curind[seqsToHandle.at(k)],seqsToHandle.at(k)).pos-prevSites.at(k).at(backind[seqsToHandle.at(k)]-1).epos-1);
                    int y;
                    y=(int)(data->getSite(curind[seqsToHandle.at(l)],seqsToHandle.at(l)).pos-prevSites.at(l).at(backind[seqsToHandle.at(l)]-1).epos-1);
                    if(x-y>=0) {
                      deltaphi=(x-y)*2*PI/nuc_per_rot;
                      deltaphi+=PI;
                      deltaphi=fmod(deltaphi,2*PI);
                      deltaphi-=PI;
                    } else {
                      deltaphi=(x-y)*2*PI/nuc_per_rot;
                      deltaphi-=PI;
                      deltaphi=fmod(deltaphi,2*PI);
                      deltaphi+=PI;
                    }
                    float curpenalty=-(mu*(x+y)/2+nu*(x-y)*(x-y)/(x+y)+xi*deltaphi*deltaphi/(x+y));
                    penalty+=curpenalty;
                  }
                }
              }

              if(penalty>0) {
                cout<<"Positive penalty, something went wrong\n";
              }

              curscore+=penalty;
            } else {
            //Otherwise, the score of the previous cell is -infinity and max of that and zero is zero
              curscore=0;
            }

            //Saving, if the score is best yet
            if(curscore>curmaxscore) {
              curmaxscore=curscore;
              for(int k=0;k<dims;k++) {
                maxbackind[k]=backind[k]+indgap[k];
              }
            }

            //Moving to the next previous cell
            //As this part is used only for gapless case, all seqs are to be handeld
            //Moving in the last dimension is done separately
            //If there is no moving in any of the seqs, all the motif combinations have been looked at and calculation is finished
            ok=0;
            for(uint k=0;k<dims-1;k++) {
              //If the seq has more previous motifs, move is made to the next one and moving is finished
              if(backind[k]<prevSites.at(k).size() && (curind[k]-backind[k]-indgap[k])>0) {
                backind[k]++;
                ok=1;
                break;
              } else {
              //Else a move is made back to the first previous motif
                backind[k]=1;
              }
            }
          }
          i--;
        }
          
      //Gaps allowed, the basic solution is used
      } else if((int)seqsToHandle.size()>=dims-gaplimit) {
        //Going through previous cells and calculating the score
        bool ok=0;
        for(int k=0;k<dims;k++) {
          if(backind[k]!=0) {
            ok=1;
          }
        }

        float curscore=0;
        while(ok) {
          //Fetching the score of the previous cell
          int place[10];
          for(int k=0;k<dims;k++) {
            place[k]=curind[k]-backind[k]-indgap[k];
          }
          curscore=dynm_Matrix->outScore(place);

          float penalty=0;
          float deltaphi;

          //If the score is larger than zero, calculation is continued
          if(curscore>0) {
            for(uint k=0;k<seqsToHandle.size();k++) {
              for(uint l=0;l<seqsToHandle.size();l++) {
                if(k<l) {
                  int x;
                  x=(int)(data->getSite(curind[seqsToHandle.at(k)],seqsToHandle.at(k)).pos-prevSites.at(k).at(backind[seqsToHandle.at(k)]-1).epos-1);
                  int y;
                  y=(int)(data->getSite(curind[seqsToHandle.at(l)],seqsToHandle.at(l)).pos-prevSites.at(l).at(backind[seqsToHandle.at(l)]-1).epos-1);
                  if(x-y>=0) {
                    deltaphi=(x-y)*2*PI/nuc_per_rot;
                    deltaphi+=PI;
                    deltaphi=fmod(deltaphi,2*PI);
                    deltaphi-=PI;
                  } else {
                    deltaphi=(x-y)*2*PI/nuc_per_rot;
                    deltaphi-=PI;
                    deltaphi=fmod(deltaphi,2*PI);
                    deltaphi+=PI;
                  }
                  float curpenalty=-(mu*(x+y)/2+nu*(x-y)*(x-y)/(x+y)+xi*deltaphi*deltaphi/(x+y));
                  penalty+=curpenalty;
                }
              }
            }

            if(penalty>0) {
              cout<<"Positive penalty, something went wrong\n";
            }

            curscore+=penalty;
          } else {
          //Else the score is -infinity and maxim of that and zero is zero
            curscore=0;
          }

          //Saving, if the score is best yet
          if(curscore>curmaxscore) {
            curmaxscore=curscore;
            for(int k=0;k<dims;k++) {
              maxbackind[k]=backind[k]+indgap[k];
            }
          }

          //Moving to the next previous cell
       
          //If there is no moving in any of the seqs to be handeld, all the motif combinations have been looked at and calculation is finished
          ok=0;
          for(uint k=0;k<seqsToHandle.size();k++) {
              //If the seq has more previous motifs, move is made to the next one and moving is finished
            if(backind[seqsToHandle.at(k)]<prevSites.at(k).size() && curind[seqsToHandle.at(k)]-backind[seqsToHandle.at(k)]-indgap[seqsToHandle.at(k)]>0) {
              backind[seqsToHandle.at(k)]++;
              ok=1;
              break;
            } else {
              //Else a move is made back to the first previous motif
              backind[seqsToHandle.at(k)]=1;
            }
          }
        }
      }

      //The part of the score that depends only on current cell
      float base=0;
      if(curmaxscore==0) {
        for(int k=0;k<dims;k++) {
          if(curSites.at(k).ID==curID) {
            base+=curSites.at(k).weight;
          }
        }
        base*=lambda*(nomIter->second-1);
      } else {
        for(int k=0;k<(int)seqsToHandle.size();k++) {
          base+=curSites.at(seqsToHandle.at(k)).weight;
        }
        base*=lambda*(seqsToHandle.size()-1);
      }

      score=curmaxscore+base;
      for(int k=0;k<dims;k++) {
        BT[k]=maxbackind[k];
      }
    }

    dynm_Matrix->inScore(score,curind);
    dynm_Matrix->inBT(BT,curind);

    if(score>0) {
      validsquares++;
    }

    //Saving best result
    if(score>bestScore) {
      bestScore=score;
      for(int k=0;k<dims;k++) {
        bestInd[k]=curind[k];
      }
    }
  }

  //Filling the matrix
  void fillm_Matrix(seqCode i) {
    int motifs=data->sequenceLens(i);

    //Recursive going through of the dimensions
    if(i<(this->dims-1)) {
      for(int j=0;j<motifs;j++) {
        curind[i]=j;
        fillm_Matrix(i+1);
      }
    } else if(i==(this->dims-1)) {
    //Last dimension, computing is done
      for(int j=0;j<motifs;j++) {
        curind[i]=j;
        fillScore();
      }

      //After the last cell the list of previous cells is emptied
      fulls.clear();
    }
  }

  void giveBestAligns(map<float, map<int, map<seqCode, vector<id_triple> > > > *aligns, map<float, map<int, vector<float> > > *scores, int k, vector<seqCode> seqs, runningInd *alignInd) {
    if(bestScore>=0) {
      vector<map<seqCode, id_triple> > revalign;
      vector<float> revScore;
      //Starting from the cell with the best score
      int ind[10];
      for(int i=0;i<dims;i++) {
        ind[i]=bestInd[i];
      }

      int ok=1;
      //Going through cells in the aligment
      while(ok) {
        //Fetching data of the motifs of current cell
        map<seqCode, id_triple> sites;

        for(int k=0;k<dims;k++) {
          id_triple site=data->getSite(ind[k],k);
          sites[seqs.at(k)]=site;
        }
        revalign.push_back(sites);
        revScore.push_back(dynm_Matrix->outScore(ind));

        uint BT[10];
        dynm_Matrix->outBT(ind, BT);

        //If the whole backtrack is 0, this is the end of the alignment
        //Otherwise the backtrack is followed
        ok=0;
        for (int k=0;k<dims;k++) {
          if(BT[k]) {
            ind[k]-=BT[k];
            ok=1;
          }
        }
      }

      map<seqCode, vector<id_triple> > align;
      vector<float> aScore;
      vector<float>::reverse_iterator rits=revScore.rbegin();
      for(int i=(revalign.size()-1);i>=0;i--) {
        map<seqCode, id_triple>::iterator it;
        for(it=revalign.at(i).begin();it!=revalign.at(i).end();it++) {
          align[it->first].push_back(it->second);
        }
        aScore.push_back(*rits);
        rits++;
      }
      int aind=alignInd->giveInd();
      (*aligns)[bestScore][aind]=align;
      (*scores)[bestScore][aind]=aScore;

      //Suboptimals
      for(int i=1;i<k;i++) {
        //The previous alignment is removes and the part it affects is recalculated
        //Fetching the list of the cells in the last alignment
        vector<vector<posCode> > oldPath=givePath();

        for(int i=oldPath.size()-1;i>=0;i--) {
          this->oldpaths.push_back(oldPath.at(i));
          fillm_MatrixPart(0,oldPath.at(i));
        }

        bestScore=-DBL_MAX;
        findBest(0);

        if(bestScore>=0) {
          revalign.clear();
          revScore.clear();
          //Starting from the cell with the best score
          int ind[10];
          for(int k=0;k<dims;k++) {
            ind[k]=bestInd[k];
          }

          int ok=1;
          //Going through the cells in the alignment
          while(ok) {
            //Fetching the data of the motifs of current cell
            map<seqCode, id_triple> sites;

            for(int k=0;k<dims;k++) {
              id_triple site=data->getSite(ind[k],k);
              sites[seqs.at(k)]=site;
            }
            revalign.push_back(sites);
            revScore.push_back(dynm_Matrix->outScore(ind));

            uint BT[10];
            dynm_Matrix->outBT(ind, BT);

            //If the whole backtrack is 0, this is the end of the alignment
            //Otherwise the backtrack is followed
            ok=0;
            for (int k=0;k<dims;k++) {
              if(BT[k]) {
                ind[k]-=BT[k];
                ok=1;
              }
            }
          }

          align.clear();
          aScore.clear();
          vector<float>::reverse_iterator rits=revScore.rbegin();
          for(int i=(revalign.size()-1);i>=0;i--) {
            map<seqCode, id_triple>::iterator it;
            for(it=revalign.at(i).begin();it!=revalign.at(i).end();it++) {
              align[it->first].push_back(it->second);
            }
            aScore.push_back(*rits);
            rits++;
          }
          int aind=alignInd->giveInd();
          (*aligns)[bestScore][aind]=align;
          (*scores)[bestScore][aind]=aScore;
        }
      }
    }
  }

  PyObject *giveBest() {
    //If there are no good cells, NULL is returned
    if(bestScore<0) {
      return NULL;
    } else {
      PyObject *goodAlign=PyList_New(0);
      if (goodAlign == NULL || PyErr_Occurred()) {
        cout<<"Creating goodalign failed\n";
        return NULL;
      }
      PyObject *seqs;
      PyObject *coords;
      PyObject *seqPos,*siteScore,*annotations;

      //Starting from the cell with the best score
      int ind[10];
      for(int k=0;k<dims;k++) {
        ind[k]=bestInd[k];
      }

      int ok=1;
      //Going through the cells in the alignment
      while(ok) {
        //Fetching the data of the motifs of current cell
        vector<id_triple> sites;
        map<motifCode,int> numberOfMotifs;
        numberOfMotifs.clear();
        char strands[10];
        siteCode motifIDs[10];

        for(int k=0;k<dims;k++) {
          id_triple site=data->getSite(ind[k],k);
          strands[k]=site.strand;
          motifIDs[k]=site.ID;
          sites.push_back(site);
          if(numberOfMotifs[sites.at(k).ID]) {
            numberOfMotifs[sites.at(k).ID]++;
          } else {
            numberOfMotifs[sites.at(k).ID]=1;
          }
        }

        //The motif most present in the cell
        map<motifCode,int>::iterator nomIter;
        for(nomIter=numberOfMotifs.begin();nomIter!=numberOfMotifs.end();nomIter++) {
          if(nomIter->second>=dims-gaplimit) {
            break;
          }
        }
        int seq=0;

        while(data->getSite(ind[seq],seq).ID!=nomIter->first) {
          seq++;
        }

        char strand=strands[seq];
        siteCode motifID=motifIDs[seq];

        seqs=PyTuple_New(dims);
        coords=PyTuple_New(dims);
        seqPos=PyTuple_New(dims);
        siteScore=PyTuple_New(dims);
        annotations=PyTuple_New(dims);

        int seqC=0;

        //Fetching backtrack
        uint BT[10];
        dynm_Matrix->outBT(ind, BT);

        int last=1;
        for(int k=0;k<dims;k++) {
          if(BT[k]) {
            last=0;
          }
        }

        for(uint i=0;i<sites.size();i++) {
          // Sequence
          PyTuple_SetItem(seqs,seqC,PyInt_FromLong(i));

          // Coordinates on the sequence
          PyObject *coord;
          coord=Py_BuildValue("(ii)",(int)sites[i].pos,(int)sites[i].epos);
          PyTuple_SetItem(coords,seqC,coord);

          if(BT[i] || last) {
            // Position on the site sequence
            PyTuple_SetItem(seqPos,seqC,PyInt_FromLong(ind[i]));

            PyTuple_SetItem(siteScore,seqC,PyFloat_FromDouble(sites[i].weight));

            // Annotation
            PyTuple_SetItem(annotations,seqC,PyString_FromString(sites[i].annot.c_str()));
          } else {
            // Position on the site sequence
            PyTuple_SetItem(seqPos,seqC,PyString_FromString("[ ]"));

            PyTuple_SetItem(siteScore,seqC,PyString_FromString(" "));

            // Annotation
            PyTuple_SetItem(annotations,seqC,PyString_FromString(sites[i].annot.c_str()));
          }

          seqC++; 
        }

        // Resizing might change the location of the tuples
        _PyTuple_Resize(&seqs,seqC);
        _PyTuple_Resize(&coords,seqC);
        _PyTuple_Resize(&seqPos,seqC);
        _PyTuple_Resize(&siteScore,seqC);
        _PyTuple_Resize(&annotations,seqC);

        float score=dynm_Matrix->outScore(ind);

        PyObject *alnRow=PyAln_New_Multi(data->factor(motifID),
          seqs,seqPos,coords,strand,score,siteScore,annotations);

        PyList_Append(goodAlign,alnRow);
        Py_DECREF(alnRow);
        alnRow=NULL;

        //If the whole backtrack is 0, this is the end of the alignment
        //Otherwise the backtrack is followed
        ok=0;
        for (int k=0;k<dims;k++) {
          if(BT[k]) {
            ind[k]-=BT[k];
            ok=1;
          }
        }
      }

      return goodAlign;
    }
  }

  vector<vector<posCode> > givePath() {
    //The path is returned in the order it is found, so [0] has the last cell of the alignment
    vector<vector<posCode> > path;

    vector<posCode> ind;
    ind.resize(dims);
      //Starting from the cell with the best score
    for(int k=0;k<dims;k++) {
      ind[k]=bestInd[k];
    }

    int ok=1;
    while(ok) {
      //Adding the cell to the path
      path.push_back(ind);

      //Fetching backtrack
      //Haetaan backtrack
      int place[10];
      for(int k=0;k<dims;k++) {
        place[k]=ind.at(k);
      }
      uint BT[10];
      dynm_Matrix->outBT(place,BT);

        //If the whole backtrack is 0, this is the end of the alignment
        //Otherwise the backtrack is followed
      ok=0;
      for (int k=0;k<dims;k++) {
        if(BT[k]) {
          ind[k]-=BT[k];
          ok=1;
        }
      }
    }
    return path;
  }

  void fillm_MatrixPart(seqCode i, vector<posCode> start) {
    int motifs=data->sequenceLens(i);
    limit[i]=data->getSite(start[i],i).epos+MAX_BP_DIST;
    //Going through dimensions recursively
    if(i<(this->dims-1)) {
      for(int j=start.at(i);j<motifs && data->getSite(j,i).pos<limit[i];j++) {
        curind[i]=j;
        fillm_MatrixPart(i+1,start);
      }
    } else if(i==(this->dims-1)) {
    //Calculating in the last dimension
      for(int j=start.at(i);j<motifs && data->getSite(j,i).pos<limit[i];j++) {
        curind[i]=j;
        //If the cell is in the list of old paths, the score is -infinity
        int out=0;
        for(int l=0;l<oldpaths.size();l++) {
          for(int k=0;k<dims;k++) {
            if(curind[k]==oldpaths.at(l).at(k)) {
              out++;
            }
          }
          if(out>=dims) {
            break;
          } else {
            out=0;
          }
        }
        if(out) {
          uint BT[10];
          int place[10];
          for(int k=0;k<dims;k++) {
            BT[k]=0;
            place[k]=start.at(k);
          }
          dynm_Matrix->inBT(BT,place);
          dynm_Matrix->inScore(-DBL_MAX,place);
        //Otherwise a new score will be calculated
        } else {
          float score=dynm_Matrix->outScore(curind);

          if(score>0) {
            fillScore();
            float newscore=dynm_Matrix->outScore(curind);
            if(newscore!=score) {
              for(int k=0;k<dims;k++) {
                limit[k]=data->getSite(curind[k],k).epos+MAX_BP_DIST;
              }
            }
          }
        }
      }

      //After the last cell, the list of previous cells is emptied
      fulls.clear();
    }
  }

  void findBest(seqCode i) {
    int motifs=data->sequenceLens(i);

    //Going through dimensions recursively
    if(i<(this->dims-1)) {
      for(int j=0;j<motifs;j++) {
        curind[i]=j;
        findBest(i+1);
      }
    } else if(i==(this->dims-1)) {
    //Calculating in the last dimension
      for(int j=0;j<motifs;j++) {
        curind[i]=j;

        float score=dynm_Matrix->outScore(curind);

        //Saving best result
        if(score>bestScore) {
          bestScore=score;
          for(int k=0;k<dims;k++) {
            bestInd[k]=curind[k];
          }
        }
      }
    }
  }

  PyObject *giveSubBest() {
    //Removing the previous alignment and recalculating the part it affects
    //Fething the cells of the previous alignment
    vector<vector<posCode> > oldPath=givePath();

    for(int i=oldPath.size()-1;i>=0;i--) {
      this->oldpaths.push_back(oldPath.at(i));
      fillm_MatrixPart(0,oldPath.at(i));
    }

    bestScore=-DBL_MAX;
    findBest(0);

    //Returning new best alignment
    return giveBest();
  }

};


typedef struct {
  PyObject_HEAD

  PyObject *bestAlignments;

  PyObject *names;

  double secs_to_align;
  double lambda;
  double xi; 
  double mu;
  double nu;
  double nuc_per_rotation;
  int gaplimit;
  int dims;

  double Rsquared;
} smalign_AlignmentObject;


static PyObject *smalign_aligndata(PyObject *self, PyObject *args) {
  cout<<"This function (smalign_aligndata) hasn't been written as it doesn't seem to get used";
  return NULL;
}

static void smalignment_dealloc(smalign_AlignmentObject *self) {
  if(self->bestAlignments) {
    Py_DECREF(self->bestAlignments);
    self->bestAlignments=NULL;
  }
  if(self->names) {
    Py_DECREF(self->names);
    self->names=NULL;
  }
  self->ob_type->tp_free((PyObject*)self);
}


extern "C" int smalignment_init(smalign_AlignmentObject *self, PyObject *args, PyObject *kwds) {

  double lambda, xi, mu, nu, nuc_per_rotation;
  string firstSeqName,secondSeqName,sequence;
  int gaplimit;
  int numofalign;
  PyObject *data;
  static char *kwlist[] = {"data","gaplimit","numofalign","lambda","xi","mu","nu","nuc_per_rotation",NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oiiddddd", kwlist, 
				  &data,&gaplimit,&numofalign,
			&lambda, &xi, &mu, &nu, &nuc_per_rotation)){
    return -1;
  }

  self->gaplimit=gaplimit;
  self->lambda=lambda;
  self->xi=xi;
  self->mu=mu;
  self->nu=nu;
  self->nuc_per_rotation=nuc_per_rotation;
  self->secs_to_align=0.0;
  self->Rsquared=0;

  tms before,after;
  long ticks_per_sec=sysconf(_SC_CLK_TCK);
  long ticks_to_align;
  times(&before);

  //Multialign
  SimpleMultiAlign *alignObject=new SimpleMultiAlign(data,gaplimit,lambda,mu,nu,xi,nuc_per_rotation);

  if(PyErr_Occurred()) {
    cout<<"Multialign failed\n";
    return -1;
  }

  //Extracting the best alignment
  self->bestAlignments=PyList_New(0);
  if(self->bestAlignments==NULL || PyErr_Occurred()) {
    cout<<"Creating bestAlignments failed\n";
    return -1;
  }
  PyList_Append(self->bestAlignments,alignObject->giveBest());

  if(PyErr_Occurred()) {
    cout<<"Fetching best alignment failed\n";
    return -1;
  }

  for(int i=1;i<numofalign;i++) {
    PyList_Append(self->bestAlignments,alignObject->giveSubBest());
    if(PyErr_Occurred()) {
      cout<<"Fetching a suboptimal alignment failed\n";
      return -1;
    }
  }

  times(&after);
  ticks_to_align=((after.tms_utime-before.tms_utime)+
    (after.tms_stime-before.tms_stime));
  self->secs_to_align+=((double)ticks_to_align)/ticks_per_sec;

  int dims=alignObject->giveDims();
  self->dims=dims;

  //Fetching sequence names
  map<string,seqCode>::iterator iter=alignObject->sequenceNames();
  self->names=PyTuple_New(dims);
  if(PyErr_Occurred()) {
    cout<<"Fetching sequence names failed\n";
    return -1;
  }
  for(int i=0;i<dims;i++,iter++) {
    PyTuple_SET_ITEM(self->names,(int)iter->second,
       PyString_FromString(iter->first.c_str()));
    if(PyErr_Occurred()) {
      cout<<"Handling sequence names failed\n";
      return -1;
    }
  }

  delete alignObject;

  return 0;
}


static PyObject *smalignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {

  smalign_AlignmentObject *self;

  self=(smalign_AlignmentObject *)type->tp_alloc(type,0);

  return (PyObject *)self;
}

static PyMethodDef smalignment_methods[] = {
  {NULL}
};

static PyMemberDef smalignment_members[] = {
  {"bestAlignments",T_OBJECT_EX, offsetof(smalign_AlignmentObject, bestAlignments), 0, "The best alignment"},
  {"names",T_OBJECT_EX, offsetof(smalign_AlignmentObject, names), 0, "The sequence names"},
  {"secs_to_align",T_DOUBLE, offsetof(smalign_AlignmentObject, secs_to_align), 0, "CPU time for Alignment"},
  {"Lambda",T_DOUBLE, offsetof(smalign_AlignmentObject,lambda), 0, "Parameter for bonus"},
  {"Xi",T_DOUBLE, offsetof(smalign_AlignmentObject,xi), 0, "Parameter for distance penalty"},
  {"Nu",T_DOUBLE, offsetof(smalign_AlignmentObject,nu), 0, "Parameter for distance difference penalty"},
  {"Mu",T_DOUBLE, offsetof(smalign_AlignmentObject,mu), 0, "Parameter for rotation penalty"},
  {"nuc_per_rotation",T_DOUBLE, offsetof(smalign_AlignmentObject,nuc_per_rotation), 0, "Parameter for nucleotides per 360 deg rotation of DNA"},
  {"gaplimit",T_INT, offsetof(smalign_AlignmentObject,gaplimit), 0, "Specifies how many sequences can have gaps at the same point."},
  {"Rsquared",T_DOUBLE, offsetof(smalign_AlignmentObject,Rsquared), 0, "This is here just to make this module compatible with the original Output.py"},
  {NULL}
};


static PyTypeObject smalign_AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "shortMultiAlign.ShortMultiAlignment",             /*tp_name*/
    sizeof(smalign_AlignmentObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)smalignment_dealloc,                         /*tp_dealloc*/
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
    "Short multiple alignment object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    smalignment_methods,             /* tp_methods */
    smalignment_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)smalignment_init,      /* tp_init */
    0,                         /* tp_alloc */
    smalignment_new,                 /* tp_new */
};


static PyMethodDef smalignMethods[] = {
  {"aligndata",  smalign_aligndata, METH_VARARGS,
   "aligns computed sequences\nArguments: data,gaplimit,numofalign,lambda,xi,mu,nu,nuc_per_rotation"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

extern "C" void initshortMultiAlign(void) {

  PyObject* m=NULL;

  smalign_AlignmentType.tp_new=PyType_GenericNew;
  if(PyType_Ready(&smalign_AlignmentType)<0)
    return;

  m=Py_InitModule("eellib.shortMultiAlign", smalignMethods);

  if(m==NULL)
    return;

  if(import_alnCols()<0)
    return;

  Py_INCREF(&smalign_AlignmentType);
  PyModule_AddObject(m, "ShortMultiAlignment", (PyObject *)&smalign_AlignmentType);
}
