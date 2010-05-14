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
#include "shortMultiAlign.cc"

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

using namespace std;


typedef struct {
  PyObject_HEAD

  PyObject *bestAlignments;
  PyObject *pwBase;

  PyObject *names;

  double secs_to_align;
  double lambda;
  double xi; 
  double mu;
  double nu;
  double nuc_per_rotation;
  int gaplimit;
  int numofalign;
  int buffer;
  int dims;

  double Rsquared;
} mfp_AlignmentObject;

static vector<vector<vector<posCode> > > getGoodplaces(PyObject *pairwise, m_Inputs *indata, mfp_AlignmentObject *self, int buffer) {

  int dims=indata->sequences();

  vector<map<posind, int> > indCounts;
  indCounts.resize(dims);

  //Every seq has it's own list in the vector, first map has the index of this seq, the second has the index of the associated, named seq
  vector<map<posind,map<seqCode,posind> > > buddies;
  buddies.resize(dims);

  //Going through the pariwise alignment file
  PyObject *pwIter=PyObject_GetIter(pairwise);
  PyObject *pwSite;

  PyObject *goodAlign=PyList_New(0);
  if(goodAlign==NULL || PyErr_Occurred()) {
    cout<<"Creating goodAlign failed\n";
    PyErr_Print();
    return (vector<vector<vector<posCode> > >) NULL;
  }

  PyObject *seqs;
  PyObject *coords;
  PyObject *seqPos,*siteScore,*annotations;

  int numOfPair=-1;
  int numOfAlign=-1;
  int numOfSite=-1;

  vector<string> seqNames;
  seqCode *seqId;
  float score=0;

  vector<map<posind, int> > thisIndCounts;
  vector<vector<posind> > pos;
  vector<map<posind,posind> > posRef;
  posRef.resize(dims);

  vector<PyObject *> pwAligns;
  //First dimension the number of pairwise alignments, second dimension the number of seqs, third dimension the number of cells in the pairwise alignment
  vector<vector<vector<posind> > > alignpos;

  //Counting occurences in pairwise alignments
  while((pwSite=PyIter_Next(pwIter))!=NULL) {
    PyObject *tmp_obj=PySequence_GetItem(pwSite,0);
    string tempstr1=string(PyString_AsString(tmp_obj));
    tmp_obj=PySequence_GetItem(pwSite,2);
    string tempstr2=string(PyString_AsString(tmp_obj));
    if(tempstr1=="###") {
      //New pair
      numOfPair++;
      numOfAlign=-1;
      score=0;
      thisIndCounts.clear();
      thisIndCounts.resize(2);
      pos.clear();
      pos.resize(2);
      pos.at(0).resize(2);
      pos.at(1).resize(2);
      seqNames.clear();
      if(seqId) {
        delete [] seqId;
        seqId=NULL;
      }
      seqId=new seqCode[2];
      for(int i=0;i<2;i++) {
        pwSite=PyIter_Next(pwIter);
        tmp_obj=PySequence_GetItem(pwSite,0);
        string tempstr=string(PyString_AsString(tmp_obj));
        seqNames.push_back(tempstr.erase(0,1));
        seqId[i]=indata->seqIdFromName(tempstr);
      }

    } else if(tempstr2=="CisModule") {
      //New alignment
      numOfAlign++;
      numOfSite=-1;
      score=0;
      //If it's not the first one, saving the last one
      //Jos linjaus ei ole eka, pistetään edellinen talteen
      if(numOfPair>0 || numOfAlign>0) {
        pwAligns.push_back(goodAlign);
        goodAlign=PyList_New(0);
        if(goodAlign==NULL || PyErr_Occurred()) {
          cout<<"Creating new goodAlign failed\n";
          PyErr_Print();
          return (vector<vector<vector<posCode> > >) NULL;
        }
      }
      alignpos.resize(alignpos.size()+1);
      alignpos.back().resize(dims);
    } else {
      //New cell
      numOfSite++;

      seqs=PyTuple_New(dims);
      coords=PyTuple_New(dims);
      seqPos=PyTuple_New(dims);
      siteScore=PyTuple_New(dims);
      annotations=PyTuple_New(dims);

      tmp_obj=PyNumber_Int(PySequence_GetItem(pwSite,3));
      pos[0][0]=(posind)PyInt_AsLong(tmp_obj);
      tmp_obj=PyNumber_Int(PySequence_GetItem(pwSite,4));
      pos[0][1]=(posind)PyInt_AsLong(tmp_obj);
      tmp_obj=PySequence_GetItem(pwSite,2);
      string motifName=string(PyString_AsString(tmp_obj));

      posRef.at(seqId[0]).insert(pair<posind,posind>(pos[0][0],pos[0][1]));
      alignpos.back().at(seqId[0]).push_back(pos[0][0]);

      // Sequence
      PyTuple_SetItem(seqs,0,PyInt_FromLong(seqId[0]));

      // Coordinates on the sequence
      PyObject *coord;
      coord=Py_BuildValue("(ii)",(int)pos[0][0],(int)pos[0][1]);
      PyTuple_SetItem(coords,0,coord);

      // Position on the site sequence
      PyTuple_SetItem(seqPos,0,PyInt_FromLong(indata->getPosCode(pos[0][0],seqId[0],motifName)));

      // Idtriple
      id_triple site=indata->getSite(indata->getPosCode(pos[0][0],seqId[0],motifName),seqId[0]);

      // Weight
      PyTuple_SetItem(siteScore,0,PyFloat_FromDouble(site.weight));

      // Annotation
      PyTuple_SetItem(annotations,0,PyString_FromString(site.annot.c_str()));


      pwSite=PyIter_Next(pwIter);

      tmp_obj=PyNumber_Int(PySequence_GetItem(pwSite,3));
      pos[1][0]=(posind)PyInt_AsLong(tmp_obj);
      tmp_obj=PyNumber_Int(PySequence_GetItem(pwSite,4));
      pos[1][1]=(posind)PyInt_AsLong(tmp_obj);
      tmp_obj=PySequence_GetItem(pwSite,2);
      motifName=string(PyString_AsString(tmp_obj));

      posRef.at(seqId[1]).insert(pair<posind,posind>(pos[1][0],pos[1][1]));
      alignpos.back().at(seqId[1]).push_back(pos[1][0]);

      // Sequence
      PyTuple_SetItem(seqs,1,PyInt_FromLong(seqId[1]));

      // Coordinates on the sequence
      coord=Py_BuildValue("(ii)",(int)pos[1][0],(int)pos[1][1]);
      PyTuple_SetItem(coords,1,coord);

      // Position on the site sequence
      PyTuple_SetItem(seqPos,1,PyInt_FromLong(indata->getPosCode(pos[1][0],seqId[1],motifName)));

      site=indata->getSite(indata->getPosCode(pos[1][0],seqId[1],motifName),seqId[1]);

      // Weight
      PyTuple_SetItem(siteScore,1,PyFloat_FromDouble(site.weight));

      // Annotation
      PyTuple_SetItem(annotations,1,PyString_FromString(site.annot.c_str()));


      score+=PyFloat_AsDouble(PyNumber_Float(PySequence_GetItem(pwSite,5)));

      char strand=site.strand;

      _PyTuple_Resize(&seqs,2);
      _PyTuple_Resize(&coords,2);
      _PyTuple_Resize(&seqPos,2);
      _PyTuple_Resize(&siteScore,2);
      _PyTuple_Resize(&annotations,2);

      PyObject *alnRow=PyAln_New_Multi(indata->factor(site.ID),
        seqs,seqPos,coords,strand,score,siteScore,annotations);
      if(PyErr_Occurred()) {
        cout<<"Creating alnRow failed\n";
        PyErr_Print();
        return (vector<vector<vector<posCode> > >) NULL;
      }

      PyList_Append(goodAlign,alnRow);
      if(PyErr_Occurred()) {
        cout<<"Appending to goodAlign failed\n";
        PyErr_Print();
        return (vector<vector<vector<posCode> > >) NULL;
      }
      Py_DECREF(alnRow);
      alnRow=NULL;

      //Counting how many times a spot has occured in alignments of this pair
      if(thisIndCounts.at(0).find(pos[0][0])!=thisIndCounts.at(0).end()) {
        thisIndCounts.at(0)[pos[0][0]]++;
      } else {
        thisIndCounts.at(0).insert(pair<posind,int>(pos[0][0],1));

        //Counting how many different pairs have had alignments where this index occurs
        if(indCounts.at(seqId[0]).find(pos[0][0])!=indCounts.at(seqId[0]).end()) {
          indCounts.at(seqId[0])[pos[0][0]]++;
        } else {
          indCounts.at(seqId[0]).insert(pair<posind,int>(pos[0][0],1));
        }
      }
      buddies.at(seqId[0])[pos[0][0]][seqId[1]]=pos[1][0];

      if(thisIndCounts.at(1).find(pos[1][0])!=thisIndCounts.at(1).end()) {
        thisIndCounts.at(1)[pos[1][0]]++;
      } else {
        thisIndCounts.at(1).insert(pair<posind,int>(pos[1][0],1));

        //Counting how many different pairs have had alignments where this index occurs
        if(indCounts.at(seqId[1]).find(pos[1][0])!=indCounts.at(seqId[1]).end()) {
          indCounts.at(seqId[1])[pos[1][0]]++;
        } else {
          indCounts.at(seqId[1]).insert(pair<posind,int>(pos[1][0],1));
        }
      }
      buddies.at(seqId[1])[pos[1][0]][seqId[0]]=pos[0][0];
    }
    Py_DECREF(tmp_obj);
    tmp_obj=NULL;
  }
  delete [] seqId;
  seqId=NULL;
  Py_DECREF(pwIter);
  pwIter=NULL;
  if(pwSite) {
    Py_DECREF(pwSite);
    pwSite=NULL;
  }
  pwAligns.push_back(goodAlign);

  vector<vector<posind> > hits;
  hits.resize(dims);
  for(int i=0;i<dims;i++) {
    map<posind,int>::iterator it;

    //Fetching seq indexes, that occur in more than one pairwise alignment
    for(it=indCounts.at(i).begin();it!=indCounts.at(i).end();it++) {
      if(it->second>1) {
        hits.at(i).push_back(it->first);
      }
    }

    sort(hits.at(i).begin(),hits.at(i).end());
  }

  //Grouping the spots, anything in MAX_BP_DIST of each other belongs in the same region
  vector<vector<vector<posind> > > posvector;
  int ind=0;
  vector<int> changeind;
  changeind.push_back(0);
  //Going through seqs
  for(seqCode i=0;i<dims;i++) {
    //Going through seq hits
    for(uint j=0;j<hits.at(i).size();j++) {
      //If starting new region, add to posvector
      if(ind>=(int)posvector.size()) {
        posvector.resize(ind+1);
        posvector.at(ind).resize(dims);
        for(int l=0;l<dims;l++) {
          posvector.at(ind).at(l).resize(2);
        }
        //Initializing places in posvector, first and last motif are the same
        for(int x=0;x<2;x++) {
          //The seq being handled is an easy thing
          posvector.at(ind).at(i).at(x)=hits.at(i).at(j);
          //For other seqs the buddylist must be checked: looking for a spot that concerns this seq and place
          map<posind,map<seqCode,posind> >::iterator posIter=buddies.at(i).find(posvector.at(ind).at(i).at(x));
          for(int k=0;k<dims;k++) {
            if(k!=i) {
              //Looking for the place in the other seq from the second part of buddylist
              map<seqCode,posind>::iterator seqIter=posIter->second.find(k);
              if(seqIter!=posIter->second.end()) {
                posvector.at(ind).at(k).at(x)=seqIter->second;
              }
            }
          }
        }
      } else {
        //If the next spot is close enough, put that as the new endpoint
        if(hits.at(i).at(j)-posRef.at(i)[posvector.at(ind).at(i).at(1)]<MAX_BP_DIST) {
          posvector.at(ind).at(i).at(1)=hits.at(i).at(j);
          //Checking the effect on other seqs
          map<posind,map<seqCode,posind> >::iterator posIter=buddies.at(i).find(posvector.at(ind).at(i).at(1));
          for(int k=0;k<dims;k++) {
            if(k!=i) {
              map<seqCode,posind>::iterator seqIter=posIter->second.find(k);
              if(seqIter!=posIter->second.end()) {
                if(seqIter->second > posvector.at(ind).at(k).at(1)) {
                  posvector.at(ind).at(k).at(1)=seqIter->second;
                }
                if(seqIter->second < posvector.at(ind).at(k).at(0) || posvector.at(ind).at(k).at(0)==0) {
                  posvector.at(ind).at(k).at(0)=seqIter->second;
                }
              }
            }
          }
        } else {
          ind++;
          posvector.resize(ind+1);
          posvector.at(ind).resize(dims);
          for(int l=0;l<dims;l++) {
            posvector.at(ind).at(l).resize(2);
          }
          for(int x=0;x<2;x++) {
            posvector.at(ind).at(i).at(x)=hits.at(i).at(j);
            map<posind,map<seqCode,posind> >::iterator posIter=buddies.at(i).find(posvector.at(ind).at(i).at(x));
            for(int k=0;k<dims;k++) {
              if(k!=i) {
                map<seqCode,posind>::iterator seqIter=posIter->second.find(k);
                if(seqIter!=posIter->second.end()) {
                  posvector.at(ind).at(k).at(x)=seqIter->second;
                }
              }
            }
          }
        }
      }
    }
    if(hits.at(i).size()>0) {
      ind++;
    }
    changeind.push_back(ind);
  }

  //Removing groups that don't have all dimensions
  vector<int> todelete;
  for(uint i=0;i<posvector.size();i++) {
    int count=0;
    for(uint j=0;j<posvector.at(i).size();j++) {
      if(posvector.at(i).at(j).at(0)) {
        count++;
      }
    }
    if(count<dims) {
      todelete.push_back(i);
    }
  }
  vector<vector<vector<posind> > >::iterator veciter=posvector.end();
  int vecind=posvector.size();
  vector<int>::iterator matchiter=todelete.end();
  while(matchiter>todelete.begin()--) {
    matchiter--;
    while(*matchiter!=vecind) {
      veciter--;
      vecind--;
    }
    posvector.erase(veciter);
  }

  //Changing the saved indexes to correspond to the thinned out vector
  for(int i=(int)todelete.size()-1;i>=0;i--) {
    for(int j=(int)changeind.size()-1;j>=0;j--) {
      if(changeind.at(j)>todelete.at(i)) {
        changeind.at(j)--;
      } else {
        break;
      }
    }
  }

  vector<vector<vector<posCode> > > retvector;

  if(posvector.size()>0) {
    //Making a list of the pairwise alignments, which the area for multiple alignment is based on
    self->pwBase=PyList_New(0);
    if(PyErr_Occurred()) {
      cout<<"Parilinjauslistan luominen epäonnistui\n";
      PyErr_Print();
      return (vector<vector<vector<posCode> > >) NULL;
    }
    //Going through seqs
    for(int i=0;i<dims;i++) {
      //Going through the part of posvector that corresponds to the current seq and is found in changeind
      for(int j=changeind.at(i);j<changeind.at(i+1);j++) {
        PyObject *pwList=PyList_New(0);
        //Going through obtained pairwise alignments
        for(int k=0;k<(int)alignpos.size();k++) {
          //Going through rows of current pairwise alignment
          for(int m=0;m<(int)alignpos.at(k).at(i).size();m++) {
            //If the current row is inside a good area, add the pairwise alignment to the list and stop going through it
            if(alignpos.at(k).at(i).at(m)>=posvector.at(j).at(i).at(0) && alignpos.at(k).at(i).at(m)<=posvector.at(j).at(i).at(1)) {
              PyList_Append(pwList,pwAligns.at(k));
              break;
            }
          }
        }
        PyList_Append(self->pwBase,pwList);
      }
    }

    //First dimension the number of potential areas, second dimension the number of seqs, third dimension 2 (start and end)
    retvector.resize(posvector.size());
    for(int i=0;i<(int)posvector.size();i++) {
      retvector.at(i).resize(dims);
      for(int j=0;j<dims;j++) {
        for(int k=0;k<2;k++) {
          retvector.at(i).at(j).push_back(indata->getPosCode(posvector[i][j][k],j));
        }
      }
    }

    //Adding a buffer to the ends
    for(int j=0;j<(int)retvector.size();j++) {
      for(int i=0;i<dims;i++) {
        int k=0;
        while(retvector.at(j).at(i).at(0)-k>=0 && indata->getSite(retvector.at(j).at(i).at(0),i).pos-indata->getSite(retvector.at(j).at(i).at(0)-k,i).epos<buffer) {
          k++;
        }
        retvector.at(j).at(i).at(0)-=k-1;
        k=0;
        while(retvector.at(j).at(i).at(1)+k<indata->sequenceLens().at(i) && indata->getSite(retvector.at(j).at(i).at(1)+k,i).pos-indata->getSite(retvector.at(j).at(i).at(1),i).epos<buffer) {
          k++;
        }
        retvector.at(j).at(i).at(1)+=k-1;
      }
    }

  } else {
    cout<<"No good places found\n";
  }
  return retvector;
}



static void mfp_dealloc(mfp_AlignmentObject *self) {
  if(self->bestAlignments) {
    Py_DECREF(self->bestAlignments);
    self->bestAlignments=NULL;
  }
  if(self->pwBase) {
    Py_DECREF(self->pwBase);
    self->pwBase=NULL;
  }
  if(self->names) {
    Py_DECREF(self->names);
    self->names=NULL;
  }
  self->ob_type->tp_free((PyObject*)self);
}


extern "C" int mfp_init(mfp_AlignmentObject *self, PyObject *args, PyObject *kwds) {

  double lambda, xi, mu, nu, nuc_per_rotation;
  int gaplimit;
  int numofalign;
  int buffer;
  PyObject *data;
  PyObject *pairwise;
  static char *kwlist[] = {"data","pairwise","gaplimit","numofalign","buffer","lambda","xi","mu","nu","nuc_per_rotation",NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOiiiddddd", kwlist,&data,&pairwise,&gaplimit,&numofalign,&buffer,&lambda, &xi, &mu, &nu, &nuc_per_rotation)){
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
  self->numofalign=numofalign;
  self->buffer=buffer;

  m_Inputs *indata=new m_Inputs(data);

  if(PyErr_Occurred()) {
    cout<<"Creating m_Inputs failed\n";
    PyErr_Print();
    return -1;
  }

  int dims=indata->sequences();
  self->dims=dims;

  //Searching for potential areas for multiple alignment from pairwise alignments
  vector<vector<vector<posCode> > > goodplaces=getGoodplaces(pairwise,indata,self,buffer);

  if(PyErr_Occurred()) {
    cout<<"Searching for good places failed\n";
    PyErr_Print();
    return -1;
  }

  if(goodplaces.size()==0) {
    return -1;
  }

  for(uint i=0;i<goodplaces.size();i++) {
    for(uint j=0;j<goodplaces.at(i).size();j++) {
      for(uint k=0;k<goodplaces.at(i).at(j).size();k++) {
        cout<<""<<goodplaces.at(i).at(j).at(k)<<" ";
      }
      cout<<"\n";
    }
    cout<<"\n";
  }

  tms before,after;
  long ticks_per_sec=sysconf(_SC_CLK_TCK);
  long ticks_to_align;
  times(&before);

  SimpleMultiAlign *alignObject;

  self->bestAlignments=PyList_New(0);
  if(self->bestAlignments==NULL || PyErr_Occurred()) {
    cout<<"Creating bestAlignments failed\n";
    PyErr_Print();
    return -1;
  }

  //Going through found potential areas
  for(uint i=0;i<goodplaces.size();i++) {
    m_Inputs *curdata=indata->getPart(goodplaces.at(i));

    //Multialign
    alignObject=new SimpleMultiAlign(curdata,gaplimit,lambda,mu,nu,xi,nuc_per_rotation);

    if(PyErr_Occurred()) {
      cout<<"Multialign failed\n";
      PyErr_Print();
      return -1;
    }

    //Extracting the best alignment
    PyList_Append(self->bestAlignments,alignObject->giveBest());
    for(int j=1;j<numofalign;j++) {
      PyList_Append(self->bestAlignments,alignObject->giveSubBest());
      if(PyErr_Occurred()) {
        cout<<"Extracting best result failed\n";
        return -1;
      }
    }

    delete alignObject;
  }

  times(&after);
  ticks_to_align=((after.tms_utime-before.tms_utime)+
    (after.tms_stime-before.tms_stime));
  self->secs_to_align+=((double)ticks_to_align)/ticks_per_sec;

  if(PyErr_Occurred()) {
    cout<<"An error spotted after timing has been stopped\n";
    PyErr_Print();
    return -1;
  }

  //Fetching sequence names
  map<string,seqCode>::iterator iter=indata->sequenceNames();
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

  return 0;
}

static PyObject *mfp_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {

  mfp_AlignmentObject *self;

  self=(mfp_AlignmentObject *)type->tp_alloc(type,0);

  return (PyObject *)self;
}

static PyMethodDef mfp_methods[] = {
  {NULL}
};

static PyMemberDef mfp_members[] = {
  {"bestAlignments",T_OBJECT_EX, offsetof(mfp_AlignmentObject, bestAlignments), 0, "The best alignments"},
  {"pwBase",T_OBJECT_EX, offsetof(mfp_AlignmentObject, pwBase), 0, "The pairwise alignments the multiple alignments are based on"},
  {"names",T_OBJECT_EX, offsetof(mfp_AlignmentObject, names), 0, "The sequence names"},
  {"secs_to_align",T_DOUBLE, offsetof(mfp_AlignmentObject, secs_to_align), 0, "CPU time for Alignment"},
  {"Lambda",T_DOUBLE, offsetof(mfp_AlignmentObject,lambda), 0, "Parameter for bonus"},
  {"Xi",T_DOUBLE, offsetof(mfp_AlignmentObject,xi), 0, "Parameter for distance penalty"},
  {"Nu",T_DOUBLE, offsetof(mfp_AlignmentObject,nu), 0, "Parameter for distance difference penalty"},
  {"Mu",T_DOUBLE, offsetof(mfp_AlignmentObject,mu), 0, "Parameter for rotation penalty"},
  {"nuc_per_rotation",T_DOUBLE, offsetof(mfp_AlignmentObject,nuc_per_rotation), 0, "Parameter for nucleotides per 360 deg rotation of DNA"},
  {"gaplimit",T_INT, offsetof(mfp_AlignmentObject,gaplimit), 0, "Specifies how many sequences can have gaps at the same point."},
  {"numofalign",T_INT, offsetof(mfp_AlignmentObject,numofalign), 0, "Tells how many aligments per likely spot are calculated."},
  {"buffer",T_INT, offsetof(mfp_AlignmentObject,buffer), 0, "The size of the buffer area on both sides of found areas to multialign."},
  {"Rsquared",T_DOUBLE, offsetof(smalign_AlignmentObject,Rsquared), 0, "This is here just to make this module compatible with the original Output.py"},
  {NULL}
};


static PyTypeObject mfp_AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "multiFromPairwise.MultiFromPairwise",             /*tp_name*/
    sizeof(mfp_AlignmentObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)mfp_dealloc,                         /*tp_dealloc*/
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
    "Multifrompairwise alignment object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    mfp_methods,             /* tp_methods */
    mfp_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)mfp_init,      /* tp_init */
    0,                         /* tp_alloc */
    mfp_new,                 /* tp_new */
};

static PyObject *mfp_aligndata(PyObject *self, PyObject *args) {
  cout<<"This function (mfp_aligndata) hasn't been written as it doesn't seem to get used";
  return NULL;
}

static PyMethodDef mfpMethods[] = {
  {"mfpaligndata",  mfp_aligndata, METH_VARARGS,
   "aligns computed sequences\nArguments: data,gaplimit,numofalign,buffer,lambda,xi,mu,nu,nuc_per_rotation"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

extern "C" void initmultiFromPairwise(void) {
  PyObject* m=NULL;

  mfp_AlignmentType.tp_new=PyType_GenericNew;
  if(PyType_Ready(&mfp_AlignmentType)<0)
    return;

  m=Py_InitModule("eellib.multiFromPairwise", mfpMethods);

  if(m==NULL)
    return;

  if(import_alnCols()<0)
    return;

  Py_INCREF(&mfp_AlignmentType);
  PyModule_AddObject(m, "MultiFromPairwise", (PyObject *)&mfp_AlignmentType);
}
