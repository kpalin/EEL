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
typedef unsigned long int ulint;

using namespace std;


class TreeNode {
  TreeNode *left;
  TreeNode *right;
  TreeNode *parent;

  //If the node is a leaf, the ID of the corresponding sequence is saved
  seqCode seqID;

  //If the node is not a leaf, a identifying number is saved so results can be asked later
  int number;

  //If the node is not a leaf, a set of local alignments will be computed for it and saved
  map<float, map<int, map<seqCode, vector<id_triple> > > > lAligns;
  map<float, map<int, vector<float> > > alignScores;

  //If representative sequences are used, they are saved here
  map<float, map<int, vector<id_triple> > > repSeqs;

  //If the results are asked, a list of the sequences is needed and will be saved here
  map<string,seqCode> names;

  public:
  TreeNode();

  TreeNode(int number) {
    this->number=number;
    parent=NULL;
    left=NULL;
    right=NULL;
    seqID=-1;
  }

  TreeNode(int number, TreeNode *parent) {
    this->parent=parent;
    this->number=number;

    left=NULL;
    right=NULL;

    seqID=-1;
  }

  TreeNode(TreeNode *parent, seqCode seqID) {
    this->parent=parent;
    number=-1;
    left=NULL;
    right=NULL;

    this->seqID=seqID;
  }

  ~TreeNode() {
    if(parent) {
      delete parent;
      parent=NULL;
    }
    if(right) {
      delete right;
      right=NULL;
    }
    if(left) {
      delete left;
      left=NULL;
    }
  }

  void addChild(TreeNode *child) {
    if(!right) {
      this->right=child;
    } else if (!left){
      this->left=child;
    } else {
      cout<<"Error: Only binary trees are accepted\n";
    }
  }

  int isLeaf() {
    return seqID>=0;
  }

  void empty() {
    lAligns.clear();
    alignScores.clear();
  }

  map<float, map<int, map<seqCode, vector<id_triple> > > > *getSet() {
    return &lAligns;
  }

  map<float, map<int, vector<float> > > getScores() {
    return alignScores;
  }

  seqCode getSeq() {
    return seqID;
  }

  map<string,seqCode>::iterator getSeqs(m_Inputs *data) {
    if (this->names.size()) {
      return this->names.begin();
    } else {
      map<float, map<int, map<seqCode, vector<id_triple> > > >::iterator it1=lAligns.begin();
      map<int, map<seqCode, vector<id_triple> > >::iterator it2=it1->second.begin();
      map<seqCode, vector<id_triple> >::iterator it3;
      for(it3=it2->second.begin();it3!=it2->second.end();it3++) {
        this->names[data->sequence(it3->first)]=it3->first;
      }
      return this->names.begin();
    }
  }

  int getSeqNum() {
    return this->names.size();
  }

  map<float, map<int, vector<id_triple> > > getRepSeqs() {
    return repSeqs;
  }

  void printTree() {
    if(seqID>=0) {
      cout<<"l"<<seqID;
    } else {
      cout<<"s"<<number;
    }
    cout<<":(";
    if(left) {
      left->printTree();
    }
    cout<<",";
    if(right) {
      right->printTree();
    }
    cout<<")";
  }

  void sumOfPairsScore(int i, int j, map<seqCode, vector<id_triple> > first, map<seqCode, vector<id_triple> > second, m_Matrix *matrices, double lambda, double mu, double nu, double xi, double nuc_per_rot) {
    map<seqCode, vector<id_triple> >::iterator it1=first.begin();
    map<seqCode, vector<id_triple> >::iterator it2=second.begin();
    if(it1->second.at(i).ID==it2->second.at(j).ID) {
      //If both seqs/alignments have the same motif, a score will be computed

      //Fetching previous motifs from the seqs
      map<seqCode, vector<id_triple> > prev1;
      int indgap1=0;
      map<seqCode, vector<id_triple> > prev2;
      int indgap2=0;

      if(i>0) {
        bool ok=1;
        int pos=i-1;
        while(ok && pos>=0) {
          if(it1->second.at(pos).epos<it1->second.at(i).pos) {
            int ok2=1;
            for(it1=first.begin();it1!=first.end();it1++) {
              if(it1->second.at(i).pos-it1->second.at(pos).epos>=MAX_BP_DIST) {
                ok2=0;
                break;
              }
            }
            if(ok2) {
              for(it1=first.begin();it1!=first.end();it1++) {
                prev1[it1->first].push_back(it1->second.at(pos));
              }
              it1=first.begin();
            } else {
              ok=0;
            }
          } else {
            indgap1++;
          }
          pos--;
        }
      }
      if(j>0) {
        bool ok=1;
        int pos=j-1;
        while(ok && pos>=0) {
          if(it2->second.at(pos).epos<it2->second.at(j).pos) {
            int ok2=1;
            for(it2=second.begin();it2!=second.end();it2++) {
              if(it2->second.at(j).pos-it2->second.at(pos).epos>=MAX_BP_DIST) {
                ok2=0;
                break;
              }
            }
            if(ok2) {
              for(it2=second.begin();it2!=second.end();it2++) {
                prev2[it2->first].push_back(it2->second.at(pos));
              }
              it2=second.begin();
            } else {
              ok=0;
            }
          } else {
            indgap2++;
          }
          pos--;
        }
      }

      //The positive part of the score, depending only on the current cell
      float base=0;
      for(it1=first.begin();it1!=first.end();it1++) {
        base+=it1->second.at(i).weight;
      }
      for(it2=second.begin();it2!=second.end();it2++) {
        base+=it2->second.at(j).weight;
      }
      base*=lambda*(first.size()+second.size()-1);

      it1=first.begin();
      it2=second.begin();

      float maxscore=0;
      vector<int> maxbacktrack;
      maxbacktrack.resize(2,0);
      //The penalty, depending on both the current and previous cells
      if(prev1.size() && prev2.size()) {
        map<seqCode, vector<id_triple> >::iterator itp1=prev1.begin();
        map<seqCode, vector<id_triple> >::iterator itp2=prev2.begin();
        for(int a=0;a<itp1->second.size();a++) {
          itp2=prev2.begin();
          for(int b=0;b<itp2->second.size();b++) {
            if(itp1->second.at(a).ID==itp2->second.at(b).ID) {
              float penalty=0;
              float deltaphi;
              for(it1=first.begin();it1!=first.end();it1++) {
                for(it2=second.begin();it2!=second.end();it2++) {
                  int x;
                  x=(int)(it1->second.at(i).pos-prev1[it1->first].at(a).epos-1);
                  int y;
                  y=(int)(it2->second.at(j).pos-prev2[it2->first].at(b).epos-1);
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
              it1=first.begin();
              it2=second.begin();

              int place[10];
              place[0]=i-a-indgap1-1;
              place[1]=j-b-indgap2-1;
              float curscore=matrices->outScore(place)+penalty;

              if(curscore>maxscore) {
                maxscore=curscore;
                maxbacktrack.at(0)=a+indgap1+1;
                maxbacktrack.at(1)=b+indgap2+1;
              }
            }
          }
        }
      }

      float score=maxscore+base;

      int place[10];
      place[0]=i;
      place[1]=j;
      matrices->inScore(score,place);
      uint BT[10];
      BT[0]=(uint)maxbacktrack.at(0);
      BT[1]=(uint)maxbacktrack.at(1);
      matrices->inBT(BT,place);
      return;
    } else {
      //If the motifs are different, score is -infinity
      int place[10];
      place[0]=i;
      place[1]=j;
      matrices->inScore(-DBL_MAX,place);
      uint BT[10];
      BT[0]=0;
      BT[1]=0;
      matrices->inBT(BT,place);
      return;
    }
  }

  void generalizedScore(int i, int j, map<seqCode, vector<id_triple> > first, map<seqCode, vector<id_triple> > second, m_Matrix *matrices, double lambda, double mu, double nu, double xi, double nuc_per_rot) {
    map<seqCode, vector<id_triple> >::iterator it1=first.begin();
    map<seqCode, vector<id_triple> >::iterator it2=second.begin();
    if(it1->second.at(i).ID==it2->second.at(j).ID) {
      //If both seqs/alignments have the same motif, a score will be computed

      //Fetching previous motifs from the seqs
      map<seqCode, vector<id_triple> > prev1;
      int indgap1=0;
      map<seqCode, vector<id_triple> > prev2;
      int indgap2=0;

      if(i>0) {
        bool ok=1;
        int pos=i-1;
        while(ok && pos>=0) {
          if(it1->second.at(pos).epos<it1->second.at(i).pos) {
            int ok2=1;
            for(it1=first.begin();it1!=first.end();it1++) {
              if(it1->second.at(i).pos-it1->second.at(pos).epos>=MAX_BP_DIST) {
                ok2=0;
                break;
              }
            }
            if(ok2) {
              for(it1=first.begin();it1!=first.end();it1++) {
                prev1[it1->first].push_back(it1->second.at(pos));
              }
              it1=first.begin();
            } else {
              ok=0;
            }
          } else {
            indgap1++;
          }
          pos--;
        }
      }
      if(j>0) {
        bool ok=1;
        int pos=j-1;
        while(ok && pos>=0) {
          if(it2->second.at(pos).epos<it2->second.at(j).pos) {
            int ok2=1;
            for(it2=second.begin();it2!=second.end();it2++) {
              if(it2->second.at(j).pos-it2->second.at(pos).epos>=MAX_BP_DIST) {
                ok2=0;
                break;
              }
            }
            if(ok2) {
              for(it2=second.begin();it2!=second.end();it2++) {
                prev2[it2->first].push_back(it2->second.at(pos));
              }
              it2=second.begin();
            } else {
              ok=0;
            }
          } else {
            indgap2++;
          }
          pos--;
        }
      }

      //The positive part of the score, depending only on the current cell
      float base=0;
      for(it1=first.begin();it1!=first.end();it1++) {
        base+=it1->second.at(i).weight;
      }
      for(it2=second.begin();it2!=second.end();it2++) {
        base+=it2->second.at(j).weight;
      }
      base*=lambda*(first.size()+second.size()-1);

      float maxscore=-DBL_MAX;
      vector<int> maxbacktrack;
      maxbacktrack.resize(2,0);

      //The penalty, depending on both the current and previous cells
      if(prev1.size() && prev2.size()) {
        map<seqCode, vector<id_triple> >::iterator itp1=prev1.begin();
        map<seqCode, vector<id_triple> >::iterator itp2=prev2.begin();
        for(int a=0;a<itp1->second.size();a++) {
          itp2=prev2.begin();
          for(int b=0;b<itp2->second.size();b++) {
            if(itp1->second.at(a).ID==itp2->second.at(b).ID) {
              //Calculating the distances
              vector<int> dists;
              for(it1=first.begin();it1!=first.end();it1++) {
                dists.push_back((int)(it1->second.at(i).pos-prev1[it1->first].at(a).epos-1));
              }
              for(it2=second.begin();it2!=second.end();it2++) {
                dists.push_back((int)(it2->second.at(j).pos-prev2[it2->first].at(b).epos-1));
              }

              int distsum=0;
              for(int c=0;c<dists.size();c++) {
                distsum+=dists.at(c);
              }
              float distmean=distsum/dists.size();

              float deltaphisqrd=0;
              for(int c=0;c<dists.size();c++) {
                float curdeltaphi;
                if(distmean-dists.at(c)>=0) {
                    curdeltaphi=(distmean-dists.at(c))*2*PI/nuc_per_rot;
                    curdeltaphi+=PI;
                    curdeltaphi=fmod(curdeltaphi,2*PI);
                    curdeltaphi-=PI;
                  } else {
                    curdeltaphi=(distmean-dists.at(c))*2*PI/nuc_per_rot;
                    curdeltaphi-=PI;
                    curdeltaphi=fmod(curdeltaphi,2*PI);
                    curdeltaphi+=PI;
                }
                deltaphisqrd+=curdeltaphi*curdeltaphi;
              }              

              float deltadistsqrd=0;
              for(int c=0;c<dists.size();c++) {
                deltadistsqrd+=(distmean-dists.at(c))*(distmean-dists.at(c));
              }

              float penalty=-(mu*distmean + nu*deltadistsqrd/distsum + xi*deltaphisqrd/distsum);

              int place[10];
              place[0]=i-a-indgap1-1;
              place[1]=j-b-indgap2-1;
              float curscore=matrices->outScore(place)+penalty;

              if(curscore>maxscore) {
                maxscore=curscore;
                maxbacktrack.at(0)=a+indgap1+1;
                maxbacktrack.at(1)=b+indgap2+1;
              }
            }
          }
        }
      }

      if(maxscore<0) {
        maxscore=base;
        maxbacktrack.at(0)=0;
        maxbacktrack.at(1)=0;
      } else {
        maxscore+=base;
      }

      int place[10];
      place[0]=i;
      place[1]=j;
      matrices->inScore(maxscore,place);
      uint BT[10];
      BT[0]=(uint)maxbacktrack.at(0);
      BT[1]=(uint)maxbacktrack.at(1);
      matrices->inBT(BT,place);
      return;
    } else {
      //If the motifs are different, score is -infinity
      int place[10];
      place[0]=i;
      place[1]=j;
      matrices->inScore(-DBL_MAX,place);
      uint BT[10];
      BT[0]=0;
      BT[1]=0;
      matrices->inBT(BT,place);

      return;
    }
  }

  void getRepAligns(map<float, map<int, vector<vector<int> > > > *aligns, map<float, map<int, vector<float> > > *alignScores,vector<id_triple> first, vector<id_triple> second, int k, int method, double lambda, double mu, double nu, double xi, double nuc_per_rot, runningInd *alignInd) {

    m_Inputs *curdata=new m_Inputs(first,second);
    int gaplimit=0;
    SimpleMultiAlign *alignObject=new SimpleMultiAlign(curdata,gaplimit,lambda,mu,nu,xi,nuc_per_rot);
    map<float, map<int, map<seqCode, vector<id_triple> > > > alignment;
    map<float, map<int, vector<float> > > scores;
    vector<seqCode> seqs;
    seqs.push_back(0);
    seqs.push_back(1);

    alignObject->giveBestAligns(&alignment, &scores, k, seqs, alignInd);

    delete alignObject;

    //The "real" alignment is changed into an aligments of indexes
    map<float, map<int, map<seqCode, vector<id_triple> > > >::iterator it;
    map<float, map<int, vector<float> > >::iterator its=scores.begin();
    vector<vector<int> > tempVec;
    for(it=alignment.begin();it!=alignment.end();it++) {
      map<int, map<seqCode, vector<id_triple> > >::iterator it2;
      map<int, vector<float> >::iterator its2=its->second.begin();
      for(it2=it->second.begin();it2!=it->second.end();it2++) {
        map<seqCode, vector<id_triple> >::iterator it3;
        vector<int> tempVec2;
        for(it3=it2->second.begin();it3!=it2->second.end();it3++) {
          int ind=0;
          for(int i=0; i<it3->second.size();i++) {
            bool found=0;
            while(!found) {
              if(it3->first==0 && ind<first.size()) {
                if(it3->second.at(i)==first.at(ind)) {
                  found=1;
                }
              }
              if(it3->first==1 && ind<second.size()) {
                if(it3->second.at(i)==second.at(ind)) {
                  found=1;
                }
              }
              if(!found) {
                ind++;
              }
              if(ind>=first.size() && ind>=second.size()) {
                cout<<"Error: id_triple wasn't found in either sequence\n";
                break;
              }
            }
            tempVec2.push_back(ind);
          }
          tempVec.push_back(tempVec2);
          tempVec2.clear();
        }
        (*aligns)[it->first][it2->first]=tempVec;
        tempVec.clear();
        (*alignScores)[its->first][its2->first]=its2->second;
        its2++;
      }
      its++;
    }
  }


  void getAligns(map<float, map<int, map<seqCode, vector<id_triple> > > > *aligns, map<float, map<int, vector<float> > > *alignScores,map<seqCode, vector<id_triple> > first, map<seqCode, vector<id_triple> > second, int k, int method, double lambda, double mu, double nu, double xi, double nuc_per_rot, vector<seqCode> seqs, runningInd *alignInd) {

    m_Matrix *matrices;
    vector<int> sizes;

    map<seqCode, vector<id_triple> >::iterator it1=first.begin();
    map<seqCode, vector<id_triple> >::iterator it2=second.begin();
    sizes.push_back(it1->second.size());
    sizes.push_back(it2->second.size());
    matrices=new m_Matrix(2, sizes);

    for(int i=0;i<it1->second.size();i++) {
      for(int j=0;j<it2->second.size();j++) {
        if(method==1) {
          sumOfPairsScore(i,j,first,second,matrices,lambda,mu,nu,xi,nuc_per_rot);
        } else if (method==2) {
          generalizedScore(i,j,first,second,matrices,lambda,mu,nu,xi,nuc_per_rot);
        }
      }
    }

    map<seqCode, vector<id_triple> > revalign;
    vector<float> revScore;

    //Finding the best score
    float bestscore=-DBL_MAX;
    vector<int> ind;
    ind.resize(2);
    int place[10];

    for(int i=0;i<sizes.at(0);i++) {
      place[0]=i;
      for(int j=0;j<sizes.at(1);j++) {
        place[1]=j;
        if(matrices->outScore(place)>bestscore) {
          bestscore=matrices->outScore(place);
          ind.at(0)=i;
          ind.at(1)=j;
        }
      }
    }

    if(bestscore>=0) {
      //Saving the best alignment
      bool ok=1;
      vector<vector<int> > alignInds;
      while(ok) {
        for(it1=first.begin();it1!=first.end();it1++) {
          revalign[it1->first].push_back(it1->second.at(ind.at(0)));
          alignInds.push_back(ind);
        }
        for(it2=second.begin();it2!=second.end();it2++) {
          revalign[it2->first].push_back(it2->second.at(ind.at(1)));
          alignInds.push_back(ind);
        }
        place[0]=ind.at(0);
        place[1]=ind.at(1);
        revScore.push_back(matrices->outScore(place));

        //Jos linjaus jatkuu, siirrytään taaksepäin
        ok=0;
        uint BT[10];
        matrices->outBT(place,BT);
        if(BT[0]) {
          ind.at(0)-=(int)BT[0];
          ok=1;
        }
        if(BT[1]) {
          ind.at(1)-=(int)BT[1];
          ok=1;
        }
      }

      map<seqCode, vector<id_triple> > align;
      vector<float> aScore;
      for(int i=(revScore.size()-1);i>=0;i--) {
        aScore.push_back(revScore.at(i));
      }
      map<seqCode, vector<id_triple> >::iterator it;
      for(it=revalign.begin();it!=revalign.end();it++) {
        for(int i=(it->second.size()-1);i>=0;i--) {
          align[it->first].push_back(it->second.at(i));
        }
      }
      int aind=alignInd->giveInd();
      (*aligns)[bestscore][aind]=align;
      (*alignScores)[bestscore][aind]=aScore;

      //Fetching suboptimal results
      for(int i=1;i<k;i++) {
        for(int j=(alignInds.size()-1);j>=0;j--) {
          int place[10];
          place[0]=alignInds.at(j).at(0);
          place[1]=alignInds.at(j).at(1);
          matrices->inScore(-DBL_MAX, place);
          uint BT[10];
          BT[0]=0;
          BT[1]=0;
          matrices->inBT(BT,place);

          vector<int> rind=alignInds.at(j);
          rind.at(1)++;
          bool ok1;
          it2=second.begin();
          if(rind.at(1)<it2->second.size()) {
            ok1=1;
          } else {
            ok1=0;
          }
          while(ok1) {
            bool ok2=1;
            while(ok2) {
              if(method==1) {
                sumOfPairsScore(rind.at(0),rind.at(1),first,second,matrices,lambda,mu,nu,xi,nuc_per_rot);
              } else if (method==2) {
                generalizedScore(rind.at(0),rind.at(1),first,second,matrices,lambda,mu,nu,xi,nuc_per_rot);
              }
            
              it2=second.begin();
              rind.at(1)++;
              if(rind.at(1)<it2->second.size()) {
                if(it2->second.at(rind.at(1)).pos-it2->second.at(alignInds.at(j).at(1)).epos>=MAX_BP_DIST) {
                  ok2=0;
                }
              } else {
                ok2=0;
              }
            }
            it1=first.begin();
            rind.at(1)=alignInds.at(j).at(1);
            rind.at(0)++;
            if(rind.at(0)<it1->second.size()) {
              if(it1->second.at(rind.at(0)).pos-it1->second.at(alignInds.at(j).at(0)).epos>=MAX_BP_DIST) {
                ok1=0;
              }
            } else {
              ok1=0;
            }
          }
        }

        //Finding new best score
        revalign.clear();
        revScore.clear();

        float bestscore=-DBL_MAX;
        vector<int> ind;
        ind.resize(2);
        int place[10];

        for(int i=0;i<sizes.at(0);i++) {
          place[0]=i;
          for(int j=0;j<sizes.at(1);j++) {
            place[1]=j;
            if(matrices->outScore(place)>bestscore) {
              bestscore=matrices->outScore(place);
              ind.at(0)=i;
              ind.at(1)=j;
            }
          }
        }

        if(bestscore>=0) {
          //Saving new best alignment
          bool ok=1;
          alignInds.clear();
          while(ok) {
            for(it1=first.begin();it1!=first.end();it1++) {
              revalign[it1->first].push_back(it1->second.at(ind.at(0)));
              alignInds.push_back(ind);
            }
            for(it2=second.begin();it2!=second.end();it2++) {
              revalign[it2->first].push_back(it2->second.at(ind.at(1)));
              alignInds.push_back(ind);
            }

            place[0]=ind.at(0);
            place[1]=ind.at(1);
            revScore.push_back(matrices->outScore(place));

            //If the alignment continues, move backwards
            ok=0;
            uint BT[10];
            BT[0]=0;
            BT[1]=0;
            matrices->outBT(place,BT);
            if(BT[0]) {
              ind.at(0)-=(int)BT[0];
              ok=1;
            }
            if(BT[1]) {
              ind.at(1)-=(int)BT[1];
              ok=1;
            }
          }

          align.clear();
          aScore.clear();
          for(int i=(revScore.size()-1);i>=0;i--) {
            aScore.push_back(revScore.at(i));
          }
          map<seqCode, vector<id_triple> >::iterator it;
          for(it=revalign.begin();it!=revalign.end();it++) {
            for(int i=(it->second.size()-1);i>=0;i--) {
              align[it->first].push_back(it->second.at(i));
            }
          }
          int aind=alignInd->giveInd();
          (*aligns)[bestscore][aind]=align;
          (*alignScores)[bestscore][aind]=aScore;
        }
      }
    }

    delete matrices;
  }

  Inputs *getPart(m_Inputs *data, seqCode seq1, seqCode seq2) {
    vector< vector<id_triple> > seq;
    map<string,seqCode> SEQ_to_id;
    vector<string> seqNames;
    map<string,motifCode> TF_to_id;
    vector<string> factorNames;

    SEQ_to_id[data->sequence(seq1)]=seq1;
    SEQ_to_id[data->sequence(seq2)]=seq2;

    seqNames.push_back(data->sequence(seq1));
    seqNames.push_back(data->sequence(seq2));

    map<string,motifCode>::iterator tfit=data->transfacNames();
    for(int i=0;i<(data->factors()/2);i++) {
      factorNames.push_back(tfit->first);
      TF_to_id[tfit->first]=tfit->second;
      tfit++;
    }

    seq.resize(2);
    for(posCode j=0;j<data->sequenceLens(seq1);j++) {
      seq.at(0).push_back(data->getSite(j,seq1));
    }
    for(posCode j=0;j<data->sequenceLens(seq2);j++) {
      seq.at(1).push_back(data->getSite(j,seq2));
    }

    return new Inputs(TF_to_id, factorNames, SEQ_to_id, seqNames, seq);
  }

  void fillTree(int method, m_Inputs *data, double lambda, double mu, double nu, double xi, double nuc_per_rot, int k, runningInd *alignInd, string pwfiles) {
    int leaves=0;

    if(right) {
      if(!right->isLeaf()) {
        //Filling the right subtree recursively
        right->fillTree(method,data,lambda,mu,nu,xi,nuc_per_rot,k,alignInd,pwfiles);
      } else {
        leaves++;
      }
    } else {
      cout<<"Error: Trying to compute in a node with no right child\n";
    }
    if(left) {
      if(!left->isLeaf()) {
        //Filling the left subtree recursively
        left->fillTree(method,data,lambda,mu,nu,xi,nuc_per_rot,k,alignInd,pwfiles);
      } else {
        leaves++;
      }
    } else {
      cout<<"Error: Trying to compute in a node with no left child\n";
    }

    //Doing the computation in this node
    if(leaves==2) {
      //Both children are leaves, aligment works the same for all scoring functions
      //Meaning, normal pairwise alignin, best alignment + k-1 suboptimal
      cout<<"Case 1\n";

      seqCode rightSeq=right->getSeq();
      seqCode leftSeq=left->getSeq();

      map<float, map<int, map<seqCode, vector<id_triple> > > > alignment;
      map<float, map<int, vector<float> > > scores;
      vector<seqCode> seqs;
      seqs.push_back(rightSeq);
      seqs.push_back(leftSeq);

      if(pwfiles.length()) {
        string filename=pwfiles+data->sequence(rightSeq)+"_"+data->sequence(leftSeq)+".align.gff";
        cout<<"Opening file "<<filename<<"\n";
        string line;
        ifstream curpwfile;
        char *s;
        s=(char*)filename.c_str();
        curpwfile.open(s);
        if(curpwfile.is_open()) {
          int firsta=1;
          float bestScore;
          map<seqCode, vector<id_triple> > align;
          vector<float>aScore;
          while(!curpwfile.eof()) {
            getline(curpwfile, line);
            
            size_t ind=line.find("CisModule");
            if(ind!=string::npos) {
              if(firsta) {
                firsta=0;
              } else {
                if(align.size()) {
                  int aind=alignInd->giveInd();
                  alignment[bestScore][aind]=align;
                  scores[bestScore][aind]=aScore;
                }
              }

              ind=line.find("|");
              for(int i=0;i<4;i++) {
                ind=line.find("|",ind+1);
              }
              size_t ind2=line.find("|",ind+1);
              istringstream mystream(line.substr(ind+1,ind2-ind-1));
              mystream >> bestScore;
              align.clear();
              aScore.clear();

            } else {
              seqCode curseq=-1;
              ind=line.find(data->sequence(rightSeq));
              if(ind!=string::npos) {
                curseq=rightSeq;
              }
              ind=line.find(data->sequence(leftSeq));
              if(ind!=string::npos) {
                curseq=leftSeq;
              }

              if(curseq>=0 && line.length()>25) {
                ind=line.find("|"); //1.
                ind=line.find("|",ind+1); //2.
                size_t ind2=line.find("|",ind+1); //3.
                string factor=line.substr(ind+1,ind2-ind-1);
                motifCode ID=data->motifIdFromName(factor);
                ind=line.find("|",ind2+1); //4.
                istringstream mystream(line.substr(ind2+1,ind-ind2-1));
                posind pos;
                mystream >> pos;
                ind2=line.find("|",ind+1); //5.
                istringstream mystream2(line.substr(ind+1,ind2-ind-1));
                posind epos;
                mystream2 >> epos;
                ind=line.find("|",ind2+1); //6.
                istringstream mystream3(line.substr(ind2+1,ind-ind2-1));
                float score;
                mystream3 >> score;
                char strand=line.at(ind+1);
                ind=line.find("|",ind+1); //7.
                ind=line.find("|",ind+1); //8.
                ind=line.find("|",ind+1); //9.
                ind2=line.find("|",ind+1); //10.
                ind=line.find("MW ");
                ind=ind+3;
                istringstream mystream4(line.substr(ind,ind2-ind));
                double weight;
                mystream4 >> weight;

                id_triple motif;
                motif.ID=ID;
                motif.weight=weight;
                motif.pos=pos;
                motif.epos=epos;
                motif.strand=strand;
              
                align[curseq].push_back(motif);
                if(curseq==rightSeq) {
                  aScore.push_back(score);
                }
              }
            }
          }
          int aind=alignInd->giveInd();
          alignment[bestScore][aind]=align;
          scores[bestScore][aind]=aScore;
          curpwfile.close();
        } else {
          cout<<"ERROR: The file "<<filename<<" couldn't be opened\n";
        }
      } else {
        //Using my own, simple multialign code with a full matrix for case 1
        m_Inputs *curdata=data->getPart(rightSeq,leftSeq);
        int gaplimit=0;
        SimpleMultiAlign *malignObject=new SimpleMultiAlign(curdata,gaplimit,lambda,mu,nu,xi,nuc_per_rot);
        malignObject->giveBestAligns(&alignment, &scores, k, seqs, alignInd);

        delete malignObject;

/*
      //Using Kimmo's simple multialign code with sparse matrix for case 1
      //This is faster and takes less memomy, but the results are different, often shorter and with smaller scores
      //The difference is probably caused by the suboptimals results, which Kimmo's multialigncode doesn't do properly
      //So either Kimmo's pairwise aligning code should be used, which looks to be tricky, or a new sparse matrix version should be implemented

      malign_AlignmentObject *alignObject=new malign_AlignmentObject();

      alignObject->CP=new struct __CPSTUF;
      alignObject->names = PyTuple_New(0);
      alignObject->bestAlignments=PyList_New(0);
      alignObject->alpha=0;
      alignObject->beta=0;
      alignObject->Rsquared=-1.0;
      alignObject->RMSE=0.0;

      alignObject->CP->indata=getPart(data,rightSeq,leftSeq);
      alignObject->CP->dynmat=new Matrix(alignObject->CP->indata);
      alignObject->lambda=lambda;
      alignObject->xi=xi; 
      alignObject->mu=mu;
      alignObject->nu=nu;
      alignObject->nuc_per_rotation=nuc_per_rot;
      alignObject->askedresults=k;
      alignObject->memSaveUsed=0;
      alignObject->secs_to_align=0;

      malignObject(alignObject);

      map<float, map<int, map<seqCode, vector<id_triple> > > > alignment;
      map<float, map<int, vector<float> > > scores;
      vector<seqCode> seqs;
      seqs.push_back(rightSeq);
      seqs.push_back(leftSeq);

      for(int i=0; i<k; i++) {
        if(alignObject->memSaveUsed==1 and alignObject->askedresults<=PyList_Size(alignObject->bestAlignments)) {
          break;
        }
        map<seqCode, vector<id_triple> > align;
        vector<float> aScore;
        findNextAlign(alignObject, &align, &aScore, seqs);
        
        if(!aScore.size()) {
          break;
        }

        int aind=alignInd->giveInd();
        float bestscore=aScore.at(aScore.size()-1);
        map<seqCode, vector<id_triple> >::iterator it;
        for(it=align.begin();it!=align.end();it++) {
          for(int j=it->second.size()-1;j>=0;j--) {
            alignment[bestscore][aind][it->first].push_back(it->second.at(j));
          }
        }
        for(int j=aScore.size()-1;j>=0;j--) {
          scores[bestscore][aind].push_back(aScore.at(j));
        }
      }
*/

      }
      map<float, map<int, map<seqCode, vector<id_triple> > > >::iterator it;
      map<float, map<int, vector<float> > >::iterator its=scores.begin();
      for(it=alignment.begin();it!=alignment.end();it++) {
        map<int, map<seqCode, vector<id_triple> > >::iterator it2;
        map<int, vector<float> >::iterator its2=its->second.begin();
        for(it2=it->second.begin();it2!=it->second.end();it2++) {
          lAligns[it->first][it2->first]=it2->second;
          alignScores[its->first][its2->first]=its2->second;
          its2++;
        }
        its++;
      }

      if(method==3) {
        //Computing the representative sequences
        map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit;
        for(rit=lAligns.rbegin();rit!=lAligns.rend();rit++) {
          map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2;
          for(rit2=rit->second.rbegin();rit2!=rit->second.rend();rit2++) {
            map<seqCode, vector<id_triple> > alignment=rit2->second;
            vector<id_triple> repSeq;
            map<seqCode, vector<id_triple> >::iterator it=alignment.begin();
            id_triple prevsite;
            vector<id_triple> prevcolumn;
            bool isprev=0;

            for(int i=0;i<it->second.size();i++) {
              vector<id_triple> column;
              for(it=alignment.begin();it!=alignment.end();it++) {
                column.push_back(it->second.at(i));
              }
              it=alignment.begin();
              id_triple repsite;

              repsite.ID=column[0].ID;
              repsite.weight=column[0].weight;
              repsite.strand=column[0].strand;

              if(isprev) {
                int dist=0;
                for(int i=0;i<column.size();i++) {
                  dist+=column[i].pos-prevcolumn[i].epos;
                }
                dist/=column.size();
                repsite.pos=prevsite.epos+dist;
                repsite.epos=repsite.pos+column[0].epos-column[0].pos;
              } else {
                repsite.pos=0;
                repsite.epos=column[0].epos-column[0].pos;
              }

              repSeq.push_back(repsite);
              prevsite=repsite;
              prevcolumn=column;
              isprev=1;
            }

            repSeqs[rit->first][rit2->first]=repSeq;
            repSeq.clear();
          }
        }
      }
    } else if(method==1 || method==2) {
      //Sum of pairs and generalized scoring work the same except for the scoring itself
      if(leaves==1) {
        //Aligning a sequence to a set of alignments
        cout<<"Case 2\n";
        map<float, map<int, map<seqCode, vector<id_triple> > > > *set;
        seqCode seqID;
        if(right->isLeaf()) {
          set=left->getSet();

          seqID=right->getSeq();
        } else {
          set=right->getSet();

          seqID=left->getSeq();
        }

        map<seqCode, vector<id_triple> > seqVec;
        seqVec[seqID]=data->getSeq(seqID);
        map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit;
        for(rit=set->rbegin();rit!=set->rend();rit++) {
          map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2;
          for(rit2=rit->second.rbegin();rit2!=rit->second.rend();rit2++) {
            map<float, map<int, map<seqCode, vector<id_triple> > > > alignment;
            map<float, map<int, vector<float> > > scores;

            vector<seqCode> seqs;
            map<seqCode, vector<id_triple> >::iterator seqit;
            for(seqit=rit2->second.begin();seqit!=rit2->second.end();seqit++) {
              seqs.push_back(seqit->first);
            }
            seqs.push_back(seqID);

            getAligns(&alignment,&scores,rit2->second,seqVec,sqrt(k),method,lambda,mu,nu,xi,nuc_per_rot,seqs,alignInd);

            map<float, map<int, map<seqCode, vector<id_triple> > > >::iterator it;
            map<float, map<int, vector<float> > >::iterator its=scores.begin();
            for(it=alignment.begin();it!=alignment.end();it++) {
              map<int, map<seqCode, vector<id_triple> > >::iterator it2;
              map<int, vector<float> >::iterator its2=its->second.begin();
              for(it2=it->second.begin();it2!=it->second.end();it2++) {
                lAligns[it->first][it2->first]=it2->second;
                alignScores[its->first][its2->first]=its2->second;
                its2++;
              }
              its++;
            }
          }
        }
      } else if(leaves==0) {
        //Aligning two sets of alignments
        cout<<"Case 3\n";
        map<float, map<int, map<seqCode, vector<id_triple> > > > *rightSet=right->getSet();
        map<float, map<int, map<seqCode, vector<id_triple> > > > *leftSet=left->getSet();

        map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator leftrit;
        map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rightrit;
        for(leftrit=leftSet->rbegin();leftrit!=leftSet->rend();leftrit++) {
          for(rightrit=rightSet->rbegin();rightrit!=rightSet->rend();rightrit++) {
            map<int, map<seqCode, vector<id_triple> > >::reverse_iterator leftrit2;
            map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rightrit2;
            for(leftrit2=leftrit->second.rbegin();leftrit2!=leftrit->second.rend();leftrit2++) {
              for(rightrit2=rightrit->second.rbegin();rightrit2!=rightrit->second.rend();rightrit2++) {
                map<float, map<int, map<seqCode, vector<id_triple> > > > alignment;
                map<float, map<int, vector<float> > > scores;

                vector<seqCode> seqs;
                map<seqCode, vector<id_triple> >::iterator seqit;
                for(seqit=rightrit2->second.begin();seqit!=rightrit2->second.end();seqit++) {
                  seqs.push_back(seqit->first);
                }
                for(seqit=leftrit2->second.begin();seqit!=leftrit2->second.end();seqit++) {
                  seqs.push_back(seqit->first);
                }

                getAligns(&alignment,&scores,rightrit2->second,leftrit2->second,1,method,lambda,mu,nu,xi,nuc_per_rot,seqs,alignInd);

                map<float, map<int, map<seqCode, vector<id_triple> > > >::iterator it;
                map<float, map<int, vector<float> > >::iterator its=scores.begin();
                for(it=alignment.begin();it!=alignment.end();it++) {
                  map<int, map<seqCode, vector<id_triple> > >::iterator it2;
                  map<int, vector<float> >::iterator its2=its->second.begin();
                  for(it2=it->second.begin();it2!=it->second.end();it2++) {
                    lAligns[it->first][it2->first]=it2->second;
                    alignScores[its->first][its2->first]=its2->second;
                    its2++;
                  }
                  its++;
                }
              }
            }
          }
        }
      }

      if(lAligns.size()>k) {
        //Saving k best alignments
        map<float, map<int, map<seqCode, vector<id_triple> > > > temp;
        map<float, map<int, vector<float> > > temps;

        map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit=lAligns.rbegin();
        map<float, map<int, vector<float> > >::reverse_iterator rits=alignScores.rbegin();

        int n=0;
        for(int i=0;n<k && i<lAligns.size();i++) {
          map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2=rit->second.rbegin();
          map<int, vector<float> >::reverse_iterator rits2=rits->second.rbegin();
          for(int j=0;n<k && j<rit->second.size();j++) {
            temp[rit->first][rit2->first]=rit2->second;
            temps[rits->first][rits2->first]=rits2->second;
            rit2++;
            rits2++;
            n++;
          }
          rit++;
          rits++;
        }
        lAligns.clear();
        alignScores.clear();

        lAligns=temp;
        alignScores=temps;
      }
    } else if(method==3) {
      //Scoring using representative sequences
      if(leaves==1) {
        cout<<"Case 2\n";
        map<float, map<int, vector<id_triple> > > repSeqs;
        seqCode seqID;
        map<float, map<int, map<seqCode, vector<id_triple> > > > *set;
        if(right->isLeaf()) {
          seqID=right->getSeq();
          repSeqs=left->getRepSeqs();
          set=left->getSet();
        } else {
          seqID=left->getSeq();
          repSeqs=right->getRepSeqs();
          set=right->getSet();
        }

        map<float, map<int, vector<id_triple> > >::reverse_iterator rit;
        vector<id_triple> seqVec=data->getSeq(seqID);
        for(rit=repSeqs.rbegin();rit!=repSeqs.rend();rit++) {
          map<int, vector<id_triple> >::reverse_iterator rit2;
          for(rit2=rit->second.rbegin();rit2!=rit->second.rend();rit2++) {
            map<float, map<int, vector<vector<int> > > > alignment;
            map<float, map<int, vector<float> > > scores;

            getRepAligns(&alignment,&scores,rit2->second,seqVec,sqrt(k),method,lambda,mu,nu,xi,nuc_per_rot,alignInd);

            //Changing the alignment of representative sequences to alignment of real sequences for printing the results
            map<float, map<int, vector<vector<int> > > >::iterator it;
            map<float, map<int, vector<float> > >::iterator its=scores.begin();
            for(it=alignment.begin();it!=alignment.end();it++) {
              map<int, vector<vector<int> > >::iterator it2;
              map<int, vector<float> >::iterator its2=its->second.begin();
              for(it2=it->second.begin();it2!=it->second.end();it2++) {
                map<seqCode, vector<id_triple> > tempMap;
                map<seqCode, vector<id_triple> >::iterator setIt;
                for(setIt=(*set)[rit->first][rit2->first].begin();setIt!=(*set)[rit->first][rit2->first].end();setIt++) {
                  vector<id_triple> tempVec;
                  for(int a=0;a<it2->second.at(0).size();a++) {
                    tempVec.push_back(setIt->second.at(it2->second.at(0).at(a)));
                  }
                  tempMap[setIt->first]=tempVec;
                  tempVec.clear();
                }

                vector<id_triple> tempVec;
                for(int a=0;a<it2->second.at(1).size();a++) {
                  tempVec.push_back(seqVec.at(it2->second.at(1).at(a)));
                }
                tempMap[seqID]=tempVec;
                tempVec.clear();

                lAligns[it->first][it2->first]=tempMap;
                alignScores[its->first][its2->first]=its2->second;
                its2++;
              }
              its++;
            }
          }
        }
      } else if(leaves==0) {
        cout<<"Case 3\n";
        map<float, map<int, vector<id_triple> > > rightRepSeqs=right->getRepSeqs();
        map<float, map<int, map<seqCode, vector<id_triple> > > > *rightSet=right->getSet();

        map<float, map<int, vector<id_triple> > > leftRepSeqs=left->getRepSeqs();
        map<float, map<int, map<seqCode, vector<id_triple> > > > *leftSet=left->getSet();

        map<float, map<int, vector<id_triple> > >::reverse_iterator leftrit;
        map<float, map<int, vector<id_triple> > >::reverse_iterator rightrit;

        for(leftrit=leftRepSeqs.rbegin();leftrit!=leftRepSeqs.rend();leftrit++) {
          for(rightrit=rightRepSeqs.rbegin();rightrit!=rightRepSeqs.rend();rightrit++) {
            map<int, vector<id_triple> >::reverse_iterator leftrit2;
            map<int, vector<id_triple> >::reverse_iterator rightrit2;
            for(leftrit2=leftrit->second.rbegin();leftrit2!=leftrit->second.rend();leftrit2++) {
              for(rightrit2=rightrit->second.rbegin();rightrit2!=rightrit->second.rend();rightrit2++) {
                map<float, map<int, vector<vector<int> > > > alignment;
                map<float, map<int, vector<float> > > scores;

                getRepAligns(&alignment,&scores,rightrit2->second,leftrit2->second,1,method,lambda,mu,nu,xi,nuc_per_rot,alignInd);

                //Changing the alignment of representative sequences to alignment of real sequences for printing the results
                map<float, map<int, vector<vector<int> > > >::iterator it;
                map<float, map<int, vector<float> > >::iterator its=scores.begin();
                for(it=alignment.begin();it!=alignment.end();it++) {
                  map<int, vector<vector<int> > >::iterator it2;
                  map<int, vector<float> >::iterator its2=its->second.begin();
                  for(it2=it->second.begin();it2!=it->second.end();it2++) {
                    map<seqCode, vector<id_triple> > tempMap;
                    map<seqCode, vector<id_triple> >::iterator setIt;
                    for(setIt=(*rightSet)[rightrit->first][rightrit2->first].begin();setIt!=(*rightSet)[rightrit->first][rightrit2->first].end();setIt++) {
                      vector<id_triple> tempVec;
                      for(int a=0;a<it2->second.at(0).size();a++) {
                        tempVec.push_back(setIt->second.at(it2->second.at(0).at(a)));
                      }
                      tempMap[setIt->first]=tempVec;
                      tempVec.clear();
                    }
                    for(setIt=(*leftSet)[leftrit->first][leftrit2->first].begin();setIt!=(*leftSet)[leftrit->first][leftrit2->first].end();setIt++) {
                      vector<id_triple> tempVec;
                      for(int a=0;a<it2->second.at(1).size();a++) {
                        tempVec.push_back(setIt->second.at(it2->second.at(1).at(a)));
                      }
                      tempMap[setIt->first]=tempVec;
                      tempVec.clear();
                    }
                    lAligns[it->first][it2->first]=tempMap;
                    alignScores[its->first][its2->first]=its2->second;
                    its2++;
                  }
                  its++;
                }
              }
            }
          }
        }
      }

      if(lAligns.size()>k) {
        //Saving k best alignment
        map<float, map<int, map<seqCode, vector<id_triple> > > > temp;
        map<float, map<int, vector<float> > > temps;

        map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit=lAligns.rbegin();
        map<float, map<int, vector<float> > >::reverse_iterator rits=alignScores.rbegin();

        int n=0;
        for(int i=0;n<k && i<lAligns.size();i++) {
          map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2=rit->second.rbegin();
          map<int, vector<float> >::reverse_iterator rits2=rits->second.rbegin();
          for(int j=0;n<k && j<rit->second.size();j++) {
            temp[rit->first][rit2->first]=rit2->second;
            temps[rits->first][rits2->first]=rits2->second;
            rit2++;
            rits2++;
            n++;
          }
          rit++;
          rits++;
        }

        lAligns.clear();
        alignScores.clear();

        lAligns=temp;
        alignScores=temps;
      }

      //Computing representative sequences
      map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit;
      for(rit=lAligns.rbegin();rit!=lAligns.rend();rit++) {
        map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2;
        for(rit2=rit->second.rbegin();rit2!=rit->second.rend();rit2++) {
          map<seqCode, vector<id_triple> > alignment=rit2->second;
          vector<id_triple> repSeq;
          map<seqCode, vector<id_triple> >::iterator it=alignment.begin();
          id_triple prevsite;
          vector<id_triple> prevcolumn;
          bool isprev=0;

          for(int i=0;i<it->second.size();i++) {
            vector<id_triple> column;
            for(it=alignment.begin();it!=alignment.end();it++) {
              column.push_back(it->second.at(i));
            }
            it=alignment.begin();
            id_triple repsite;

            repsite.ID=column[0].ID;
            repsite.weight=column[0].weight;
            repsite.strand=column[0].strand;

            if(isprev) {
              int dist=0;
              for(int i=0;i<column.size();i++) {
                dist+=column[i].pos-prevcolumn[i].epos;
              }
              dist/=column.size();
              repsite.pos=prevsite.epos+dist;
              repsite.epos=repsite.pos+column[0].epos-column[0].pos;
            } else {
              repsite.pos=0;
              repsite.epos=column[0].epos-column[0].pos;
            }

            repSeq.push_back(repsite);
            prevsite=repsite;
            prevcolumn=column;
            isprev=1;
          }

          repSeqs[rit->first][rit2->first]=repSeq;
        }
      }
    }
  }
};

class TreeMultiAlign {
  double lambda;
  double xi;
  double mu;
  double nu;
  double nuc_per_rot;

  int dims;

  m_Inputs *data;

  vector<TreeNode *> treeVec;

  runningInd *nodeInd;
  runningInd *alignInd;

  int method;
  int k;
  string pwfiles;

public:
  TreeMultiAlign() { return; }
  TreeMultiAlign(PyObject* inpSeq, PyObject* tree, double lambda, double mu, double nu, double xi, double nuc_per_rot, int k, int method, string pwfiles) {

    //Saving parameters
    this->lambda=lambda;
    this->mu=mu;
    this->nu=nu;
    this->xi=xi;
    this->nuc_per_rot=nuc_per_rot;
    this->k=k;
    this->method=method;
    this->pwfiles=pwfiles;

    //Handling input file
    this->data=new m_Inputs(inpSeq);

    this->dims=data->sequences();

    nodeInd=new runningInd();

    makeTree(tree, nodeInd, &treeVec);

    alignInd=new runningInd();

    treeVec.at(0)->fillTree(method,data,lambda,mu,nu,xi,nuc_per_rot,k,alignInd,pwfiles);
  }

  ~TreeMultiAlign() {
    if(data) {
      delete data;
      data=NULL;
    }
  }

  map<string,seqCode>::iterator sequenceNames() {
    return data->sequenceNames();
  }

  map<string,seqCode>::iterator sequenceNames(int number) {
    return treeVec.at(number)->getSeqs(this->data);
  }

  int sequenceNum(int number) {
    return treeVec.at(number)->getSeqNum();
  }

  int giveDims() {
    return dims;
  }

  TreeNode *makeTree(PyObject* tree, runningInd *nodeInd, vector<TreeNode *> *treeVec) {

    TreeNode *root=new TreeNode(nodeInd->giveInd());
    treeVec->push_back(root);

    string treeString=string(PyString_AsString(tree));

    bool dist=(treeString.find_first_of(':')!=treeString.npos);

    size_t iter=treeString.find_last_of(')');
    while(iter>0 && (treeString.at(iter)==',' || treeString.at(iter)==')')) {
      iter--;

      if(dist) {
        iter=treeString.find_last_of(':',iter);
        iter--;
      }

      iter=addNode(root, treeString, iter, dist, nodeInd, treeVec);
    }

    root->printTree();
    cout<<"\n";
    return root;
  }

  size_t addNode(TreeNode *parent, string treeString, size_t iter, bool dist, runningInd *nodeInd, vector<TreeNode *> *treeVec) {
    size_t newiter;
    TreeNode *node;
    if(treeString.at(iter)==')') {
      //The node has children
      node=new TreeNode(nodeInd->giveInd(), parent);
      treeVec->push_back(node);

      while(iter>0 && (treeString.at(iter)==',' || treeString.at(iter)==')')) {
        iter--;

        if(dist) {
          iter=treeString.find_last_of(':',iter);
          iter--;
        }

        iter=addNode(node, treeString, iter, dist, nodeInd, treeVec);
      }
      newiter=iter-1;
    } else {
      //The node is a leaf
      //Reading the sequence name
      newiter=treeString.find_last_of(",(",iter);
      string seqName=treeString.substr(newiter+1,iter-newiter);

      //Fetching the sequence code
      seqCode seqID=data->findSeq(seqName);

      //Creating new node
      node=new TreeNode(parent,seqID);
    }
    //Adding new node as the child
    parent->addChild(node);
    //Returning the new iterator
    return newiter;
  }

  PyObject *giveBest(int number) {
    map<float, map<int, map<seqCode, vector<id_triple> > > > *aligns=treeVec.at(number)->getSet();
    map<float, map<int, vector<float> > > scores=treeVec.at(number)->getScores();
    PyObject *goodAlign=PyList_New(0);
    if (goodAlign == NULL || PyErr_Occurred()) {
      cout<<"Creating goodalign failed\n";
      return NULL;
    }
    if(aligns->size()) {
      map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit=aligns->rbegin();
      map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2=rit->second.rbegin();
      map<seqCode, vector<id_triple> > align=rit2->second;
      map<float, map<int, vector<float> > >::reverse_iterator rits=scores.rbegin();
      map<int, vector<float> >::reverse_iterator rits2=rits->second.rbegin();
      vector<float> aScore=rits2->second;

      PyObject *seqs;
      PyObject *coords;
      PyObject *seqPos,*siteScore,*annotations;

      map<seqCode, vector<id_triple> >::iterator it=align.begin();
      for(int j=0;j<it->second.size();j++) {
        char strand=it->second.at(j).strand;
        siteCode motifID=it->second.at(j).ID;

        seqs=PyTuple_New(align.size());
        coords=PyTuple_New(align.size());
        seqPos=PyTuple_New(align.size());
        siteScore=PyTuple_New(align.size());
        annotations=PyTuple_New(align.size());

        int seqC=0;

        for(int i=0;i<align.size();i++) {
          // Sequence
          PyTuple_SetItem(seqs,seqC,PyInt_FromLong(i));

          // Coordinates on the sequence
          PyObject *coord;
          coord=Py_BuildValue("(ii)",(int)it->second.at(j).pos,(int)it->second.at(j).epos);
          PyTuple_SetItem(coords,seqC,coord);

          // Position on the site sequence
          PyTuple_SetItem(seqPos,seqC,PyInt_FromLong(j));

          PyTuple_SetItem(siteScore,seqC,PyFloat_FromDouble(it->second.at(j).weight));

          // Annotation
          PyTuple_SetItem(annotations,seqC,PyString_FromString(it->second.at(j).annot.c_str()));

          seqC++;   
          it++;
        }
        it=align.begin();
        // Resizing might change the location of the tuples
        _PyTuple_Resize(&seqs,seqC);
        _PyTuple_Resize(&coords,seqC);
        _PyTuple_Resize(&seqPos,seqC);
        _PyTuple_Resize(&siteScore,seqC);
        _PyTuple_Resize(&annotations,seqC);

        float score=aScore.at(j);

        PyObject *alnRow=PyAln_New_Multi(data->factor(motifID),
          seqs,seqPos,coords,strand,score,siteScore,annotations);

        PyList_Append(goodAlign,alnRow);
        Py_DECREF(alnRow);
        alnRow=NULL;

      }
    }
    return goodAlign;
  }

  //Gives the x best alignment, if x<=k
  PyObject *giveSubBest(int number, int x) {
    map<float, map<int, map<seqCode, vector<id_triple> > > > *aligns=treeVec.at(number)->getSet();
    map<float, map<int, vector<float> > > scores=treeVec.at(number)->getScores();
    PyObject *goodAlign=PyList_New(0);
    if (goodAlign == NULL || PyErr_Occurred()) {
      cout<<"Creating goodalign failed\n";
      return NULL;
    }
    if(x>k) {
      cout<<"Can't give "<<x<<"th best alignment, as only "<<k<<" were calculated.\n";
      return NULL;
    } else {
      int n=0;
      map<float, map<int, map<seqCode, vector<id_triple> > > >::reverse_iterator rit=aligns->rbegin();
      map<int, map<seqCode, vector<id_triple> > >::reverse_iterator rit2=rit->second.rbegin();
      while(n<x && rit!=aligns->rend()) {
        rit2=rit->second.rbegin();
        while(n<x && rit2!=rit->second.rend()) {
          n++;
          if(n<x) {
            rit2++;
          }
        }
        if(n<x) {
          rit++;
        }
      }

      if(n==x) {
        map<seqCode, vector<id_triple> > align=rit2->second;
        vector<float> aScore=scores[rit->first][rit2->first];

        PyObject *seqs;
        PyObject *coords;
        PyObject *seqPos,*siteScore,*annotations;

        map<seqCode, vector<id_triple> >::iterator it=align.begin();
        for(int j=0;j<it->second.size();j++) {
          char strand=it->second.at(j).strand;
          siteCode motifID=it->second.at(j).ID;

          seqs=PyTuple_New(align.size());
          coords=PyTuple_New(align.size());
          seqPos=PyTuple_New(align.size());
          siteScore=PyTuple_New(align.size());
          annotations=PyTuple_New(align.size());

          int seqC=0;

          for(int i=0;i<align.size();i++) {
            // Sequence
            PyTuple_SetItem(seqs,seqC,PyInt_FromLong(i));

            // Coordinates on the sequence
            PyObject *coord;
            coord=Py_BuildValue("(ii)",(int)it->second.at(j).pos,(int)it->second.at(j).epos);
            PyTuple_SetItem(coords,seqC,coord);

            // Position on the site sequence
            PyTuple_SetItem(seqPos,seqC,PyInt_FromLong(j));

            PyTuple_SetItem(siteScore,seqC,PyFloat_FromDouble(it->second.at(j).weight));

            // Annotation
            PyTuple_SetItem(annotations,seqC,PyString_FromString(it->second.at(j).annot.c_str()));

            seqC++;     
            it++;      
          }
          it=align.begin();

          // Resizing might change the location of the tuples
          _PyTuple_Resize(&seqs,seqC);
          _PyTuple_Resize(&coords,seqC);
          _PyTuple_Resize(&seqPos,seqC);
          _PyTuple_Resize(&siteScore,seqC);
          _PyTuple_Resize(&annotations,seqC);

          float score=aScore.at(j);

          PyObject *alnRow=PyAln_New_Multi(data->factor(motifID),
            seqs,seqPos,coords,strand,score,siteScore,annotations);

          PyList_Append(goodAlign,alnRow);
          Py_DECREF(alnRow);
          alnRow=NULL;
        }
        return goodAlign;
      } else {
        cout<<"Can't give "<<x<<"th best alignment, as only "<<aligns->size()<<" were found.\n";
        return NULL;
      }
    }
  }
};


typedef struct {
  PyObject_HEAD

  TreeMultiAlign *alignObject;

  PyObject *bestAlignments;

  PyObject *names;

  double secs_to_align;
  double lambda;
  double xi; 
  double mu;
  double nu;
  double nuc_per_rotation;
  int k;
  int method;
  int dims;
  int numofalign;

  char *pwfiles;

  double Rsquared;
} tmalign_AlignmentObject;


static PyObject *tmalign_aligndata(PyObject *self, PyObject *args) {
  cout<<"This function (tmalign_aligndata) hasn't been written as it doesn't seem to get used";
  return NULL;
}


static void tmalignment_dealloc(tmalign_AlignmentObject *self) {
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

static void tmalignment_getNodeAlignment(tmalign_AlignmentObject *self, PyObject *args) {
  int number;
  int numofalign;
  PyArg_ParseTuple(args, "ii", &number, &numofalign);
  if(PyErr_Occurred()) {
    cout<<"Extracting parameters failed\n";
  }

  PyObject *bestAlignments=PyList_New(0);
  if(bestAlignments==NULL || PyErr_Occurred()) {
    cout<<"Creating bestAlignments failed\n";
  }
  PyList_Append(bestAlignments,self->alignObject->giveBest(number));

  if(PyErr_Occurred()) {
    cout<<"Extracting the best alignment failed\n";
  }

  for(int i=2;i<=numofalign;i++) {
    if(self->alignObject->giveSubBest(number,i)) {
      PyList_Append(bestAlignments,self->alignObject->giveSubBest(number,i));
      if(PyErr_Occurred()) {
        cout<<"Extracting a suboptimal alignment failed\n";
      }
    } else {
      break;
    }
  }

  Py_DECREF(self->bestAlignments);
  if(PyErr_Occurred()) {
    cout<<"Decreffing old bestAlignments failed\n";
  }
  self->bestAlignments=bestAlignments;
  if(self->bestAlignments==NULL || PyErr_Occurred()) {
    cout<<"Putting new bestAlignments in the object failed\n";
  }

  Py_DECREF(self->names);
  if(PyErr_Occurred()) {
    cout<<"Decreffing old names failed\n";
  }
  map<string,seqCode>::iterator iter=self->alignObject->sequenceNames(number);
  int size=self->alignObject->sequenceNum(number);
  self->names=PyTuple_New(size);
  if(PyErr_Occurred()) {
    cout<<"Extracting sequence names failed\n";
  }
  for(int i=0;i<size;i++,iter++) {
    PyTuple_SET_ITEM(self->names,i,
       PyString_FromString(iter->first.c_str()));
    if(PyErr_Occurred()) {
      cout<<"Handlign sequence names failed\n";
    }
  }
}

extern "C" int tmalignment_init(tmalign_AlignmentObject *self, PyObject *args, PyObject *kwds) {

  double lambda, xi, mu, nu, nuc_per_rotation;
  string firstSeqName,secondSeqName,sequence;
  int numofalign,k,method;
  PyObject *pwfiles;
  PyObject *data;
  PyObject *tree;
  static char *kwlist[] = {"data","tree","k","numofalign","method","lambda","xi","mu","nu","nuc_per_rotation","pwfiles",NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOiiidddddO", kwlist, 
				  &data, &tree, &k, &numofalign, &method,
			&lambda, &xi, &mu, &nu, &nuc_per_rotation, &pwfiles)){
    return -1;
  }

  self->lambda=lambda;
  self->k=k;
  self->method=method;
  self->numofalign=numofalign;
  self->xi=xi;
  self->mu=mu;
  self->nu=nu;
  self->nuc_per_rotation=nuc_per_rotation;
  self->secs_to_align=0.0;
  self->Rsquared=0;
  self->pwfiles=PyString_AsString(pwfiles);
  string pwstring=string(self->pwfiles);

  tms before,after;
  long ticks_per_sec=sysconf(_SC_CLK_TCK);
  long ticks_to_align;
  times(&before);

  self->alignObject=new TreeMultiAlign(data,tree,lambda,mu,nu,xi,nuc_per_rotation,k,method,pwstring);

  if(PyErr_Occurred()) {
    cout<<"Aligning failed\n";
    return -1;
  }

  //Fetching the best result from the alignment
  self->bestAlignments=PyList_New(0);
  if(self->bestAlignments==NULL || PyErr_Occurred()) {
    cout<<"Creating bestAlignments failed\n";
    return -1;
  }
  PyList_Append(self->bestAlignments,self->alignObject->giveBest(0));

  if(PyErr_Occurred()) {
    cout<<"Extracting best alignment failed\n";
    return -1;
  }

  for(int i=2;i<=numofalign;i++) {
    if(self->alignObject->giveSubBest(0,i)) {
      PyList_Append(self->bestAlignments,self->alignObject->giveSubBest(0,i));
      if(PyErr_Occurred()) {
        cout<<"Extracting suboptimal alignment failed\n";
        return -1;
      }
    }
  }

  times(&after);
  ticks_to_align=((after.tms_utime-before.tms_utime)+
    (after.tms_stime-before.tms_stime));
  self->secs_to_align+=((double)ticks_to_align)/ticks_per_sec;

  int dims=self->alignObject->giveDims();
  self->dims=dims;

  map<string,seqCode>::iterator iter=self->alignObject->sequenceNames();
  self->names=PyTuple_New(dims);
  if(PyErr_Occurred()) {
    cout<<"Extracting sequence names failed\n";
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

static PyObject *tmalignment_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {

  tmalign_AlignmentObject *self;

  self=(tmalign_AlignmentObject *)type->tp_alloc(type,0);

  return (PyObject *)self;
}

static PyMethodDef tmalignment_methods[] = {
  {"getNodeAlignments", (PyCFunction)tmalignment_getNodeAlignment, METH_VARARGS,
   "Fetch the list of alignments from a node in the tree."},
  {NULL}
};

static PyMemberDef tmalignment_members[] = {
  {"bestAlignments",T_OBJECT_EX, offsetof(tmalign_AlignmentObject, bestAlignments), 0, "The best alignment"},
  {"names",T_OBJECT_EX, offsetof(tmalign_AlignmentObject, names), 0, "The sequence names"},
  {"secs_to_align",T_DOUBLE, offsetof(tmalign_AlignmentObject, secs_to_align), 0, "CPU time for Alignment"},
  {"Lambda",T_DOUBLE, offsetof(tmalign_AlignmentObject,lambda), 0, "Parameter for bonus"},
  {"Xi",T_DOUBLE, offsetof(tmalign_AlignmentObject,xi), 0, "Parameter for distance penalty"},
  {"Nu",T_DOUBLE, offsetof(tmalign_AlignmentObject,nu), 0, "Parameter for distance difference penalty"},
  {"Mu",T_DOUBLE, offsetof(tmalign_AlignmentObject,mu), 0, "Parameter for rotation penalty"},
  {"nuc_per_rotation",T_DOUBLE, offsetof(tmalign_AlignmentObject,nuc_per_rotation), 0, "Parameter for nucleotides per 360 deg rotation of DNA"},
  {"k",T_INT, offsetof(tmalign_AlignmentObject,k), 0, "Specifies how many alignments are saved in the sublevels."},
  {"method",T_INT, offsetof(tmalign_AlignmentObject,method), 0, "Specifies which method of scoring is used. 1=Sum of pairs, 2=Score generalized from pairwise, 3=Score based on representative sequences."},
  {"Rsquared",T_DOUBLE, offsetof(tmalign_AlignmentObject,Rsquared), 0, "This is here just to make this module compatible with the original Output.py"},
  {NULL}
};


static PyTypeObject tmalign_AlignmentType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "treeMultiAlign.TreeMultiAlignment",             /*tp_name*/
    sizeof(tmalign_AlignmentObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)tmalignment_dealloc,                         /*tp_dealloc*/
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
    "Tree multiple alignment object",           /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    tmalignment_methods,             /* tp_methods */
    tmalignment_members,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)tmalignment_init,      /* tp_init */
    0,                         /* tp_alloc */
    tmalignment_new,                 /* tp_new */
};


static PyMethodDef tmalignMethods[] = {
  {"aligndata",  tmalign_aligndata, METH_VARARGS,
   "aligns computed sequences\nArguments: data,tree,numofalign,k,method,lambda,xi,mu,nu,nuc_per_rotation,pwfiles"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

extern "C" void inittreeMultiAlign(void) {

  PyObject* m=NULL;

  tmalign_AlignmentType.tp_new=PyType_GenericNew;
  if(PyType_Ready(&tmalign_AlignmentType)<0)
    return;

  m=Py_InitModule("eellib.treeMultiAlign", tmalignMethods);

  if(m==NULL)
    return;

  if(import_alnCols()<0)
    return;

  Py_INCREF(&tmalign_AlignmentType);
  PyModule_AddObject(m, "TreeMultiAlignment", (PyObject *)&tmalign_AlignmentType);
}
