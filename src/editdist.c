#include <Python.h>
#include <string.h> 

#include <limits.h>

int min2(int a,int b) 
{
  return (a<b?a:b);
}

int min3(int a,int b,int c) 
{
  return min2(a,min2(b,c));
}


static PyObject *
spam_editDistance(self, args)
    PyObject *self;
    PyObject *args;
{
  int const delcost=1,chgcost=1;
  int i,j;
  unsigned int smaxlen,sminlen,arrayLen,editDist,countEnds=1;
  int *prev,*curr,*tmp;
  char *smin,*smax;

  if (!PyArg_ParseTuple(args, "ss|d",&smin,&smax,&countEnds)){
    Py_INCREF(Py_None);
    return Py_None;
  }

  sminlen=strlen(smin);
  smaxlen=strlen(smax);

  if(smaxlen<sminlen) {
    i=smaxlen;
    smaxlen=sminlen;
    sminlen=i;

    tmp=smin;
    smin=smax;
    smax=tmp;
  }

  arrayLen=sminlen+1;

  prev=malloc(arrayLen*sizeof(int));
  curr=malloc(arrayLen*sizeof(int));

  if(!prev || !curr) {
    Py_INCREF(Py_None);
    return Py_None;
  }


  /* Set initial values. Penalty for shorter string. */
  for(i=0;i<arrayLen;i++) {
    prev[i]=i;
  }


  editDist=sminlen+smaxlen;
  curr[0]=0;
  for(i=0;i<smaxlen;i++) {
    if(countEnds) {
      curr[0]=i+1;
    }
    for(j=0;j<sminlen;j++) {
      if(smax[i]==smin[j]) {
	curr[j+1]=min3(prev[j],curr[j]+delcost,prev[j+1]+chgcost);
      } else {
	curr[j+1]=min3(prev[j]+chgcost,curr[j]+delcost,prev[j+1]+delcost);
      }
    }
    editDist=min2(editDist,curr[sminlen]);
    /* Done for one row. Now get a new curr. */
    tmp=curr;
    curr=prev;
    prev=tmp;
  }

  if(countEnds) {
    editDist=prev[sminlen];
  }
  free(prev);
  free(curr);

  return Py_BuildValue("i",editDist);
}


int IUPACmatch(char a,char b)
{
  char *IUPACclass[]={"RGA",
		     "YTC",
		     "MAC",
		     "KGT",
		     "SGC",
		     "WAT",
		     "HACT",
		     "BGTC",
		     "VGCA",
		     "DGAT",
		     "NACGT",
		      NULL};
  int i,j,match_a=0,match_b=0;

  if(a==b) return 1;

  for(i=0;IUPACclass[i];i++) {
    match_a=0;
    match_b=0;
    for(j=0;j<strlen(IUPACclass[i]);j++) {
      match_a|= (IUPACclass[i][j]==a);
      match_b|= (IUPACclass[i][j]==b);
    }
    if(match_a && match_b) {
      return 1;
    }
  }
  return 0;
}

static PyObject *
spam_alignSeq(self, args)
    PyObject *self;
    PyObject *args;
{
  int const delcost=1,chgcost=1;
  int flipFlag=0;
  int i,j,aln_p;
  unsigned int smaxlen,sminlen,arrayLen,editDist,countEnds=1;
  unsigned int useIUPAC=0;
  char *smin,*smax,*stmp,nothing,*backTrace;
  int *curr,*prev,*tmp;
  char *minaln,*maxaln;
  PyObject *ret_val;

  if (!PyArg_ParseTuple(args, "ss|dd",&smin,&smax,&countEnds,&useIUPAC)){
    Py_INCREF(Py_None);
    return Py_None;
  }

  sminlen=strlen(smin);
  smaxlen=strlen(smax);

  if(smaxlen<sminlen) {
    flipFlag=1;
    i=smaxlen;
    smaxlen=sminlen;
    sminlen=i;
    stmp=smin;
    smin=smax;
    smax=stmp;
  }

  arrayLen=sminlen+1;

  backTrace=malloc((sminlen+1)*(smaxlen+1)*sizeof(char));

#define BT_MATCH 1
#define BT_MININS 2
#define BT_MAXINS 3
#define BT_POS(mat,x,y) (char)mat[(x)*(sminlen)+(y)]

  prev=malloc(arrayLen*sizeof(int));
  curr=malloc(arrayLen*sizeof(int));

  if(!prev || !curr) {
    Py_INCREF(Py_None);
    return Py_None;
  }


  /* Set initial values. Penalty for shorter string. */
  for(i=0;i<arrayLen;i++) {
    prev[i]=i;
  }


  editDist=sminlen+smaxlen;
  curr[0]=0;
  for(i=0;i<smaxlen;i++) {
    if(countEnds) {
      curr[0]=i+1;
    }
    for(j=0;j<sminlen;j++) {
      /* Look for equality */
      if((useIUPAC&&IUPACmatch(smax[i],smin[j])) || smax[i]==smin[j]) {
	if (prev[j]<curr[j]+delcost && prev[j]<prev[j+1]+delcost) {
	  BT_POS(backTrace,i,j)=BT_MATCH;
	  curr[j+1]=prev[j];
	} else if (prev[j]>curr[j]+delcost && prev[j+1]+delcost>curr[j]+delcost) {
	  BT_POS(backTrace,i,j)=BT_MAXINS;
	  curr[j+1]=curr[j]+delcost;
	} else {
	  BT_POS(backTrace,i,j)=BT_MININS;
	  curr[j+1]=prev[j+1]+delcost;
	}
      } else {
	if (prev[j]+chgcost<curr[j]+delcost && prev[j]+chgcost<prev[j+1]+delcost) {
	  BT_POS(backTrace,i,j)=BT_MATCH;
	  curr[j+1]=prev[j]+chgcost;
	} else if (prev[j]+chgcost>curr[j]+delcost && prev[j+1]+delcost>curr[j]+delcost) {
	  BT_POS(backTrace,i,j)=BT_MAXINS;
	  curr[j+1]=curr[j]+delcost;
	} else {
	  BT_POS(backTrace,i,j)=BT_MININS;
	  curr[j+1]=prev[j+1]+delcost;
	}
      }
    }
    editDist=min2(editDist,curr[sminlen]);
    /* Done for one row. Now get a new curr. */
    tmp=curr;
    curr=prev;
    prev=tmp;
  }
		   

		   
  i=smaxlen-1;
  j=sminlen-1;
  minaln=(char*)malloc((sminlen+smaxlen+1)*sizeof(char));
  maxaln=(char*)malloc((sminlen+smaxlen+1)*sizeof(char));
  memset(minaln,0,sminlen+smaxlen+1);
  memset(maxaln,0,sminlen+smaxlen+1);
  aln_p=sminlen+smaxlen-1;

  while (i>=0 && j>=0) {
    switch(BT_POS(backTrace,i,j)) {
    case BT_MATCH:
      minaln[aln_p]=smin[j];
      maxaln[aln_p]=smax[i];
      i--;
      j--;
      break;
    case BT_MAXINS:
      minaln[aln_p]=smin[j];
      maxaln[aln_p]='-';
      j--;
      break;
    case BT_MININS:
      minaln[aln_p]='-';
      maxaln[aln_p]=smax[i];
      i--;
      break;

    }
    aln_p--;
  }

  while(i>=0) {
    minaln[aln_p]='-';
    maxaln[aln_p]=smax[i];
    i--;
    aln_p--;
  }
  while(j>=0) {
    minaln[aln_p]=smin[j];
    maxaln[aln_p]='-';
    j--;
    aln_p--;
  }


  if(countEnds) {
    editDist=prev[sminlen];
  }
  if(flipFlag) {
    stmp=minaln;
    minaln=maxaln;
    maxaln=stmp;
  }
  ret_val=Py_BuildValue("iss",editDist,minaln+aln_p+1,maxaln+aln_p+1);


  free(backTrace);
  free(prev);
  free(curr);
  free(minaln);
  free(maxaln);

  return ret_val;
}






static PyMethodDef EditDistMethods[] = {
    {"editDistance",  spam_editDistance, METH_VARARGS,
     "Computes the edit distance of two strings. If third argument is non zero, discard the ends of the longer sequence."},
    {"alignSeq",  spam_alignSeq, METH_VARARGS,
     "Computes the edit distance of two strings. If third argument is non zero, discard the ends of the longer sequence."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};



void
initeditdist(void)
{
    (void) Py_InitModule("editdist", EditDistMethods);
}

