//
// $Log$
// Revision 1.6  2006/08/10 11:36:10  kpalin
// Another try on 64bit port.
//
// Revision 1.5  2006/08/10 11:38:44  kpalin
// Port to 64bit. Changed custom bit32 type to uint32_t.
//
// Revision 1.4  2005/05/19 07:49:42  kpalin
// Merged Waterman-Eggert style suboptimal alignments and
// SNP matching.
//
// Revision 1.3.2.3  2005/05/12 11:07:12  kpalin
// Made matrix to compile.
//
// Revision 1.3.2.2  2005/05/09 07:38:02  kpalin
// Fixed some half ready work.
//
// Revision 1.3.2.1  2005/05/09 07:29:10  kpalin
// Reflecting few cleanups in _c_matrix.cc
//
// Revision 1.3  2005/03/03 09:01:59  kpalin
// Presumably working with proper output.
//
// Revision 1.2  2005/02/25 09:28:35  kpalin
// Few additions for the less inefficient scanner version.
//
// Revision 1.1  2005/02/23 13:40:25  kpalin
// Initial Clean up.
//
// 


#ifndef LARGE_AFFY_DELTA
#define LARGE_AFFY_DELTA 1.0
#endif

#include<stdint.h>

// Range for rounding the matrix values for p-value thresholing
// Increasing this gives more accurate results slowly.
#define ROUNDING_RANGE 1000.0


struct __BGdataCPP {
  vector<unsigned long int> counts;
  vector<unsigned long int> shortCounts;
};

typedef struct {
  PyObject_HEAD

  PyObject *bgSample;
  //unsigned int sampleLen;

  PyObject *strBuf,*py_read,*py_readParam;
  unsigned int bytes_read,buf_p;

  unsigned int order;
  unsigned int qgram;

  uint32_t shiftMask;

  int totalCounts;

  uint32_t streamHistory;

  // Number of nucleotides added to the stream.
  unsigned int streamCount;
  struct __BGdataCPP *CP;
} matrix_bgObject;



class TFBSscan {
  void nextACGTsingle(char chr,int snpCode);
  void doubleHistory();
  void nextACGT(char const chr,int fc=0,int tc=-1);
public:
  PyObject *py_matrix;  // Pointer to the matrix itself
  double bound;   // Cutoff
  int mat_length;   //Matrix length
  vector<vector<double> > M;   // Parsed matrix for easy access
  deque<deque<double> > history;
  deque<deque<double> >  compl_history;
  
  //deque<class SNPdat> SNPs;

  TFBSscan(PyObject *mat,double cutoff);
  
  double matItem(int i,char nucl,int crick=0);
  static int ACGTtoCode(char nucl);
  double setSNPscoreDif(class SNPdat &snp,int crick=0);
  void nextChar(char chr);
  void halfHistories();
  vector<double> WatsonScore();
  vector<double> CrickScore();
  unsigned int allelCount() {return this->history.size(); }
  int const length() const { return this->mat_length; }
};



class SNPdat{
public:
  char ambig;
  char allele;
  int pos;
  double scoreDif;

  SNPdat(char amb,char all,int p) {ambig=amb;allele=all;pos=p;scoreDif=0.0;}
  int operator==(SNPdat &other);
  int diffAllele(SNPdat &other);
  PyObject *buildPySNP(int refPos);
  char *alleles();
};



class TFBShit {
private:
  TFBShit() {return; };
public:
  TFBSscan *mat;
  unsigned int pos;
  double score;  //Max score
  double minScore;
  char strand;
  vector<class SNPdat> sigGenotype; // Significant genotype

  TFBShit(TFBSscan* mat,unsigned int seqPos,char strand);
  void addHit(double Score,vector<class SNPdat> &genotype);

  PyObject *buildPySNPs();
};



class TFBShelper {
  int haveBG;
  deque<matrix_bgObject> bg;
  deque<deque<double> > probBuffer;
  int maxLen;
  vector<TFBSscan*> &matricies;
  unsigned int seqCount;


  void doubleBackground();
  
  void nextACGT(char chr,unsigned int startFrom=0,unsigned int upTo=0);
  void removeScannerHistories();
public:
  unsigned int SNPcount() {return this->SNPs.size()/2; }
  unsigned int bgOrder() { return (this->haveBG?this->bg[0].order:0); }
  deque<class SNPdat> SNPs;
  unsigned int matrixCount() {return this->matricies.size(); }
  unsigned int allelCount() {return 1<<this->SNPcount(); }
  TFBShelper(matrix_bgObject *bg,vector<TFBSscan*> &matricies);
  void nextChar(char chr);
  double getBGprob(int matInd,int snpCode);
  vector<double> getBGprobs(int matInd);
  unsigned int seqPos() {return seqCount; }

  vector<TFBShit*> getMatches();
  vector<SNPdat> getSNPs(int snpCode,int matInd=-1,int crick=0);
  ~TFBShelper() { for(unsigned int i=0;i<this->matricies.size();i++) delete this->matricies[i];}
};

