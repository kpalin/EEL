//
// $Log$
// 


typedef unsigned long int bit32;

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

  bit32 shiftMask;

  int totalCounts;

  bit32 streamHistory;

  // Number of nucleotides added to the stream.
  unsigned int streamCount;
  struct __BGdataCPP *CP;
} matrix_bgObject;


class TFBSscan {
  void nextACGTsingle(char chr,int snpCode);
  void nextACGT(char chr,int fromCode=0);
public:
  PyObject *py_matrix;  // Pointer to the matrix itself
  double bound;   // Cutoff
  int mat_length;   //Matrix length
  vector<vector<double> > M;   // Parsed matrix for easy access
  deque<deque<double> > history;
  deque<deque<double> >  compl_history;
  
  deque<class SNPdat> SNPs;

  TFBSscan(PyObject *mat,double cutoff);
  void nextChar(char chr);
  void halfHistories();
  vector<double> WatsonScore();
  vector<double> CrickScore();
  int length() { return this->mat_length; }
};



class SNPdat{
public:
  char ambig;
  char allele;
  int pos;
  SNPdat(char amb,char all,int p) {ambig=amb;allele=all;pos=p;}
  int operator==(SNPdat &other);
  int diffAllele(SNPdat &other);
  PyObject *buildPySNP(int refPos);
};



class TFBShit {
private:
  TFBShit() {return; };
public:
  TFBSscan *mat;
  unsigned int pos;
  double score;
  char strand;
  vector<class SNPdat> sigGenotype; // Significant genotype

  TFBShit(TFBSscan* mat,unsigned int seqPos,char strand);
  void addHit(double Score,vector<class SNPdat> genotype);

  PyObject *buildPySNPs();
};



class TFBShelper {
  int haveBG;
  deque<matrix_bgObject> bg;
  deque<deque<double> > probBuffer;
  int maxLen;
  vector<TFBSscan*> &matricies;
  unsigned int seqCount;


  void nextACGT(char chr,unsigned int startFrom=0);
public:
  deque<class SNPdat> SNPs;
  unsigned int matrixCount() {return this->matricies.size(); }
  unsigned int allelCount() {return this->probBuffer.size(); }
  TFBShelper(matrix_bgObject *bg,vector<TFBSscan*> &matricies);
  void nextChar(char chr);
  double getBGprob(int matInd,int snpCode);
  vector<double> getBGprobs(int matInd);
  unsigned int seqPos() {return seqCount; }

  vector<TFBShit*> getMatches();
  vector<SNPdat> getSNPs(int snpCode);
};

