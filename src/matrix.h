//
// $Log$
// Revision 1.11  2008/06/06 10:56:35  jazkorho
// Finalized new TFBS searching code. Zero-order background version now uses an Aho-Corasick based filter, but the original implementation is still used with Markov background.
//
// Revision 1.10  2008/05/23 12:10:10  jazkorho
// Switched to a simpler TFBS scanning algorithm due to problems in tests on Murska cluster
//
// Revision 1.9  2008/05/19 08:14:40  jazkorho
// Rewrote TFBS search code.
//
// Revision 1.8  2008/02/29 09:10:09  kpalin
// Report all snps
//
// Revision 1.7.2.1  2008/01/21 12:51:02  kpalin
// Now report all, even weak, SNPs hitting binding sites.
//
// Revision 1.7  2006/11/13 12:37:09  kpalin
// Added a constant for rounding range. Linear dependency on the speed of
// the p-value threshold computation. Unknown, probably mostly positive,
// dependency on the accuracy of the thresholding.
//
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
#define LARGE_AFFY_DELTA 0.0
#endif

#ifndef SEQ_BUFFER_SIZE
#define SEQ_BUFFER_SIZE 15000000
#endif

#ifndef SCANNING_WINDOW_SIZE
#define SCANNING_WINDOW_SIZE 9
#endif


#include<stdint.h>

// Range for rounding the matrix values for p-value thresholing
// Increasing this gives more accurate results slowly.
#define ROUNDING_RANGE 1000.0

typedef vector<int> intArray;
typedef vector<double> doubleArray;
typedef vector<intArray> intMatrix;
typedef vector<doubleArray> doubleMatrix;
typedef vector<unsigned char> charArray;

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


// TFBS search stuff

// *********************
// Zero-order BG version


// Struct for comparing matrix rows with sort()
struct compareRows
{
    const doubleArray *goodness;
    bool operator() (int i, int j)
    {
        return ( (*goodness)[i] > (*goodness)[j] );
    }
};

// Output list element for mm algorithms
struct OutputListElementMulti
{
    double score;
    int matrix;
    bool full;
};

// Struct used in construction of multiple matrix AC automaton
struct ConstructionACStateMulti
{
    int transition[4];
    list<OutputListElementMulti> output;
};

// Struct used in final multiple matrix AC automaton
struct FinalACStateMulti
{
    FinalACStateMulti * transition[4];
    list<OutputListElementMulti> output;
};

// Misc. data structure used in AC automaton construction
struct ConstructionQueueElementMulti
{
    int prev;
    int i;
    list<OutputListElementMulti> scores;
};

// Algorithms

static PyObject * matrix_getAllTFBSzeroOrderBG(PyObject *self, PyObject *args);

doubleArray expectedDifferences(const doubleMatrix &mat, const doubleArray &bg);

void multipleMatrixAhoCorasickLookaheadFiltration(const charArray &s, const intArray &start_pos, const intArray &end_pos, const vector<doubleMatrix> &matrices, const doubleArray &bg, const doubleArray &tol, PyObject *ret_dict, vector<PyObject*> py_matrices, const charArray &strands);
void getHitsWithSNPs(const charArray &s, const intArray &snp_pos, const doubleMatrix &p, const doubleArray &bg, const double cutoff, PyObject *ret_dict, PyObject *py_matrix, const char strand);

// *****************
// Markov BG version


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

static PyObject * matrix_getAllTFBSMarkovBG(PyObject *self, PyObject *args);

