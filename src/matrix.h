//
// $Log$
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

// Zero-order BG version

static PyObject * matrix_getAllTFBSzeroOrderBG(PyObject *self, PyObject *args);

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



doubleArray expectedDifferences(const doubleMatrix &mat, const doubleArray &bg);

void multipleMatrixAhoCorasickLookaheadFiltration(const charArray &s, const intArray &start_pos, const intArray &end_pos, const vector<doubleMatrix> &matrices, const doubleArray &bg, const doubleArray &tol, PyObject *ret_dict, vector<PyObject*> py_matrices, const charArray &strands);
void getHitsWithSNPs(const charArray &s, const intArray &snp_pos, const doubleMatrix &p, const doubleArray &bg, const double cutoff, PyObject *ret_dict, PyObject *py_matrix, const char strand);

// Markov BG version

static PyObject * matrix_getAllTFBSMarkovBG(PyObject *self, PyObject *args);
void getHitsWithMarkovBG(const charArray &sequence, const doubleMatrix &pssm,  matrix_bgObject * bg, const doubleArray &bgProps, const double cutoff, PyObject *ret_dict, PyObject *py_matrix, const char strand);
void getHitsWithMarkovBG_v2(const charArray &sequence, const doubleMatrix &pssm,  matrix_bgObject * bg, const doubleArray &bgProps, const double cutoff, PyObject *ret_dict, PyObject *py_matrix, const char strand);




