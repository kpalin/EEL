#include <iostream>
#include <ext/hash_map>
#include <vector>

using std::cout;
using std::cerr;
using std::vector;
using __gnu_cxx::hash_map;

typedef vector<int> intArray;
typedef vector<double> doubleArray;
typedef hash_map<int, double> myHashMap;
typedef vector< intArray > intMatrix;
typedef vector<unsigned char> charArray;


charArray randomString(const int n, const int numA) 
{
    charArray ret;
    for (int i = 0; i < n; ++i) {
	double r = rand() / (1.0 + RAND_MAX);
	ret.push_back((char) (numA * r));
    }
    return ret;
}

//this one uses table
int tresholdFromP(const intMatrix &mat, const double &p)
{
    int numA = mat.size();
    int n = mat[0].size();

    int maxT = 0;
    int minV = INT_MAX;

    for (int i = 0; i < n; ++i) {
	int max = mat[0][i];
	int min = max;
	for (int j = 1; j < numA; ++j) {
	    int v = mat[j][i];
	    if (max < v)
		max = v;
	    else if (min > v)
		min = v;
	}
	maxT += max;
	if (minV > min)
	    minV = min;
    }
    int R = maxT - n * minV;

    //cout << "Size " << R << "\n" ;

    // use indexes, map, own implementation...
    doubleArray table0(R+1, 0.0);
    doubleArray table1(R+1, 0.0);

    for (int j = 0; j < numA; ++j)
	table0[mat[j][0] - minV] += (1.0 / numA); // change this to use own background model

    for (int i = 1; i < n; ++i) {
	for (int j = 0; j < numA; ++j) {
	    int s = mat[j][i] - minV;
	    for (int r = s; r < R+1; ++r)
		table1[r] += (1.0 / numA) * table0[r - s]; // change this to use own background model
	}
	for (int r = 0; r < R+1; ++r) {
	    table0[r] = table1[r];
	    table1[r] = 0.0;
	}
    }

    double sum = 0.0;

    for (int r = R; r >= 0; --r) {
	sum += table0[r];
	//cout << "sum = " << sum << "\n" ;
	if (sum > p) {
	    //cout << "tol = " << r << "\n" ;
	    return (r + n * minV + 1);
	}
    }
    cerr << "Error: No treshold found!";
    return 0;
    //cout << "sum = " << sum << "\n";
}

// same as above, uses hashmap instead of arrays...
// somewhat slower...
int tresholdFromP2(const intMatrix &mat, const double &p)
{
    int numA = mat.size();
    int n = mat[0].size();

    int maxT = 0;
    int minV = INT_MAX;

    for (int i = 0; i < n; ++i) {
	int max = mat[0][i];
	int min = max;
	for (int j = 1; j < numA; ++j) {
	    int v = mat[j][i];
	    if (max < v)
		max = v;
	    else if (min > v)
		min = v;
	}
	maxT += max;
	if (minV > min)
	    minV = min;
    }

    myHashMap table0;
    myHashMap table1;

    for (int j = 0; j < numA; ++j)
	table0[ mat[j][0] ] += (1.0 / numA); // change this to use own background model

    for (int i = 1; i < n; ++i) {
	//cout << "Size = " << table0.size() << "\n";
	for (int j = 0; j < numA; ++j) {
	    int s = mat[j][i];
	    for (myHashMap::iterator it = table0.begin(); 
		 it != table0.end(); ++it) {
		table1[ it->first + s ] += (1.0 / numA) * (it->second); // change this to use own background model
	    }
	}
	table0 = table1;
	table1.clear();
    }

    //cout << "maxT = " << maxT << " minV = " << minV << "\n";
    double sum = 0.0;
    for (int r = maxT; r >= n * minV; --r) {
	sum += table0[r];
	//cout << "sum = " << sum << "\n" ;
	if (sum > p) {
	    //cout << "tol = " << r << "\n" ;
	    return (r + 1);
	}
    }
    cerr << "Error: No treshold found!";
    return 0;
}


//just a driver to test the p-value computation
int main()
{
    int numA = 4;
    int m = 5;
    int n = 10000000;

    charArray sequence( randomString(n, numA) );

    intMatrix mat;
    /*for (int i = 0; i < numA; ++i) {
	mat.push_back( intArray() );
	for (int j = 0; j < m; ++j) {
	    double r = rand() / (1.0 + RAND_MAX);
	    mat[i].push_back( (int) (1000 * r) );
	}
	}*/
    mat.push_back(intArray());
    mat[0].push_back(16);
    mat[0].push_back(16);
    mat[0].push_back(16);
    mat[0].push_back(0);
    mat[0].push_back(1);
    mat.push_back(intArray());
    mat[1].push_back(0);
    mat[1].push_back(0);
    mat[1].push_back(0);
    mat[1].push_back(0);
    mat[1].push_back(9);
    mat.push_back(intArray());
    mat[2].push_back(0);
    mat[2].push_back(0);
    mat[2].push_back(0);
    mat[2].push_back(16);
    mat[2].push_back(1);
    mat.push_back(intArray());
    mat[3].push_back(0);
    mat[3].push_back(0);
    mat[3].push_back(0);
    mat[3].push_back(0);
    mat[3].push_back(5);

    for (double p = 0; p < 0.01; p+=0.001) {
	int T = tresholdFromP2(mat, p);
	cout << "T2 = " << T << " T1 = " << tresholdFromP(mat, p) << "\n";
	int hits = 0;
	for (int i = 0; i <= n - m; ++i) {
	    int score = 0;
	    for (int j = 0; j < m; ++j)
		score += mat[sequence[i + j]][j];
	    if (score >= T)
		++hits;
	    
	}
	cout << "p = " << p << " hits " << (1.0 * hits / (n - m + 1)) << "\n";
	cout << "E = " << p*(n - m + 1) << " hits " << (1.0 * hits ) << "\n";
	
    }
    return 0;
}
