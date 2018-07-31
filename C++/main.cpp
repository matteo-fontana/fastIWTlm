#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "GetPot"
#include <ctime>
#include <iomanip>

#define INFO true	// Flag for printing information about the code flow
#define SHOW true	// Flag for printing information about the partial results
#define TIME true	// Flag for printing information about the elapsed time

using std::cout;
using std::endl;
using std::max;
using std::setw;

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
typedef Eigen::Array<double, Eigen::Dynamic, 1, Eigen::ColMajor> VectorType;
typedef std::string AlterType;

// Print a message on standard output
void print_mex (const std::string &);

// Generate random matrix with elements jittered around two mean values
MatrixType random_matrix (int, int, double, double);

// Generate constant vector with elements jittered around two mean values
VectorType constant_vector (int, double, double);

// Check inputs
bool check_data (const MatrixType &, const MatrixType &);
bool check_mu (const VectorType &, const MatrixType &);
bool check_alt (const AlterType &);
bool check_B (int);
bool check_maxrow (int, int);
bool check_THREADS (int);

// Compute T^2 statistic given the mean differences, the domain cardinality, and the test type
VectorType compute_T2 (const VectorType &, int, const AlterType &);

int main(int argc, char* argv[]) {
	
	clock_t tempo1 = clock();
	
	// Read parameters from command line
	GetPot commandLine(argc,argv);
	int B = commandLine("B", 1000);
	int maxrow = commandLine("maxrow", 0);
	int THREADS = commandLine("THREADS", 1);
	bool paired = commandLine("paired", false);
	bool recycle = commandLine("recycle", false);
	AlterType alt = commandLine("alt", "two.sided");
	
	// Read dimensions from command line
	int n1 = commandLine("n1", 50);
	int n2 = commandLine("n2", 50);
	int p = commandLine("p", 10);
	
	// Generate data1, data2, and mu automatically
	MatrixType data1(random_matrix(n1,p,0,2));
	MatrixType data2(random_matrix(n2,p,0,0));
	VectorType mu(constant_vector(p,0,0));
	
	/*
	// Read dimensions from text file
	int n1, n2, p;
	std::ifstream Param("Param.txt", std::ifstream::in);
	Param >> n1 >> n2 >> p;
	Param.close();
	
	// Read data1 from text file
	MatrixType data1(MatrixType::Zero(n1,p));
	std::ifstream Data1("Data1.txt", std::ifstream::in);
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < p; j++)
			Data1 >> data1(i,j);
	Data1.close();
	
	// Read data2 from text file
	MatrixType data2(MatrixType::Zero(n2,p));
	std::ifstream Data2("Data2.txt", std::ifstream::in);
	for (int i = 0; i < n2; i++)
		for (int j = 0; j < p; j++)
			Data2 >> data2(i,j);
	Data2.close();
	
	// Read mu from text file
	VectorType mu(VectorType::Zero(p));
	std::ifstream Mean0("Mean0.txt", std::ifstream::in);
	for (int j = 0; j < p;  j++)
		Mean0 >> mu(j);
	Mean0.close();
	*/
	
	// Adjust data1 for tilde test
	for (int i = 0; i < n1; i++) data1.row(i) -= mu;
	
	// Check inputs
	if (INFO) {
		print_mex("Inputs check");
		cout << "data1 and data2: " << check_data(data1,data2) << endl;
		cout << "mu: " << check_mu(mu,data1) << endl;
		cout << "alternative: " << check_alt(alt) << endl;
		cout << "B: " << check_B(B) << endl;
		cout << "maxrow: " << check_maxrow(maxrow,p) << endl;
		cout << "THREADS: " << check_THREADS(THREADS) << endl;
		if (paired) cout << "n1==n2: " << int(n1==n2) << endl;
	}
	
	clock_t tempo3 = clock();
	
	if (INFO) print_mex("T0");
	
	// Initialize vector of column-by-column mean differences
	VectorType delta0(VectorType::Zero(p));
	
	// Compute vector of column-by-column mean differences
	for (int j = 0; j < p; j++)
		delta0(j) = data1.col(j).mean() - data2.col(j).mean();
	
	// Initialize and compute vector of T0 statistics
	VectorType T0 (compute_T2(delta0, p, alt));
	
	if (SHOW) cout << round(T0.transpose()*1000)/1000 << endl;
	
	clock_t tempo4 = clock();
	
	if (INFO) print_mex("Pointwise p-values");
	
	// Initialize counter vector
	VectorType count(VectorType::Zero(p));
	
	// Initialize matrix of T^2 statistics of the permutations
	MatrixType T_perm(MatrixType::Zero(B,p));
	
	// Perform B permutations for computing pointwise p-values
	#pragma omp parallel for schedule(static) num_threads(THREADS)
	for (int b = 0; b < B; b++) {
		
		// Declare temporary variables
		double temp1, temp2;
		
		// Initialize vector of mean differences
		VectorType delta(VectorType::Zero(p));
		
		// Paired test
		if (paired) {
			
			// Declare vector of permutation indices
			std::vector<int> indices(n1);
			
			// Generate binary random sequence
			for (int i = 0; i < n1; i++)
				indices[i] = rand() % 2;
			
			// Compute column-by-column mean differences
			for (int j = 0; j < p; j++) {
				temp1 = 0;
				temp2 = 0;
				for (int i = 0; i < n1; i++) {
					if (indices[i]==1) {
						temp1 += data1(i,j);
						temp2 += data2(i,j);
					}
					else {
						temp1 += data2(i,j);
						temp2 += data1(i,j);
					}
				}
				delta(j) = (temp1-temp2) / n1;
			}
		}
		
		// Unpaired test
		else {
			
			// Declare vector of permuation indices
			std::vector<int> indices(n1+n2);
			
			// Declare auxiliary variables
			bool ok;
			int k;
			
			// Sample integers from 1 to n1+n2
			for (int i = 0; i < (n1+n2); i++) {
				indices[i] = rand() % (n1+n2);
				if (i > 0) {
					ok = false;
					while (ok == false) {
						k = 0;
						while (indices[i] != indices[k])
							k++;
						if (k < i)
							indices[i] = rand() % (n1+n2);
						else
							ok = true;
					}
				}
			}
			
			// Compute column-by-column mean differences
			for (int j = 0; j < p; j++) {
				temp1 = 0;
				temp2 = 0;
				for (int i = 0; i < n1; i++) {
					if (indices[i] < n1)
						temp1 += data1(indices[i],j);
					else
						temp1 += data2(indices[i]-n1,j);
				}
				for (int i = n1; i < (n1+n2); i++) {
					if (indices[i] < n1)
						temp2 += data1(indices[i],j);
					else
						temp2 += data2(indices[i]-n1,j);
				}
				delta(j) = temp1/n1 - temp2/n2;
			}
		}
		
		// Compute T^2 statistic of the current permutation
		T_perm.row(b) = compute_T2(delta, p, alt);
		
		// Possibly increase counter
		for (int j = 0; j < p; j++) 
			if (T_perm(b,j) >= T0(j))
				count(j) += 1;
	}
	
	// Initialize and compute pointwise p-values vector
	VectorType pvalue_point(count/B);
	
	if (SHOW) cout << round(pvalue_point.transpose()*1000)/1000 << endl;
	
	clock_t tempo5 = clock();
	
	if (INFO) print_mex("Interval-wise p-values");
	
	// Initialize intervalwise p-values matrix
	MatrixType pvalue_inter(MatrixType::Zero(p,p));
	
	// With recycle
	if (recycle==true) {
		
		// Cycle on subintervals length
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int i = p - 2; i >= maxrow; i--) {
			
			// Initialize subinterval length
			int len(p - i);
			
			// Declare counter variable
			unsigned cont;
			
			// Declare temporary variables
			double T0_temp, T_temp;
			
			// Cycle on domain points
			for (int j = 0; j < p; j++) {
				
				// Sum T0 statistics involved
				if (j + len > p)
					T0_temp = T0.tail(p-j).sum() + T0.head(len-(p-j)).sum();
				else
					T0_temp = T0.segment(j,len).sum();
				
				// Initialize counter variable
				cont = 0;
				
				// Cycle on all B permutations
				for (int b = 0; b < B; b++) {
					
					// Sum permutation T^2 statistics involved
					if (j + len > p)
						T_temp = T_perm.row(b).tail(p-j).sum() + T_perm.row(b).head(len-(p-j)).sum();
					else
						T_temp = T_perm.row(b).segment(j,len).sum();
					
					// Possibly increase counter
					if (T_temp >= T0_temp)
						cont++;
				}
				
				// Compute intervalwise p-value
				pvalue_inter(i,j) = cont / double(B);
			}
		}
	}
	
	// Without recycle
	else {
		
		// Cycle on subintervals length
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int i = p - 2; i >= maxrow; i--) {
			
			// Initialize subinterval length
			int len(p - i);
			
			// Declare counter variable
			unsigned cont;
			
			// Declare temporary variables
			double T0_temp, T_temp;
			
			// Cycle on domain points
			for (int j = 0; j <= i; j++) {
				
				// Sum T0 statistics involved
				T0_temp = T0.segment(j,len).sum();
				
				// Initialize conter variable
				cont = 0;
				
				// Cycle on all B permutations
				for (int b = 0; b < B; b++) {
					
					// Sum permutation T^2 statistics involved
					T_temp = T_perm.block(b,j,1,len).sum();
					
					// Possibly increase counter
					if (T_temp >= T0_temp)
						cont++;
				}
				
				// Compute intervalwise p-value
				pvalue_inter(i,j) = cont / double(B);
			}
		}
	}
	
	// Subintervals of length 1 are just points
	pvalue_inter.row(p-1) = pvalue_point;
	
	if (SHOW) cout << round(pvalue_inter*1000)/1000 << endl;
	
	clock_t tempo6 = clock();
	
	if (INFO) print_mex("Corrections");
	
	// Initialize corrected pointwise p-values vector
	VectorType pvalue_corr(VectorType::Zero(p));
	
	// Declare temporary variable
	double temp_last;
	
	// No corrections needed
	if (maxrow==p-1)
		pvalue_corr = pvalue_inter.row(p-1);
	
	// Corrections needed
	else {
		
		// Cycle on pvalue_inter rows
		for (int i = maxrow; i < p; i++) {
			
			// With recycle
			if (recycle==true) {
				
				// Cycle on domain points
				for (int j = p - 1; j >= 0; j--) {
					
					// Store useful value in temporary variable
					if (j==p-1) temp_last = pvalue_corr(p-1);
					
					// First row
					if (i==maxrow)
						pvalue_corr(j) = pvalue_inter(i,j);
					
					// First column
					else if (j==0)
						pvalue_corr(j) = max(pvalue_inter(i,j), max(temp_last, pvalue_corr(j)));
					
					// Center
					else
						pvalue_corr(j) = max(pvalue_inter(i,j), max(pvalue_corr(j-1), pvalue_corr(j)));
				}
			}
			
			// Without recycle
			else {
				
				// Cycle on domain points
				for (int j = i; j >= 0; j--) {
					
					// First column
					if (j==0)
						pvalue_corr(j) = max(pvalue_inter(i,j), pvalue_corr(j));
					
					// Diagonal
					else if (i==j)
						pvalue_corr(j) = max(pvalue_inter(i,j), pvalue_corr(j-1));
					
					// Center
					else
						pvalue_corr(j) = max(pvalue_inter(i,j), max(pvalue_corr(j-1), pvalue_corr(j)));
				}
			}
		}
	}
	
	clock_t tempo7 = clock();
	
	if (SHOW) cout << round(pvalue_corr.transpose()*1000)/1000 << endl;
	
	// Compute elapsed times
	if (TIME) {
		double FACT(CLOCKS_PER_SEC/1000);
		int SW(6);
		print_mex("Elapsed time");
		cout << "Read:  " << setw(SW) << round(double(tempo3-tempo1)/FACT) << " ms" << endl;
		cout << "T0:    " << setw(SW) << round(double(tempo4-tempo3)/FACT) << " ms" << endl;
		cout << "Point: " << setw(SW) << round(double(tempo5-tempo4)/FACT) << " ms" << endl;
		cout << "Inter: " << setw(SW) << round(double(tempo6-tempo5)/FACT) << " ms" << endl;
		cout << "Corr:  " << setw(SW) << round(double(tempo7-tempo6)/FACT) << " ms" << endl;
		cout << "TOTAL: " << setw(SW) << round(double(tempo7-tempo1)/CLOCKS_PER_SEC) << " s" << endl << endl;
	}
	
	return 0;
}

MatrixType random_matrix (int n, int p, double m1, double m2) {
	
	MatrixType data(MatrixType::Random(n,p));
	
	data.block(0, 0, n, floor(p/2.)) += m1;
	data.block(0, floor(p/2.), n, ceil(p/2.)) += m2;
	
	return data;
}

VectorType constant_vector (int p, double m1, double m2) {
	
	VectorType mu(VectorType::Zero(p));
	
	mu << VectorType::Constant(floor(p/2.),m1), VectorType::Constant(ceil(p/2.),m2);
	
	return mu;
}

void print_mex (const std::string & s){
	cout << endl << endl << "*** " << s << " ***" << endl << endl;
}

bool check_data (const MatrixType & d1, const MatrixType & d2) {
	return (d1.cols() == d2.cols());
}

bool check_mu (const VectorType & m, const MatrixType & d) {
	return (m.size() == d.cols());
}

bool check_alt (const AlterType & a) {
	return !(a != "two.sided" && a != "greater" && a != "less");
}

bool check_B (int B) {
	return (B > 0);
}

bool check_maxrow (int m, int p) {
	return (m >= 0 && m < p);
}

bool check_THREADS (int T) {
	return (T > 0);
}

VectorType compute_T2 (const VectorType & delta, int p, const AlterType & alt) {
	
	VectorType T(delta);
	
	if 		(!alt.compare("greater")) {
		for (int j = 0; j < p; j++)
			// Ignore negative values
			if (delta(j) < 0) T(j) = 0;
	}
	else if (!alt.compare("less")) {
		for (int j = 0; j < p; j++)
			// Ignore positive values
			if (delta(j) > 0) T(j) = 0;
	}
	
	return T * T;
}
