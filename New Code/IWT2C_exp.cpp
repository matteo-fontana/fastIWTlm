#include <RcppEigen.h>
#include <omp.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <ctime>
#include <iomanip>

#define INFO true	// Code flow
#define SHOW false	// Partial esults
#define TIME true	// Elapsed times

using std::endl;
using std::max;
using namespace Rcpp;
using namespace RcppEigen;

typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
typedef Eigen::Array<double, Eigen::Dynamic, 1, Eigen::ColMajor> VectorType;
typedef std::string AlterType;

void print_mex (const std::string &);

bool check_data (const MatrixType &, const MatrixType &);
bool check_alt (const AlterType &);
bool check_B (int);
bool check_maxrow (int, int);
bool check_THREADS (int);

VectorType compute_T2 (const VectorType &, int, const AlterType &);

//[[Rcpp::export]]
List IWT2C_exp(NumericMatrix Data1, NumericMatrix Data2 , NumericVector Mu, int B = 1000, AlterType alt = "two.sided", int maxrow = 0,
             bool paired = false, bool recycle = false, int THREADS = 1) {
	
	clock_t tempo10 = clock();
  
  
  //input Matrix Features
  int n1 = Data1.nrow();
  int n2 = Data2.nrow();
  int p = Data1.ncol();
  //Save Data inputs to MatrixType or VectorType Objects
  
  MatrixType data1 = MatrixType::Map(Data1.begin(),n1,p);
  MatrixType data2 = MatrixType::Map(Data2.begin(),n2,p); 
  VectorType mu = VectorType::Map(Mu.begin(),p);
  
  
	
	// Read dimensions from text file
	//int n1, n2, p;
	//std::ifstream Param("Param.txt", std::ifstream::in);
	//Param >> n1 >> n2 >> p;
	//Param.close();
	
	// Read data1 from text file
	//MatrixType data1(MatrixType::Zero(n1,p));
	//std::ifstream Data1("Data1.txt", std::ifstream::in);
	//for (int i = 0; i < n1; i++)
	//	for (int j = 0; j < p; j++)
	//		Data1(i,j) = data1(i,j);
	//Data1.close();
	
	// Read data2 from text file
	//MatrixType data2(MatrixType::Zero(n2,p));
	//std::ifstream Data2("Data2.txt", std::ifstream::in);
	//for (int i = 0; i < n2; i++)
	//	for (int j = 0; j < p; j++)
	//		Data2(i,j) = data2(i,j);
	//Data2.close();
	
	// Read mu from text file
	//VectorType mu(VectorType::Zero(p));
	//std::ifstream Mean0("Mean0.txt", std::ifstream::in);
	//for (int j = 0; j < p;  j++)
	//	Mean0 >> mu(j);
	//Mean0.close();
	
	// Adjust data1 for tilde test
	for (int i = 0; i < n1; i++) data1.row(i) -= mu;
	
	clock_t tempo11 = clock();
	
	// Check input parameters
	if (INFO) {
		print_mex("Inputs check");
		Rcout << "Data matrices\t\t" << check_data(data1,data2) << endl;
		Rcout << "Permutations B\t\t" << check_B(B) << endl;
		Rcout << "Alternative\t\t" << check_alt(alt) << endl;
		Rcout << "Truncation maxrow\t" << check_maxrow(maxrow,p) << endl;
		Rcout << "Number of threads\t" << check_THREADS(THREADS) << endl;
		if (paired) Rcout << "Paired test allowed\t" << int(n1==n2) << endl;
	}
	
	// T0
	if (INFO) print_mex("Start computing T0");
	
	VectorType delta0(VectorType::Zero(p));
	
	clock_t tempo20 = clock();
	
	for (int j = 0; j < p; j++)
		delta0(j) = data1.col(j).mean() - data2.col(j).mean();
	
	VectorType T0 (compute_T2(delta0, p, alt));
	
	clock_t tempo21 = clock();
	
	if (SHOW) Rcout << "T0 = [ " << T0.transpose() << "]" << endl;
	
	if (INFO) print_mex("End computing T0");
	
	// Pointwise p-values
	if (INFO) print_mex("Start computing point-wise pvalues");
	
	VectorType count(VectorType::Zero(p));
	MatrixType T_perm(MatrixType::Zero(B,p));
	
	clock_t tempo30 = clock();
	
    #pragma omp parallel for schedule(static) num_threads(THREADS)
	for (int b = 0; b < B; b++) {
		
		double temp1, temp2;
		VectorType delta(VectorType::Zero(p));
		
		if (paired) {
			
			std::vector<int> indices(n1);
			
			for (int i = 0; i < n1; i++)
				indices[i] = rand() % 2;
			
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
		else {
			
			std::vector<int> indices(n1+n2);
			bool ok;
			int k;
			
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
		
		T_perm.row(b) = compute_T2(delta, p, alt);
		for (int j = 0; j < p; j++) 
			if (T_perm(b,j) >= T0(j))
				count(j) += 1;
		
	}
	
	clock_t tempo31 = clock();
	
	VectorType pvalue_point(count/B);
	if (SHOW) Rcout << "Pointwise p-values:" << pvalue_point.transpose() << endl;
	
	if (INFO) print_mex("End computing pointwise p-values");
	
	// Intervalwise p-values
	if (INFO) print_mex("Start computing intervalwise p-values");
	
	MatrixType pvalue_inter(MatrixType::Zero(p,p) * R_NaN);
	
	clock_t tempo40 = clock();
	
	if (recycle==true) {
		
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int i = p - 2; i >= maxrow; i--) {
			
			int len(p - i);
			unsigned cont;
			double T0_temp, T_temp;
			
			for (int j = 0; j < p; j++) {
				
				if (j + len > p)
					T0_temp = T0.tail(p-j).sum() + T0.head(len-(p-j)).sum();
				else
					T0_temp = T0.segment(j,len).sum();
				
				cont = 0;
				
				for (int b = 0; b < B; b++) {
					
					if (j + len > p)
						T_temp = T_perm.row(b).tail(p-j).sum() + T_perm.row(b).head(len-(p-j)).sum();
					else
						T_temp = T_perm.row(b).segment(j,len).sum();
					
					if (T_temp >= T0_temp)
						cont++;
				}
				pvalue_inter(i,j) = cont / double(B);
			}
		}
	}
	else {
	    
		#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
		for (int i = p - 2; i >= maxrow; i--) {
			
			int len(p - i);
			unsigned cont;
			double T0_temp, T_temp;
			
			for (int j = 0; j <= i; j++) {
				
				T0_temp = T0.segment(j,len).sum();
				cont = 0;
				
				for (int b = 0; b < B; b++) {
					
					T_temp = T_perm.block(b,j,1,len).sum();
					
					if (T_temp >= T0_temp)
						cont++;
				}
	            pvalue_inter(i,j) = cont / double(B);
			}
		}
	}
	
	clock_t tempo41 = clock();
	
	pvalue_inter.row(p-1) = pvalue_point;
	
	if (SHOW) Rcout << "Intervalwise p-values:" << endl << pvalue_inter << endl;
	
	if (INFO) print_mex("End computing intervalwise p-values");
	
	// Corrections
	if (INFO) print_mex("Start computing corrections");
	
	clock_t tempo50 = clock();
	
	VectorType pvalue_corr(VectorType::Zero(p));
	double temp_last;
	
	if (maxrow==p-1)
		pvalue_corr = pvalue_inter.row(p-1);
	else {
		for (int i = maxrow; i < p; i++) {
			
			if (recycle==true) {
				for (int j = p - 1; j >= 0; j--) {
					
					if (j==p-1) temp_last = pvalue_corr(p-1);
					
					if (i==maxrow)
						pvalue_corr(j) = pvalue_inter(i,j);
					else if (j==0)
						pvalue_corr(j) = max(pvalue_inter(i,j), max(temp_last, pvalue_corr(j)));
					else
						pvalue_corr(j) = max(pvalue_inter(i,j), max(pvalue_corr(j-1), pvalue_corr(j)));
				}
			}
			else {
				for (int j = i; j >= 0; j--) {
					if (j==0)
						pvalue_corr(j) = max(pvalue_inter(i,j), pvalue_corr(j));
					else if (i==j)
						pvalue_corr(j) = max(pvalue_inter(i,j), pvalue_corr(j-1));
					else
						pvalue_corr(j) = max(pvalue_inter(i,j), max(pvalue_corr(j-1), pvalue_corr(j)));
				}
			}
		}
	}
	
	clock_t tempo51 = clock();
	
	if (SHOW) Rcout << "Corrected p-values:" << pvalue_corr.transpose() << endl;
	
	if (INFO) print_mex("End computing corrections");
	
	// Elapsed times
	if (TIME) {
		int TOT_TIME(tempo51-tempo10);
		Rcout << "TOTAL\tRead\tT0\tPoint\tInter\tCorr" << endl;
		Rcout << round(TOT_TIME/CLOCKS_PER_SEC) << "\t"
			  << round(double(tempo11-tempo10)/CLOCKS_PER_SEC) << "\t"
			  << round(double(tempo21-tempo20)/CLOCKS_PER_SEC) << "\t"
			  << round(double(tempo31-tempo30)/CLOCKS_PER_SEC) << "\t"
			  << round(double(tempo41-tempo40)/CLOCKS_PER_SEC) << "\t"
			  << round(double(tempo51-tempo50)/CLOCKS_PER_SEC) << endl;
	}
	
	List RES;
	RES["T0"]    = T0;
	RES["point"] = pvalue_point;
	RES["inter"] = pvalue_inter;
	RES["corre"] = pvalue_corr;
	return RES;
}

void print_mex (const std::string & s){
	Rcout << endl << "*** " << s << " ***" << endl << endl;
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
			if (delta(j) < 0) T(j) = 0;
	}
	else if (!alt.compare("less")) {
		for (int j = 0; j < p; j++)
			if (delta(j) > 0) T(j) = 0;
	}
	
	return T * T;
}
