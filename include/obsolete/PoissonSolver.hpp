#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

//#include "math_utils.h"
//#include "exon_junction_extractor.h"
#include <vector>
#include <iostream>

class isoform_opt {

const double epsilon = 1e-10;
const int max_iteration = 10000;
bool do_EM = false;

/*double sum_matrix(const std::vector<std::vector<double> > &matrix) {
	double result = 0;
	for (int i = 0; i < (int)matrix.size(); i++)
		for (int j = 0; j < (int)matrix[i].size(); j++)
			result += matrix[i][j];
	return result;
}*/

/*std::vector<std::vector<double> > inv_matrix(const std::vector<std::vector<double> > &matrix) {
	int m = (int)matrix.size(), n = (int)matrix[0].size();
	MATRIX input, *output;
	input.init((short)m, (short)n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			input.put(i, j, matrix[i][j]);
	output = input.inverse();
	std::vector<std::vector<double> > result;
	result.resize(m);
	for (int i = 0; i < m; i++) {
		result[i].resize(m);
		for (int j = 0; j < m; j++)
			result[i][j] = output->get(i,j);
	}
	delete output;
	return result;
}*/

public:
	bool verbose;
	int M; //number of categories
	int K; //number of isoforms
// for parameters
//	std::vector<double> N; //category read counts
//	std::vector<std::vector<double> > A; //category sampling rate matrix; A[i][j] is the sampling rate for isoform i and category j

	isoform_opt() {verbose = false;}

	double log_likelihood(const std::vector<std::vector<double> > &A, const std::vector<double> &N, const std::vector<double> &X) {
		double result = 0;
		double temp1 = 0;	
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j]*A[j][i];
			}
			temp1 += temp2;
		}
		result -= temp1;
		temp1 = 0;
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1 += log(temp2 + epsilon) * N[i];
		}
		result += temp1;
		return result;
	}

	double grad_log_likelihood(const std::vector<std::vector<double> > &A, const std::vector<double> &N, const std::vector<double> &X, int k) {
		double result = 0;
		double temp1 = 0;	
		for (int i = 0; i < M; i++) {
			temp1 += A[k][i];
		}
		result -= temp1;
		temp1 = 0;
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1 += N[i] * A[k][i] / (temp2 + epsilon);
		}
		result += temp1;
		return result;
	}

	std::vector<double> grad_log_likelihood(const std::vector<std::vector<double> > &A, const std::vector<double> &N, const std::vector<double> &X) {
		std::vector<double> grad;
		grad.resize(K);
		for (int i = 0; i < K; i++) {
			grad[i] = grad_log_likelihood(A, N, X, i);
		}
		return grad;
	}

	std::vector<std::vector<double> > hessian_log_likelihood(const std::vector<std::vector<double> > &A, const std::vector<double> &N, const std::vector<double> &X) {
		std::vector<std::vector<double> > hessian;
		hessian.resize(K);
		for (int i = 0; i < K; i++) hessian[i].resize(K);
		std::vector<double> temp1;
		temp1.resize(M);
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1[i] = N[i] / (temp2 + epsilon) / (temp2 + epsilon);
		}
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < K; l++) {
				hessian[k][l] = 0;
				for (int i = 0; i < M; i++) {
					hessian[k][l] -= temp1[i] * A[k][i] * A[l][i];
				}
				if (k == l) hessian[k][l] -= 1e-6;
			}
		}
		return hessian;
	}

	std::vector<std::vector<double> > obs_fisher(const std::vector<std::vector<double> > &A, const std::vector<double> &N, const std::vector<double> &X) {
		std::vector<std::vector<double> > temp = hessian_log_likelihood(A, N, X);
		for (int i = 0; i < K; i++)
			for (int j = 0; j < K; j++)
				temp[i][j] = -temp[i][j];
		return temp;
	}

/*	std::vector<std::vector<double> > inv_obs_fisher(const std::vector<double> &X) {
		return inv_matrix(obs_fisher(X));
	}*/
	
	std::vector<std::vector<double> > fisher(const std::vector<std::vector<double> > &A, const std::vector<double> &N, const std::vector<double> &X) {
		std::vector<std::vector<double> > fisher;
		fisher.resize(K);
		for (int i = 0; i < K; i++) fisher[i].resize(K);
		std::vector<double> temp1;
		temp1.resize(M);
		for (int i = 0; i < M; i++) {
			double temp2 = 0;
			for (int j = 0; j < K; j++) {
				temp2 += X[j] * A[j][i];
			}
			temp1[i] = 1 / (temp2 + epsilon);
		}
		for (int k = 0; k < K; k++) {
			for (int l = 0; l < K; l++) {
				fisher[k][l] = 0;
				for (int i = 0; i < M; i++) {
					fisher[k][l] += temp1[i] * A[k][i] * A[l][i];
				}
				if (k == l) fisher[k][l] += 1e-6;
			}
		}
		return fisher;
	}

/*	std::vector<std::vector<double> > inv_fisher(const std::vector<double> &X) {
		return inv_matrix(fisher(X));
	}*/

	double fmax_coord(const std::vector<std::vector<double> > &A, const std::vector<double> &N, std::vector<double> &X, int k) {
		double f_value = log_likelihood(A, N, X);
		double grad = grad_log_likelihood(A, N, X, k);
		if (grad > 0) grad = 1; else if (grad < 0) grad = -1; else return f_value;
		double x_old = X[k];
		double step = 1e6;
		int iteration = 0;
		while (true) {
			X[k] = x_old + step * grad;
			if (X[k] < 0) X[k] = 0;
			double new_f_value = log_likelihood(A, N, X);
			if (new_f_value > f_value) return new_f_value;
			if (step > epsilon) step /= 2; else break;
			if (iteration < max_iteration) iteration++; else break;
		}
		//if (verbose && step <= epsilon) cout << "warning: step <= epsilon\n";
		if (iteration >= max_iteration) std::cout << "warning: iteration >= max_iteration\n";
		X[k] = x_old;
		return f_value;
	}

	double fmax_coord(const std::vector<std::vector<double> > &A, const std::vector<double> &N, std::vector<double> &X, bool use_EM) {
		double f_value = log_likelihood(A, N, X);
		int iteration = 0;
		while (true) {
			if (verbose) {
				std::cout << "iteration: " << iteration << " value: " << f_value << " X: ";
				for (int i = 0; i < K; i++) std::cout << X[i] << ",";
				std::cout << std::endl;
			}
			if (use_EM) {
				std::vector<double> temp1(M, 0);
				for (int j = 0; j < M; j++) {
					for (int i = 0; i < K; i++) {
						temp1[j] += A[i][j] * X[i];
					}
					if (temp1[j] < epsilon) temp1[j] = epsilon;
				}
				std::vector<double> X_old = X;
				for (int i = 0; i < K; i++) {
					double temp2 = 0, temp3 = 0;
					for (int j = 0; j < M; j++) {
						temp2 += A[i][j] * N[j] / temp1[j];
						temp3 += A[i][j];
					}
					if (temp3 < epsilon) temp3 = epsilon;
					X[i] = X_old[i] * temp2 / temp3;
					if (X[i] < 1e-4) X[i] = 0;
				}
			} else { //coordinate wise descent
				for (int i = 0; i < K; i++) {
					fmax_coord(A, N, X, i);
				}
			}
			double new_f_value = log_likelihood(A, N, X);
			if (fabs(f_value - new_f_value) < epsilon) {
				if (verbose) std::cout << "number of iteration: " << iteration << std::endl;
				return new_f_value;
			}
			f_value = new_f_value;
			if (iteration < max_iteration) iteration++; else break;
		}
		if (iteration >= max_iteration) std::cout << "warning: iteration >= max_iteration\n";
		return f_value;
	}
};


inline void test_solve_likelihood() {
	isoform_opt iso_opt;
	bool do_EM = false;
	int M = 3;
	int K = 2;
	iso_opt.M = M;
	iso_opt.K = K;
	std::vector<double> N;
	N.resize(M);
	N[0] = 100;
	N[1] = 700;
	N[2] = 500;
	std::vector<std::vector<double> > A;
	A.resize(K);
	A[0].resize(M);
	A[0][0] = 1;
	A[0][1] = 1;
	A[0][2] = 0;
	A[1].resize(M);
	A[1][0] = 0;
	A[1][1] = 1;
	A[1][2] = 1;
	iso_opt.verbose = true;
	std::vector<double> X;
	X.resize(2);
	X[0] = 1;
	X[1] = 1;
	iso_opt.fmax_coord(A, N, X, do_EM);
	std::cout << "X: " << X[0] << "\t" << X[1] << std::endl;
	std::vector<std::vector<double> > H = iso_opt.hessian_log_likelihood(A, N, X);
	std::cout << "H: ";
	std::cout << "\t" << H[0][0] << "\t" << H[0][1] << std::endl;
	std::cout << "\t" << H[1][0] << "\t" << H[1][1] << std::endl;	
	std::vector<std::vector<double> > OF = iso_opt.obs_fisher(A, N, X);
	std::cout << "OF: ";
	std::cout << "\t" << OF[0][0] << "\t" << OF[0][1] << std::endl;
	std::cout << "\t" << OF[1][0] << "\t" << OF[1][1] << std::endl;	
	std::vector<std::vector<double> > F = iso_opt.fisher(A, N, X);
	std::cout << "F: ";
	std::cout << "\t" << F[0][0] << "\t" << F[0][1] << std::endl;
	std::cout << "\t" << F[1][0] << "\t" << F[1][1] << std::endl;	
/*		std::vector<std::vector<double> > IOF = inv_obs_fisher(X);
	cout << "IOF: ";
	cout << "\t" << IOF[0][0] << "\t" << IOF[0][1] << std::endl;
	cout << "\t" << IOF[1][0] << "\t" << IOF[1][1] << std::endl;	
	std::vector<std::vector<double> > IF = inv_fisher(X);
	cout << "IF: ";
	cout << "\t" << IF[0][0] << "\t" << IF[0][1] << std::endl;
	cout << "\t" << IF[1][0] << "\t" << IF[1][1] << std::endl;	*/
}

inline void solve_likelihood(const std::vector<std::vector<double> > &rates, const std::vector<double> &counts, std::vector<double> &exp, bool use_EM = false) {
	isoform_opt iso_opt;
	iso_opt.M = (int)counts.size();
	iso_opt.K = (int)rates.size();
	exp.resize(iso_opt.K);
	for (int i = 0; i < iso_opt.K; i++) {
		exp[i] = 1;
	}
	iso_opt.fmax_coord(rates, counts, exp, use_EM);
}

inline void solve_likelihood_t(const std::vector<std::vector<double> > &rates, const std::vector<double> &counts, std::vector<double> &exp, bool use_EM = false) {
	isoform_opt iso_opt;
	iso_opt.M = (int)counts.size();
	//__ASSERT(iso_opt.M > 0, "internal error: empty rates.\n");
	iso_opt.K = (int)rates[0].size();
	std::vector<std::vector<double> > A;
	A.resize(iso_opt.K);
	for (int i = 0; i < iso_opt.K; i++) A[i].resize(iso_opt.M);
	for (int i = 0; i < iso_opt.M; i++) {
		for (int j = 0; j < iso_opt.K; j++) {
			A[j][i] = rates[i][j];
		}
	}
	exp.resize(iso_opt.K);
	for (int i = 0; i < iso_opt.K; i++) {
		exp[i] = 1;
	}
	iso_opt.fmax_coord(A, counts, exp, use_EM);
}


#endif // POISSON_SOLVER_HPP