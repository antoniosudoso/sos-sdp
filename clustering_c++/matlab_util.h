#ifndef MATLAB_CALLER_MATLAB_UTIL_H
#define MATLAB_CALLER_MATLAB_UTIL_H


#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include <armadillo>

std::unique_ptr<matlab::engine::MATLABEngine> init_matlab(std::string &session_name, const char *sdp_solver_folder);
std::unique_ptr<matlab::engine::MATLABEngine> start_matlab(const char *sdp_solver_folder);
matlab::data::TypedArray<double> arma_to_matlab_vector(matlab::data::ArrayFactory &factory, arma::vec &v);
matlab::data::TypedArray<double> arma_to_matlab_matrix(matlab::data::ArrayFactory &factory, arma::mat &X);
matlab::data::SparseArray<double> arma_to_matlab_sparse(matlab::data::ArrayFactory &factory, arma::sp_mat &X);
matlab::data::CellArray arma_to_matlab_cell(matlab::data::ArrayFactory &factory, arma::sp_mat &A);
matlab::data::CellArray arma_to_matlab_cell(matlab::data::ArrayFactory &factory, std::vector<arma::sp_mat> &B_vector);
arma::mat matlab_to_arma_matrix(matlab::data::TypedArray<double> &X_matlab);
arma::vec matlab_to_arma_vector(matlab::data::TypedArray<double> &v_matlab);
std::vector<arma::sp_mat> matlab_to_arma_sp_mat_vector(matlab::data::CellArray &B_matlab);
arma::sp_mat matlab_to_arma_sparse(matlab::data::SparseArray<double> &X_matlab);


#endif //MATLAB_CALLER_MATLAB_UTIL_H
