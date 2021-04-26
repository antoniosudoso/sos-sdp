#include "sdp_solver_util.h"

arma::sp_mat build_A(int m, int n) {
    arma::sp_mat A(m, n * n);
    // Ze = e
    for (int i = 0; i < n; i++) {
        arma::sp_mat temp_A(n, n);
        temp_A.row(i) = 0.5 * arma::ones(n).t();
        temp_A.col(i) = 0.5 * arma::ones(n);
        temp_A(i, i) = 1;
        A.row(i) = temp_A.as_col().t();
    }
    // trace(Z) = k
    A.row(m - 1) = arma::speye(n, n).as_col().t();
    return A;
}

arma::vec build_b(int m, int k) {
    // Ze = e
    arma::vec b_sum_one = arma::ones(m - 1);
    // trace(Z) = k
    arma::vec b_trace(1);
    b_trace(0) = k;
    return arma::join_cols(b_sum_one, b_trace);
}


// BUILD DATA FOR THE SDP BRANCH AND BOUND

arma::sp_mat build_A_must_link(int m, arma::sp_mat &TTt) {
    arma::sp_mat A(m, (m - 1) * (m - 1));
    arma::vec e_l = arma::vec(arma::diagvec(TTt));
    // Ze_l = e
    for (int k = 0; k < m - 1; k++) {
        arma::sp_mat temp_A(m - 1, m - 1);
        temp_A.row(k) = 0.5 * e_l.t();
        temp_A.col(k) = 0.5 * e_l;
        temp_A(k, k) = e_l(k);
        A.row(k) = temp_A.as_col().t();
    }
    // <TTt, Z> = k
    A.row(m - 1) = TTt.as_col().t();
    return A;
}

arma::vec build_b_must_link(int m, int k) {
    return build_b(m, k);
}

arma::sp_mat build_A_cannot_link(arma::sp_mat &old_A, int i, int j) {
    if (i >= j) {
        std::fprintf(stderr, "build_A_cannot_link(): branching node must be sorted!");
        exit(EXIT_FAILURE);
    }
    int m = old_A.n_rows;
    int n = old_A.n_cols;
    arma::sp_mat A = arma::resize(old_A, m + 1, n);
    // Z_ij = 0
    int sqrt_n = std::sqrt(n);
    arma::sp_mat diff(sqrt_n, sqrt_n);
    diff(i, j) = 0.5;
    diff(j, i) = 0.5;
    A.row(m) = diff.as_col().t();
    return A;
}

arma::vec build_b_cannot_link(arma::vec &old_b) {
    int m = old_b.size();
    arma::vec b(old_b);
    b.resize(m + 1);
    // Z_ij = 0
    b(m) = 0;
    return b;
}