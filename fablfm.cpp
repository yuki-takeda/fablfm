#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "EigenAdditional/Csv"


using namespace std;
using namespace Eigen;
using namespace eigen_csv_util;

class fablfm {

    int N; // number of data
    int K; // number of latent features
    int D; // dimension of data
    double delta_; // threshold
    double mu_threshold_; // threshold of shrinkage acceleration
    double fic_lower_bound_;
    bool debug_;
    MatrixXd w_; // D * K matrix
    MatrixXd new_w_; // D * K matrix
    MatrixXd mu_; // K * N matrix
    MatrixXd new_mu_; // K * N matrix
    MatrixXd mu_pow_; // K * N matrix
    MatrixXd sum_eq_zn_znt_; // K * K matrix
    MatrixXd x_; // D * N matrix
    VectorXd b_; // D vector
    VectorXd new_b_; // D vector
    VectorXd lambda_; // D vector
    VectorXd new_lambda_; // D vector
    VectorXd rho_; // K vector
    VectorXd new_rho_; // K vector

    void initialize(MatrixXd input);
    MatrixXd update_w(MatrixXd const& x, VectorXd const& b, MatrixXd const& mu, MatrixXd const& sum_eq_zn_znt);
    VectorXd update_lambda(MatrixXd const& x, VectorXd const& b, MatrixXd const& w, MatrixXd const& mu, MatrixXd const& sum_eq_zn_znt);
    VectorXd update_b(MatrixXd const& x, MatrixXd const& w, MatrixXd const& mu);
    MatrixXd update_mu(MatrixXd const& x, MatrixXd const& w, VectorXd const& b, VectorXd const& lambda, VectorXd const& rho, MatrixXd const& current_mu);
    VectorXd update_rho(MatrixXd const& mu);
    void update_parameters();
    double c(int n, int k, MatrixXd const& x, MatrixXd const& w, VectorXd const& b, VectorXd const& lambda, MatrixXd const& mu);
    double eta(double pi_k);
    VectorXd eta(VectorXd const& pi);
    double g(double arg);

    double calc_lower_bound(MatrixXd const& x, MatrixXd const& w, VectorXd const& b, VectorXd const& lambda, VectorXd const& rho, MatrixXd const& mu, MatrixXd const& sum_eq_zn_znt);
    bool all_one_int(VectorXi v);
    bool all_one(VectorXd v);
    bool all_nearly_one(VectorXd v);
    bool all_zero(VectorXd v);
    bool nearly_all_zero(VectorXd v);
    double calc_mu_factor(int row, int col, int dim_z, VectorXd x, MatrixXd w, VectorXd b, VectorXd lambda, VectorXd rho);
    double calc_entropy(MatrixXd const& mu);
    MatrixXd calc_sum_eq_zn_znt(MatrixXd const& mu, MatrixXd const& mu_pow);
    MatrixXd calc_mu_pow(MatrixXd const& mu);

public:
    fablfm(bool debug = false);
    void solve(MatrixXd input);
    void set_parameters(MatrixXd w, VectorXd b, VectorXd lambda, VectorXd rho);
    void set_threshold(double threshold, double mu_threshold);
    double get_fic_lower_bound();
    void print_parameters(bool show_mu = false);
    void dump_parameters();
    void dump_mu(int loop_counter, int mu_loop_counter);
    MatrixXd discriminate(MatrixXd input);
};

fablfm::fablfm(bool debug) {

    debug_ = debug;
}

void fablfm::solve(MatrixXd input) {

    if (debug_) cout << "*****init*****" << endl;
    initialize(input);
    if (debug_) print_parameters();
    if (debug_) cout << "calc_sum_mu_pow" << endl;
    mu_pow_ = calc_mu_pow(mu_);
    if (debug_) cout << "calc_sum_eq_zn_znt" << endl;
    sum_eq_zn_znt_ = calc_sum_eq_zn_znt(mu_, mu_pow_);
    fic_lower_bound_ = calc_lower_bound(x_, w_, b_, lambda_, rho_, mu_, sum_eq_zn_znt_);
    double old_flb = fic_lower_bound_;
    if (debug_) {
        cout << "FIC lower bound = " << fic_lower_bound_ << endl;
        cout << "init end" << endl;
    }
    int loop_counter = 1;
    bool convergence = false;
    while (!convergence) {
        if (debug_) cout << "*****loop:" << loop_counter << "*****" << endl;

        if (debug_) cout << "update_rho" << endl;
        new_rho_ = update_rho(mu_);
        if (debug_) cout << "update_w" << endl;
        new_w_ = update_w(x_, b_, mu_, sum_eq_zn_znt_);
        if (debug_) cout << "update_lambda" << endl;
        new_lambda_ = update_lambda(x_, b_, w_, mu_, sum_eq_zn_znt_);
        if (debug_) cout << "update_b" << endl;
        new_b_ = update_b(x_, w_, mu_);

        update_parameters();

        if (debug_) cout << "update_mu" << endl;
        mu_ = update_mu(x_, w_, b_, lambda_, rho_, mu_);

        // HACK: 汚い!!!
        // shrinkage
        vector<VectorXd> mu_temp;
        vector<VectorXd> w_temp;
        for (int k = 0; k < mu_.rows(); k++) {
            if (nearly_all_zero(mu_.row(k))) {
                if (debug_) cout << "shrinkage!!" << endl;
            } else if (all_nearly_one(mu_.row(k))) {
                if (debug_) cout << "shrinkage!!" << endl;
                b_ += w_.col(k);
            } else {
                mu_temp.push_back(mu_.row(k));
                w_temp.push_back(w_.col(k));
            }
        }
        mu_.resize(mu_temp.size(), N);
        w_.resize(D, w_temp.size());
        for (int k = 0; k < mu_temp.size(); k++) {
            for (int n = 0; n < mu_.cols(); n++) {
                mu_(k, n) = mu_temp[k](n);
            }
            for (int d = 0; d < D; ++d) {
                w_(d, k) = w_temp[k](d);
            }
        }
        K = mu_.rows();
        rho_ = update_rho(mu_);

        //**************debug************************
        if (debug_) print_parameters();
        //**************debug************************

        old_flb = fic_lower_bound_;
        if (debug_) cout << "calc_sum_mu_pow" << endl;
        mu_pow_ = calc_mu_pow(mu_);
        if (debug_) cout << "calc_sum_eq_zn_znt" << endl;
        sum_eq_zn_znt_ = calc_sum_eq_zn_znt(mu_, mu_pow_);
        fic_lower_bound_ = calc_lower_bound(x_, w_, b_, lambda_, rho_, mu_, sum_eq_zn_znt_);
        if (debug_) {
            cout << "FIC lower bound = " << fic_lower_bound_ << endl;
            cout << "old FIC lower bound = " << old_flb << endl;
            cout << "K = " << K << endl;
        }
        if (abs(fic_lower_bound_ - old_flb) < delta_) {
            convergence = true;
        }
        loop_counter++;
    }
}

void fablfm::initialize(MatrixXd input) {

    x_ = input;

    N = x_.cols();
    D = x_.rows();
    K = w_.cols();

    mu_.resize(K, N);
    new_mu_.resize(K, N);
    for (int k = 0; k < mu_.rows(); k++) {
        for (int n = 0; n < mu_.cols(); n++) {
            mu_(k, n) = rho_(k);
        }
    }
}

// パラメータの初期値を設定する
void fablfm::set_parameters(MatrixXd w, VectorXd b, VectorXd lambda, VectorXd rho) {

    w_ = w;
    b_ = b;
    lambda_ = lambda;
    rho_ = rho;
}

// 閾値を設定する
void fablfm::set_threshold(double threshold, double mu_threshold) {

    delta_ = threshold;
    mu_threshold_ = mu_threshold;
}

MatrixXd fablfm::update_w(MatrixXd const& x, VectorXd const& b, MatrixXd const& mu, MatrixXd const& sum_eq_zn_znt) {

    MatrixXd sum_x_bar_mu_t = MatrixXd::Zero(D, K);
    for (int n = 0; n < mu.cols(); n++) {
        sum_x_bar_mu_t += (x.col(n) - b) * mu.col(n).transpose();
    }

    return sum_x_bar_mu_t * sum_eq_zn_znt.inverse();
}

VectorXd fablfm::update_lambda(MatrixXd const& x, VectorXd const& b, MatrixXd const& w, MatrixXd const& mu, MatrixXd const& sum_eq_zn_znt) {

    VectorXd new_lambda = VectorXd::Zero(x.rows());
    for (int d = 0; d < x.rows(); d++) {
        double sum_temp = 0.0;
        double x_bar_n_d = 0.0;
        double mu_n_w_d = 0.0;
        for (int n = 0; n < x.cols(); n++) {
            x_bar_n_d = x(d, n) - b(d);
            mu_n_w_d = (double) (mu.col(n).transpose() * w.row(d).transpose());
            sum_temp -= 2 * x_bar_n_d * mu_n_w_d;
            sum_temp += x_bar_n_d * x_bar_n_d;
        }
        sum_temp += w.row(d) * sum_eq_zn_znt * w.row(d).transpose();
        new_lambda(d) = N / sum_temp;
    }
    return new_lambda;
}

VectorXd fablfm::update_b(MatrixXd const& x, MatrixXd const& w, MatrixXd const& mu) {

    VectorXd new_b = VectorXd::Zero(D);
    for (int n = 0; n < x.cols(); n++) {
        new_b += x.col(n) - w * mu.col(n);
    }
    new_b /= N;
    return new_b;
}

// パラメータを新しいものにおきかえる
void fablfm::update_parameters() {

    w_ = new_w_;
    lambda_ = new_lambda_;
    b_ = new_b_;
    rho_ = new_rho_;
}

// ベクトルの要素が全て1であるかどうかを判定する
bool fablfm::all_one_int(VectorXi v) {

    for (int i = 0; i < v.size(); i++) {
        if (v(i) == 0) {
            return false;
        }
    }
    return true;
}

// ベクトルの要素が全てほぼ1であるかどうかを判定する
bool fablfm::all_nearly_one(VectorXd v) {

    double sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += v(i);
    }
    // return (sum / v.size() == 1.0);
    return (sum / v.size() >= 1.0 - mu_threshold_);
}

// ベクトルの要素がほぼ全て0であるかどうかを判定する
bool fablfm::nearly_all_zero(VectorXd v) {

    double sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += v(i);
    }
    return (sum / v.size() <= mu_threshold_);
}

// ベルヌーイ分布のパラメータを更新する
VectorXd fablfm::update_rho(MatrixXd const& mu) {

    VectorXd new_rho;
    new_rho.resize(K);
    for (int k = 0; k < new_rho.size(); k++) {
        double rho_k = 0;
        for (int n = 0; n < mu.cols(); n++) {
            rho_k += mu(k, n);
        }
        new_rho(k) = rho_k / N;
    }
    return new_rho;
}

MatrixXd fablfm::update_mu(MatrixXd const& x, MatrixXd const& w, VectorXd const& b, VectorXd const& lambda, VectorXd const& rho, MatrixXd const& current_mu) {

    MatrixXd mu = current_mu;
    MatrixXd new_mu = current_mu;
    VectorXd rho_tilde = rho;

    for (int k = 0; k < rho.size(); k++) {

        VectorXd c_k;
        c_k.resize(x.cols());
        for (int n = 0; n < mu.cols(); n++) {
            double c_nk = c(n, k, x, w, b, lambda, new_mu);
            c_k(n) = c_nk;
        }

        bool mu_convergence = false;
        while (!mu_convergence) {
            // update new_mu
            for (int n = 0; n < mu.cols(); n++) {
                new_mu(k, n) = g(c_k(n) + eta(rho(k)) - D / (2 * N * rho_tilde(k)));
            }
            // convergence check
            double mu_diff = 0.0;
            for (int n = 0; n < mu.cols(); n++) {
                mu_diff += abs(mu(k, n) - new_mu(k, n));
            }
            if ((mu_diff / N) < 0.0001) {
                mu_convergence = true;
            }
            // update rho_tilde
            double sum_rho_k = 0;
            for (int n = 0; n < new_mu.cols(); n++) {
                sum_rho_k += new_mu(k, n);
            }
            rho_tilde(k) = sum_rho_k / N;
            // update mu
            for (int n = 0; n < mu.cols(); n++) {
                mu(k, n) = new_mu(k, n);
            }
        }
    }
    return mu;
}

double fablfm::c(int n, int k, MatrixXd const& x, MatrixXd const& w, VectorXd const& b, VectorXd const& lambda, MatrixXd const& mu) {

    MatrixXd sum_temp;
    sum_temp = x.col(n) - b;
    for (int l = 0; l < w.cols(); l++) {
        if (l != k) {
            sum_temp -= mu(l, n) * w.col(l);
        }
    }
    sum_temp -= 0.5 * w.col(k);
    return (w.col(k).transpose() * lambda.asDiagonal() * sum_temp)(0,0);
}

double fablfm::eta(double pi_k) {

    return log(pi_k / (1.0 - pi_k));
}

VectorXd fablfm::eta(VectorXd const& pi) {

    VectorXd eta_pi(pi.size());
    for (int k = 0; k < pi.size(); k++) {
        eta_pi(k) = eta(pi(k));
    }
    return eta_pi;
}

double fablfm::g(double arg) {

    return 1.0 / (1.0 + exp(- arg));
}

double fablfm::calc_lower_bound(MatrixXd const& x, MatrixXd const& w, VectorXd const& b, VectorXd const& lambda, VectorXd const& rho, MatrixXd const& mu, MatrixXd const& sum_eq_zn_znt) {

    int K = mu.rows();
    int D = x.rows();
    int N = x.cols();
    MatrixXd sum_mun_xbarnt = MatrixXd::Zero(K, D);
    MatrixXd sum_xbarn_xbarnt = MatrixXd::Zero(D, D);
    VectorXd sum_mun = VectorXd::Zero(K);
    double sum_log_1_exp_eta_rhok = 0.0;
    double sum_log_lambdad = 0.0;
    double sum_penalty = 0.0;
    for (int n = 0; n < N; n++) {
        sum_mun_xbarnt += mu.col(n) * (x.col(n) - b).transpose();
        sum_xbarn_xbarnt += (x.col(n) - b) * (x.col(n) - b).transpose();
        sum_mun += mu.col(n);
    }
    for (int k = 0; k < K; k++) {
        sum_log_1_exp_eta_rhok += log(1 + exp(eta(rho(k))));
        double sum_munk = 0.0;
        for (int n = 0; n < N; n++) {
            sum_munk += mu(k, n);
        }
        sum_penalty += log(N * rho(k)) + (sum_munk / N - rho(k)) / rho(k);
    }
    for (int d = 0; d < D; d++) {
        sum_log_lambdad += log(lambda(d));
    }

    double eq_log_p_x_mid_z_theta = -0.5 * (
            lambda.asDiagonal() * w * sum_eq_zn_znt * w.transpose()
            - 2 * lambda.asDiagonal() * w * sum_mun_xbarnt
            + lambda.asDiagonal() * sum_xbarn_xbarnt
        ).trace() +  0.5 * N * sum_log_lambdad - 0.5 * D * N * log(2 * M_PI); 
    double eq_log_p_z_mid_rho = eta(rho).transpose() * sum_mun - N * sum_log_1_exp_eta_rhok;

    double penalty = -0.5 * D * sum_penalty - (2 * D + K) * 0.5 * log(N);
    double entropy = calc_entropy(mu);

    return eq_log_p_x_mid_z_theta + eq_log_p_z_mid_rho + penalty + entropy;
}

double fablfm::calc_entropy(MatrixXd const& mu) {

    double entropy = 0.0;
    for (int k = 0; k < mu.rows(); k++) {
        for (int n = 0; n < mu.cols(); n++) {
            double munk  = mu(k, n);
            // log(0)によるnanを回避
            if (munk != 0.0 && munk != 1.0) {
                entropy += (- munk * log(munk) - (1 - munk) * log(1 - munk));
            }
        }
    }
    return entropy;
}

MatrixXd fablfm::calc_sum_eq_zn_znt(MatrixXd const& mu, MatrixXd const& mu_pow) {

    int K = mu.rows();
    MatrixXd sum_eq_zn_znt = MatrixXd::Zero(K, K);
    MatrixXd diag;
    for (int n = 0; n < N; n++) {
        MatrixXd diag =  (mu.col(n) - mu_pow.col(n)).asDiagonal();
        sum_eq_zn_znt += mu.col(n) * mu.col(n).transpose() + diag;
    }
    return sum_eq_zn_znt;
}

MatrixXd fablfm::calc_mu_pow(MatrixXd const& mu) {

    MatrixXd mu_pow = mu;
    for (int k = 0; k < mu_pow.rows(); k++) {
        for (int n = 0; n < mu_pow.cols(); n++) {
            mu_pow(k, n) = mu_pow(k, n) * mu_pow(k, n);
        }
    }
    return mu_pow;
}

// パラメータの値をCSVファイルにダンプする
void fablfm::dump_parameters() {

    stringstream ss;
    ss << "dump/w.csv";
    string w_dump_path = ss.str();

    ss.str("");
    ss.clear();
    ss << "dump/b.csv";
    string b_dump_path = ss.str();

    ss.str("");
    ss.clear();
    ss << "dump/lambda.csv";
    string lambda_dump_path = ss.str();

    ss.str("");
    ss.clear();
    ss << "dump/rho.csv";
    string rho_dump_path = ss.str();

    ss.str("");
    ss.clear();
    ss << "dump/mu.csv";
    string mu_dump_path = ss.str();

    matrix_to_csv(w_, w_dump_path);
    vector_to_csv(b_, b_dump_path);
    vector_to_csv(lambda_, lambda_dump_path);
    vector_to_csv(rho_, rho_dump_path);
    matrix_to_csv(mu_, mu_dump_path);
}

// パラメータの値を出力する
void fablfm::print_parameters(bool show_mu) {

    cout << "rho = " << endl;
    cout << rho_ << endl;
    cout << "w = " << endl;
    cout << w_ << endl;
    cout << "b = " << endl;
    cout << b_ << endl;
    cout << "lambda = " << endl;
    cout << lambda_ << endl;
    if (show_mu) {
        cout << "mu = " << endl;
        cout << mu_ << endl;
    }
}

double fablfm::get_fic_lower_bound() {

    return fic_lower_bound_;
}

MatrixXd fablfm::discriminate(MatrixXd input) {

    // N = input.cols();
    initialize(input);
    return update_mu(input, w_, b_, lambda_, rho_, mu_);
}
