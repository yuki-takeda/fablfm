#include <iostream>
#include <string>
#include <vector>

namespace eigen_csv_util {

    Eigen::VectorXd csv_to_vector(string path) {
        csv_util cu;
        std::vector<double> stlvec = cu.csv2dvec(path);
        Eigen::VectorXd eigen_vector(stlvec.size());
        for (int i = 0; i < stlvec.size(); i++) {
            eigen_vector(i) = stlvec[i];
        }
        return eigen_vector;
    }

    Eigen::MatrixXd csv_to_matrix(string path) {
        csv_util cu;
        std::vector< std::vector<double> > stlmat = cu.csv2dmat(path);
        Eigen::MatrixXd eigen_matrix(stlmat.size(), stlmat[0].size());
        for (int i = 0; i < stlmat.size(); i++) {
            for (int j = 0; j < stlmat[0].size(); j++) {
                eigen_matrix(i, j) = stlmat[i][j];
            }
        }
        return eigen_matrix;
    }

    void vector_to_csv(Eigen::VectorXd eigen_vector, string path) {
        csv_util cu;
        std::vector<double> stlvec;
        for (int i = 0; i < eigen_vector.size(); i++) {
            stlvec.push_back(eigen_vector(i));
        }
        cu.dvec2csv(stlvec, path);

    }

    void matrix_to_csv(Eigen::MatrixXd eigen_matrix, string path) {
        csv_util cu;
        std::vector< std::vector<double> > stlmat(eigen_matrix.rows());
        for (int i = 0; i < eigen_matrix.rows(); i++) {
            for (int j = 0; j < eigen_matrix.cols(); j++) {
                stlmat[i].push_back(eigen_matrix(i, j));
            }
        }
        cu.dmat2csv(stlmat, path);
    }
}

