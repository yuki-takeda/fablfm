#include <iostream>
#include <chrono>
#include "Eigen/Core"
#include "Eigen/LU"
#include "Eigen/Dense"
#include "EigenAdditional/Csv"
#include "fablfm.cpp"


using namespace std;
using namespace Eigen;
using namespace eigen_csv_util;


int main(int argc,char *argv[]){
    
    // コマンドライン引数の処理
    int n_components;
    bool debug;
    string input_path;
    for (int i = 1; i < argc; i++) {
        // options.push_back(argv[i]);
        string option = argv[i];
        // cout << option << endl;
        if (option.find("--n_components=") != string::npos) {
            // コンポーネント数
            n_components = stoi(option.substr(option.find("=") + 1));
        }else if (option.find("--debug=true") != string::npos){
            // デバッグモード
            debug = true; 
        } else {
            // 入力データのファイルパス
            input_path = option;
        }
    }

    std::chrono::system_clock::time_point  start, end; // 型は auto で可
    start = std::chrono::system_clock::now(); // 計測開始時間

    MatrixXd input = csv_to_matrix(input_path);

    int d = input.rows();
    int k = n_components;
    int hoge;
    std::srand((unsigned int) time(0) + reinterpret_cast<uintptr_t>(&hoge));
    MatrixXd w = MatrixXd::Random(d, k) * 5.0;
    VectorXd b = VectorXd::Zero(d);
    VectorXd lambda = VectorXd::Random(d) + VectorXd::Constant(d, 1.0);
    VectorXd rho = (VectorXd::Random(k) + VectorXd::Constant(k, 1.0)) * 0.5;

    fablfm fl(debug);
    fl.set_parameters(w, b, lambda, rho);
    fl.set_threshold(0.001, 0.01);
    fl.solve(input);

    end = std::chrono::system_clock::now();  // 計測終了時間
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count(); //処理に要した時間をミリ秒に変換
    cout << "elapsed time: " << elapsed << "ms" << endl;
    
    return 0;
}
