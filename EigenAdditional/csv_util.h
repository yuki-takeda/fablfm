#ifndef INCLUDED_Sample_h_
#define INCLUDED_Sample_h_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;

class csv_util {
public:
    // convert CSV to double vector
    vector<double> csv2dvec(string path);
    // convert CSV to double matrix
    vector< vector<double> > csv2dmat(string path);
    // convert CSV to integer vector
    vector<int> csv2ivec(string path);
    // convert CSV to integer matrix
    vector< vector<int> > csv2imat(string path);
    // convert CSV to string vector
    vector<string> csv2svec(string path);
    // convert CSV to string matrix
    vector< vector<string> > csv2smat(string path);
    // convert CSV line to string svector
    vector<string> line2svec(string str);
    // convert CSV cell (string data) to double
    double cell2double(string cell);
    // count CSV column
    int count_column(string path);
    // count CSV cells per line
    int count_cell(string str);
    // count CSV row(line)
    int count_row(string path);
    // write double vector to CSV
    void dvec2csv(vector<double> dvec, string path);
    // write matrix to CSV
    void dmat2csv(vector< vector<double> > dmat, string path);
};

#endif
