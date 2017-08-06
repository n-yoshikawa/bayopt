#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>
#include <utility>

#include "GaussianProcess.h"

double CDF(double x) {
    return 0.5 + 0.5 * erf(x * M_SQRT1_2);
}

int main(int argc, char** argv) {
    // load train data
    std::vector<std::vector<double>> X;
    std::vector<double> t;

    std::ifstream ifs("experiment2-train.csv");
    std::string line;
    std::getline(ifs, line); // ignore the first line
    while(std::getline(ifs, line)) {
        std::stringstream ss(line);
        std::string elem;
        int row = 0;
        std::vector<double> row_x;
        double row_t;
        while(std::getline(ss, elem, ',')) {
            std::stringstream ss2(elem);
            row++;
            if(row==1) {
                std::string row_name;
                ss2 >> row_name;
            } else {
                double e;
                ss2 >> e;
                if(row <= 21) row_x.push_back(e);
                else row_t = e;
            }
        }
        X.push_back(row_x);
        t.push_back(row_t);
    }
    ifs.close();

    GaussianProcess<std::vector<double>> GP(X, t);

    // load and predict
    std::vector<std::pair<double, std::string>> pred;
    ifs.open("experiment2-pred.csv");
    std::getline(ifs, line); // ignore the first line
    while(std::getline(ifs, line)) {
        std::stringstream ss(line);
        std::string elem;
        int row = 0;
        std::vector<double> row_x;
        std::string row_name;
        while(std::getline(ss, elem, ',')) {
            std::stringstream ss2(elem);
            row++;
            if(row==1) {
                ss2 >> row_name;
            } else {
                double e;
                ss2 >> e;
                if(row <= 21) row_x.push_back(e);
                else break;
            }
        }
        double t_pred = GP.predict(row_x).mu;
        pred.push_back(std::make_pair(t_pred, row_name));
    }

    // show result
    std::sort(pred.begin(), pred.end(), std::greater<std::pair<double, std::string>>());

    auto it = std::find_if(pred.begin(), pred.end(),
    [](const std::pair<double, std::string>& e){ return e.second == "GAYF";} );
    std::cout << "GAYF found at " << std::distance(pred.begin(), it) << std::endl;
}
