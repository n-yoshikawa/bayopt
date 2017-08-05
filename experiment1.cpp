#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <limits>
#include <algorithm>

#include "GaussianProcess.h"

double CDF(double x) {
    return 0.5 + 0.5 * erf(x * M_SQRT1_2);
}

int main(int argc, char** argv) {
    // parse input
    if(argc != 2) {
        std::cerr << "no csv is selected" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::vector<std::vector<double>> X;
    std::vector<double> t;
    std::vector<std::string> name;

    std::ifstream ifs(argv[1]);
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
                name.push_back(row_name);
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

    // show experiment number
    for(int ex=0; ex<155; ++ex) {
        std::cout << ex+1 << ",";
    }
    std::cout << std::endl;

    // random experiment
    /*for(int seed=0; seed<1000; ++seed) {
        std::vector<int> seq(X.size());
        for(size_t i=0; i<X.size(); i++) seq[i] = i;

        std::mt19937 mt(seed);

        std::shuffle(seq.begin(), seq.end(), mt);
        double current_best = t[seq[0]];

        for(int ex=1; ex<155; ++ex) {
            current_best = std::max(current_best, t[seq[ex]]);
            std::cout << current_best << ",";
        }
        std::cout << std::endl;
    }*/ 

    // bayesian optimization
    for(int seed=0; seed<1000; ++seed) {
        std::vector<int> isUsed(X.size(), 0);
        std::vector<std::vector<double>> X_done;
        std::vector<double> t_done;

        std::mt19937 mt(seed);
        std::uniform_int_distribution<int> sample(2, X.size()-1);  

        // initial experiment
        double current_best = -std::numeric_limits<double>::max();
        int a = sample(mt);
        X_done.push_back(X[a]);
        t_done.push_back(t[a]);
        current_best = t[a];
        isUsed[a] = 1;

        GaussianProcess<std::vector<double>> GP(X_done, t_done);

        for(int ex=0; ex<154; ++ex) {
            int bestIndex = -1;
            double best_candidate = -std::numeric_limits<double>::max();
            for(size_t i=0; i<X.size(); ++i) {
                if(isUsed[i] == 0) {
                    muAndSigma pred = GP.predict(X[i]);
                    double pi = CDF((pred.mu - current_best)/pred.sigma);
                    if(pi > best_candidate) {
                        best_candidate = pi;
                        bestIndex = i;
                    }
                }
            }
            isUsed[bestIndex] = 1;
            GP.addData(X[bestIndex], t[bestIndex]);
            if(t[bestIndex] > current_best) current_best = t[bestIndex];
            std::cout << current_best << ",";
        }
        std::cout << std::endl;
    }
}
