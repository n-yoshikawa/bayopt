#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>

#include "GaussianProcess.h"

int main(int argc, char** argv) {
    /* parse input */
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

    std::vector<int> isUsed(X.size(), 0);

    /* do bayopt */
    std::vector<std::vector<double>> X_done;
    std::vector<double> t_done;

    std::random_device rnd;
    std::mt19937 mt(rnd());
    std::uniform_int_distribution<int> sample(2, X.size()-1);  

    double current_best = -1e50;
    for(int i=0; i<20; ++i) {
        int a = sample(mt);
        if(isUsed[a] == 0) {
            //std::cout << "done experiment: " << name[a] << std::endl;
            X_done.push_back(X[a]);
            t_done.push_back(t[a]);
            current_best = std::max(t[a], current_best);
            isUsed[a] = 1;
        }
    }

    GaussianProcess<std::vector<double>> GP(X_done, t_done);

    for(int ex=0; ex<160; ++ex) {
        if(isUsed[0] == 1) {
            std::cout << "YFP found at " << ex << std::endl;
            std::exit(0);
        }

        double bestPI = -100;
        int bestIndex = -1;
        double best_candidate = -1e50;
        for(size_t i=0; i<X.size(); ++i) {
            if(isUsed[i] == 0) {
                muAndSigma pred = GP.predict(X[i]);
                if(pred.mu > best_candidate) {
                    best_candidate = pred.mu;
                    bestIndex = i;
                }
            }
        }
        isUsed[bestIndex] = 1;
        GP.addData(X[bestIndex], t[bestIndex]);
    }
}
