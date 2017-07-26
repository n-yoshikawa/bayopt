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


    /* do bayopt */
    std::vector<double> log;
    for(int repeat=0; repeat<100; ++repeat) {
        std::vector<int> isUsed(X.size(), 0);
        std::vector<std::vector<double>> X_done;
        std::vector<double> t_done;

        std::random_device rnd;
        std::mt19937 mt(rnd());
        std::uniform_int_distribution<int> sample(2, X.size()-1);  

        double current_best = -std::numeric_limits<double>::max();
        for(int i=0; i<20; ++i) {
            int a = sample(mt);
            if(isUsed[a] == 0) {
                X_done.push_back(X[a]);
                t_done.push_back(t[a]);
                current_best = std::max(t[a], current_best);
                isUsed[a] = 1;
            }
        }

        GaussianProcess<std::vector<double>> GP(X_done, t_done);

        for(int ex=0; ex<160; ++ex) {
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
            if(isUsed[0] == 1) {
                std::cout << repeat << ": YFP found at " << ex << std::endl;
                log.push_back(ex);
                break;
            }
            if(bestIndex == -1) {
                std::cout << "Abort." << std::endl;
                break;
            }
            GP.addData(X[bestIndex], t[bestIndex]);
            if(t[bestIndex] > current_best) current_best = t[bestIndex];
        }
    }
    int min = *std::min_element(log.begin(), log.end());
    int max = *std::max_element(log.begin(), log.end());
    double mean = std::accumulate(log.begin(), log.end(), 0) / static_cast<double>(log.size());
    double sq_sum = std::inner_product(log.begin(), log.end(), log.begin(), 0.0);
    double sd = std::sqrt(sq_sum / log.size() - mean * mean);
    std::cout << "mean: " << mean << " sd: " << sd << " min: " << min << " max: " << max << std::endl;
}