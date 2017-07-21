#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <chrono>

#include "Eigen/Dense"
#include "Eigen/Cholesky"

#include "GaussianProcess.h"


int main(void) {
    std::vector<double> X{
        0.000000,
        0.111111,
        0.222222,
        0.333333,
        0.444444,
        0.555556,
        0.666667
    };
    std::vector<double> t{
        0.349486,
        0.830839,
        1.007332,
        0.971507,
        0.133066,
        0.166823,
        -0.848307
    };
    GaussianProcess<double> GP(X, t);

    std::vector<double> x_plot, y_pred_plot, y_true_plot, s_plot;
    for (double p = 0; p < 1000; p++) {
        double x_pred = p / 1000.0;
        x_plot.push_back(x_pred);
        auto pred = GP.predict(x_pred);
        y_pred_plot.push_back(pred.mu);
        s_plot.push_back(pred.sigma);
        double y_true = sin(2 * M_PI * x_pred);
        y_true_plot.push_back(y_true);
    }

    FILE *gp = popen("gnuplot -persist", "w");
    fprintf(gp, "set multiplot\n");
    fprintf(gp, "set nokey\n");
    fprintf(gp, "set xrange [0:1]\n");
    fprintf(gp, "set yrange [-1.4:1.4]\n");
    fprintf(gp, "plot '-' with filledcurves lc rgb \"pink\"\n");
    for(size_t i=0; i<x_plot.size(); ++i) {
        fprintf(gp, "%f\t%f\t%f\n", x_plot[i], y_pred_plot[i]-s_plot[i]*2,
                                               y_pred_plot[i]+s_plot[i]*2);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "plot '-' with lines lt 1 lc rgb \"red\"\n");
    for(size_t i=0; i<x_plot.size(); ++i) {
        fprintf(gp, "%f\t%f\n", x_plot[i], y_pred_plot[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "plot '-' with lines lt 1 lc rgb \"dark-green\"\n");
    for(size_t i=0; i<x_plot.size(); ++i) {
        fprintf(gp, "%f\t%f\n", x_plot[i], y_true_plot[i]);
    }
    fprintf(gp, "e\n");

    fprintf(gp, "plot '-' with points pt 7 lc rgb \"blue\"\n");
    for (size_t i=0; i<X.size(); ++i) {
        fprintf(gp, "%f\t%f\n", X[i], t[i]);
    }
    fprintf(gp, "e\n");
    fprintf(gp, "set nomultiplot\n");
    fflush(gp);
    getchar();
    fprintf(gp, "exit\n");
    pclose(gp);
}