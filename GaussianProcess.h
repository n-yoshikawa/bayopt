#include <vector>
#include <iostream>

#include "Eigen/Dense"

struct muAndSigma {
    double mu, sigma;
};

template <typename T>
class GaussianProcess
{
    std::vector<T> x_data;
    Eigen::MatrixXd Cn;
    Eigen::VectorXd t_data;
    Eigen::LLT<Eigen::MatrixXd> cholCn;
    Eigen::MatrixXd L; 
    Eigen::VectorXd alpha;

    size_t N;
    const double beta = 30.0;

    double kernel(T x1, T x2, double t = 3.0) {
        double r = pow(x1 - x2, 2);
        return exp(-t * r);
    }

  public:
    GaussianProcess(std::vector<T> X, std::vector<double> t)
        : x_data{X}, N{X.size()}
    {
        assert(X.size() == t.size());
        N = X.size();
        Cn.resize(N, N);
        t_data.resize(N);
        for (int i = 0; i < N; ++i) {
            t_data(i) = t[i];
            for (int j = 0; j < N; ++j) {
                Cn(i, j) = kernel(x_data[i], x_data[j]) + (i == j ? 1.0 / beta : 0);
            }
        }
        cholCn = Eigen::LLT<Eigen::MatrixXd>(Cn);
        L = cholCn.matrixL();
        alpha = cholCn.solve(t_data);
    }

    void addData(T x, double t) {
        N++;
        x_data.push_back(x);
        t_data.conservativeResize(N);
        t_data(N-1) = t;
        Cn.conservativeResize(N, N);
        for (int i = N - 1; i < N; ++i) {
            double e = kernel(x_data[i], x) + (i == N-1 ? 1.0 / beta : 0);
            Cn(i, N - 1) = e;
            Cn(N-1, N - 1) = e;
        }
        cholCn = Eigen::LLT<Eigen::MatrixXd>(Cn);
    }

    muAndSigma predict(T x_pred)
    {
        Eigen::VectorXd k(N);
        for (int i = 0; i < N; ++i)
        {
            k(i) = kernel(x_data[i], x_pred);
        }
        auto v = L.triangularView<Eigen::Lower>().solve(k);
        double m = k.transpose() * alpha;
        double s = sqrt(kernel(x_pred, x_pred) + 1.0/beta - v.transpose() * v);
        return {m, s};
    }
};