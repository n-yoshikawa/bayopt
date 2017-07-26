#include <vector>
#include <iostream>
#include <cassert>

#include "Eigen/Dense"

struct muAndSigma {
    double mu, sigma;
};

template <typename T>
class GaussianProcess
{
    std::vector<T> x_data;
    std::vector<double> t_store;
    Eigen::MatrixXd Cn;
    Eigen::VectorXd t_data;
    Eigen::LLT<Eigen::MatrixXd> cholCn;
    Eigen::MatrixXd L; 
    Eigen::VectorXd alpha;

    size_t N;
    const double beta = 30.0;

    double kernel(const T& x1, const T& x2, double t = 1e-3);

  public:
    GaussianProcess(std::vector<T> X, std::vector<double> t)
        : x_data{X}, t_store{t}, N{X.size()}
    {
        assert(X.size() == t.size());
        N = X.size();
        Cn.resize(N, N);
        t_data.resize(N);
        for (int i = 0; i < N; ++i) {
            t_data(i) = t_store[i];
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
        t_store.push_back(t);
        Cn.resize(N, N);
        t_data.resize(N);
        for (int i = 0; i < N; ++i) {
            t_data(i) = t_store[i];
            for (int j = 0; j < N; ++j) {
                Cn(i, j) = kernel(x_data[i], x_data[j]) + (i == j ? 1.0 / beta : 0);
            }
        }
        cholCn = Eigen::LLT<Eigen::MatrixXd>(Cn);
        L = cholCn.matrixL();
        alpha = cholCn.solve(t_data);
    }

    muAndSigma predict(T x_pred)
    {
        Eigen::VectorXd k(N);
        for (int i = 0; i < N; ++i)
        {
            k(i) = kernel(x_data[i], x_pred);
        }

        //Eigen::VectorXd v = cholCn.matrixL().solve(k);
        Eigen::VectorXd v = L.triangularView<Eigen::Lower>().solve(k);
        double m = k.transpose() * alpha;
        double s2 = kernel(x_pred, x_pred) + 1.0/beta - v.transpose() * v;
        assert(s2 > 0);
        double s = sqrt(s2);
        return {m, s};
    }
};

template <typename T>
double GaussianProcess<T>::kernel(const T& x1, const T& x2, double t) {
    double r = pow(x1 - x2, 2);
    return exp(-t * r);
}

template <>
double GaussianProcess<std::vector<double>>::kernel(const std::vector<double>& x1, const std::vector<double>& x2, double t) {
    assert(x1.size() == x2.size());
    double r = 0;
    for(size_t i=0; i<x1.size(); ++i) {
        r += pow(x1[i] - x2[i], 2);
    }
    return exp(-t * r);
}
