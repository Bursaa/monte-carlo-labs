#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <functional>

using namespace std;

const double M = 10;
double rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dist(0.0, 1.0); // Poprawiony typ

    return dist(gen);
}

int randomInt(int N)
{
    return static_cast<int>(rnd() * N);
}

double C1(double x)
{
    return 1 + tanh(x);
}
double C2(double x)
{
    return 1 / (1 + x * x);
}
double C3(double x)
{
    return pow(cos(M_PI * x), 10.0);
}
tuple<double, double> elementary_method(int N, double a, double b, int M, function<double(double)> fun, vector<double> &interval)
{
    double g1 = 0;
    double g2 = 0;

    double step = (b - a) / M;

    for (int i = 0; i < N; i++)
    {
        double u = rnd();
        double xi = a + (b - a) * u;

        double fxi = fun(xi);
        g1 += (b - a) * fxi;
        g2 += pow((b - a) * fxi, 2.0);

        int k = static_cast<int>((xi - a) / step);
        if (k == M)
            k = M - 1; // Na wypadek gdy xi == b
        interval[k] += 1.0 / N;
    }
    double sigma = (g2 / N - (g1 / N) * (g1 / N)) / N;

    return {g1 / N, sigma};
}

tuple<double, double, vector<double>> systematical_draw_method(int N, double a, double b, int M, function<double(double)> fun, vector<double> &interval)
{
    vector<double> g1_m(M);
    vector<double> g2_m(M);
    vector<double> p_i(M);
    vector<double> sigma_m(M);
    double g1 = 0;
    double sigma = 0;
    p_i.assign(M, 1.0 / M);
    double step = (b - a) / M;
    int Nm = N / M;

    for (int m = 0; m < M; m++)
    {
        double xm = a + m * step;
        double xm1 = xm + step;
        for (int i = 0; i < Nm; i++)
        {
            double u = rnd();
            double xi = xm + step * u;
            double fxi = fun(xi);
            g1_m[m] += (b - a) * fxi / Nm;
            g2_m[m] += pow((b - a) * fxi, 2.0) / Nm;
            interval[m] += 1.0 / N;
        }
        g1 += p_i[m] * g1_m[m];
        sigma += p_i[m] * p_i[m] * (g2_m[m] - g1_m[m] * g1_m[m]) / Nm;
        sigma_m[m] = sqrt(g2_m[m] - g1_m[m] * g1_m[m]);
    }

    return {g1, sigma, sigma_m};
}

tuple<double, double> layer_draw_method(int N, double a, double b, int M, function<double(double)> fun, vector<double> &interval)
{
    vector<double> g1_m(M);
    vector<double> g2_m(M);
    vector<double> p_i(M);
    double g1 = 0;
    double sigma = 0;
    p_i.assign(M, 1.0 / M);
    double step = (b - a) / M;

    int Ns = N >= 1000 ? 1000 : 100;
    auto [g, sig, sigma_m] = systematical_draw_method(Ns, a, b, M, fun, interval);
    double sum_sigma = 0;
    for (int m = 0; m < M; m++)
    {
        sum_sigma += sigma_m[m] * p_i[m];
    }

    interval.assign(M, 0);
    for (int m = 0; m < M; m++)
    {
        double xm = a + m * step;
        double xm1 = xm + step;
        double Nm = p_i[m] * sigma_m[m] / sum_sigma * N;
        for (int i = 0; i < Nm; i++)
        {
            double u = rnd();
            double xi = xm + (xm1 - xm) * u;
            double fxi = fun(xi);
            g1_m[m] += (b - a) * fxi / Nm;
            g2_m[m] += pow((b - a) * fxi, 2.0) / Nm;
            interval[m] += 1.0 / N;
        }
        g1 += p_i[m] * g1_m[m];
        sigma += p_i[m] * p_i[m] * (g2_m[m] - g1_m[m] * g1_m[m]) / Nm;
    }

    return {g1, sigma};
}

int main()
{
    vector<int> k_values = {2, 3, 4, 5, 6, 7};
    vector<function<double(double)>> functions = {C1, C2, C3};
    vector<string> names = {"C1", "C2", "C3"};
    vector<double> a = {-3, 0, 0};
    vector<double> b = {3, 10, 1};
    for (int i = 0; i <= functions.size(); i++)
    {
        for (auto k : k_values)
        {
            int N = pow(10, k);
            vector<double> interval1(M);
            vector<double> interval2(M);
            vector<double> interval3(M);
            interval1.assign(M, 0.0);
            interval2.assign(M, 0.0);
            interval3.assign(M, 0.0);

            auto [g1, sig1] = elementary_method(N, a[i], b[i], M, functions[i], interval1);
            auto [g2, sig2, sigm] = systematical_draw_method(N, a[i], b[i], M, functions[i], interval2);
            auto [g3, sig3] = layer_draw_method(N, a[i], b[i], M, functions[i], interval3);

            string filename1 = names[i] + "_intervals_k=" + to_string(k) + ".txt";
            std::string filename2 = names[i] + "_values_k=" + to_string(k) + ".txt";
            ofstream file1("data/" + filename1);
            ofstream file2("data/" + filename2);

            for (int j = 0; j < M; j++)
            {
                file1 << interval1[j] << " " << interval2[j] << " " << interval3[j] << endl;
            }
            file2 << g1 << " " << sig1 << endl;
            file2 << g2 << " " << sig2 << endl;
            file2 << g3 << " " << sig3 << endl;

            file1.close();
            file2.close();
        }
    }
    return 0;
}