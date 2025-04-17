#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
#include "particle_translation.h"
#include <sstream>
#include <iomanip>

using namespace std;

const double D = 1.0;
const double dt = 0.1;
const double sigma_dt = sqrt(2 * D * dt);

double norm_rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::normal_distribution<> normal_dist(0, sigma_dt);

    return normal_dist(gen);
}

double sum(const std::vector<double> &x, std::vector<double> y = {})
{
    if (y.empty())
    {
        y = std::vector<double>(x.size(), 1.0);
    }

    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        result += x[i] * y[i];
    }

    return result / x.size();
}

void zad1()
{
    const int t_max = 100;
    vector<int> k_values = {2, 3, 4, 5};
    const int time_steps = (int)t_max / dt;

    for (auto k : k_values)
    {
        int N = pow(10, k);
        vector<vector<double>> x(time_steps, vector<double>(N, 0.0));
        vector<vector<double>> y(time_steps, vector<double>(N, 0.0));
        string filename1 = "Zad1_x(i,t)_k=" + to_string(k) + ".txt";
        string filename2 = "Zad1_y(i,t)_k=" + to_string(k) + ".txt";
        string filename3 = "D_values_k=" + to_string(k) + ".txt";
        ofstream x_file("data/" + filename1);
        ofstream y_file("data/" + filename2);
        ofstream D_file("data/" + filename3);

        for (int t_index = 0; t_index < time_steps - 1; t_index++)
        {
            for (int i = 0; i < N; i++)
            {
                double dx = norm_rnd();
                double dy = norm_rnd();

                x[t_index + 1][i] = x[t_index][i] + dx;
                y[t_index + 1][i] = y[t_index][i] + dy;

                x_file << x[t_index + 1][i] << " ";
                y_file << y[t_index + 1][i] << " ";
            }
            double x_sr = sum(x[t_index + 1]);
            double y_sr = sum(y[t_index + 1]);
            double x_sr2 = sum(x[t_index + 1], x[t_index + 1]);
            double y_sr2 = sum(y[t_index + 1], y[t_index + 1]);
            double xy_sr = sum(x[t_index + 1], y[t_index + 1]);
            double Dxx = (x_sr2 - x_sr * x_sr) / (2 * (t_index + 1) * dt);
            double Dyy = (y_sr2 - y_sr * y_sr) / (2 * (t_index + 1) * dt);
            double Dxy = (xy_sr - x_sr * y_sr) / (2 * (t_index + 1) * dt);
            D_file << Dxx << " " << Dyy << " " << Dxy << endl;

            x_file << endl;
            y_file << endl;
        }
        x_file.close();
        y_file.close();
        D_file.close();
    }
}

void zad2(double Ra, double omega)
{
    const int Nmax = 10000;
    const double D = 1.0;
    const double tmax = 1000.0;
    const int N = tmax / dt;
    const double xr = 0.0, yr = 0.0, Rr = 5.0; // główny obszar
    const double xa = 3.0, ya = 0.0;           // absorber
    const double xs = -4.5, ys = 0.0;          // źródło
    const double dn = omega * dt;
    vector<bool> theta(Nmax, 0);
    vector<vector<double>> x(N, vector<double>(Nmax, xs));
    vector<vector<double>> y(N, vector<double>(Nmax, ys));
    std::ostringstream suffix;
    suffix << "_omega=" << omega << "_Ra=" << std::fixed << std::setprecision(1) << Ra;

    string filename1 = "Zad2_x(i,t)" + suffix.str() + ".txt";
    string filename2 = "Zad2_y(i,t)" + suffix.str() + ".txt";
    string filename3 = "Zad2_n(t)" + suffix.str() + ".txt";
    ofstream x_file("data/" + filename1);
    ofstream y_file("data/" + filename2);
    ofstream n_file("data/" + filename3);

    n_file << 0 << " ";
    for (int it = 0; it < N - 1; it++)
    {
        int n = 0;
        int n_new = 0;
        for (int i = 0; i < Nmax; i++)
        {
            if (theta[i] == 0 && n_new < dn)
            {
                theta[i] = 1;
                x[it][i] = xs;
                y[it][i] = ys;
                n_new++;
            }
            if (theta[i] == 1)
            {
                double x1 = x[it][i];
                double y1 = y[it][i];
                bool theta_i = theta[i];
                double length = 2 * Rr;
                double dx = norm_rnd();
                double dy = norm_rnd();

                double x2 = x1 + dx;
                double y2 = y1 + dy;

                do
                {
                    particle_translation(x1, y1, x2, y2, xr, yr, Rr, xa, ya, Ra, theta_i, length);
                } while (length > pow(10, -6));
                x[it + 1][i] = x1;
                y[it + 1][i] = y1;
                theta[i] = theta_i;
                if (theta_i == 1)
                {
                    n++;
                }
            }
            x_file << x[it + 1][i] << " ";
            y_file << y[it + 1][i] << " ";
        }
        n_file << n << " ";
        x_file << endl;
        y_file << endl;
    }
    x_file.close();
    y_file.close();
    n_file.close();
}

int main()
{
    zad1();
    vector<double> Ra_arr = {0.1, 0.5};
    vector<double> omega_arr = {10, 50, 100};
    for (auto omega : omega_arr)
    {
        for (auto Ra : Ra_arr)
        {
            zad2(Ra, omega);
        }
    }
    return 0;
}