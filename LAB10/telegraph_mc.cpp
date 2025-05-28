#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <fstream>
#include <random>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

// Parametry linii transmisyjnej
const double L = 0.25e-6; // H/m
const double C = 100e-12; // F/m
const double R = 12.5;    // Ohm/m
const double G = 0.5e-3;  // S/m
const double l = 2.0;     // m
const double Rg = 75.0;   // Ohm
const double Rl = 12.5;   // Ohm

// Parametry źródła
const double nu = 1e9;        // Hz
const double t0 = 7.5e-9;     // s
const double sigma = 0.75e-9; // s

// Wyliczane parametry
const double c = 1.0 / sqrt(L * C);
const double R0 = sqrt(L / C);
const double mu = G / C;
const double lambda = 0.5 * (R / L - G / C);
const double zeta = R0 / (R0 + Rg);
const double Gamma_g = (Rg - R0) / (Rg + R0);
const double Gamma_l = (Rl - R0) / (Rl + R0);

extern double u_xt_exact(double x, double t, double t0, double freq, double sigma, double R, double G, double L, double C, double Rg, double Rl, double length, int number_nodes, int n_sum_terms);

double rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(0, 1);

    return dist(gen);
}
string format_filename(int npaths)
{
    ostringstream oss;
    oss << "data/u(x,t)_npaths=" << npaths << ".txt";
    return oss.str();
}

// Potencjał źródła
double Vg(double t)
{
    return sin(2 * M_PI * nu * t) * exp(-pow(t - t0, 2) / (2 * sigma * sigma));
}

// Główna funkcja do wyznaczania f(x, t) lub b(x, t)
double simulate_fb(double xstart, double tstart, int npaths, int sign)
{
    double sum = 0.0;

    for (int n = 0; n < npaths; ++n)
    {
        double x = xstart;
        double t = tstart;
        double eta = 1.0;
        int current_sign = sign;

        while (t > 0)
        {
            double U1 = rnd();
            double s = -log(U1) / (lambda + mu);

            if (current_sign == -1)
            { // f
                if ((x - c * s) > 0)
                {
                    eta *= (lambda / (lambda + mu));
                }
                else
                {
                    s = x / c;
                    sum += eta * zeta * Vg(t - s);
                    eta *= Gamma_g;
                }
                x = x - c * s;
            }
            else if (current_sign == 1)
            { // b
                if ((x + c * s) < l)
                {
                    eta *= (lambda / (lambda + mu));
                }
                else
                {
                    s = (l - x) / c;
                    eta *= Gamma_l;
                }
                x = x + c * s;
            }

            t -= s;
            current_sign = -current_sign;
        }
    }

    return sum / npaths;
}

int main()
{
    // Parametry symulacji
    double t_max = 50e-9;                                // czas [s]
    double dt = 1e-9;                                    // krok w czasie [s]
    vector<int> npaths = {(int)1e3, (int)1e4, (int)1e5}; // liczba ścieżek MC
    double dx = 0.01;                                    // krok w przestrzeni [m]

    for (auto npath : npaths)
    {
        string filename = format_filename(npath);
        ofstream file(filename);
        for (double t = 0.0; t <= t_max + 1e-12; t += dt)
        {
            file << t;
            for (double x = 0.0; x <= l + 1e-12; x += dx)
            {
                double f_xt = simulate_fb(x, t, npath, -1);
                double b_xt = simulate_fb(x, t, npath, 1);
                double u_xt = f_xt + b_xt;
                file << " " << u_xt;
            }
            file << "\n";
            cout << "Zapisano t = " << t * 1e9 << " ns " << "npaths= " << npath << endl;
        }
        file.close();
    }
    cout << "Zapisuje u(x,t) exact" << endl;
    ofstream file_exact("data/u(x,t)_exact.txt");
    for (double t = 0.0; t <= t_max + 1e-12; t += dt)
    {
        file_exact << t;
        for (double x = 0.0; x <= l + 1e-12; x += dx)
        {
            double u_xt = u_xt_exact(x, t, t0, nu, sigma, R, G, L, C, Rg, Rl, l, 1000, 100);
            file_exact << " " << u_xt;
        }
        file_exact << "\n";
        cout << "Zapisano t = " << t * 1e9 << " ns dla exact" << endl;
    }
    file_exact.close();

    return 0;
}