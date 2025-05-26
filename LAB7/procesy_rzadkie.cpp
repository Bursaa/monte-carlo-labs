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
#include <sstream>
#include <iomanip>

using namespace std;

const double k1 = 1.0;
const double k2 = 1.0;
const double k3 = 0.001;
const double k4 = 0.01;
const int N = 50;
const int Pmax = 5;
const double tmax = 200.0;
const int x1_0 = 120;
const int x2_0 = 80;
const int x3_0 = 1;
const double delta_t = tmax / N;

double rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(0, 1);

    return dist(gen);
}

double Gamma(int gamma_number, int x1, int x2, int x3)
{
    if (gamma_number == 1)
    {
        return k1;
    }
    if (gamma_number == 2)
    {
        return k2;
    }
    if (gamma_number == 3)
    {
        return k3 * x1 * x2;
    }
    if (gamma_number == 4)
    {
        return k4 * x3;
    }
    return 0.0;
}

tuple<int, int, int> Gamma_affection(int m_number, int x1, int x2, int x3)
{
    if (m_number == 1)
    {
        return {x1 + 1, x2, x3};
    }
    if (m_number == 2)
    {
        return {x1, x2 + 1, x3};
    }
    if (m_number == 3)
    {
        return {x1 - 1, x2 - 1, x3 + 1};
    }
    if (m_number == 4)
    {
        return {x1, x2, x3 - 1};
    }
    return {x1, x2, x3};
}

int get_m(double Gamma_max, int x1, int x2, int x3)
{
    double U2 = rnd();
    for (int n = 1; n <= 4; n++)
    {
        double licz = 0;
        for (int i = 1; i <= n; i++)
        {
            licz += Gamma(i, x1, x2, x3);
        }
        if (U2 <= licz / Gamma_max)
        {
            return n;
        }
    }
    return 4;
}

void zad1()
{
    vector<double> h0(N, 0.0);
    vector<double> h1(N, 0.0);
    vector<double> h2(N, 0.0);
    vector<int> ncount(N, 0);
    for (int p = 0; p < Pmax; p++)
    {
        double t = 0.0;
        int x1 = x1_0;
        int x2 = x2_0;
        int x3 = x3_0;
        string filename1 = "x(t)_p=" + to_string(p + 1) + ".txt";
        ofstream x_file("data/" + filename1);
        x_file << t << " " << x1 << " " << x2 << " " << x3 << endl;
        while (t <= tmax)
        {
            double Gamma_max = Gamma(1, x1, x2, x3) + Gamma(2, x1, x2, x3) + Gamma(3, x1, x2, x3) + Gamma(4, x1, x2, x3);
            double dt = -1 / Gamma_max * log(rnd());
            int m = get_m(Gamma_max, x1, x2, x3);
            tie(x1, x2, x3) = Gamma_affection(m, x1, x2, x3);
            t += dt;
            x_file << t << " " << x1 << " " << x2 << " " << x3 << endl;
        }
        x_file.close();
    }
}

int main()
{
    zad1();
    return 0;
}