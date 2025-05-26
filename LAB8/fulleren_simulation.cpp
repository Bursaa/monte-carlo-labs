#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <chrono>
#include "vector_operator.h"

using namespace std;
using namespace std::chrono;

const double R0 = 1.315;  // A
const double R1 = 1.7;    // A
const double R2 = 2.0;    // A
const double D_e = 6.325; // eV
const double S = 1.29;
const double lambda = 1.5; // A^(-1)
const double delta = 0.80469;
const double a_0 = 0.011304;
const double c_0 = 19.0;
const double d_0 = 2.5;

double rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(0, 1);

    return dist(gen);
}
double norm(const vector<double> &r)
{
    double result = 0.0;
    for (int i = 0; i < r.size(); i++)
    {
        result += r[i] * r[i];
    }
    return sqrt(result);
}

vector<double> cart2spher(const vector<double> &r)
{
    double r_norm = norm(r);
    double phi = atan2(r[1], r[0]);
    double theta = acos(r[2] / r_norm);
    return {r_norm, phi, theta};
}

vector<double> spher2cart(const vector<double> &r)
{
    double x = r[0] * sin(r[2]) * cos(r[1]);
    double y = r[0] * sin(r[2]) * sin(r[1]);
    double z = r[0] * cos(r[2]);
    return {x, y, z};
}
int found_peak(const vector<double> &pcf)
{
    int n = pcf.size();
    if (n == 0)
        return -1; // brak danych

    int max_index = 0;
    double max_value = pcf[0];

    for (int i = 1; i < n; ++i)
    {
        if (pcf[i] > max_value)
        {
            max_value = pcf[i];
            max_index = i;
        }
    }

    return max_index;
}

double f_cut(double &r)
{
    if (r <= R1)
    {
        return 1;
    }
    else if (r > R2)
    {
        return 0;
    }
    else
    {
        return 1 / 2 * (1 + cos((r - R1) / (R2 - R1) * M_PI));
    }
}

double V_R(double &r)
{
    return D_e / (S - 1) * exp(-sqrt(2.0 * S) * lambda * (r - R0));
}

double V_A(double &r)
{
    return D_e * S / (S - 1) * exp(-sqrt(2.0 / S) * lambda * (r - R0));
}

double g(double &cos_theta)
{
    return a_0 * (1 + pow(c_0 / d_0, 2.0) - (c_0 * c_0 / (pow(d_0, 2.0) + pow(1 + cos_theta, 2.0))));
}

double B(const vector<vector<double>> &r, int &i, int &j, const vector<double> &r_ij_vec, double &norm_rij, bool &updated_B)
{
    double zeta_ij = 0.0;

    for (int k = 0; k < r.size(); ++k)
    {
        if (k == i || k == j)
            continue;

        vector<double> r_ik_vec = r[k] - r[i];
        double norm_rik = norm(r_ik_vec);
        if (norm_rik >= R2)
            continue;

        double cos_theta = dot_product(r_ij_vec, r_ik_vec) / (norm_rij * norm_rik);
        if (cos_theta > 0 && updated_B)
            zeta_ij = 10;
        zeta_ij += f_cut(norm_rik) * g(cos_theta);
    }

    return 1.0 / pow(1.0 + zeta_ij, delta);
}

double calculate_V(const vector<vector<double>> &r, int &i, bool &updated_B)
{
    double V_i = 0.0;
    for (int j = i + 1; j < r.size(); ++j)
    {
        vector<double> r_ij_vec = r[j] - r[i];
        double norm_rij = norm(r_ij_vec);
        if (norm_rij > R2)
            continue;

        vector<double> reverser = {-1.0, -1.0, -1.0};
        vector<double> r_ji_vec = reverser * r_ij_vec;
        double B1 = B(r, i, j, r_ij_vec, norm_rij, updated_B);
        double B2 = B(r, j, i, r_ji_vec, norm_rij, updated_B);
        double B_dash = 0.5 * (B1 + B2);
        double fcut_val = f_cut(norm_rij);
        double VR_val = V_R(norm_rij);
        double VA_val = V_A(norm_rij);
        V_i += fcut_val * (VR_val - B_dash * VA_val);
    }

    return V_i;
}
double calculate_V_tot(vector<vector<double>> &r, bool updated_B = false)
{
    double V_tot = 0.0;
    for (int i = 0; i < r.size(); i++)
    {
        V_tot += calculate_V(r, i, updated_B);
    }
    return V_tot;
}

vector<double> compute_pcf(const vector<vector<double>> &r, int M)
{

    vector<double> pcf(M, 0.0);
    int n = r.size();
    double r_sr = 0.0;
    for (int i = 0; i < n; ++i)
    {
        r_sr += norm(r[i]);
    }
    r_sr /= n;

    double Omega = 4.0 * M_PI * pow(r_sr, 2);
    double r_max = 2.5 * r_sr;
    double delta_r = r_max / M;

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double r_ij = norm(r[j] - r[i]);
            int m = floor(r_ij / delta_r);
            if (m < M)
            {
                double dOmega = 2 * M_PI * r_ij * delta_r;
                pcf[m] += 2.0 * Omega / (n * n * dOmega);
            }
        }
    }
    return pcf;
}

vector<vector<double>> loadDataFromFile(const string &filename)
{
    ifstream file(filename);
    vector<vector<double>> data;

    if (!file.is_open())
    {
        cerr << "Nie można otworzyć pliku: " << filename << endl;
        return data;
    }

    string line;
    while (getline(file, line))
    {
        vector<double> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, ','))
        {
            row.push_back(stod(value));
        }
        data.push_back(row);
    }

    file.close();
    return data;
}

void mc_move_single_atom(vector<vector<double>> &r, int &i, double &beta, double &omega_r, double &omega_phi, double &omega_theta, bool &updated_B)
{
    const vector<double> r_old = r[i];               // zachowaj starą pozycję
    const vector<double> r_spher = cart2spher(r[i]); // konwersja tylko raz

    double u1 = rnd();
    double u2 = rnd();
    double u3 = rnd();
    double u4 = rnd();

    double r_new = r_spher[0] + r_spher[0] * (2 * u1 - 1) * omega_r;
    double phi_new = r_spher[1] + r_spher[1] * (2 * u2 - 1) * omega_phi;
    double theta_new = r_spher[2] + r_spher[2] * (2 * u3 - 1) * omega_theta;

    // Napraw zakresy
    if (phi_new < 0)
        phi_new += 2 * M_PI;
    if (phi_new > 2 * M_PI)
        phi_new -= 2 * M_PI;
    if (theta_new < 0 || theta_new > M_PI)
        theta_new = r_spher[2];

    vector<double> r_new_cart = spher2cart({r_new, phi_new, theta_new});

    double V_old = calculate_V(r, i, updated_B);
    r[i] = r_new_cart;
    double V_new = calculate_V(r, i, updated_B);

    if (exp(-beta * (V_new - V_old)) < u4)
    {
        r[i] = r_old; // odrzuć krok
    }
}

void mc_change_r(vector<vector<double>> &r, double &beta, double &W_all, bool &updated_B)
{
    double u1 = rnd();
    double u2 = rnd();
    vector<vector<double>> r_new = r;
    double scaler = (1 + W_all * (2 * u1 - 1));
    vector<double> vec_scaler = {scaler, scaler, scaler};
    for (int i = 0; i < r.size(); i++)
    {
        r_new[i] = r[i] * vec_scaler;
    }
    double V_tot_old = calculate_V_tot(r, updated_B);
    double V_tot_new = calculate_V_tot(r_new, updated_B);
    double p_acc = min(1.0, exp(-beta * (V_tot_new - V_tot_old)));
    if (p_acc >= u2)
    {
        r = r_new;
    }
}

void SA_algorithm(int n, double r_init, double it_max, int M, double omega_r, double omega_theta, double omega_phi, double W_all, double beta_min, double beta_max, double p, bool updated_B, string filename)
{
    ofstream values("data/values_" + filename + ".txt");
    vector<vector<double>> r(n, vector<double>(3, 0.0));
    for (int i = 0; i < n; i++)
    {
        double r_mag = r_init;
        double phi = 2 * M_PI * rnd();
        double theta = M_PI * rnd();
        r[i] = spher2cart({r_mag, phi, theta});
    }
    for (double it = 0.0; it < it_max; it += 1.0)
    {
        double beta = beta_min + pow(it / it_max, p) * (beta_max - beta_min);
        for (int i = 0; i < n; i++)
        {
            mc_move_single_atom(r, i, beta, omega_r, omega_phi, omega_theta, updated_B);
        }
        mc_change_r(r, beta, W_all, updated_B);
        if ((static_cast<int>(it) + 1) % 100 == 0)
        {
            double V_t = calculate_V_tot(r, updated_B);
            double r_sr = 0.0;
            for (int i = 0; i < n; i++)
            {
                r_sr += norm(r[i]);
            }
            r_sr /= n;
            values << it << " " << V_t << " " << r_sr << " " << beta << endl;
        }
        if ((static_cast<int>(it) + 1) % 1000 == 0)
        {
            cout << it << endl;
        }
    }
    string fname2 = "data/pcf_" + filename + ".txt";
    ofstream pcf_file(fname2);
    vector<double> pcf = compute_pcf(r, M);
    for (auto p : pcf)
    {
        pcf_file << p << endl;
    }
    string fname = "data/polygons_" + filename + ".txt";
    int m = found_peak(pcf);
    double r_sr = 0.0;
    for (int i = 0; i < n; i++)
    {

        r_sr += norm(r[i]);
    }
    r_sr /= n;
    double r_max = 2.5 * r_sr;
    double delta_r = r_max / M;
    cout << "DONE2" << endl;
    write_polygons_from_atoms(1.5, n, r, fname.c_str());
    values.close();
}

void zad1()
{
    const int M = 100;
    cout << "ZADANIE 1" << endl;
    string filename = "data/data1.txt";
    vector<vector<double>> data = loadDataFromFile(filename);
    double V_tot = calculate_V_tot(data);
    cout << "V_tot: " << V_tot << endl;
    vector<double> pcf = compute_pcf(data, M);
    string fname = "data/polygons_zad1.txt";
    string fname2 = "data/pcf_zad1.txt";
    ofstream file(fname2);
    for (int i = 0; i < M; i++)
    {
        file << pcf[i] << endl;
    }
    write_polygons_from_atoms(1.6, 60, data, fname);
    file.close();
}

void zad2()
{
    cout << "ZADANIE 2" << endl;
    const int n = 60;
    const double omega_r = pow(10.0, -4.0);
    const double omega_theta = 0.05;
    const double omega_phi = 0.05;
    const double W_all = pow(10.0, -4.0);
    const int M = 100;
    const double beta_min = 1.0;
    const double beta_max = 100.0;
    const double p = 2.0;
    const double it_max = 100000.0;
    const double r_init = 3.5;
    const bool updated_B = false;

    SA_algorithm(n, r_init, it_max, M, omega_r, omega_theta, omega_phi, W_all, beta_min, beta_max, p, updated_B, "zad2");
}
void zad3()
{
    cout << "ZADANIE 3" << endl;
    const int n = 60;
    const double omega_r = pow(10.0, -4.0);
    const double omega_theta = 0.05;
    const double omega_phi = 0.05;
    const double W_all = pow(10.0, -5.0);
    const int M = 100;
    const double beta_min = 1.0;
    const double beta_max = 100.0;
    const double p = 2.0;
    const double it_max = 100000.0;
    const double r_init = 3.5;
    const bool updated_B = true;

    SA_algorithm(n, r_init, it_max, M, omega_r, omega_theta, omega_phi, W_all, beta_min, beta_max, p, updated_B, "zad3");
}

void zad4()
{
    cout << "ZADANIE 4" << endl;
    const int n = 60;
    const double omega_r = pow(10.0, -5.0);
    const double omega_theta = 0.05;
    const double omega_phi = 0.05;
    const double W_all = pow(10.0, -5.0);
    const int M = 100;
    const double beta_min = 1.0;
    const double beta_max = 100.0;
    const double p = 2.0;
    const double it_max = 100000.0;
    const double r_init = 2.5;
    const bool updated_B = true;

    SA_algorithm(n, r_init, it_max, M, omega_r, omega_theta, omega_phi, W_all, beta_min, beta_max, p, updated_B, "zad4");
}

void zad5()
{
    cout << "ZADANIE 5" << endl;
    const int n = 60;
    const double omega_r = pow(10.0, -5.0);
    const double omega_theta = 0.05;
    const double omega_phi = 0.05;
    const double W_all = pow(10.0, -5.0);
    const int M = 100;
    const double beta_min = 1.0;
    const double beta_max = 100.0;
    const double p = 1.5;
    const double it_max = 100000.0;
    const double r_init = 2.5;
    const bool updated_B = true;

    SA_algorithm(n, r_init, it_max, M, omega_r, omega_theta, omega_phi, W_all, beta_min, beta_max, p, updated_B, "zad5");
}

void zad6()
{
    cout << "ZADANIE 6" << endl;
    const double omega_r = pow(10.0, -4.0);
    const double omega_theta = 0.05;
    const double omega_phi = 0.05;
    const double W_all = pow(10.0, -4.0);
    const int M = 100;
    const double beta_min = 1.0;
    const double beta_max = 100.0;
    const double p = 2.0;
    const double it_max = 100000.0;
    const double r_init = 2.5;
    const bool updated_B = true;

    for (int n = 30; n <= 40; n++)
    {
        string filename = "zad6_n=" + to_string(n);
        SA_algorithm(n, r_init, it_max, M, omega_r, omega_theta, omega_phi, W_all, beta_min, beta_max, p, updated_B, filename);
    }
}

int main()
{
    zad1();
    zad2();
    zad3();
    zad4();
    zad5();
    return 0;
}