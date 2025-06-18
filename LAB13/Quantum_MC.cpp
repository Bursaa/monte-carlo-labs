#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <iomanip>
#include <filesystem>
#include <sstream>
#include <unordered_set>

using namespace std;

const int M = 200;        // liczba podprzedziałów histogramu
const double rmax = 12.0; // maksymalny promień
const double dr_hist = rmax / M;
double rnd()
{
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_real_distribution<> dist(0, 1);

    return dist(gen);
}
unordered_set<int> generate_log_intervals(int max_n)
{
    unordered_set<int> frames;
    for (int exp = 0;; ++exp)
    {
        int step = pow(10, exp);
        int end = step * 100;
        if (step > max_n)
            break;

        for (int val = step; val < end; val += step)
        {
            if (val <= max_n)
                frames.insert(val);
        }
    }
    return frames;
}

double PsiT(double r, double a, double c)
{
    return (1.0 + c * r) * exp(-a * r);
}

double epsilon_local(double r, double a, double c)
{
    double numerator = -a * a * c * r * r + (-a * a + 4 * a * c - 2 * c) * r + 2 * a - 2 * c - 2;
    double denominator = 2 * c * r * r + 2 * r;
    return numerator / denominator;
}

double p_density(double r, double a, double c)
{
    double psi = PsiT(r, a, c);
    return r * r * psi * psi;
}

void metropolis_mc(double a, double c, double delta_r, int N, double &energy, double &variance, const unordered_set<int> *frame_steps = nullptr, vector<double> *hist = nullptr)
{
    double r = rnd();
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < N; ++i)
    {
        double r_new = r + delta_r * (2.0 * rnd() - 1.0);

        double p_old = p_density(r, a, c);
        double p_new = p_density(r_new, a, c);

        double accept_prob = min(1.0, p_new / p_old);
        if (rnd() <= accept_prob && r_new > 0)
            r = r_new;

        double e_loc = epsilon_local(r, a, c);
        sum += e_loc;
        sum2 += e_loc * e_loc;
        if (hist && r < rmax)
        {
            int bin = static_cast<int>(floor(r / dr_hist));
            (*hist)[bin] += 1.0 / (N * dr_hist);
        }
        if (frame_steps && frame_steps->count(i + 1) && hist)
        {

            // zapis histogramu
            ostringstream fname_hist;
            fname_hist << "data/hist_a=" << a
                       << "_c=" << c
                       << "_" << i + 1 << ".dat";
            ofstream fout_hist(fname_hist.str());
            for (int k = 0; k < M; ++k)
            {
                double r_val = (k + 0.5) * dr_hist;
                double exact = (a == 1.0 && c == 0.0) ? 2 * exp(-1 * r_val) : 1 / (2 * sqrt(2)) * (2 - r_val) * exp(-r_val / 2);
                fout_hist << r_val << " " << (*hist)[k] << " " << r_val * r_val * exact * exact << "\n";
            }
        }
    }

    energy = sum / N;
    variance = sum2 / N - energy * energy;
}

void generate_histogram_sequence(double a, double c, double delta_r)
{
    int N = 1e8;
    unordered_set<int> hist_frames = generate_log_intervals(N);

    vector<double> hist(M, 0.0);
    double E, var;
    metropolis_mc(a, c, delta_r, N, E, var, &hist_frames, &hist);
}

void generate_energy_map_sequence(double delta_r)
{
    int N = 1e6;
    ostringstream fname;
    fname << "data/Energies" << ".dat";
    ofstream fout(fname.str());
    for (double a = 0.3; a <= 1.2; a += 0.02)
    {
        for (double c = -0.7; c <= 0.3; c += 0.02)
        {
            double E, var;
            metropolis_mc(a, c, delta_r, N, E, var);
            fout << fixed << setprecision(5)
                 << a << " " << c << " " << E << " " << var << " " << sqrt(var) << "\n";
        }
    }
}

int main()
{
    const double delta_r = 0.1;

    cout << "[1] Generowanie mapy energii...\n";
    generate_energy_map_sequence(delta_r);

    cout << "[2] Generowanie histogramu w krokach animacyjnych dla a=1, c=0...\n";
    generate_histogram_sequence(1.0, 0.0, delta_r);

    cout << "[3] Generowanie histogramu w krokach animacyjnych dla a=0.5, c=-0.5...\n";
    generate_histogram_sequence(0.5, -0.5, delta_r);

    cout << "Zakonczono. Dane klatkowe zapisane w folderze 'data'.\n";
    return 0;
}
