#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <iomanip>
#include <filesystem>
#include <sstream>

using namespace std;

const int M = 200;       // liczba podprzedziałów histogramu
const double rmax = 8.0; // maksymalny promień
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
    for (int exp = 2; exp <= log10(max_n); ++exp)
    {
        int base = pow(10, exp);
        for (int i = 1; i <= 100; ++i)
        {
            int val = i * base / 100;
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

void metropolis_mc(double a, double c, double delta_r, int N, double &energy, double &variance, vector<double> *hist = nullptr, const unordered_set<int> *frame_steps = nullptr, bool save_hist = false)
{
    double r = 1.0; // startowa wartość
    double sum = 0.0, sum2 = 0.0;

    int burn_in = N / 10;
    int count = 0;

    if (hist)
    {
        hist->assign(M, 0.0);
    }

    for (int i = 0; i < N + burn_in; ++i)
    {
        double r_new = r + delta_r * (2.0 * uniform(rng) - 1.0);
        if (r_new <= 0.0)
            continue;

        double p_old = p_density(r, a, c);
        double p_new = p_density(r_new, a, c);

        double accept_prob = min(1.0, p_new / p_old);
        if (uniform(rng) < accept_prob)
            r = r_new;

        if (i >= burn_in)
        {
            double e_loc = epsilon_local(r, a, c);
            sum += e_loc;
            sum2 += e_loc * e_loc;
            ++count;

            if (hist && r < rmax)
            {
                int bin = static_cast<int>(floor(r / dr_hist));
                (*hist)[bin] += 1.0 / (N * dr_hist);

                if (frame_steps && frame_steps->count(count) && save_hist)
                {
                    // zapis histogramu
                    ostringstream fname_hist;
                    fname_hist << "data/hist_" << setfill('0') << setw(4) << count << ".dat";
                    ofstream fout_hist(fname_hist.str());
                    for (int k = 0; k < M; ++k)
                    {
                        double r_val = (k + 0.5) * dr_hist;
                        double exact = (a == 1.0 && c == 0.0) ? 4 * r_val * r_val * exp(-2 * r_val) : 0.0;
                        fout_hist << r_val << " " << (*hist)[k] << " " << exact << "\n";
                    }
                }
            }
        }
    }

    energy = sum / count;
    variance = sum2 / count - energy * energy;
}

void generate_histogram_sequence(double a, double c, double delta_r)
{
    int N = 1e6;
    unordered_set<int> hist_frames = generate_log_intervals(N);

    vector<double> hist(M, 0.0);
    double E, var;
    metropolis_mc(a, c, delta_r, N, E, var, &hist, &hist_frames, true);
}

void generate_energy_map_sequence(double delta_r)
{
    int N = 1e6;
    unordered_set<int> frames = generate_log_intervals(N);

    for (int frame = 0; frame < steps.size(); ++frame)
    {
        int N = steps[frame];
        ostringstream fname;
        fname << "data/en_" << setfill('0') << setw(3) << frame << ".dat";
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
}

int main()
{
    const double delta_r = 0.1;

    cout << "[1] Generowanie mapy energii w krokach animacyjnych...\n";
    generate_energy_map_sequence(delta_r);

    cout << "[2] Generowanie histogramu w krokach animacyjnych dla a=1, c=0...\n";
    generate_histogram_sequence(1.0, 0.0, delta_r);

    cout << "Zakończono. Dane klatkowe zapisane w folderze 'data'.\n";
    return 0;
}
