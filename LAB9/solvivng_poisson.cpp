#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <random>
#include <iomanip>

using namespace std;

const int nx = 30;
const int ny = 30;
const double delta = 0.1;
const double epsilon = 1.0;
const double VL = 1.0;
const double VB = -1.0;
const double VT = -1.0;
const double xmax = delta * nx;
const double ymax = delta * ny;
const double rho_max = 1.0;
const double sigma_rho = xmax / 10.0;

double rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist(0, 1);

    return dist(gen);
}

vector<vector<double>> init_rho()
{
    vector<vector<double>> rho(nx + 1, vector<double>(ny + 1));
    double x, y;
    for (int i = 0; i <= nx; ++i)
    {
        x = i * delta;
        for (int j = 0; j <= ny; ++j)
        {
            y = j * delta;
            double dx = x - xmax / 2.0;
            double dy = y - ymax / 2.0;
            rho[i][j] = rho_max * exp(-(dx * dx + dy * dy) / (2 * sigma_rho * sigma_rho));
        }
    }
    return rho;
}

void write_to_file(const string &filename, const vector<vector<double>> &mat)
{
    ofstream file("data/" + filename);
    for (int i = 0; i <= nx; ++i)
    {
        for (int j = 0; j <= ny; ++j)
        {
            file << setw(12) << mat[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

void relaxacja(vector<vector<double>> &Vrel, const vector<vector<double>> &rho)
{
    const double omega = 1.8;
    const double tol = 1e-6;
    const int itmax = 10000;
    double Fold = 0, Fnew = 0;

    // Warunki brzegowe Dirichleta
    for (int j = 0; j <= ny; ++j)
        Vrel[0][j] = VL * sin(M_PI * j * delta / ymax);

    for (int i = 0; i <= nx; ++i)
    {
        Vrel[i][0] = VB * sin(M_PI * i * delta / xmax);
        Vrel[i][ny] = VT * sin(M_PI * i * delta / xmax);
    }

    for (int it = 1; it < itmax; ++it)
    {
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                Vrel[i][j] = (1.0 - omega) * Vrel[i][j] + omega / 4.0 * (Vrel[i + 1][j] + Vrel[i - 1][j] + Vrel[i][j + 1] + Vrel[i][j - 1] + delta * delta * rho[i][j] / epsilon);
            }
        }

        for (int j = 1; j < ny; ++j) // Neumann
            Vrel[nx][j] = Vrel[nx - 1][j];

        Fold = Fnew;
        Fnew = 0.0;
        for (int i = 1; i < nx; ++i)
        {
            for (int j = 1; j < ny; ++j)
            {
                double Ex = (Vrel[i + 1][j] - Vrel[i - 1][j]) / (2.0 * delta);
                double Ey = (Vrel[i][j + 1] - Vrel[i][j - 1]) / (2.0 * delta);
                Fnew += (Ex * Ex + Ey * Ey) / 2.0 - rho[i][j] * Vrel[i][j];
            }
        }

        if (fabs((Fnew - Fold) / Fnew) < tol)
            break;
    }
}

void mc_method(vector<vector<double>> &Vmc, vector<vector<double>> &sigmaV, vector<vector<double>> &absorbed,
               const vector<vector<double>> &rho, int Nchains, int nlength, bool block_after_estimate)
{
    Vmc.assign(nx + 1, vector<double>(ny + 1, 0.0));
    sigmaV.assign(nx + 1, vector<double>(ny + 1, 0.0));
    absorbed.assign(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<int>> B(nx + 1, vector<int>(ny + 1, 0)); // Dirichlet mask

    // Warunki Dirichleta
    for (int j = 0; j <= ny; ++j)
    {
        Vmc[0][j] = VL * sin(M_PI * j * delta / ymax);
        B[0][j] = 1;
    }
    for (int i = 0; i <= nx; ++i)
    {
        Vmc[i][0] = VB * sin(M_PI * i * delta / xmax);
        Vmc[i][ny] = VT * sin(M_PI * i * delta / xmax);
        B[i][0] = B[i][ny] = 1;
    }

    for (int i0 = 1; i0 < nx; ++i0)
    {
        for (int j0 = 1; j0 < ny; ++j0)
        {
            double sum1 = 0, sum2 = 0;
            int kchains = 0;

            for (int n = 0; n < Nchains; ++n)
            {
                int i = i0, j = j0;
                double g = 0;

                for (int s = 0; s < nlength; ++s)
                {
                    int m = int(rnd() * 4);
                    if (m == 0)
                        i--;
                    else if (m == 1)
                        i++;
                    else if (m == 2)
                        j--;
                    else if (m == 3)
                        j++;

                    if (i == nx + 1)
                        i = nx - 1;

                    if (i < 0 || j < 0 || i > nx || j > ny)
                        break;

                    if (B[i][j])
                    {
                        double dV = Vmc[i][j] + g;
                        sum1 += dV;
                        sum2 += dV * dV;
                        kchains++;
                        break;
                    }

                    g += delta * delta * rho[i][j] / (4.0 * epsilon);
                }
            }

            if (kchains > 0)
            {
                double V1 = sum1 / kchains;
                double V2 = sum2 / kchains;
                Vmc[i0][j0] = V1;
                sigmaV[i0][j0] = sqrt((V2 - V1 * V1) / kchains);
                absorbed[i0][j0] = double(kchains) / Nchains;
            }

            if (block_after_estimate)
                B[i0][j0] = 1;
        }
    }
}

int main()
{
    auto rho = init_rho();

    // Nadrelaksacja
    vector<vector<double>> Vrel(nx + 1, vector<double>(ny + 1, 0.0));
    relaxacja(Vrel, rho);
    write_to_file("Vrel.txt", Vrel);

    // Metoda MC (można zmieniać parametry)
    vector<vector<double>> Vmc, sigmaV, absorbed;
    mc_method(Vmc, sigmaV, absorbed, rho, 100, 100, false);
    write_to_file("Vmc_N=100_B=0.txt", Vmc);
    write_to_file("SigmaV_N=100_B=0.txt", sigmaV);
    write_to_file("Absorbed_N=100_B=0.txt", absorbed);

    mc_method(Vmc, sigmaV, absorbed, rho, 100, 100, true);
    write_to_file("Vmc_N=100_B=1.txt", Vmc);
    write_to_file("SigmaV_N=100_B=1.txt", sigmaV);
    write_to_file("Absorbed_N=100_B=1.txt", absorbed);

    mc_method(Vmc, sigmaV, absorbed, rho, 300, 300, false);
    write_to_file("Vmc_N=300_B=0.txt", Vmc);
    write_to_file("SigmaV_N=300_B=0.txt", sigmaV);
    write_to_file("Absorbed_N=300_B=0.txt", absorbed);

    mc_method(Vmc, sigmaV, absorbed, rho, 300, 300, true);
    write_to_file("Vmc_N=300_B=1.txt", Vmc);
    write_to_file("SigmaV_N=300_B=1.txt", sigmaV);
    write_to_file("Absorbed_N=300_B=1.txt", absorbed);

    cout << "Obliczenia zakonczone. Wyniki zapisane do plikow." << endl;
    return 0;
}
