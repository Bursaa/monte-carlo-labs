#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

const double RA = 1;
const double RB = sqrt(2) * RA;
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

double calculate_2D_BM()
{
    double u1 = rnd();
    double u2 = rnd();
    return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}

double normalize_XY(double x, double y)
{
    double norm = sqrt(x * x + y * y);
    return x / norm;
}

double isInCircle(vector<double> r, double R, vector<double> XY)
{
    double x = XY[0];
    double y = XY[1];
    return sqrt((x - r[0]) * (x - r[0]) + (y - r[1]) * (y - r[1])) <= R;
}

vector<vector<double>> generate_N_uniform_circle_points(vector<double> r, double R, double N, string filename = "")
{
    vector<vector<double>> XY_points(N, vector<double>(2));
    ofstream file("data/" + filename);
    for (int i = 0; i < N; i++)
    {

        double X_norm = calculate_2D_BM();
        double Y_norm = calculate_2D_BM();

        double X_jedn = normalize_XY(X_norm, Y_norm);
        double Y_jedn = normalize_XY(Y_norm, X_norm);
        double u3 = rnd();
        XY_points[i][0] = sqrt(u3) * R * X_jedn + r[0];
        XY_points[i][1] = sqrt(u3) * R * Y_jedn + r[1];

        if (file && filename != "")
        {
            file << XY_points[i][0] << " " << XY_points[i][1] << endl;
        }
    }
    file.close();
    return XY_points;
};

int main()
{
    int N = pow(10, 4);
    string alpha = "A";
    vector<vector<double>> XY_points_KA(N, vector<double>(2));
    vector<vector<double>> XY_points_KB(N, vector<double>(2));
    vector<double> rA_beg = {RA + RB, 0};
    vector<double> rB_beg = {0, 0};
    XY_points_KA = generate_N_uniform_circle_points(rA_beg, RA, N, "KA_points.txt");
    XY_points_KB = generate_N_uniform_circle_points(rB_beg, RB, N, "KB_points.txt");

    N = pow(10, 6);
    vector<double> rA = {0, 0};
    vector<double> rB = {0, 0};
    string rAString = rA[0] == 0 ? "=" : "!=";
    string filename = "mu_points_" + alpha + "_xA" + rAString + "0.txt";
    ofstream file("data/" + filename);
    cout << filename << endl;
    for (int k = 2; k <= 6; k++)
    {
        int n = pow(10, k);
        XY_points_KA = generate_N_uniform_circle_points(rA, RA, N);
        XY_points_KB = generate_N_uniform_circle_points(rB, RB, N);

        int sum = 0.0;
        for (int i = 0; i < n; i++)
        {
            int j = randomInt(N);
            vector<double> point = alpha == "A" ? XY_points_KA[j] : XY_points_KB[j];
            sum += alpha == "A" ? isInCircle(rB, RB, point) : isInCircle(rA, RA, point);
        }
        double mu_1 = 0;
        double mu_2 = 0;
        if (alpha == "A")
        {
            mu_1 = M_PI * pow(RA, 2) * sum / n;
            mu_2 = M_PI * pow(RA, 2) * mu_1;
        }
        else
        {
            mu_1 = M_PI * pow(RB, 2) * sum / n;
            mu_2 = M_PI * pow(RB, 2) * mu_1;
        }
        double sigma = sqrt((mu_2 - mu_1 * mu_1) / n);
        file << mu_1 << " " << sigma << "\n";
    }
    file.close();
    return 0;
}