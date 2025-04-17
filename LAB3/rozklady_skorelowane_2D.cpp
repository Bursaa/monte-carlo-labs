#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

double rnd()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dist(0.0, 1.0); // Poprawiony typ

    return dist(gen);
}

double calculate_2D_BM(double u1, double u2)
{
    return sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
}

double normalize_XY(double x, double y, double u)
{
    double norm = sqrt(x * x + y * y);
    return sqrt(u) * x / norm;
}

std::vector<double> afinical_transform_to_elipse(double x, double y, double b1, double b2, double alpha)
{

    double alpha_rad = alpha * M_PI / 180;
    double x_rotated = b1 * x * cos(alpha_rad) - b2 * y * sin(alpha_rad);
    double y_rotated = b1 * x * sin(alpha_rad) + b2 * y * cos(alpha_rad);

    return {x_rotated, y_rotated};
}
int main()
{
    int N = std::pow(10, 4);
    std::ofstream file1("data/points.txt");
    std::ofstream file2("data/macierz_jedn.txt");
    std::ofstream file3("data/macierz_norm.txt");

    double alpha = 45;
    double b1 = 1;
    double b2 = 0.2;

    vector<vector<double>> XY_elipse_jedn(N, vector<double>(2));
    vector<vector<double>> XY_elipse_norm(N, vector<double>(2));

    for (int i = 0; i < N; i++)
    {
        double u1 = rnd();
        double u2 = rnd();

        double X_norm = calculate_2D_BM(u1, u2);

        u1 = rnd();
        u2 = rnd();
        double Y_norm = calculate_2D_BM(u1, u2);

        u1 = rnd();
        double X_jedn = sqrt(u1) * X_norm / sqrt(X_norm * X_norm + Y_norm * Y_norm);
        double Y_jedn = sqrt(u1) * Y_norm / sqrt(X_norm * X_norm + Y_norm * Y_norm);
        vector<double> xy_elipse_jedn = afinical_transform_to_elipse(X_jedn, Y_jedn, b1, b2, alpha);
        vector<double> xy_elipse_norm = afinical_transform_to_elipse(X_norm, Y_norm, b1, b2, alpha);

        XY_elipse_jedn[i] = xy_elipse_jedn;
        XY_elipse_norm[i] = xy_elipse_norm;

        file1 << X_norm << " " << Y_norm << " " << X_jedn << " " << Y_jedn << " " << xy_elipse_jedn[0] << " " << xy_elipse_jedn[1] << " " << xy_elipse_norm[0] << " " << xy_elipse_norm[1] << "\n";
    }

    double sumx_jedn = 0;
    double sumy_jedn = 0;
    double sumx2_jedn = 0;
    double sumy2_jedn = 0;
    double sumxy_jedn = 0;
    double sumx_norm = 0;
    double sumy_norm = 0;
    double sumx2_norm = 0;
    double sumy2_norm = 0;
    double sumxy_norm = 0;

    for (int i = 0; i < N; i++)
    {
        double x_jedn = XY_elipse_jedn[i][0];
        double y_jedn = XY_elipse_jedn[i][1];
        sumx_jedn += x_jedn;
        sumy_jedn += y_jedn;
        sumx2_jedn += x_jedn * x_jedn;
        sumy2_jedn += y_jedn * y_jedn;
        sumxy_jedn += x_jedn * y_jedn;

        double x_norm = XY_elipse_norm[i][0];
        double y_norm = XY_elipse_norm[i][1];
        sumx_norm += x_norm;
        sumy_norm += y_norm;
        sumx2_norm += x_norm * x_norm;
        sumy2_norm += y_norm * y_norm;
        sumxy_norm += x_norm * y_norm;
    }

    file2 << sumx2_jedn - sumx_jedn * sumx_jedn << " " << sumxy_jedn - sumx_jedn * sumx_jedn << "\n"
          << sumxy_jedn - sumx_jedn * sumx_jedn << " " << sumy2_jedn - sumy_jedn * sumy_jedn << "\n";

    file3 << sumx2_norm - sumx_norm * sumx_norm << " " << sumxy_norm - sumx_norm * sumx_norm << "\n"
          << sumxy_norm - sumx_norm * sumx_norm << " " << sumy2_norm - sumy_norm * sumy_norm << "\n";

    file1.close();
    file2.close();
    file3.close();

    return 0;
}