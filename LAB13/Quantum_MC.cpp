#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "photon_diffusion.h"

using namespace std;
// --- Parametry ogólne ---
const int NX = 100;
const int NY = 100;
const int N = 1e6; // Liczba ścieżek do symulacji

void simulate_case(int case_id, double rx, double ry, vector<Layer> &layers)
{
    PHOTON_DIFFUSION_2D ob;

    ob.xmax = 0.2;
    ob.x_source = 0.1;
    ob.dx_source = 0.0;
    ob.x_detect = 0.15;
    ob.dx_detect = 0.01;
    ob.nx = NX;
    ob.ny = NY;
    ob.rx0 = rx;
    ob.ry0 = ry;
    ob.nlayers = layers.size() - 1;

    for (int i = 1; i <= ob.nlayers; ++i)
    {
        ob.layers_data[i][0] = layers[i].mua; // µa
        ob.layers_data[i][1] = layers[i].mus; // µs
        ob.layers_data[i][2] = layers[i].d;   // d
        ob.layers_data[i][3] = layers[i].g;   // g
        ob.layers_data[i][4] = layers[i].n;   // n
    }

    ob.write_all_paths = 0;
    ob.write_source_detection_paths = 0;

    ob.init();

    for (int k = 0; k < N; ++k)
    {
        ob.single_path();
    }

    string filename;

    // Zapis absorption
    filename = "data/absorption_" + to_string(case_id) + ".dat";
    ofstream f_abs(filename);
    f_abs << scientific << setprecision(6);
    for (int i = 0; i <= ob.nx; ++i)
    {
        for (int j = 0; j <= ob.ny; ++j)
            f_abs << ob.absorption[i][j] << " ";
        f_abs << "\n";
    }
    f_abs.close();

    // Zapis reflectance
    filename = "data/reflectance_" + to_string(case_id) + ".dat";
    ofstream f_ref(filename);
    f_ref << scientific << setprecision(6);
    for (int i = 0; i <= ob.nx; ++i)
        f_ref << ob.reflectance[i] << "\n";
    f_ref.close();

    // Zapis transmittance
    filename = "data/transmittance_" + to_string(case_id) + ".dat";
    ofstream f_trans(filename);
    f_trans << scientific << setprecision(6);
    for (int i = 0; i <= ob.nx; ++i)
        f_trans << ob.transmittance[i] << "\n";
    f_trans.close();

    cout << "Zakonczono symulacje przypadku #" << case_id << endl;
}

int main()
{
    int n_layers = 3; // Liczba warstw
    vector<Layer> layers(n_layers + 1);
    layers[1] = Layer(1.0, 10.0, 0.02, 0.75, 1.3);   // Warstwa 1
    layers[2] = Layer(1.0, 190.0, 0.02, 0.075, 1.5); // Warstwa 2
    layers[3] = Layer(10.0, 90.0, 0.02, 0.95, 1.0);  // Warstwa 3

    // --- ZADANIE 2: Wiązka ukośna, wewnętrzne odbicia ---
    simulate_case(1, 0.8, 0.6, layers);

    layers[2].n = 2.5; // n2
    simulate_case(2, 0.8, 0.6, layers);

    layers[1].n = 1.0; // n1
    layers[2].n = 1.5; // n2
    simulate_case(3, 0.8, 0.6, layers);

    layers[1].n = 1.0;    // n1
    layers[2].n = 1.5;    // n2
    layers[2].mus = 10.0; // mus2
    simulate_case(4, 0.8, 0.6, layers);

    // --- ZADANIE 3: Wiązka prostopadła ---

    layers[1].n = 1.3; // n1
    layers[2].n = 1.0; // n2
    simulate_case(5, 0.0, 1.0, layers);

    layers[1].n = 1.0;     // n1
    layers[2].n = 1.5;     // n2
    layers[2].mua = 10.0;  // µa2
    layers[2].mus = 210.0; // µs2
    simulate_case(6, 0.0, 1.0, layers);

    layers[1].n = 1.0;     // n1
    layers[2].n = 1.5;     // n2
    layers[2].mua = 1.0;   // µa2
    layers[2].mus = 210.0; // µs2
    simulate_case(7, 0.0, 1.0, layers);

    layers[1].n = 1.0;     // n1
    layers[2].n = 1.5;     // n2
    layers[2].mua = 10.0;  // µa2
    layers[2].mus = 210.0; // µs2
    layers[2].g = 0.75;    // g2
    simulate_case(8, 0.0, 1.0, layers);

    return 0;
}
