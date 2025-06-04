#ifndef PHOTON_DIFFUSION_H
#define PHOTON_DIFFUSION_H

#include <vector>
#include <string>

using std::vector;

class Layer
{
public:
    double mua; // Współczynnik absorpcji (μa)
    double mus; // Współczynnik rozpraszania (μs)
    double d;   // Grubość warstwy (d)
    double g;   // Anizotropia (g)
    double n;   // Współczynnik załamania (n)

    // Konstruktor domyślny
    Layer()
        : mua(1.0), mus(10.0), d(0.02), g(0.75), n(1.0) {}

    // Konstruktor z parametrami
    Layer(double mua, double mus, double d, double g, double n)
        : mua(mua), mus(mus), d(d), g(g), n(n) {}
};

class PHOTON_DIFFUSION_2D
{
private:
    class BEAM
    {
    public:
        double w;
        double x, y;         // położenia
        double x_new, y_new; // nowe proponowane położenia
        double rx, ry;       // kierunek wiązki
        int layer;           // numer warstwy
        bool alive = true;   // czy wiązka jeszcze żyje
        int length;          // flaga
        vector<double> path; // rejestr ścieżki
    };

public:
    // zmienne konfiguracyjne
    int nlayers;
    int nx, ny;
    double xmax, ymax, dx, dy;
    double x_source, dx_source;
    double x_detect, dx_detect;
    double p_min;
    double w_min;
    double rx0, ry0;

    int write_all_paths;
    int write_source_detection_paths;

    // wiązka
    BEAM beam;

    double abs_specular;

    // tablice
    vector<vector<double>> absorption;
    vector<double> reflectance;
    vector<double> transmittance;
    vector<vector<double>> layers_data;

    // funkcje
    PHOTON_DIFFUSION_2D();
    void init();
    double uniform();
    void single_path();
    void roulette();
    void calculate_new_position();
    void scatter_in_layer();
    void scatter_up_down_boundary();
    double sign(double);
    void segment_intersection(double x1, double y1, double x2, double y2,
                              double x3, double y3, double x4, double y4,
                              double &x_cross, double &y_cross, int &icross);
    void write_paths_to_file();
};

#endif // PHOTON_DIFFUSION_H
