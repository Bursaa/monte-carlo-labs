#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include "dsmc_2d.h"

using namespace std;
void merge_files(const string &file1, const string &file2, const string &output)
{
    ifstream in1(file1), in2(file2);
    ofstream out(output);

    if (!in1 || !in2 || !out)
    {
        cerr << "Błąd otwarcia plików.\n";
        return;
    }

    out << in1.rdbuf();
    out << in2.rdbuf();

    in1.close();
    in2.close();
    out.close();
}
int main()
{
    for (int i = 1; i <= 5; i++)
    {
        cout << "Przetwarzanie danych z pliku dane" << i << ".dat" << endl;
        DSMC_2D ob;
        string data_file = "data/dane" + to_string(i) + ".dat";
        string output_catalog = "data/wyniki" + to_string(i);
        ob.read(data_file.c_str()); // wczytujemy dane zpliku wejściowego
        ob.init();
        ob.nthreads = 8;                                 // obliczenia na 8 rdzeniach
        ob.icol = 1;                                     // cząstki zderzają się
        ob.evolution(0.0, 2500, output_catalog.c_str()); // wykonujemy 5 tysięcy kroków ( tmax - nieznany )
    }

    DSMC_2D ob;
    ob.read("data/dane_left.dat");
    ob.init();
    ob.write_position_velocity("data/rv_left.dat");

    ob.read("data/dane_right.dat");
    ob.init();
    ob.write_position_velocity("data/rv_right.dat");

    merge_files("data/rv_left.dat", "data/rv_right.dat", "data/rv_all.dat");
    ob.read("data/dane_6.dat"); // wczytujemy dane zpliku wejściowego
    ob.init();
    ob.nthreads = 8; // obliczenia na 8 rdzeniach
    ob.icol = 1;     // cząstki zderzają się
    string output_catalog = "data/wyniki6";
    ob.evolution(0.0, 2500, output_catalog.c_str()); // wykonujemy 5 tysięcy kroków ( tmax - nieznany )
}
