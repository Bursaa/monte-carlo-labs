#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

const int N = 1000000;

double f(double x)
{
    return 4.0 / 5.0 * (1 + x - pow(x, 3));
}

double complex_distribution()
{
    double g1 = 4.0 / 5.0;
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;

    if (u1 <= g1)
    {
        return u2;
    }
    else
    {
        return sqrt(1 - sqrt(1 - u2));
    }
}

double markov_chain(double x_i, double delta)
{
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double x_new = x_i + (2 * u1 - 1) * delta;
    double p_acc = f(x_new) / f(x_i);

    if (x_new >= 0 && x_new <= 1 && u2 <= p_acc)
    {
        return x_new;
    }
    else
    {
        return x_i;
    }
}

double elimination_method_function()
{
    while (true)
    {
        double u1 = (double)rand() / RAND_MAX;
        double u2 = 1.15 * ((double)rand() / RAND_MAX);
        if (u2 <= f(u1))
        {
            return u1;
        }
    }
}

int main()
{
    srand(time(0));

    std::ofstream file("data/sampling_results.txt");
    file << "function,Complex Distribution,Elimination Method,Markov Chain (delta=0.5),Markov Chain (delta=0.05)\n";

    double markov_chain_d1 = (double)rand() / RAND_MAX;
    double markov_chain_d2 = (double)rand() / RAND_MAX;

    for (int i = 0; i < N; i++)
    {
        double function_value = f((double)i / N);
        double complex_value = complex_distribution();
        double elimination_value = elimination_method_function();
        if (i > 0)
        {
            markov_chain_d1 = markov_chain(markov_chain_d1, 0.5);
            markov_chain_d2 = markov_chain(markov_chain_d2, 0.05);
        }
        file << function_value << "," << complex_value << "," << elimination_value << ","
             << markov_chain_d1 << "," << markov_chain_d2 << "\n";
    }

    file.close();
    return 0;
}
