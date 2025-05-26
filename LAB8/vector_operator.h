#include <vector>
#include <cassert>
#include <iostream>
#include <numeric>
template <class T>
std::vector<T> operator-(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    assert(vec1.size() == vec2.size());
    std::vector<T> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); i++)
        result[i] = vec1[i] - vec2[i];
    return result;
}

// Dodawanie dwóch wektorów
template <class T>
std::vector<T> operator+(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    assert(vec1.size() == vec2.size());
    std::vector<T> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); i++)
        result[i] = vec1[i] + vec2[i];
    return result;
}

// Iloczyn elementów dwóch wektorów (nie iloczyn skalarny!)
template <class T>
std::vector<T> operator*(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    assert(vec1.size() == vec2.size());
    std::vector<T> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); i++)
        result[i] = vec1[i] * vec2[i];
    return result;
}

// Mnożenie wektora przez liczbę (z lewej strony)
template <class T>
std::vector<T> operator*(double number, const std::vector<T> &vec)
{
    std::vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); i++)
        result[i] = number * vec[i]; // <-- tu był błąd: `vector[i]` -> `vec[i]`
    return result;
}

// Iloczyn skalarny (dot product)
double dot_product(const std::vector<double> &a, const std::vector<double> &b)
{
    assert(a.size() == b.size());
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

void write_polygons_from_atoms(double rmax, int n, std::vector<std::vector<double>> atom, const std::string plik);