#include <vector>
#include <cassert>
#include <iostream>
#include <numeric>
template <class T>
std::vector<T> operator-(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    assert(vec1.size() == vec2.size());
    const size_t size = vec1.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; i++)
    {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
};

template <class T>
std::vector<T> operator+(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    assert(vec1.size() == vec2.size());
    const size_t size = vec1.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; i++)
    {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
};

template <class T>
std::vector<T> operator*(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    assert(vec1.size() == vec2.size());
    const size_t size = vec1.size();
    std::vector<T> result(size);
    for (size_t i = 0; i < size; i++)
    {
        result[i] = vec1[i] * vec2[i];
    }
    return result;
};

double dot_product(const std::vector<double> &a, const std::vector<double> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

void write_polygons_from_atoms(double rmax, int n, std::vector<std::vector<double>> atom, const std::string plik);