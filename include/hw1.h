#ifndef AP_HW1_H
#define AP_HW1_H
#include <vector>
#include <stdint.h>
using Matrix = std::vector<std::vector<double>>;
namespace algebra {
Matrix zeros(std::size_t n, std::size_t m);
Matrix ones(std::size_t n, std::size_t m);
Matrix random(std::size_t n, std::size_t m, double min, double max);
void show(const Matrix& matrix);
Matrix multiply(const Matrix& matrix, double c);
Matrix multiply(const Matrix& matrix1, const Matrix& matrix2);
Matrix sum(const Matrix& matrix, double c);
Matrix sum(const Matrix& matrix1, const Matrix& matrix2);
Matrix transpose(const Matrix& matrix);
Matrix minor(const Matrix& matrix, std::size_t n, std::size_t m);
double determinant(const Matrix& matrix);
Matrix inverse(const Matrix& matrix);
Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis=0);
Matrix ero_swap(const Matrix& matrix, std::size_t r1, std::size_t r2);
Matrix ero_multiply(const Matrix& matrix, std::size_t r, double c); 

}

#endif //AP_HW1_H
