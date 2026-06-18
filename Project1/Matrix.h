#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"

class Matrix {
private:
    int mNumRows;
    int mNumCols;
    double** mData;

    void Allocate();
    void Deallocate();

public:
    Matrix(int numRows, int numCols);
    Matrix(const Matrix& other);
    ~Matrix();

    int GetNumRows() const;
    int GetNumCols() const;

    double& operator()(int i, int j);
    double  operator()(int i, int j) const;

    Matrix& operator=(const Matrix& other);
    Matrix  operator-() const;
    Matrix  operator+(const Matrix& other) const;
    Matrix  operator-(const Matrix& other) const;
    Matrix  operator*(double scalar) const;
    Matrix  operator*(const Matrix& other) const;
    Vector  operator*(const Vector& v) const;
    friend Matrix operator*(double scalar, const Matrix& m);

    Matrix Transpose() const;
    double Determinant() const;
    Matrix Inverse() const;
    Matrix PseudoInverse() const;
};

#endif
