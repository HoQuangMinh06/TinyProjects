#include "Matrix.h"
#include <cassert>
#include <cmath>
#include <algorithm>

void Matrix::Allocate() {
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; i++)
        mData[i] = new double[mNumCols]();
}

void Matrix::Deallocate() {
    for (int i = 0; i < mNumRows; i++)
        delete[] mData[i];
    delete[] mData;
}

Matrix::Matrix(int numRows, int numCols) : mNumRows(numRows), mNumCols(numCols) {
    assert(numRows > 0 && numCols > 0);
    Allocate();
}

Matrix::Matrix(const Matrix& other) : mNumRows(other.mNumRows), mNumCols(other.mNumCols) {
    Allocate();
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            mData[i][j] = other.mData[i][j];
}

Matrix::~Matrix() {
    Deallocate();
}

int Matrix::GetNumRows() const { return mNumRows; }
int Matrix::GetNumCols() const { return mNumCols; }

double& Matrix::operator()(int i, int j) {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i - 1][j - 1];
}

double Matrix::operator()(int i, int j) const {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i - 1][j - 1];
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
        for (int i = 0; i < mNumRows; i++)
            for (int j = 0; j < mNumCols; j++)
                mData[i][j] = other.mData[i][j];
    }
    return *this;
}

Matrix Matrix::operator-() const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            result.mData[i][j] = -mData[i][j];
    return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            result.mData[i][j] = mData[i][j] + other.mData[i][j];
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            result.mData[i][j] = mData[i][j] - other.mData[i][j];
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            result.mData[i][j] = mData[i][j] * scalar;
    return result;
}

Matrix operator*(double scalar, const Matrix& m) {
    return m * scalar;
}

Matrix Matrix::operator*(const Matrix& other) const {
    assert(mNumCols == other.mNumRows);
    Matrix result(mNumRows, other.mNumCols);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < other.mNumCols; j++)
            for (int k = 0; k < mNumCols; k++)
                result.mData[i][j] += mData[i][k] * other.mData[k][j];
    return result;
}

Vector Matrix::operator*(const Vector& v) const {
    assert(mNumCols == v.GetSize());
    Vector result(mNumRows);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            result[i] += mData[i][j] * v[j];
    return result;
}

Matrix Matrix::Transpose() const {
    Matrix result(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; i++)
        for (int j = 0; j < mNumCols; j++)
            result.mData[j][i] = mData[i][j];
    return result;
}

double Matrix::Determinant() const {
    assert(mNumRows == mNumCols);
    int n = mNumRows;
    Matrix temp(*this);
    double det = 1.0;

    for (int k = 0; k < n; k++) {
        int pivotRow = k;
        double maxVal = std::abs(temp.mData[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (std::abs(temp.mData[i][k]) > maxVal) {
                maxVal = std::abs(temp.mData[i][k]);
                pivotRow = i;
            }
        }
        if (maxVal < 1e-14) return 0.0;
        if (pivotRow != k) {
            std::swap(temp.mData[pivotRow], temp.mData[k]);
            det *= -1.0;
        }
        det *= temp.mData[k][k];
        for (int i = k + 1; i < n; i++) {
            double factor = temp.mData[i][k] / temp.mData[k][k];
            for (int j = k; j < n; j++)
                temp.mData[i][j] -= factor * temp.mData[k][j];
        }
    }
    return det;
}

Matrix Matrix::Inverse() const {
    assert(mNumRows == mNumCols);
    int n = mNumRows;

    Matrix aug(n, 2 * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug.mData[i][j] = mData[i][j];
        aug.mData[i][i + n] = 1.0;
    }

    for (int k = 0; k < n; k++) {
        int pivotRow = k;
        double maxVal = std::abs(aug.mData[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (std::abs(aug.mData[i][k]) > maxVal) {
                maxVal = std::abs(aug.mData[i][k]);
                pivotRow = i;
            }
        }
        assert(maxVal > 1e-14);
        std::swap(aug.mData[pivotRow], aug.mData[k]);

        double pivot = aug.mData[k][k];
        for (int j = 0; j < 2 * n; j++)
            aug.mData[k][j] /= pivot;

        for (int i = 0; i < n; i++) {
            if (i == k) continue;
            double factor = aug.mData[i][k];
            for (int j = 0; j < 2 * n; j++)
                aug.mData[i][j] -= factor * aug.mData[k][j];
        }
    }

    Matrix inv(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv.mData[i][j] = aug.mData[i][j + n];
    return inv;
}

Matrix Matrix::PseudoInverse() const {
    Matrix At = Transpose();
    if (mNumRows >= mNumCols) {
        Matrix AtA = At * (*this);
        return AtA.Inverse() * At;
    } else {
        Matrix AAt = (*this) * At;
        return At * AAt.Inverse();
    }
}
