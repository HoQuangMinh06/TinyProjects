#include "LinearSystem.h"
#include <cassert>
#include <cmath>
#include <algorithm>

LinearSystem::LinearSystem(const Matrix& A, const Vector& b) {
    assert(A.GetNumRows() == A.GetNumCols());
    assert(A.GetNumRows() == b.GetSize());
    mSize = b.GetSize();
    mpA = new Matrix(A);
    mpb = new Vector(b);
}

LinearSystem::~LinearSystem() {
    delete mpA;
    delete mpb;
}

Vector LinearSystem::Solve() {
    Matrix A(*mpA);
    Vector b(*mpb);
    int n = mSize;

    for (int k = 1; k <= n; k++) {
        int pivotRow = k;
        double maxVal = std::abs(A(k, k));
        for (int i = k + 1; i <= n; i++) {
            if (std::abs(A(i, k)) > maxVal) {
                maxVal = std::abs(A(i, k));
                pivotRow = i;
            }
        }
        assert(maxVal > 1e-14);


        if (pivotRow != k) {
            for (int j = 1; j <= n; j++) {
                double tmp = A(k, j);
                A(k, j) = A(pivotRow, j);
                A(pivotRow, j) = tmp;
            }
            double tmp = b(k);
            b(k) = b(pivotRow);
            b(pivotRow) = tmp;
        }

        for (int i = k + 1; i <= n; i++) {
            double factor = A(i, k) / A(k, k);
            for (int j = k; j <= n; j++)
                A(i, j) -= factor * A(k, j);
            b(i) -= factor * b(k);
        }
    }

    Vector x(n);
    for (int i = n; i >= 1; i--) {
        double sum = b(i);
        for (int j = i + 1; j <= n; j++)
            sum -= A(i, j) * x(j);
        x(i) = sum / A(i, i);
    }
    return x;
}

OverDeterminedSystem::OverDeterminedSystem(const Matrix& A, const Vector& b) {
    mpA = new Matrix(A);
    mpb = new Vector(b);
}

OverDeterminedSystem::~OverDeterminedSystem() {
    delete mpA;
    delete mpb;
}

Vector OverDeterminedSystem::SolvePseudoInverse() {
    Matrix Aplus = mpA->PseudoInverse();
    return Aplus * (*mpb);
}
Vector OverDeterminedSystem::SolveTikhonov(double lambda) {
    Matrix& A = *mpA;
    Matrix At = A.Transpose();
    Matrix AtA = At * A;
    int n = AtA.GetNumRows();
    for (int i = 1; i <= n; i++)
        AtA(i, i) += lambda;
    return AtA.Inverse() * (At * (*mpb));
}


PosSymLinSystem::PosSymLinSystem(const Matrix& A, const Vector& b)
    : LinearSystem(A, b) {
    int n = A.GetNumRows();
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            assert(std::abs(A(i, j) - A(j, i)) < 1e-10);
}

Vector PosSymLinSystem::Solve() {
    int n = mSize;
    Matrix& A = *mpA;
    Vector& b = *mpb;

    Vector x(n);
    Vector r(b);
    Vector p(r);
    double rsold = 0.0;
    for (int i = 1; i <= n; i++)
        rsold += r(i) * r(i);

    for (int iter = 0; iter < n * 1000; iter++) {
        Vector Ap = A * p;

        double pAp = 0.0;
        for (int i = 1; i <= n; i++)
            pAp += p(i) * Ap(i);

        if (std::abs(pAp) < 1e-14) break;

        double alpha = rsold / pAp;

        for (int i = 1; i <= n; i++) {
            x(i) += alpha * p(i);
            r(i) -= alpha * Ap(i);
        }

        double rsnew = 0.0;
        for (int i = 1; i <= n; i++)
            rsnew += r(i) * r(i);

        if (std::sqrt(rsnew) < 1e-10) break;

        double beta = rsnew / rsold;
        for (int i = 1; i <= n; i++)
            p(i) = r(i) + beta * p(i);

        rsold = rsnew;
    }
    return x;
}
