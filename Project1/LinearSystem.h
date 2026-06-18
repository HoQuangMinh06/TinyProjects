#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Matrix.h"
#include "Vector.h"

class LinearSystem {
protected:
    int mSize;
    Matrix* mpA;
    Vector* mpb;

private:
    LinearSystem(const LinearSystem& other);

public:
    LinearSystem(const Matrix& A, const Vector& b);
    ~LinearSystem();

    virtual Vector Solve();
};

class PosSymLinSystem : public LinearSystem {
public:
    PosSymLinSystem(const Matrix& A, const Vector& b);

    Vector Solve();
};

class OverDeterminedSystem {
private:
    Matrix* mpA;
    Vector* mpb;

    OverDeterminedSystem(const OverDeterminedSystem& other);

public:
    OverDeterminedSystem(const Matrix& A, const Vector& b);
    ~OverDeterminedSystem();

    Vector SolvePseudoInverse();
    Vector SolveTikhonov(double lambda);
};

#endif
