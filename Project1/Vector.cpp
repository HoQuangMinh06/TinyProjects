#include "Vector.h"
#include <cassert>
#include <cmath>

Vector::Vector(int size) : mSize(size) {
    assert(size > 0);
    mData = new double[mSize]();
}

Vector::Vector(const Vector& other) : mSize(other.mSize) {
    mData = new double[mSize];
    for (int i = 0; i < mSize; i++)
        mData[i] = other.mData[i];
}

Vector::~Vector() {
    delete[] mData;
}

int Vector::GetSize() const {
    return mSize;
}

Vector& Vector::operator=(const Vector& other) {
    if (this != &other) {
        assert(mSize == other.mSize);
        for (int i = 0; i < mSize; i++)
            mData[i] = other.mData[i];
    }
    return *this;
}

Vector Vector::operator-() const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++)
        result.mData[i] = -mData[i];
    return result;
}

Vector Vector::operator+(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++)
        result.mData[i] = mData[i] + other.mData[i];
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++)
        result.mData[i] = mData[i] - other.mData[i];
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++)
        result.mData[i] = mData[i] * scalar;
    return result;
}

Vector operator*(double scalar, const Vector& v) {
    return v * scalar;
}

double& Vector::operator[](int index) {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

double Vector::operator[](int index) const {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

double& Vector::operator()(int index) {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

double Vector::operator()(int index) const {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}
