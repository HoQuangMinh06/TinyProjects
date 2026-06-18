#ifndef VECTOR_H
#define VECTOR_H

class Vector {
private:
    int mSize;
    double* mData;

public:
    explicit Vector(int size);
    Vector(const Vector& other);
    ~Vector();

    int GetSize() const;

    Vector& operator=(const Vector& other);
    Vector  operator-() const;
    Vector  operator+(const Vector& other) const;
    Vector  operator-(const Vector& other) const;
    Vector  operator*(double scalar) const;
    friend Vector operator*(double scalar, const Vector& v);

    double& operator[](int index);
    double  operator[](int index) const;
    double& operator()(int index);
    double  operator()(int index) const;
};

#endif
