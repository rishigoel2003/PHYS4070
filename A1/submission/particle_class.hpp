#ifndef PARTICLE_CLASS_HPP
#define PARTICLE_CLASS_HPP

#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>

using namespace std;


// This .hpp file has the vector class and the overloads we use in all the questions
// The idea of this file is to define a struct (not sure why its not a class) for the bodies that we consider in this question





using Vector = array<double, 2>; // Alias for a 2D vector

// Defines a 'particle' object: position, velocity, mass
struct Particle {
    Vector r;
    Vector v;
    double mass;
};


// overloading vector addition for our 2d vectors
Vector operator+(const Vector &a, const Vector &b) {
    return {a[0] + b[0], a[1] + b[1]};
}

// Overload for vector subtraction
Vector operator-(const Vector &a, const Vector &b) {
    return {a[0] - b[0], a[1] - b[1]};
}

// Dot product overload, I guess its not really an overload as its actually just a function
double dot(const Vector &a, const Vector &b) {
    return a[0] * b[0] + a[1] * b[1];
}

// overload for scalar multiplication for vectors
Vector operator*(double scalar, const Vector &v) {
    return {scalar * v[0], scalar * v[1]};
}

#endif // PARTICLE_CLASS_HPP
