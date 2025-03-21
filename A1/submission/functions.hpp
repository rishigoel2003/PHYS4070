#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <iostream>
#include <array>
#include <cmath>
#define M_PI 3.14159265358979323846
#include <fstream>

//including the particle class as we use it to find the acceleration and RK4 steps etc
#include "particle_class.hpp"


//This file has the functions we use for the gravitational acceleration between 2 bodies as well as the RK4 steps for 2 bodies
// and then the same for the 3 body systems which we use in B2



// Gravitational acceleration of particle 1 due to particle 2
Vector a_gravity(const Particle &p1, const Particle &p2) {
    Vector r12 = p2.r - p1.r; // Displacement vector
    double r12_mag = sqrt(dot(r12, r12)); // Distance magnitude

    if (r12_mag == 0) return {0, 0}; // Avoid division by zero (because we dont stop the objects if they collide)

    return (p2.mass / pow(r12_mag,3)) * r12; // Newton's law of gravity (we set G=1)
}



// RK4 step for motion
void rk4_step(Particle &p1, const Particle &p2, double dt) {
    // Compute k1 (eulers method to estimate derivative at the start point)
    Vector k1_r = p1.v;
    Vector k1_v = a_gravity(p1, p2);

    // Compute k2 (derivatives at the point calculated after using k1 approximate derivative)
    Vector k2_r = p1.v + (dt / 2) * k1_v;
    Vector k2_v = a_gravity({p1.r + (dt / 2) * k1_r, p1.v, p1.mass}, p2);

    // Compute k3 (approximate derivative using k2)
    Vector k3_r = p1.v + (dt / 2) * k2_v;
    Vector k3_v = a_gravity({p1.r + (dt / 2) * k2_r, p1.v, p1.mass}, p2);

    // Compute k4 (using k3 and double step size)
    Vector k4_r = p1.v + dt * k3_v;
    Vector k4_v = a_gravity({p1.r + dt * k3_r, p1.v, p1.mass}, p2);

    // Update particle 1s position and velocity via the weighted average to get agreement with error dt^4 (quartic agreement)
    p1.r = p1.r + (dt / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    p1.v = p1.v + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
}




// Gravitational acceleration of particle 1 due to particle 2 and particle 3, note this returns a vector of the (x,y) acceleration
Vector a_gravity_3body(const Particle &p1, const Particle &p2, const Particle &p3) {
    Vector r12 = p2.r - p1.r; // Displacement vector
    Vector r13 = p3.r - p1.r; // same for the 1 and 3 particle

    double r12_mag = sqrt(dot(r12, r12)); // Distance magnitude
    double r13_mag = sqrt(dot(r13, r13)); // Distance magnitude


    if (r12_mag == 0 || r13_mag == 0) return {0, 0}; // Avoid division by zero

    return (p2.mass / pow(r12_mag,3)) * r12 + (p3.mass / pow(r13_mag,3)) * r13; // Newton's law of gravity (ignoring G) for multiple bodies
}




// RK4 step for motion with 3 bodies (basically just using acceleration gravity 3 body and updating simultaneously)
void rk4_step_3body(Particle &p1, Particle &p2, const Particle &p3, double dt) {
    // Compute k1 for particle 1 and 2
    Vector k1_r = p1.v;
    Vector k1_v = a_gravity_3body(p1, p2, p3);

    Vector p2_k1_r = p2.v;
    Vector p2_k1_v = a_gravity_3body(p2, p1, p3);

    // Compute k2 for particle 1 and 2 using the k1 from both previous ones
    Vector k2_r = p1.v + (dt / 2) * k1_v;
    Vector k2_v = a_gravity_3body({p1.r + (dt / 2) * k1_r, p1.v, p1.mass}, {p2.r + (dt / 2) * p2_k1_r, p2.v, p2.mass},p3);

    Vector p2_k2_r = p2.v + (dt / 2) * p2_k1_v;
    Vector p2_k2_v = a_gravity_3body({p2.r + (dt / 2) * p2_k1_r, p2.v, p2.mass},{p1.r + (dt / 2) * k1_r, p1.v, p1.mass},p3);


    // Compute k3 for both using k2
    Vector k3_r = p1.v + (dt / 2) * k2_v;
    Vector k3_v = a_gravity_3body({p1.r + (dt / 2) * k2_r, p1.v, p1.mass}, {p2.r + (dt / 2) * p2_k2_r, p2.v, p2.mass} ,p3);

    Vector p2_k3_r = p2.v + (dt / 2) * p2_k2_v;
    Vector p2_k3_v = a_gravity_3body({p2.r + (dt / 2) * p2_k2_r, p2.v, p2.mass},{p1.r + (dt / 2) * k2_r, p1.v, p1.mass},p3);

    // Compute k4 for both using k3 
    Vector k4_r = p1.v + dt * k3_v;
    Vector k4_v = a_gravity_3body({p1.r + dt * k3_r, p1.v, p1.mass}, {p2.r + dt * p2_k3_r,p1.v,p2.mass},p3);

    Vector p2_k4_r = p2.v + dt * p2_k3_v;
    Vector p2_k4_v = a_gravity_3body({p2.r + dt * p2_k3_r,p2.v,p2.mass},{p1.r + dt * k3_r, p1.v, p1.mass},p3);

    
    // Update particle 1s position and velocity
    p1.r = p1.r + (dt / 6) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    p1.v = p1.v + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);

    
    // Update particle 2s position and velocity
    p2.r = p2.r + (dt / 6) * (p2_k1_r + 2 * p2_k2_r + 2 * p2_k3_r + p2_k4_r);
    p2.v = p2.v + (dt / 6) * (p2_k1_v + 2 * p2_k2_v + 2 * p2_k3_v + p2_k4_v);
}



#endif // FUNCTIONS_HPP
