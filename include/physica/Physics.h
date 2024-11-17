#ifndef PHYSICA_PHYSICS_H
#define PHYSICA_PHYSICS_H

#include "Vector2D.h"
#include <cmath>

namespace physica {

    class Physics {
    public:
        // Calculate force using Newton's second law: F = m * a
        static Vector2D calculateForce(float mass, const Vector2D& acceleration) {
            return acceleration * mass;
        }

        // Calculate acceleration: a = F / m
        static Vector2D calculateAcceleration(float mass, const Vector2D& force) {
            return force * (1.0f / mass);
        }

        // Calculate velocity from initial velocity and acceleration over time
        static Vector2D calculateVelocity(const Vector2D& initialVelocity, const Vector2D& acceleration, float time) {
            return initialVelocity + acceleration * time;
        }

        // Calculate displacement from initial position, velocity, and time (using constant acceleration)
        static Vector2D calculateDisplacement(const Vector2D& initialPosition, const Vector2D& initialVelocity, const Vector2D& acceleration, float time) {
            return initialPosition + initialVelocity * time + acceleration * (0.5f * time * time);
        }

        // Calculate kinetic energy (KE = 0.5 * m * v^2)
        static float calculateKineticEnergy(float mass, const Vector2D& velocity) {
            return 0.5f * mass * velocity.magnitude() * velocity.magnitude();
        }

        // Calculate potential energy (PE = m * g * h)
        static float calculatePotentialEnergy(float mass, float height, float gravity = 9.81f) {
            return mass * gravity * height;
        }

        // Calculate work done by a force (W = F * d * cos(θ))
        static float calculateWork(const Vector2D& force, const Vector2D& displacement) {
            return force.dot(displacement);
        }

        // Calculate gravitational force between two masses (F = G * (m1 * m2) / r^2)
        static float calculateGravitationalForce(float mass1, float mass2, float distance, float G = 6.67430e-11f) {
            return G * (mass1 * mass2) / (distance * distance);
        }

        // Calculate momentum (p = m * v)
        static Vector2D calculateMomentum(float mass, const Vector2D& velocity) {
            return velocity * mass;
        }

        // Calculate impulse (J = F * Δt)
        static Vector2D calculateImpulse(const Vector2D& force, float deltaTime) {
            return force * deltaTime;
        }

        // Calculate the acceleration due to gravity (g = G * m / r^2)
        static float calculateGravity(float mass, float distanceFromEarthCenter, float G = 6.67430e-11f) {
            return (G * mass) / (distanceFromEarthCenter * distanceFromEarthCenter);
        }

        // Check if two objects collide based on their positions and radii
        static bool checkCollision(const Vector2D& position1, float radius1, const Vector2D& position2, float radius2) {
            float distance = (position1 - position2).magnitude();
            return distance < (radius1 + radius2);
        }

        // Calculate elastic collision velocity (v1' = ((m1 - m2) * v1 + 2 * m2 * v2) / (m1 + m2))
        static Vector2D calculateElasticCollisionVelocity(float m1, const Vector2D& v1, float m2, const Vector2D& v2) {
            return ((v1 * (m1 - m2)) + (v2 * 2.0f * m2)) / (m1 + m2);
        }

        // Calculate the centripetal force (Fc = m * v^2 / r)
        static float calculateCentripetalForce(float mass, const Vector2D& velocity, float radius) {
            return mass * velocity.magnitude() * velocity.magnitude() / radius;
        }
        
        // Calculate the pressure exerted by a fluid (P = F / A)
        static float calculatePressure(float force, float area) {
            return force / area;
        }

        // Calculate torque (τ = r × F = r * F * sin(θ))
        static float calculateTorque(const Vector2D& force, const Vector2D& leverArm) {
            return leverArm.cross(force);
        }

        // Calculate angular velocity (ω = v / r)
        static float calculateAngularVelocity(const Vector2D& velocity, float radius) {
            return velocity.magnitude() / radius;
        }

        // Calculate angular acceleration (α = τ / I)
        static float calculateAngularAcceleration(float torque, float momentOfInertia) {
            return torque / momentOfInertia;
        }

        // Calculate moment of inertia for a solid sphere (I = 2/5 * m * r^2)
        static float calculateMomentOfInertiaSphere(float mass, float radius) {
            return (2.0f / 5.0f) * mass * radius * radius;
        }

        // Calculate moment of inertia for a solid cylinder (I = 1/2 * m * r^2)
        static float calculateMomentOfInertiaCylinder(float mass, float radius) {
            return 0.5f * mass * radius * radius;
        }

        // Calculate drag force (F_d = 1/2 * C_d * ρ * A * v^2)
        static float calculateDragForce(float dragCoefficient, float airDensity, float crossSectionalArea, const Vector2D& velocity) {
            return 0.5f * dragCoefficient * airDensity * crossSectionalArea * velocity.magnitude() * velocity.magnitude();
        }

        // Calculate buoyant force (F_b = ρ * V * g)
        static float calculateBuoyantForce(float fluidDensity, float volume, float gravity = 9.81f) {
            return fluidDensity * volume * gravity;
        }

        // Calculate angular momentum (L = I * ω)
        static float calculateAngularMomentum(float momentOfInertia, float angularVelocity) {
            return momentOfInertia * angularVelocity;
        }

        // Calculate harmonic oscillator displacement (x = A * cos(ωt + φ))
        static float calculateHarmonicOscillatorDisplacement(float amplitude, float angularFrequency, float time, float phase = 0.0f) {
            return amplitude * std::cos(angularFrequency * time + phase);
        }

        // Calculate spring force (F = -k * x)
        static float calculateSpringForce(float springConstant, float displacement) {
            return -springConstant * displacement;
        }

        // Calculate wave speed (v = f * λ)
        static float calculateWaveSpeed(float frequency, float wavelength) {
            return frequency * wavelength;
        }

        // Calculate the Reynolds number (Re = ρ * v * L / μ)
        static float calculateReynoldsNumber(float fluidDensity, float velocity, float characteristicLength, float dynamicViscosity) {
            return (fluidDensity * velocity * characteristicLength) / dynamicViscosity;
        }

        // Calculate heat transfer (Q = m * c * ΔT)
        static float calculateHeatTransfer(float mass, float specificHeat, float temperatureChange) {
            return mass * specificHeat * temperatureChange;
        }

        // Calculate efficiency (Efficiency = (Useful Energy Output / Energy Input) * 100)
        static float calculateEfficiency(float usefulEnergyOutput, float energyInput) {
            return (usefulEnergyOutput / energyInput) * 100.0f;
        }

        // Calculate Doppler Effect for sound (f' = f * (v + vr) / (v - vs))
        static float calculateDopplerEffect(float sourceFrequency, float speedOfSound, float velocityObserver, float velocitySource) {
            return sourceFrequency * (speedOfSound + velocityObserver) / (speedOfSound - velocitySource);
        }

        // Calculate centripetal acceleration (a = v^2 / r)
        static float calculateCentripetalAcceleration(const Vector2D& velocity, float radius) {
            return velocity.magnitude() * velocity.magnitude() / radius;
        }
        
        // Calculate escape velocity (v = sqrt(2 * G * M / r))
        static float calculateEscapeVelocity(float mass, float radius, float G = 6.67430e-11f) {
            return std::sqrt(2.0f * G * mass / radius);
        }

        // Calculate mechanical advantage of a lever (MA = load arm / effort arm)
        static float calculateMechanicalAdvantage(float loadArm, float effortArm) {
            return loadArm / effortArm;
        }
    };
}

#endif // PHYSICA_PHYSICS_H
