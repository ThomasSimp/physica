#ifndef PHYSICA_VECTOR2D_H
#define PHYSICA_VECTOR2D_H

#include <iostream>
#include <cmath>

namespace physica {

    class Vector2D {
    public:
        float x, y;

        Vector2D() : x(0), y(0) {}
        Vector2D(float x, float y) : x(x), y(y) {}

        // Vector subtraction (operator-)
        Vector2D operator-(const Vector2D& other) const {
            return Vector2D(x - other.x, y - other.y);
        }

        // Vector scalar multiplication (operator*)
        Vector2D operator*(float scalar) const {
            return Vector2D(x * scalar, y * scalar);
        }

        // Vector scalar division (operator/)
        Vector2D operator/(float scalar) const {
            return Vector2D(x / scalar, y / scalar);
        }

        // Add vectors (operator+)
        Vector2D operator+(const Vector2D& other) const {
            return Vector2D(x + other.x, y + other.y);
        }

        // Dot product of two vectors
        float dot(const Vector2D& other) const {
            return x * other.x + y * other.y;
        }

        // Cross product of two vectors (2D scalar equivalent)
        float cross(const Vector2D& other) const {
            return x * other.y - y * other.x;
        }

        // Magnitude (length) of the vector
        float magnitude() const {
            return std::sqrt(x * x + y * y);
        }

        // Normalize the vector (make it unit length)
        Vector2D normalize() const {
            float mag = magnitude();
            return Vector2D(x / mag, y / mag);
        }

        // Print the vector (x, y)
        void print() const {
            std::cout << "(" << x << ", " << y << ")";
        }
    };
}

#endif // PHYSICA_VECTOR2D_H
