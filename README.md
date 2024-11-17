# Physica - A C++ Physics and Vector Library

`Physica` is a lightweight C++ library designed for physics simulations and vector calculations. It provides a set of useful functions for vector arithmetic, physical formulas, and utility methods to assist with 2D physics simulations.

## Features

- **Vector2D Operations**: Includes standard vector operations like addition, subtraction, scalar multiplication, and dot product.
- **Physics Calculations**: Provides functions for calculating force, velocity, displacement, kinetic energy, potential energy, and more.
- **Collision Detection**: Simple collision detection using bounding circles.
- **Elastic Collision**: Calculate the new velocity of two objects after an elastic collision.
- **Gravitational Force**: Calculate the gravitational force between two objects based on their masses and distance.
- **Kinetic and Potential Energy**: Methods for calculating kinetic and potential energy in physics simulations.
- **Impulse and Momentum**: Functions for calculating impulse and momentum based on force and velocity.
- **Centripetal Force**: Calculates the centripetal force required for circular motion.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Functions](#functions)
  - [Vector2D Class](#vector2d-class)
  - [Physics Class](#physics-class)
- [License](#license)

## Installation

### Clone the Repository

To use `Physica` in your C++ projects, simply clone the repository:

```bash
git clone https://github.com/yourusername/physica.git
```

### Build the Library

1. Navigate to the project directory.
2. Compile the files using your preferred C++ compiler (e.g., `g++`).

For example:

```bash
g++ main.cpp physica/Vector2D.cpp physica/Physics.cpp -o testPhysica -std=c++17
```

## Usage

Once you have built the project, you can use the library in your C++ code by including the relevant header files:

```cpp
#include "physica/Vector2D.h"
#include "physica/Physics.h"
```

## Examples

### Vector2D Example

```cpp
#include <iostream>
#include "physica/Vector2D.h"

int main() {
    physica::Vector2D v1(3.0f, 4.0f);
    physica::Vector2D v2(1.0f, 2.0f);

    // Adding two vectors
    physica::Vector2D result = v1 + v2;
    result.print();  // Output: (4, 6)

    return 0;
}
```

### Physics Example

```cpp
#include <iostream>
#include "physica/Vector2D.h"
#include "physica/Physics.h"

int main() {
    physica::Vector2D velocity(5.0f, 3.0f);
    physica::Vector2D acceleration(0.5f, 0.2f);
    float mass = 10.0f;

    // Calculate force
    physica::Vector2D force = physica::Physics::calculateForce(mass, acceleration);
    force.print();  // Output: (5, 2)

    return 0;
}
```

## Functions

### `Vector2D Class`

#### Constructor

- `Vector2D(float x, float y)`: Creates a vector with the specified `x` and `y` values.

#### Methods

- `Vector2D operator+(const Vector2D& other) const`: Adds two vectors.
- `Vector2D operator-(const Vector2D& other) const`: Subtracts two vectors.
- `Vector2D operator*(float scalar) const`: Multiplies the vector by a scalar.
- `Vector2D operator/(float scalar) const`: Divides the vector by a scalar.
- `float dot(const Vector2D& other) const`: Calculates the dot product of two vectors.
- `float magnitude() const`: Returns the magnitude (length) of the vector.
- `Vector2D normalize() const`: Normalizes the vector to a unit vector.
- `void print() const`: Prints the vector in the format `(x, y)`.

### `Physics Class`

#### Functions

- `static Vector2D calculateForce(float mass, const Vector2D& acceleration)`: Calculates the force based on mass and acceleration.
- `static Vector2D calculateDisplacement(const Vector2D& initialPosition, const Vector2D& velocity, const Vector2D& acceleration, float time)`: Calculates the displacement using kinematic equations.
- `static float calculateKineticEnergy(float mass, const Vector2D& velocity)`: Calculates the kinetic energy.
- `static float calculatePotentialEnergy(float mass, float height, float gravity)`: Calculates the potential energy.
- `static float calculateWork(const Vector2D& force, const Vector2D& displacement)`: Calculates the work done based on force and displacement.
- `static float calculateGravitationalForce(float mass1, float mass2, float distance)`: Calculates the gravitational force between two masses.
- `static Vector2D calculateMomentum(float mass, const Vector2D& velocity)`: Calculates the momentum of an object.
- `static Vector2D calculateImpulse(const Vector2D& force, float time)`: Calculates the impulse given force and time.
- `static Vector2D calculateElasticCollisionVelocity(float m1, const Vector2D& v1, float m2, const Vector2D& v2)`: Calculates the velocity after an elastic collision.
- `static float calculateCentripetalForce(float mass, const Vector2D& velocity, float radius)`: Calculates the centripetal force required for circular motion.
- `static bool checkCollision(const Vector2D& position1, float radius1, const Vector2D& position2, float radius2)`: Checks if two objects are colliding.

## License

`Physica` is licensed under the MIT License. See the LICENSE file for more details.
