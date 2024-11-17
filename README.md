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
git clone https://github.com/ThomasSimp/physica.git
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

## License

`Physica` is licensed under the MIT License. See the LICENSE file for more details.
