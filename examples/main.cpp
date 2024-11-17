#include <iostream>
#include "physica/Vector2D.h"
#include "physica/Physics.h"

int main() {
    using namespace physica;

    // Create vectors for testing
    Vector2D velocity(5.0f, 3.0f);
    Vector2D acceleration(0.5f, 0.2f);
    Vector2D initialPosition(0.0f, 0.0f);
    Vector2D force = Physics::calculateForce(10.0f, acceleration);
    Vector2D displacement = Physics::calculateDisplacement(initialPosition, velocity, acceleration, 2.0f);
    float mass = 10.0f;
    float height = 5.0f;
    float gravity = 9.81f;

    // Print initial vectors
    std::cout << "Initial Velocity: ";
    velocity.print();
    std::cout << "\n";

    std::cout << "Acceleration: ";
    acceleration.print();
    std::cout << "\n";

    std::cout << "Force: ";
    force.print();
    std::cout << "\n";

    std::cout << "Displacement after 2 seconds: ";
    displacement.print();
    std::cout << "\n";

    // Test Physics functions
    float kineticEnergy = Physics::calculateKineticEnergy(mass, velocity);
    std::cout << "Kinetic Energy: " << kineticEnergy << " J\n";

    float potentialEnergy = Physics::calculatePotentialEnergy(mass, height, gravity);
    std::cout << "Potential Energy: " << potentialEnergy << " J\n";

    float workDone = Physics::calculateWork(force, displacement);
    std::cout << "Work Done: " << workDone << " J\n";

    // Gravitational force between two masses
    float gravitationalForce = Physics::calculateGravitationalForce(mass, 5.0f, 10.0f);
    std::cout << "Gravitational Force: " << gravitationalForce << " N\n";

    // Test momentum and impulse
    Vector2D momentum = Physics::calculateMomentum(mass, velocity);
    std::cout << "Momentum: ";
    momentum.print();
    std::cout << "\n";

    Vector2D impulse = Physics::calculateImpulse(force, 2.0f);
    std::cout << "Impulse: ";
    impulse.print();
    std::cout << "\n";

    // Test elastic collision velocity
    Vector2D v1(5.0f, 3.0f);
    Vector2D v2(3.0f, 4.0f);
    float m1 = 2.0f;
    float m2 = 3.0f;
    Vector2D elasticCollisionVelocity = Physics::calculateElasticCollisionVelocity(m1, v1, m2, v2);
    std::cout << "Elastic Collision Velocity: ";
    elasticCollisionVelocity.print();
    std::cout << "\n";

    // Test centripetal force
    float radius = 10.0f;
    float centripetalForce = Physics::calculateCentripetalForce(mass, velocity, radius);
    std::cout << "Centripetal Force: " << centripetalForce << " N\n";

    // Test gravity calculation
    float distanceFromEarthCenter = 6371000.0f; // Approximate distance in meters
    float gravityAtSurface = Physics::calculateGravity(mass, distanceFromEarthCenter, gravity);
    std::cout << "Gravity at Surface: " << gravityAtSurface << " m/s^2\n";

    // Test collision detection
    Vector2D position1(0.0f, 0.0f);
    Vector2D position2(1.0f, 1.0f);
    float radius1 = 1.0f;
    float radius2 = 1.0f;
    bool isColliding = Physics::checkCollision(position1, radius1, position2, radius2);
    std::cout << "Collision Detected: " << (isColliding ? "Yes" : "No") << "\n";

    // Calculate new velocity after 2 seconds
    Vector2D newVelocity = Physics::calculateVelocity(velocity, acceleration, 2.0f);
    std::cout << "New Velocity after 2 seconds: ";
    newVelocity.print();
    std::cout << "\n";

    return 0;
}
