#ifndef GRAVITY_SIMULATOR_H
#define GRAVITY_SIMULATOR_H

#include <cmath>
#include <vector>

#include "Body.h"
#include "InitialCondition.h"

class GravitySimulator {
public:
    explicit GravitySimulator(const SimulationConfig& config)
        : config_(config) {
        generateBodies();
    }

    int bodyCount() const {
        return static_cast<int>(bodies_.size());
    }

    double minMass() const {
        return config_.minMass;
    }

    double maxMass() const {
        return config_.maxMass;
    }

    double gravity() const {
        return config_.gravity;
    }

    double linearAttractionGravityMultiplier() const {
        return config_.linearAttractionGravityMultiplier;
    }

    double softening() const {
        return config_.softening;
    }

    double maxInitialRadius() const {
        return config_.maxInitialRadius;
    }

    bool isMergingEnabled() const {
        return config_.enableMerging;
    }

    double mergeDistance() const {
        return config_.mergeDistance;
    }

    int gasParticleCount() const {
        int count = 0;
        for (const Body& body : bodies_) {
            if (body.isGas()) {
                count++;
            }
        }
        return count;
    }

    int initialGasParticleCount() const {
        return config_.gasParticleCount;
    }

    int gasDampingInterval() const {
        return config_.gasDampingInterval;
    }

    double gasVelocityDampingGamma() const {
        return config_.gasVelocityDampingGamma;
    }

    unsigned int seed() const {
        return config_.seed;
    }

    const std::vector<Body>& bodies() const {
        return bodies_;
    }

    Vec2 centerOfMass() const {
        Vec2 weightedPosition(0.0, 0.0);
        double totalMass = 0.0;

        for (const Body& body : bodies_) {
            weightedPosition += body.position() * body.mass();
            totalMass += body.mass();
        }

        if (totalMass <= 0.0) {
            return Vec2(0.0, 0.0);
        }

        return weightedPosition / totalMass;
    }

    Vec2 coreCenter() const {
        Vec2 center = centerOfMass();
        center = centerOfMassWithinRadius(center, 50.0);
        center = centerOfMassWithinRadius(center, 35.0);
        center = centerOfMassWithinRadius(center, 25.0);
        center = centerOfMassWithinRadius(center, 15.0);
        center = centerOfMassWithinRadius(center, 8.0);
        return center;
    }

    Vec2 angularMomentum() const {
        Vec2 total(0.0, 0.0, 0.0);

        for (const Body& body : bodies_) {
            const Vec2& position = body.position();
            const Vec2& velocity = body.velocity();
            double mass = body.mass();

            total.x += mass * (position.y * velocity.z - position.z * velocity.y);
            total.y += mass * (position.z * velocity.x - position.x * velocity.z);
            total.z += mass * (position.x * velocity.y - position.y * velocity.x);
        }

        return total;
    }

    double angularMomentumZ() const {
        return angularMomentum().z;
    }

    double totalAngularMomentum() const {
        Vec2 value = angularMomentum();
        return std::sqrt(value.normSquared());
    }

    void computeForces() {
        for (Body& body : bodies_) {
            body.resetForce();
        }

        for (int i = 0; i < bodyCount(); i++) {
            for (int j = i + 1; j < bodyCount(); j++) {
                Vec2 r = bodies_[j].position() - bodies_[i].position();
                double distSquared = r.normSquared() + config_.softening * config_.softening;
                double dist = std::sqrt(distSquared);
                Vec2 direction = r / dist;

                double forceMagnitude =
                    config_.gravity * bodies_[i].mass() * bodies_[j].mass() / distSquared;

                Vec2 force = direction * forceMagnitude;
                bodies_[i].addForce(force);
                bodies_[j].addForce(force * (-1.0));
            }
        }

        applyLinearAttraction();
    }

    void step(double dt) {
        for (Body& body : bodies_) {
            body.applyHalfKick(dt);
        }

        for (Body& body : bodies_) {
            body.drift(dt);
        }

        computeForces();

        for (Body& body : bodies_) {
            body.applyHalfKick(dt);
        }
    }

    void applyMergePolicy() {
        if (config_.enableMerging) {
            mergeCloseBodies(config_.mergeDistance);
        }
    }

    void applyGasDamping() {
        for (Body& body : bodies_) {
            if (body.isGas()) {
                body.dampVelocity(config_.gasVelocityDampingGamma);
            }
        }
    }

    void mergeCloseBodies(double mergeDistance) {
        for (int i = 0; i < bodyCount(); i++) {
            for (int j = i + 1; j < bodyCount(); j++) {
                Vec2 r = bodies_[j].position() - bodies_[i].position();

                if (r.normSquared() < mergeDistance * mergeDistance) {
                    double m1 = bodies_[i].mass();
                    double m2 = bodies_[j].mass();
                    double newMass = m1 + m2;

                    Vec2 newPosition =
                        (bodies_[i].position() * m1 + bodies_[j].position() * m2) / newMass;
                    Vec2 newVelocity =
                        (bodies_[i].velocity() * m1 + bodies_[j].velocity() * m2) / newMass;
                    bool isGas = bodies_[i].isGas() || bodies_[j].isGas();

                    bodies_[i].resetAfterMerge(newPosition, newVelocity, newMass);
                    bodies_[i].setGas(isGas);
                    bodies_.erase(bodies_.begin() + j);
                    j--;
                }
            }
        }
    }

private:
    SimulationConfig config_;
    std::vector<Body> bodies_;

    Vec2 centerOfMassWithinRadius(const Vec2& center, double radius) const {
        Vec2 weightedPosition(0.0, 0.0);
        double totalMass = 0.0;
        int selectedCount = 0;
        double radiusSquared = radius * radius;

        for (const Body& body : bodies_) {
            Vec2 relativePosition = body.position() - center;
            if (relativePosition.normSquared() > radiusSquared) {
                continue;
            }

            weightedPosition += body.position() * body.mass();
            totalMass += body.mass();
            selectedCount++;
        }

        if (totalMass <= 0.0 || selectedCount < minimumCoreParticleCount()) {
            return center;
        }

        return weightedPosition / totalMass;
    }

    static int minimumCoreParticleCount() {
        return 20;
    }

    void applyLinearAttraction() {
        double densityCoefficient = config_.linearAttractionGravityMultiplier;
        if (densityCoefficient == 0.0) {
            return;
        }

        for (Body& body : bodies_) {
            double radius = std::sqrt(body.position().normSquared());
            if (radius <= 1e-9) {
                continue;
            }

            double darkMatterEnclosedMass =
                densityCoefficient * darkMatterEnclosedMassShape(radius);
            double forceScale =
                -config_.gravity * body.mass() * darkMatterEnclosedMass
                / (radius * radius * radius);

            body.addForce(body.position() * forceScale);
        }
    }

    static double darkMatterEnclosedMassShape(double radius) {
        double scaleRadius = 30.0;
        double x = radius / scaleRadius;
        return 4.0 * std::acos(-1.0)
            * scaleRadius * scaleRadius * scaleRadius
            * (x - std::atan(x));
    }

    void generateBodies() {
        InitialConditionGenerator generator(config_);

        bodies_.reserve(config_.bodyCount);

        for (int i = 0; i < config_.bodyCount; i++) {
            bodies_.push_back(generator.nextBody(i < config_.gasParticleCount));
        }
    }
};

#endif
