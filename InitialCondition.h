#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <algorithm>
#include <cmath>
#include <random>

#include "Body.h"

namespace SimulationDefaults {
    constexpr int bodyCount = 1000;
    constexpr double gasVelocityDampingGamma = 0.99;
    constexpr double darkMatterDensityAlpha = 0.0238;
}

struct SimulationConfig {
    int bodyCount = SimulationDefaults::bodyCount;
    double gravity = 0.05;
    double linearAttractionGravityMultiplier = SimulationDefaults::darkMatterDensityAlpha;
    double softening = 0.1;
    double meanMass = 0.25;
    double massStd = 0.05;
    double minMass = 0.05;
    double maxMass = 2.0;
    double maxInitialRadius = 48.0;
    double tangentVelocityMean = 1.0;
    double tangentVelocityStd = 0.05;
    double radialVelocityMean = 0.0;
    double radialVelocityStd = 0.02;
    double velocityRadiusCap = 50.0;
    bool enableMerging = false;
    double mergeDistance = 0.01;
    int gasParticleCount = 750;
    int gasDampingInterval = 10;
    double gasVelocityDampingGamma = SimulationDefaults::gasVelocityDampingGamma;
    unsigned int seed = 42;
};

class InitialConditionGenerator {
public:
    explicit InitialConditionGenerator(const SimulationConfig& config)
        : config_(config),
          gen_(config.seed),
          radiusDist_(0.0, 1.0),
          thetaDist_(0.0, 2.0 * pi()),
          tangentVelocityDist_(config.tangentVelocityMean, config.tangentVelocityStd),
          radialVelocityDist_(config.radialVelocityMean, config.radialVelocityStd),
          massDist_(config.meanMass, config.massStd) {}

    Body nextBody(bool isGas) {
        double mass = sampleMass();
        Vec2 position = samplePosition();
        Vec2 velocity = sampleVelocity(position);

        return Body(mass, position, velocity, isGas);
    }

private:
    const SimulationConfig& config_;
    std::mt19937 gen_;
    std::uniform_real_distribution<double> radiusDist_;
    std::uniform_real_distribution<double> thetaDist_;
    std::normal_distribution<double> tangentVelocityDist_;
    std::normal_distribution<double> radialVelocityDist_;
    std::normal_distribution<double> massDist_;

    static double pi() {
        return 3.141592653589793;
    }

    double sampleMass() {
        return std::clamp(massDist_(gen_), config_.minMass, config_.maxMass);
    }

    Vec2 samplePosition() {
        double radius = sampleRadius();
        double theta = thetaDist_(gen_);

        return Vec2(radius * std::cos(theta), radius * std::sin(theta), 0.0);
    }

    double sampleRadius() {
        return config_.maxInitialRadius * std::sqrt(radiusDist_(gen_));
    }

    Vec2 sampleVelocity(const Vec2& position) {
        Vec2 radial = unitRadialVector(position);
        Vec2 tangent(-radial.y, radial.x, 0.0);

        double radius = std::sqrt(position.normSquaredXY());
        double circularVelocity = approximateCircularVelocity(radius);
        double tangentVelocity = tangentVelocityDist_(gen_) * circularVelocity;
        double radialVelocity = radialVelocityDist_(gen_) * circularVelocity;

        return tangent * tangentVelocity + radial * radialVelocity;
    }

    double approximateCircularVelocity(double radius) const {
        if (radius <= 1e-9) {
            return 0.0;
        }

        double enclosedMass = expectedEnclosedMass(radius);
        double softenedRadiusSquared = radius * radius + config_.softening * config_.softening;
        double softenedRadius = std::sqrt(softenedRadiusSquared);
        double darkMatterEnclosedMass =
            config_.linearAttractionGravityMultiplier * darkMatterEnclosedMassShape(radius);

        double velocitySquared =
            config_.gravity * enclosedMass * radius * radius
            / (softenedRadiusSquared * softenedRadius)
            + config_.gravity * darkMatterEnclosedMass / radius;

        return std::sqrt(std::max(0.0, velocitySquared));
    }

    double expectedEnclosedMass(double radius) const {
        double cappedRadius = std::clamp(radius, 0.0, config_.maxInitialRadius);
        double enclosedFraction =
            (cappedRadius * cappedRadius) / (config_.maxInitialRadius * config_.maxInitialRadius);
        return config_.bodyCount * config_.meanMass * enclosedFraction;
    }

    static double darkMatterEnclosedMassShape(double radius) {
        double scaleRadius = 30.0;
        double x = radius / scaleRadius;
        return 4.0 * pi() * scaleRadius * scaleRadius * scaleRadius * (x - std::atan(x));
    }

    static Vec2 unitRadialVector(const Vec2& position) {
        double radius = std::sqrt(position.normSquaredXY()) + 1e-9;
        return Vec2(position.x / radius, position.y / radius, 0.0);
    }
};

#endif
