#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <algorithm>
#include <cmath>
#include <random>

#include "Body.h"

enum class VelocityInitializationMode {
    PositionIndependentGaussian,
    CircularReferenceGaussian
};

namespace SimulationDefaults {
    constexpr int bodyCount = 1000;
    constexpr double gasVelocityDampingGamma = 0.99;
    constexpr double alphaR48DarkMatter90 = 0.01128;
    constexpr double alphaR120DarkMatter90 = 0.00248;
    constexpr double maxInitialRadius = 48.0;
}

struct SimulationConfig {
    int bodyCount = SimulationDefaults::bodyCount;
    double gravity = 0.05;
    double linearAttractionGravityMultiplier = SimulationDefaults::alphaR48DarkMatter90;
    double softening = 0.1;
    double meanMass = 0.25;
    double massStd = 0.05;
    double minMass = 0.05;
    double maxMass = 2.0;
    double maxInitialRadius = SimulationDefaults::maxInitialRadius;
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
    int configurationId = 1;
    VelocityInitializationMode velocityMode =
        VelocityInitializationMode::PositionIndependentGaussian;
};

inline int normalizedConfigurationId(int id) {
    return id >= 1 && id <= 4 ? id : 1;
}

inline const char* configurationName(int id) {
    switch (normalizedConfigurationId(id)) {
        case 1:
            return "R=48, position-independent Gaussian velocity";
        case 2:
            return "R=48, circular-reference Gaussian multiplier";
        case 3:
            return "R=120, position-independent Gaussian velocity";
        case 4:
            return "R=120, circular-reference Gaussian multiplier";
        default:
            return "R=48, position-independent Gaussian velocity";
    }
}

inline const char* velocityDistributionName(VelocityInitializationMode mode) {
    switch (mode) {
        case VelocityInitializationMode::CircularReferenceGaussian:
            return "normal_circular_reference_multiplier";
        case VelocityInitializationMode::PositionIndependentGaussian:
        default:
            return "normal_position_independent_velocity";
    }
}

inline const char* velocityRDependencyText(VelocityInitializationMode mode) {
    switch (mode) {
        case VelocityInitializationMode::CircularReferenceGaussian:
            return "multiplied_by_circular_reference_velocity";
        case VelocityInitializationMode::PositionIndependentGaussian:
        default:
            return "none";
    }
}

inline const char* velocityFolderToken(VelocityInitializationMode mode) {
    switch (mode) {
        case VelocityInitializationMode::CircularReferenceGaussian:
            return "velN_vcRef";
        case VelocityInitializationMode::PositionIndependentGaussian:
        default:
            return "velN_Rind";
    }
}

inline const char* logRootDirectoryForConfiguration(int id) {
    switch (normalizedConfigurationId(id)) {
        case 1:
            return "spiralSimulLogs1";
        case 2:
            return "spiralSimulLogs2";
        case 3:
            return "spiralSimulLogs3";
        case 4:
            return "spiralSimulLogs4";
        default:
            return "spiralSimulLogs1";
    }
}

inline double viewScaleForConfiguration(int id) {
    switch (normalizedConfigurationId(id)) {
        case 3:
        case 4:
            return 260.0;
        case 1:
        case 2:
        default:
            return 60.0;
    }
}

inline void applyConfiguration(SimulationConfig& config, int id) {
    config.configurationId = normalizedConfigurationId(id);

    switch (config.configurationId) {
        case 1:
            config.maxInitialRadius = 48.0;
            config.linearAttractionGravityMultiplier =
                SimulationDefaults::alphaR48DarkMatter90;
            config.velocityMode =
                VelocityInitializationMode::PositionIndependentGaussian;
            break;
        case 2:
            config.maxInitialRadius = 48.0;
            config.linearAttractionGravityMultiplier =
                SimulationDefaults::alphaR48DarkMatter90;
            config.velocityMode =
                VelocityInitializationMode::CircularReferenceGaussian;
            break;
        case 3:
            config.maxInitialRadius = 120.0;
            config.linearAttractionGravityMultiplier =
                SimulationDefaults::alphaR120DarkMatter90;
            config.velocityMode =
                VelocityInitializationMode::PositionIndependentGaussian;
            break;
        case 4:
            config.maxInitialRadius = 120.0;
            config.linearAttractionGravityMultiplier =
                SimulationDefaults::alphaR120DarkMatter90;
            config.velocityMode =
                VelocityInitializationMode::CircularReferenceGaussian;
            break;
    }
}

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

        double tangentVelocity = tangentVelocityDist_(gen_);
        double radialVelocity = radialVelocityDist_(gen_);

        if (config_.velocityMode == VelocityInitializationMode::CircularReferenceGaussian) {
            double radius = std::sqrt(position.normSquaredXY());
            double circularVelocity = approximateCircularVelocity(radius);
            tangentVelocity *= circularVelocity;
            radialVelocity *= circularVelocity;
        }

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
