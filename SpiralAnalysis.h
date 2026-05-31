#ifndef SPIRAL_ANALYSIS_H
#define SPIRAL_ANALYSIS_H

#include <SFML/Graphics.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "DistributionGrid.h"
#include "GravitySimulator.h"

struct FourierSnapshot {
    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double a4 = 0.0;
};

struct FourierModeValue {
    double amplitude = 0.0;
    double phase = 0.0;
    double mass = 0.0;
};

struct RadialA2Bin {
    double rMin = 0.0;
    double rMax = 0.0;
    double a2 = 0.0;
    double phase = 0.0;
    double contrast = 0.0;
    double mass = 0.0;
};

struct LogSpiralTransformResult {
    bool valid = false;
    int mode = 2;
    double score = 0.0;
    double p = 0.0;
    double phase = 0.0;
    double pitchAngleDegrees = 0.0;
    double mass = 0.0;
};

struct SpiralMetrics {
    FourierSnapshot snapshot;
    double phase2 = 0.0;
    double innerA2 = 0.0;
    double innerPhase = 0.0;
    double outerA2 = 0.0;
    double outerPhase = 0.0;
    double armContrast = 0.0;
    LogSpiralTransformResult m2LogSpiral;
    std::vector<RadialA2Bin> radialBins;
};

class SpiralAnalyzer {
public:
    explicit SpiralAnalyzer(const GravitySimulator& simulator)
        : simulator_(simulator) {}

    FourierSnapshot computeSnapshot() const {
        FourierSnapshot snapshot;
        snapshot.a1 = computeFourierMode(1);
        snapshot.a2 = computeFourierMode(2);
        snapshot.a3 = computeFourierMode(3);
        snapshot.a4 = computeFourierMode(4);
        return snapshot;
    }

    SpiralMetrics computeMetrics() const {
        SpiralMetrics metrics;
        metrics.snapshot = computeSnapshot();

        FourierModeValue globalMode2 = computeFourierModeInRange(2, 0.0, analysisRadius());
        FourierModeValue innerMode2 = computeFourierModeInRange(2, 0.0, innerOuterSplitRadius());
        FourierModeValue outerMode2 = computeFourierModeInRange(
            2,
            innerOuterSplitRadius(),
            analysisRadius()
        );

        metrics.phase2 = globalMode2.phase;
        metrics.innerA2 = innerMode2.amplitude;
        metrics.innerPhase = innerMode2.phase;
        metrics.outerA2 = outerMode2.amplitude;
        metrics.outerPhase = outerMode2.phase;
        metrics.radialBins = computeRadialA2Bins();
        metrics.armContrast = summarizeContrast(metrics.radialBins);
        metrics.m2LogSpiral = computeBestLogSpiralTransform(2);

        return metrics;
    }

    double computeFourierMode(int mode) const {
        return computeFourierModeInRange(mode, 0.0, analysisRadius()).amplitude;
    }

    FourierModeValue computeFourierModeInRange(int mode, double rMin, double rMax) const {
        std::complex<double> sum(0.0, 0.0);
        double totalMass = 0.0;
        Vec2 coreCenter = simulator_.coreCenter();

        for (const Body& body : simulator_.bodies()) {
            Vec2 relativePosition = body.position() - coreCenter;
            double radius = std::sqrt(relativePosition.normSquaredXY());
            if (radius < rMin || radius >= rMax) {
                continue;
            }

            double theta = std::atan2(relativePosition.y, relativePosition.x);
            std::complex<double> phase = std::exp(std::complex<double>(0.0, mode * theta));

            sum += body.mass() * phase;
            totalMass += body.mass();
        }

        if (totalMass <= 0.0) {
            return FourierModeValue{};
        }

        FourierModeValue value;
        value.amplitude = std::abs(sum) / totalMass;
        value.phase = std::atan2(sum.imag(), sum.real()) / static_cast<double>(mode);
        value.mass = totalMass;
        return value;
    }

    std::vector<RadialA2Bin> computeRadialA2Bins() const {
        std::vector<RadialA2Bin> bins;
        double width = radialBinWidth();
        int binCount = static_cast<int>(std::ceil(analysisRadius() / width));

        for (int i = 0; i < binCount; i++) {
            RadialA2Bin bin;
            bin.rMin = i * width;
            bin.rMax = std::min(analysisRadius(), (i + 1) * width);

            FourierModeValue mode2 = computeFourierModeInRange(2, bin.rMin, bin.rMax);
            bin.a2 = mode2.amplitude;
            bin.phase = mode2.phase;
            bin.mass = mode2.mass;
            bin.contrast = computeArmContrast(bin.rMin, bin.rMax);

            bins.push_back(bin);
        }

        return bins;
    }

    std::string classifyTwoArmSpiral(const std::vector<double>& a2History) const {
        double mean;
        double sigma;
        double relativeSigma;
        summarize(a2History, mean, sigma, relativeSigma);

        if (mean > 0.15 && relativeSigma < 0.15) {
            return "Stable two-arm spiral structure detected.";
        }

        if (mean > 0.15) {
            return "Transient or fluctuating two-arm spiral structure.";
        }

        return "No strong two-arm spiral structure.";
    }

    void summarize(
        const std::vector<double>& values,
        double& mean,
        double& sigma,
        double& relativeSigma
    ) const {
        mean = 0.0;
        sigma = 0.0;
        relativeSigma = 0.0;

        if (values.empty()) {
            return;
        }

        for (double value : values) {
            mean += value;
        }
        mean /= values.size();

        double variance = 0.0;
        for (double value : values) {
            variance += (value - mean) * (value - mean);
        }
        variance /= values.size();

        sigma = std::sqrt(variance);

        if (mean != 0.0) {
            relativeSigma = sigma / mean;
        }
    }

    static double analysisRadius() {
        return 50.0;
    }

    static double innerOuterSplitRadius() {
        return 20.0;
    }

    static double radialBinWidth() {
        return 10.0;
    }

    static double logSpiralMinRadius() {
        return 2.0;
    }

    static double logSpiralMaxRadius() {
        return 30.0;
    }

private:
    const GravitySimulator& simulator_;

    LogSpiralTransformResult computeBestLogSpiralTransform(int mode) const {
        LogSpiralTransformResult best;
        best.mode = mode;

        const double pMin = -25.0;
        const double pMax = 25.0;
        const double pStep = 0.25;

        for (double p = pMin; p <= pMax + 1e-9; p += pStep) {
            LogSpiralTransformResult candidate = computeLogSpiralTransform(mode, p);
            if (candidate.valid && candidate.score > best.score) {
                best = candidate;
            }
        }

        return best;
    }

    LogSpiralTransformResult computeLogSpiralTransform(int mode, double p) const {
        std::complex<double> sum(0.0, 0.0);
        double totalMass = 0.0;
        Vec2 coreCenter = simulator_.coreCenter();

        for (const Body& body : simulator_.bodies()) {
            Vec2 relativePosition = body.position() - coreCenter;
            double radius = std::sqrt(relativePosition.normSquaredXY());
            if (radius < logSpiralMinRadius() || radius >= logSpiralMaxRadius()) {
                continue;
            }

            double theta = std::atan2(relativePosition.y, relativePosition.x);
            double phase = mode * theta + p * std::log(radius);
            sum += body.mass() * std::exp(std::complex<double>(0.0, -phase));
            totalMass += body.mass();
        }

        if (totalMass <= 0.0) {
            return LogSpiralTransformResult{};
        }

        LogSpiralTransformResult result;
        result.valid = true;
        result.mode = mode;
        result.score = std::abs(sum) / totalMass;
        result.p = p;
        result.phase = std::atan2(sum.imag(), sum.real());
        result.pitchAngleDegrees =
            std::atan2(static_cast<double>(mode), std::max(std::abs(p), 1e-9))
            * 180.0 / std::acos(-1.0);
        result.mass = totalMass;
        return result;
    }

    double computeArmContrast(double rMin, double rMax) const {
        const int thetaBinCount = 36;
        std::vector<double> thetaMass(thetaBinCount, 0.0);
        Vec2 coreCenter = simulator_.coreCenter();

        for (const Body& body : simulator_.bodies()) {
            Vec2 relativePosition = body.position() - coreCenter;
            double radius = std::sqrt(relativePosition.normSquaredXY());
            if (radius < rMin || radius >= rMax) {
                continue;
            }

            double theta = std::atan2(relativePosition.y, relativePosition.x);
            double fullCircle = 2.0 * std::acos(-1.0);
            if (theta < 0.0) {
                theta += fullCircle;
            }

            int bin = static_cast<int>(theta / fullCircle * thetaBinCount);
            bin = std::clamp(bin, 0, thetaBinCount - 1);
            thetaMass[bin] += body.mass();
        }

        std::sort(thetaMass.begin(), thetaMass.end());
        int groupSize = std::max(1, thetaBinCount / 5);
        double low = 0.0;
        double high = 0.0;

        for (int i = 0; i < groupSize; i++) {
            low += thetaMass[i];
            high += thetaMass[thetaBinCount - 1 - i];
        }

        low /= groupSize;
        high /= groupSize;

        if (high <= 0.0 && low <= 0.0) {
            return 0.0;
        }

        return high / std::max(low, 1e-9);
    }

    static double summarizeContrast(const std::vector<RadialA2Bin>& bins) {
        double weightedContrast = 0.0;
        double totalMass = 0.0;

        for (const RadialA2Bin& bin : bins) {
            if (bin.mass <= 0.0) {
                continue;
            }

            weightedContrast += bin.contrast * bin.mass;
            totalMass += bin.mass;
        }

        if (totalMass <= 0.0) {
            return 0.0;
        }

        return weightedContrast / totalMass;
    }
};

class DistributionLogger {
public:
    DistributionLogger(
        const DistributionGrid& grid,
        int bodyCount,
        double tvelMean,
        double tvelStd,
        double rvelMean,
        double rvelStd
    )
        : grid_(grid),
          bodyCount_(bodyCount),
          tvelMean_(tvelMean),
          tvelStd_(tvelStd),
          rvelMean_(rvelMean),
          rvelStd_(rvelStd) {}

    const std::string& path() const {
        return summaryFilePath_;
    }

    const std::string& runDirectory() const {
        return runDirectory_;
    }

    void prepare(
        const GravitySimulator& simulator,
        int steps,
        double dt,
        int ignoreSteps,
        int analyzeInterval
    ) {
        std::filesystem::create_directories(logRootDirectory(simulator));
        runDirectory_ = nextRunDirectory(simulator);
        std::filesystem::create_directories(runDirectory_);

        summaryFilePath_ = runDirectory_ + "/summary.txt";

        std::ofstream fout(summaryFilePath_, std::ios::trunc);
        writeInitialConditionBlock(fout, simulator, steps, dt, ignoreSteps, analyzeInterval);
        fout << '\n';
        fout << "========== Fourier Analysis ==========" << '\n';
        fout << "Status                 = running" << '\n';
        fout << '\n';
        fout << "Snapshots are saved every 200 steps as PNG images." << '\n';
    }

    void writeSnapshot(
        const GravitySimulator& simulator,
        int step,
        const FourierSnapshot& snapshot,
        double dt
    ) const {
        const unsigned int imageWidth = 1000;
        const unsigned int imageHeight = 750;
        sf::Image image({imageWidth, imageHeight}, sf::Color(5, 7, 12));

        for (const Body& body : simulator.bodies()) {
            int x;
            int y;
            if (!projectToImage(body.position(), imageWidth, imageHeight, x, y)) {
                continue;
            }

            double normalizedMass = body.mass() / std::max(simulator.maxMass(), 1e-9);
            int radius = static_cast<int>(1 + 3 * std::sqrt(std::clamp(normalizedMass, 0.0, 1.0)));
            std::uint8_t brightness = static_cast<std::uint8_t>(
                120 + 135 * std::sqrt(std::clamp(normalizedMass, 0.0, 1.0))
            );
            drawDisk(image, x, y, radius, sf::Color(brightness, brightness, 255));
        }

        Vec2 core = simulator.coreCenter();
        int coreX;
        int coreY;
        if (projectToImage(core, imageWidth, imageHeight, coreX, coreY)) {
            drawCross(image, coreX, coreY, 9, sf::Color(255, 70, 70));
        }

        std::ostringstream fileName;
        fileName << runDirectory_ << "/snapshot_step_"
                 << std::setw(6) << std::setfill('0') << step << ".png";
        bool saved = image.saveToFile(fileName.str());

        (void)dt;
        (void)snapshot;
        (void)saved;
    }

    void writeSummary(
        const GravitySimulator& simulator,
        double mean,
        double sigma,
        double relativeSigma,
        double meanM2LogSpiral,
        double sigmaM2LogSpiral,
        double relativeSigmaM2LogSpiral,
        const LogSpiralTransformResult& latestM2LogSpiral,
        int analyzedSamples,
        int steps,
        double dt,
        int ignoreSteps,
        int analyzeInterval,
        double initialAngularMomentum,
        double finalAngularMomentum,
        const std::string& result
    ) const {
        std::ofstream fout(summaryFilePath_, std::ios::trunc);
        writeInitialConditionBlock(fout, simulator, steps, dt, ignoreSteps, analyzeInterval);

        fout << '\n';
        fout << "========== Fourier Analysis ==========" << '\n';
        fout << "Fourier center         = CoreCenter" << '\n';
        fout << "Core radius            = 20" << '\n';
        fout << "Fourier radius         = " << SpiralAnalyzer::analysisRadius() << '\n';
        fout << "Mean A2                = " << mean << '\n';
        fout << "Std(A2)                = " << sigma << '\n';
        fout << "Relative Std(A2)       = " << relativeSigma << '\n';
        fout << "Mean m=2 log score     = " << meanM2LogSpiral << '\n';
        fout << "Std m=2 log score      = " << sigmaM2LogSpiral << '\n';
        fout << "Relative Std m=2 log   = " << relativeSigmaM2LogSpiral << '\n';
        fout << "Latest m=2 log score   = " << latestM2LogSpiral.score << '\n';
        fout << "Latest m=2 p           = " << latestM2LogSpiral.p << '\n';
        fout << "Latest m=2 pitch deg   = " << latestM2LogSpiral.pitchAngleDegrees << '\n';
        fout << "Latest m=2 phase       = " << latestM2LogSpiral.phase << '\n';
        fout << "m=2 log radius range   = "
             << SpiralAnalyzer::logSpiralMinRadius()
             << " to " << SpiralAnalyzer::logSpiralMaxRadius() << '\n';
        fout << "Analyzed samples       = " << analyzedSamples << '\n';
        fout << "Result                 = " << result << '\n';
        fout << '\n';
        fout << "========== Conservation ==========" << '\n';
        fout << "Initial angular Lz     = " << initialAngularMomentum << '\n';
        fout << "Final angular Lz       = " << finalAngularMomentum << '\n';
        fout << "Relative Lz drift      = "
             << relativeAngularMomentumDrift(initialAngularMomentum, finalAngularMomentum) << '\n';
        fout << '\n';
        fout << "========== Output ==========" << '\n';
        fout << "Snapshot interval      = 200 steps" << '\n';
        fout << "Snapshot format        = PNG image" << '\n';

        appendAnalysisIndex(
            simulator,
            mean,
            sigma,
            meanM2LogSpiral,
            sigmaM2LogSpiral,
            result
        );
    }

private:
    DistributionGrid grid_;
    int bodyCount_;
    double tvelMean_;
    double tvelStd_;
    double rvelMean_;
    double rvelStd_;
    std::string runDirectory_;
    std::string summaryFilePath_;

    std::string initialConditionFolderName(const GravitySimulator& simulator) const {
        std::ostringstream oss;
        oss << "R" << numberForFolder(simulator.maxInitialRadius())
            << "_" << velocityFolderToken(simulator.velocityMode())
            << "_vtM" << numberForFolder(tvelMean_)
            << "_vtS" << numberForFolder(tvelStd_)
            << "_vrM" << numberForFolder(rvelMean_)
            << "_vrS" << numberForFolder(rvelStd_);
        return oss.str();
    }

    std::string nextRunDirectory(const GravitySimulator& simulator) const {
        std::ostringstream base;
        base << logRootDirectory(simulator) << "/" << initialConditionFolderName(simulator)
             << "_seed" << seedForFolder(simulator.seed());

        for (int run = 1; run < 10000; run++) {
            std::ostringstream candidate;
            candidate << base.str()
                      << "_run" << std::setw(3) << std::setfill('0') << run;

            if (!std::filesystem::exists(candidate.str())) {
                return candidate.str();
            }
        }

        return base.str() + "_run9999";
    }

    void writeInitialConditionBlock(
        std::ofstream& fout,
        const GravitySimulator& simulator,
        int steps,
        double dt,
        int ignoreSteps,
        int analyzeInterval
    ) const {
        fout << "========== Initial Conditions ==========" << '\n';
        fout << "Configuration id       = " << simulator.configurationId() << '\n';
        fout << "Configuration          = "
             << configurationName(simulator.configurationId()) << '\n';
        fout << "Number of bodies       = " << bodyCount_ << '\n';
        fout << "Velocity distribution  = "
             << velocityDistributionName(simulator.velocityMode()) << '\n';
        fout << "Velocity R dependency  = "
             << velocityRDependencyText(simulator.velocityMode()) << '\n';
        fout << "Tangential input mean  = " << tvelMean_ << '\n';
        fout << "Tangential input std   = " << tvelStd_ << '\n';
        fout << "Radial input mean      = " << rvelMean_ << '\n';
        fout << "Radial input std       = " << rvelStd_ << '\n';
        if (simulator.velocityMode() == VelocityInitializationMode::CircularReferenceGaussian) {
            fout << "Velocity formula       = v = Normal(vt_mean, vt_std)*vc(R)*e_theta"
                 << " + Normal(vr_mean, vr_std)*vc(R)*e_r" << '\n';
        } else {
            fout << "Velocity formula       = v = Normal(vt_mean, vt_std)*e_theta"
                 << " + Normal(vr_mean, vr_std)*e_r" << '\n';
        }
        fout << "Random seed            = " << simulator.seed() << '\n';
        fout << "G                      = " << simulator.gravity() << '\n';
        fout << "Dark matter alpha      = "
             << simulator.linearAttractionGravityMultiplier() << '\n';
        fout << "Dark matter density    = alpha / (1 + (r/30)^2)" << '\n';
        fout << "Dark matter M(<r)      = 4*pi*alpha*30^3*(r/30 - atan(r/30))" << '\n';
        fout << "softening              = " << simulator.softening() << '\n';
        fout << "Initial xy distribution= uniform disk" << '\n';
        fout << "Initial disk radius    = " << simulator.maxInitialRadius() << '\n';
        fout << "Initial z              = 0" << '\n';
        fout << "Merging enabled        = " << boolText(simulator.isMergingEnabled()) << '\n';
        fout << "Merge distance         = " << simulator.mergeDistance() << '\n';
        fout << "Initial gas particles  = " << simulator.initialGasParticleCount() << '\n';
        fout << "Current gas particles  = " << simulator.gasParticleCount() << '\n';
        fout << "Gas damping interval   = " << simulator.gasDampingInterval() << '\n';
        fout << "Gas damping gamma      = " << simulator.gasVelocityDampingGamma() << '\n';
        fout << "Total steps            = " << steps << '\n';
        fout << "dt                     = " << dt << '\n';
        fout << "Ignored initial steps  = " << ignoreSteps << '\n';
        fout << "Analyze interval       = " << analyzeInterval << '\n';
        fout << "View scale             = " << grid_.view().scale << '\n';
    }

    void appendAnalysisIndex(
        const GravitySimulator& simulator,
        double meanA2,
        double stdA2,
        double meanM2LogSpiral,
        double stdM2LogSpiral,
        const std::string& result
    ) const {
        const std::string indexPath = logRootDirectory(simulator) + "/analysis_index.csv";
        bool needsHeader = !std::filesystem::exists(indexPath)
            || std::filesystem::file_size(indexPath) == 0;

        std::ofstream fout(indexPath, std::ios::app);
        if (needsHeader) {
            fout << "run_directory,"
                 << "velocity_distribution,"
                 << "velocity_R_dependency,"
                 << "seed,"
                 << "vt_mean,"
                 << "vt_std,"
                 << "vr_mean,"
                 << "vr_std,"
                 << "gravity_formula,"
                 << "G,"
                 << "dark_matter_alpha,"
                 << "gas_gamma,"
                 << "mean_A2,"
                 << "std_A2,"
                 << "mean_m2_log_score,"
                 << "std_m2_log_score,"
                 << "result\n";
        }

        fout << csvEscape(runDirectory_) << ','
             << csvEscape(velocityDistributionName(simulator.velocityMode())) << ','
             << csvEscape(velocityRDependencyText(simulator.velocityMode())) << ','
             << simulator.seed() << ','
             << tvelMean_ << ','
             << tvelStd_ << ','
             << rvelMean_ << ','
             << rvelStd_ << ','
             << csvEscape(gravityFormulaText()) << ','
             << simulator.gravity() << ','
             << simulator.linearAttractionGravityMultiplier() << ','
             << simulator.gasVelocityDampingGamma() << ','
             << meanA2 << ','
             << stdA2 << ','
             << meanM2LogSpiral << ','
             << stdM2LogSpiral << ','
             << csvEscape(result) << '\n';
    }

    bool projectToImage(
        const Vec2& position,
        unsigned int imageWidth,
        unsigned int imageHeight,
        int& x,
        int& y
    ) const {
        const ViewConfig& view = grid_.view();
        double worldWidth = view.scale;
        double worldHeight = view.scale * static_cast<double>(imageHeight) / imageWidth;

        double left = view.viewCenter.x - 0.5 * worldWidth;
        double top = view.viewCenter.y + 0.5 * worldHeight;

        x = static_cast<int>((position.x - left) / worldWidth * imageWidth);
        y = static_cast<int>((top - position.y) / worldHeight * imageHeight);

        return x >= 0 && y >= 0
            && x < static_cast<int>(imageWidth)
            && y < static_cast<int>(imageHeight);
    }

    static void drawDisk(
        sf::Image& image,
        int cx,
        int cy,
        int radius,
        const sf::Color& color
    ) {
        sf::Vector2u size = image.getSize();
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                if (dx * dx + dy * dy > radius * radius) {
                    continue;
                }

                int x = cx + dx;
                int y = cy + dy;
                if (x < 0 || y < 0 ||
                    x >= static_cast<int>(size.x) || y >= static_cast<int>(size.y)) {
                    continue;
                }

                sf::Color current = image.getPixel({static_cast<unsigned int>(x), static_cast<unsigned int>(y)});
                image.setPixel(
                    {static_cast<unsigned int>(x), static_cast<unsigned int>(y)},
                    sf::Color(
                        std::max(current.r, color.r),
                        std::max(current.g, color.g),
                        std::max(current.b, color.b)
                    )
                );
            }
        }
    }

    static void drawCross(
        sf::Image& image,
        int cx,
        int cy,
        int radius,
        const sf::Color& color
    ) {
        for (int d = -radius; d <= radius; d++) {
            setPixel(image, cx + d, cy, color);
            setPixel(image, cx, cy + d, color);
        }
    }

    static void setPixel(
        sf::Image& image,
        int x,
        int y,
        const sf::Color& color
    ) {
        sf::Vector2u size = image.getSize();
        if (x < 0 || y < 0 ||
            x >= static_cast<int>(size.x) || y >= static_cast<int>(size.y)) {
            return;
        }

        image.setPixel({static_cast<unsigned int>(x), static_cast<unsigned int>(y)}, color);
    }

    static std::string numberForFolder(double value) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4) << value;
        std::string s = oss.str();

        for (char& c : s) {
            if (c == '-') {
                c = 'm';
            } else if (c == '.') {
                c = 'p';
            }
        }

        return s;
    }

    static std::string seedForFolder(unsigned int seed) {
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << seed;
        return oss.str();
    }

    static std::string logRootDirectory(const GravitySimulator& simulator) {
        return logRootDirectoryForConfiguration(simulator.configurationId());
    }

    static std::string gravityFormulaText() {
        return "pairwise softened Newtonian; dark matter rho=alpha/(1+(r/30)^2), "
               "M(<r)=4*pi*alpha*30^3*(r/30-atan(r/30)), F=-G*m*M(<r)*r/r^3";
    }

    static std::string csvEscape(const std::string& value) {
        bool needsQuotes = false;
        for (char c : value) {
            if (c == ',' || c == '"' || c == '\n' || c == '\r') {
                needsQuotes = true;
                break;
            }
        }

        if (!needsQuotes) {
            return value;
        }

        std::string escaped = "\"";
        for (char c : value) {
            if (c == '"') {
                escaped += "\"\"";
            } else {
                escaped += c;
            }
        }
        escaped += '"';
        return escaped;
    }

    static double relativeAngularMomentumDrift(
        double initialAngularMomentum,
        double finalAngularMomentum
    ) {
        if (initialAngularMomentum == 0.0) {
            return finalAngularMomentum - initialAngularMomentum;
        }

        return (finalAngularMomentum - initialAngularMomentum) / initialAngularMomentum;
    }

    static const char* boolText(bool value) {
        return value ? "true" : "false";
    }
};

#endif
