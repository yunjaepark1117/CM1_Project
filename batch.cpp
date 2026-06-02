#include <iomanip>
#include <iostream>
#include <chrono>
#include <sstream>
#include <string>
#include <vector>

#include "DistributionGrid.h"
#include "GravitySimulator.h"
#include "SpiralAnalysis.h"

using namespace std;

struct BatchRequest {
    double tvelMean;
    double tvelStd;
    double rvelMean;
    double rvelStd;
    unsigned int seed;
};

static const int batchConfigurationChoice = 2;

static const vector<BatchRequest> requests = {
    // For configurations 1/3: {tangential velocity mean, std, radial velocity mean, std, seed}
    // For configurations 2/4: {tangential circular multiplier mean, std,
    //                          radial circular multiplier mean, std, seed}
    

    {0.8, 0.0001, 0, 0.0001, 51},
    {0.8, 0.0001, 0, 0.0001, 52},
    {0.8, 0.0001, 0, 0.0001, 53},

    {0.8, 0.0002, 0, 0.0002, 51},
    {0.8, 0.0002, 0, 0.0002, 52},
    {0.8, 0.0002, 0, 0.0002, 53},

    {0.8, 0.0005, 0, 0.0005, 51},
    {0.8, 0.0005, 0, 0.0005, 52},
    {0.8, 0.0005, 0, 0.0005, 53},

    {0.8, 0.001, 0, 0.001, 51},
    {0.8, 0.001, 0, 0.001, 52},
    {0.8, 0.001, 0, 0.001, 53},

    {0.8, 0.002, 0, 0.0002, 51},
    {0.8, 0.002, 0, 0.0002, 52},
    {0.8, 0.002, 0, 0.002, 53},

    {0.8, 0.005, 0, 0.005, 51},
    {0.8, 0.005, 0, 0.005, 52},
    {0.8, 0.005, 0, 0.005, 53},

    

    
    

    
};

static string formatDuration(double seconds) {
    long long totalSeconds = static_cast<long long>(seconds + 0.5);
    long long hours = totalSeconds / 3600;
    long long minutes = (totalSeconds % 3600) / 60;
    long long secs = totalSeconds % 60;

    ostringstream oss;
    if (hours > 0) {
        oss << hours << "h ";
    }
    if (hours > 0 || minutes > 0) {
        oss << minutes << "m ";
    }
    oss << secs << "s";
    return oss.str();
}

static void printProgressLine(
    int requestIndex,
    int requestCount,
    int savedSnapshots,
    int step,
    int steps,
    chrono::steady_clock::time_point batchStart
) {
    using Clock = chrono::steady_clock;

    int completedSteps = (requestIndex - 1) * steps + min(step + 1, steps);
    int totalSteps = requestCount * steps;
    double progress = totalSteps > 0
        ? static_cast<double>(completedSteps) / static_cast<double>(totalSteps)
        : 0.0;

    double elapsedSeconds =
        chrono::duration<double>(Clock::now() - batchStart).count();
    double remainingSeconds = progress > 0.0
        ? elapsedSeconds * (1.0 / progress - 1.0)
        : 0.0;

    cout << '\r'
         << "request " << requestIndex << " / " << requestCount
         << " | snapshot " << savedSnapshots
         << " saved | step " << setw(4) << step
         << " | elapsed " << formatDuration(elapsedSeconds)
         << " | ETA " << formatDuration(remainingSeconds)
         << " | progress " << fixed << setprecision(1) << progress * 100.0 << "%"
         << string(8, ' ') << flush;
}

static void runSingleRequest(
    const BatchRequest& request,
    int requestIndex,
    int requestCount,
    chrono::steady_clock::time_point batchStart
) {
    const int N = SimulationDefaults::bodyCount;
    const int steps = 2000;
    const double dt = 0.2;
    const int ignoreSteps = 1000;
    const int analyzeInterval = 1;
    const int snapshotInterval = 200;
    const int progressInterval = 20;

    SimulationConfig config;
    config.bodyCount = N;
    applyConfiguration(config, batchConfigurationChoice);
    config.tangentVelocityMean = request.tvelMean;
    config.tangentVelocityStd = request.tvelStd;
    config.radialVelocityMean = request.rvelMean;
    config.radialVelocityStd = request.rvelStd;
    config.seed = request.seed;

    GravitySimulator simulator(config);
    ViewConfig view;
    view.scale = viewScaleForConfiguration(config.configurationId);
    DistributionGrid grid(view, simulator.minMass(), simulator.maxMass());
    SpiralAnalyzer analyzer(simulator);
    DistributionLogger logger(
        grid,
        N,
        request.tvelMean,
        request.tvelStd,
        request.rvelMean,
        request.rvelStd
    );

    vector<double> a2History;
    vector<double> m2LogSpiralHistory;
    SpiralMetrics latestMetrics;
    double initialAngularMomentum = simulator.angularMomentumZ();
    int savedSnapshots = 0;

    cout << "Starting request " << requestIndex << " / " << requestCount
         << " | config " << config.configurationId
         << " | vt input mean/std = " << request.tvelMean << " / " << request.tvelStd
         << " | vr input mean/std = " << request.rvelMean << " / " << request.rvelStd
         << " | seed = " << request.seed << endl;

    logger.prepare(simulator, steps, dt, ignoreSteps, analyzeInterval);
    latestMetrics = analyzer.computeMetrics();
    simulator.computeForces();
    printProgressLine(
        requestIndex,
        requestCount,
        savedSnapshots,
        0,
        steps,
        batchStart
    );

    for (int t = 0; t < steps; t++) {
        simulator.step(dt);
        if (simulator.gasDampingInterval() > 0 &&
            (t + 1) % simulator.gasDampingInterval() == 0) {
            simulator.applyGasDamping();
        }
        simulator.applyMergePolicy();

        latestMetrics = analyzer.computeMetrics();
        FourierSnapshot snapshot = latestMetrics.snapshot;

        if (t >= ignoreSteps && t % analyzeInterval == 0) {
            a2History.push_back(snapshot.a2);
            m2LogSpiralHistory.push_back(latestMetrics.m2LogSpiral.score);
        }

        if (t % snapshotInterval == 0) {
            logger.writeSnapshot(simulator, t, snapshot, dt);
            savedSnapshots++;
            printProgressLine(
                requestIndex,
                requestCount,
                savedSnapshots,
                t,
                steps,
                batchStart
            );
        } else if (t % progressInterval == 0) {
            printProgressLine(
                requestIndex,
                requestCount,
                savedSnapshots,
                t,
                steps,
                batchStart
            );
        }
    }

    double meanA2;
    double stdA2;
    double relativeStdA2;
    double meanM2LogSpiral;
    double stdM2LogSpiral;
    double relativeStdM2LogSpiral;
    analyzer.summarize(a2History, meanA2, stdA2, relativeStdA2);
    analyzer.summarize(
        m2LogSpiralHistory,
        meanM2LogSpiral,
        stdM2LogSpiral,
        relativeStdM2LogSpiral
    );

    string result = analyzer.classifyTwoArmSpiral(a2History);
    double finalAngularMomentum = simulator.angularMomentumZ();

    logger.writeSummary(
        simulator,
        meanA2,
        stdA2,
        relativeStdA2,
        meanM2LogSpiral,
        stdM2LogSpiral,
        relativeStdM2LogSpiral,
        latestMetrics.m2LogSpiral,
        static_cast<int>(a2History.size()),
        steps,
        dt,
        ignoreSteps,
        analyzeInterval,
        initialAngularMomentum,
        finalAngularMomentum,
        result
    );

    printProgressLine(
        requestIndex,
        requestCount,
        savedSnapshots,
        steps - 1,
        steps,
        batchStart
    );
    cout << " | done | folder " << logger.runDirectory()
         << string(12, ' ') << endl;
}

int main() {
    cout.setf(ios::unitbuf);

    if (requests.empty()) {
        cout << "No batch requests are configured." << endl;
        return 0;
    }

    cout << "Batch requests = " << requests.size() << endl;
    cout << "Batch configuration = " << normalizedConfigurationId(batchConfigurationChoice)
         << " (" << configurationName(batchConfigurationChoice) << ")" << endl;
    cout << "Fixed gas gamma = " << SimulationDefaults::gasVelocityDampingGamma << endl;
    SimulationConfig displayConfig;
    applyConfiguration(displayConfig, batchConfigurationChoice);
    cout << "Dark matter alpha = "
         << displayConfig.linearAttractionGravityMultiplier << endl;

    chrono::steady_clock::time_point batchStart = chrono::steady_clock::now();
    for (int i = 0; i < static_cast<int>(requests.size()); i++) {
        runSingleRequest(
            requests[i],
            i + 1,
            static_cast<int>(requests.size()),
            batchStart
        );
    }

    cout << "All batch simulations finished." << endl;
    return 0;
}
