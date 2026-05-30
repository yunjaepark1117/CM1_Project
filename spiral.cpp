#include <iostream>
#include <string>
#include <vector>

#include "DistributionGrid.h"
#include "GravitySimulator.h"
#include "SpiralAnalysis.h"

using namespace std;

int main() {
    const int N = SimulationDefaults::bodyCount;
    double tvelMean;
    double tvelStd;
    double rvelMean;
    double rvelStd;

    cout << "Mean tangential circular multiplier (vt/vc): ";
    cin >> tvelMean;

    cout << "Std of tangential circular multiplier: ";
    cin >> tvelStd;

    cout << "Mean radial circular fraction (vr/vc): ";
    cin >> rvelMean;

    cout << "Std of radial circular fraction: ";
    cin >> rvelStd;

    SimulationConfig config;
    config.bodyCount = N;
    config.tangentVelocityMean = tvelMean;
    config.tangentVelocityStd = tvelStd;
    config.radialVelocityMean = rvelMean;
    config.radialVelocityStd = rvelStd;

    GravitySimulator simulator(config);
    ViewConfig view;
    view.scale = 60.0;
    DistributionGrid grid(view, simulator.minMass(), simulator.maxMass());
    SpiralAnalyzer analyzer(simulator);
    DistributionLogger logger(grid, N, tvelMean, tvelStd, rvelMean, rvelStd);

    int steps = 2000;
    double dt = 0.2;
    int ignoreSteps = 1000;
    int analyzeInterval = 1;
    int snapshotInterval = 200;
    double initialAngularMomentum = simulator.angularMomentumZ();
    vector<double> a2History;
    vector<double> m2LogSpiralHistory;
    SpiralMetrics latestMetrics;

    logger.prepare(simulator, steps, dt, ignoreSteps, analyzeInterval);
    latestMetrics = analyzer.computeMetrics();

    cout << "Number of bodies = " << N << endl;
    cout << "Using dt = " << dt << endl;
    cout << "Using gas gamma = " << config.gasVelocityDampingGamma << endl;
    cout << "Using dark matter alpha = "
         << config.linearAttractionGravityMultiplier << endl;

    for (int t = 0; t < steps; t++) {
        simulator.computeForces();
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
            cout << "step " << t
                 << " | A1 = " << snapshot.a1
                 << " | A2 = " << snapshot.a2
                 << " | A3 = " << snapshot.a3
                 << " | A4 = " << snapshot.a4;

            if (t < ignoreSteps) {
                cout << "   ignored";
            }

            cout << endl;
            logger.writeSnapshot(simulator, t, snapshot, dt);
        }
    }

    if (a2History.empty()) {
        cout << "No data to analyze. steps must be larger than ignoreSteps." << endl;
        return 0;
    }

    double mean;
    double sigma;
    double relativeSigma;
    analyzer.summarize(a2History, mean, sigma, relativeSigma);
    double meanM2LogSpiral;
    double sigmaM2LogSpiral;
    double relativeSigmaM2LogSpiral;
    analyzer.summarize(
        m2LogSpiralHistory,
        meanM2LogSpiral,
        sigmaM2LogSpiral,
        relativeSigmaM2LogSpiral
    );
    string result = analyzer.classifyTwoArmSpiral(a2History);
    double finalAngularMomentum = simulator.angularMomentumZ();

    logger.writeSummary(
        simulator,
        mean,
        sigma,
        relativeSigma,
        meanM2LogSpiral,
        sigmaM2LogSpiral,
        relativeSigmaM2LogSpiral,
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

    cout << endl;
    cout << "========== Spiral Analysis ==========" << endl;
    cout << "Mean A2                = " << mean << endl;
    cout << "Std(A2)                = " << sigma << endl;
    cout << "Relative Std(A2)       = " << relativeSigma << endl;
    cout << "Analyzed samples       = " << a2History.size() << endl;
    cout << "Result                 = " << result << endl;
    cout << "Result folder          = " << logger.runDirectory() << endl;
    cout << "Summary file           = " << logger.path() << endl;

    return 0;
}
