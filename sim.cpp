#include <SFML/Graphics.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <deque>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "DistributionGrid.h"
#include "GravitySimulator.h"
#include "SpiralAnalysis.h"

using namespace std;

class SimulationWindow {
public:
    SimulationWindow(
        double tvelMean,
        double tvelStd,
        double rvelMean,
        double rvelStd,
        int configurationId
    )
        : window_(sf::VideoMode({1100, 800}), "N-body Spiral Simulator"),
          view_(initialView(configurationId)),
          tvelMean_(tvelMean),
          tvelStd_(tvelStd),
          rvelMean_(rvelMean),
          rvelStd_(rvelStd) {
        window_.setFramerateLimit(60);
        window_.setView(view_);
        fontLoaded_ = font_.openFromFile("/System/Library/Fonts/Menlo.ttc");
    }

    bool isOpen() const {
        return window_.isOpen();
    }

    void handleEvents() {
        while (const optional<sf::Event> event = window_.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window_.close();
            }

            if (const auto* key = event->getIf<sf::Event::KeyPressed>()) {
                if (key->code == sf::Keyboard::Key::Space) {
                    paused_ = !paused_;
                } else if (key->code == sf::Keyboard::Key::Equal) {
                    zoom(0.9f);
                } else if (key->code == sf::Keyboard::Key::Hyphen) {
                    zoom(1.1f);
                }
            }

            if (const auto* pressed = event->getIf<sf::Event::MouseButtonPressed>()) {
                if (pressed->button == sf::Mouse::Button::Left) {
                    dragging_ = true;
                    lastMouseWorld_ = window_.mapPixelToCoords(pressed->position, view_);
                }
            }

            if (const auto* released = event->getIf<sf::Event::MouseButtonReleased>()) {
                if (released->button == sf::Mouse::Button::Left) {
                    dragging_ = false;
                }
            }

            if (const auto* moved = event->getIf<sf::Event::MouseMoved>()) {
                if (dragging_) {
                    sf::Vector2f current = window_.mapPixelToCoords(moved->position, view_);
                    sf::Vector2f delta = lastMouseWorld_ - current;
                    view_.move(delta);
                    window_.setView(view_);
                    lastMouseWorld_ = window_.mapPixelToCoords(moved->position, view_);
                }
            }
        }
    }

    void draw(
        const GravitySimulator& simulator,
        int step,
        int steps,
        double dt,
        const SpiralMetrics& metrics
    ) {
        (void)metrics;
        window_.clear(sf::Color(5, 7, 12));
        window_.setView(view_);

        drawTrails();
        Vec2 coreCenter = simulator.coreCenter();
        drawBodies(simulator);
        drawCoreCenter(coreCenter);
        drawHud(simulator, step, steps, dt, metrics);
        drawPolarDensityMinimap(simulator, coreCenter);

        window_.display();
        updateTitle(simulator, step, steps, dt);
    }

    void updateTrail(const GravitySimulator& simulator) {
        if (trails_.size() != simulator.bodies().size()) {
            resetTrails(simulator);
        }

        const vector<Body>& bodies = simulator.bodies();
        for (size_t i = 0; i < bodies.size(); i++) {
            trails_[i].push_back(toScreenPoint(bodies[i].position()));
            while (trails_[i].size() > trailLength_) {
                trails_[i].pop_front();
            }
        }
    }

    bool isPaused() const {
        return paused_;
    }

    void resetTrails(const GravitySimulator& simulator) {
        trails_.assign(simulator.bodies().size(), deque<sf::Vector2f>());

        const vector<Body>& bodies = simulator.bodies();
        for (size_t i = 0; i < bodies.size(); i++) {
            trails_[i].push_back(toScreenPoint(bodies[i].position()));
        }
    }

private:
    struct LogSpiralFit {
        bool valid = false;
        int startRadialBin = 0;
        double slope = 0.0;
        double intercept = 0.0;
        double rSquared = 0.0;
        double pValue = 1.0;
        double adjustedPValue = 1.0;
        double pitchAngleDegrees = 0.0;
        vector<int> thetaBins;
    };

    sf::RenderWindow window_;
    sf::View view_;
    sf::Font font_;
    vector<deque<sf::Vector2f> > trails_;
    bool dragging_ = false;
    bool paused_ = false;
    bool fontLoaded_ = false;
    sf::Vector2f lastMouseWorld_;
    double tvelMean_;
    double tvelStd_;
    double rvelMean_;
    double rvelStd_;
    const size_t trailLength_ = 80;

    static sf::FloatRect initialView(int configurationId) {
        if (normalizedConfigurationId(configurationId) >= 3) {
            return sf::FloatRect({-140.0f, -105.0f}, {280.0f, 210.0f});
        }

        return sf::FloatRect({-70.0f, -52.0f}, {140.0f, 104.0f});
    }

    void drawTrails() {
        for (const deque<sf::Vector2f>& trail : trails_) {
            if (trail.size() < 2) {
                continue;
            }

            vector<sf::Vertex> vertices;
            vertices.reserve(trail.size());

            for (size_t i = 0; i < trail.size(); i++) {
                float age = static_cast<float>(i + 1) / static_cast<float>(trail.size());
                sf::Color color(70, 130, 230, static_cast<uint8_t>(6 + 32 * age));
                vertices.push_back(sf::Vertex{trail[i], color});
            }

            window_.draw(vertices.data(), vertices.size(), sf::PrimitiveType::LineStrip);
        }
    }

    void drawBodies(const GravitySimulator& simulator) {
        for (const Body& body : simulator.bodies()) {
            float radius = particleRadius(body.mass(), simulator.maxMass());
            sf::CircleShape circle(radius);
            circle.setOrigin({radius, radius});
            circle.setPosition(toScreenPoint(body.position()));
            circle.setFillColor(colorByMass(body.mass(), simulator.maxMass(), body.isGas()));
            window_.draw(circle);
        }
    }

    void drawCoreCenter(const Vec2& coreCenter) {
        sf::Vector2f center = toScreenPoint(coreCenter);
        float size = 0.55f;

        sf::Vertex horizontal[] = {
            sf::Vertex{{center.x - size, center.y}, sf::Color(255, 80, 80)},
            sf::Vertex{{center.x + size, center.y}, sf::Color(255, 80, 80)}
        };
        sf::Vertex vertical[] = {
            sf::Vertex{{center.x, center.y - size}, sf::Color(255, 80, 80)},
            sf::Vertex{{center.x, center.y + size}, sf::Color(255, 80, 80)}
        };

        window_.draw(horizontal, 2, sf::PrimitiveType::Lines);
        window_.draw(vertical, 2, sf::PrimitiveType::Lines);
    }

    void drawHud(
        const GravitySimulator& simulator,
        int step,
        int steps,
        double dt,
        const SpiralMetrics& metrics
    ) {
        if (!fontLoaded_) {
            return;
        }

        sf::View oldView = window_.getView();
        window_.setView(window_.getDefaultView());

        ostringstream hud;
        hud << fixed << setprecision(3)
            << "step " << step << " / " << steps << '\n'
            << "N=" << simulator.bodyCount()
            << "  gas=" << simulator.gasParticleCount()
            << "  dt=" << dt << '\n'
            << "gas gamma/interval = " << simulator.gasVelocityDampingGamma()
            << " / " << simulator.gasDampingInterval() << '\n'
            << "m2 log score/pitch = " << metrics.m2LogSpiral.score
            << " / " << metrics.m2LogSpiral.pitchAngleDegrees << '\n'
            << "vt velocity mean/std = " << tvelMean_ << " / " << tvelStd_ << '\n'
            << "vr velocity mean/std = " << rvelMean_ << " / " << rvelStd_ << '\n'
            << (paused_ ? "PAUSED" : "RUNNING")
            << "  space: pause  +/-: zoom  drag: pan";

        sf::Text text(font_, hud.str(), 16);
        text.setPosition({14.0f, 12.0f});
        text.setFillColor(sf::Color(230, 238, 245));

        sf::FloatRect bounds = text.getGlobalBounds();
        sf::RectangleShape background({bounds.size.x + 18.0f, bounds.size.y + 16.0f});
        background.setPosition({8.0f, 8.0f});
        background.setFillColor(sf::Color(0, 0, 0, 150));

        window_.draw(background);
        window_.draw(text);
        window_.setView(oldView);
    }

    void updateTitle(
        const GravitySimulator& simulator,
        int step,
        int steps,
        double dt
    ) {
        ostringstream title;
        title << fixed << setprecision(3)
              << "N-body | step " << step << "/" << steps
              << " | N=" << simulator.bodyCount()
              << " | dt=" << dt
              << " | tv=(" << tvelMean_ << "," << tvelStd_ << ")"
              << " | rv=(" << rvelMean_ << "," << rvelStd_ << ")"
              << (paused_ ? " | paused" : "")
              << " | +/-: zoom, drag: pan";
        window_.setTitle(title.str());
    }

    void drawPolarDensityMinimap(const GravitySimulator& simulator, const Vec2& center) {
        const int radialBins = 6;
        const int thetaBins = 15;
        const double maxRadius = 30.0;
        const float padding = 10.0f;
        const float plotRadius = 82.0f;
        const float panelWidth = padding * 2.0f + 2.0f * plotRadius;
        const float panelHeight = panelWidth + 48.0f;
        const double fullCircle = 2.0 * std::acos(-1.0);

        vector<vector<double> > massGrid(
            radialBins,
            vector<double>(thetaBins, 0.0)
        );

        for (const Body& body : simulator.bodies()) {
            Vec2 relativePosition = body.position() - center;
            double radius = std::sqrt(relativePosition.normSquaredXY());
            if (radius < 0.0 || radius >= maxRadius) {
                continue;
            }

            double theta = std::atan2(relativePosition.y, relativePosition.x);
            if (theta < 0.0) {
                theta += fullCircle;
            }

            int rBin = static_cast<int>(radius / maxRadius * radialBins);
            int tBin = static_cast<int>(theta / fullCircle * thetaBins);
            rBin = std::clamp(rBin, 0, radialBins - 1);
            tBin = std::clamp(tBin, 0, thetaBins - 1);
            massGrid[rBin][tBin] += body.mass();
        }

        LogSpiralFit bestFit = computeBestLogSpiralFit(massGrid, maxRadius, fullCircle);

        sf::View oldView = window_.getView();
        window_.setView(window_.getDefaultView());

        sf::Vector2u windowSize = window_.getSize();
        sf::Vector2f origin(
            static_cast<float>(windowSize.x) - panelWidth - 12.0f,
            12.0f
        );

        sf::RectangleShape background({panelWidth, panelHeight});
        background.setPosition(origin);
        background.setFillColor(sf::Color(0, 0, 0, 155));
        background.setOutlineThickness(1.0f);
        background.setOutlineColor(sf::Color(210, 220, 235, 70));
        window_.draw(background);

        sf::Vector2f plotCenter(
            origin.x + padding + plotRadius,
            origin.y + padding + plotRadius
        );

        for (int r = 0; r < radialBins; r++) {
            vector<int> thetaOrder(thetaBins);
            for (int theta = 0; theta < thetaBins; theta++) {
                thetaOrder[theta] = theta;
            }

            std::sort(
                thetaOrder.begin(),
                thetaOrder.end(),
                [&](int a, int b) {
                    return massGrid[r][a] < massGrid[r][b];
                }
            );

            vector<int> rankByTheta(thetaBins, 0);
            for (int rank = 0; rank < thetaBins; rank++) {
                rankByTheta[thetaOrder[rank]] = rank;
            }

            for (int theta = 0; theta < thetaBins; theta++) {
                sf::Color color = densityRankColor(rankByTheta[theta]);

                float innerRadius = plotRadius * static_cast<float>(r) / radialBins;
                float outerRadius = plotRadius * static_cast<float>(r + 1) / radialBins;
                double startAngle = fullCircle * static_cast<double>(theta) / thetaBins;
                double endAngle = fullCircle * static_cast<double>(theta + 1) / thetaBins;

                sf::ConvexShape sector = makeAnnularSector(
                    plotCenter,
                    innerRadius,
                    outerRadius,
                    startAngle,
                    endAngle
                );
                sector.setFillColor(color);
                sector.setOutlineThickness(0.35f);
                sector.setOutlineColor(sf::Color(4, 6, 10, 120));
                window_.draw(sector);
            }
        }

        sf::CircleShape outline(plotRadius);
        outline.setOrigin({plotRadius, plotRadius});
        outline.setPosition(plotCenter);
        outline.setFillColor(sf::Color::Transparent);
        outline.setOutlineThickness(1.0f);
        outline.setOutlineColor(sf::Color(215, 225, 245, 85));
        window_.draw(outline);

        drawLogSpiralFitOverlay(bestFit, plotCenter, plotRadius, radialBins, thetaBins, fullCircle);
        drawLogSpiralFitText(bestFit, origin + sf::Vector2f(10.0f, panelWidth + 8.0f));

        window_.setView(oldView);
    }

    LogSpiralFit computeBestLogSpiralFit(
        const vector<vector<double> >& massGrid,
        double maxRadius,
        double fullCircle
    ) const {
        const int radialBins = static_cast<int>(massGrid.size());
        const int thetaBins = radialBins > 0 ? static_cast<int>(massGrid[0].size()) : 0;
        const int fitLength = 5;
        const int peakChoices = 2;
        const int combinations = 1 << fitLength;
        int testsEvaluated = 0;
        LogSpiralFit best;
        best.rSquared = -1.0;

        if (radialBins < fitLength || thetaBins < peakChoices) {
            return LogSpiralFit{};
        }

        vector<vector<int> > topTheta(radialBins, vector<int>(peakChoices, 0));
        vector<vector<double> > topMass(radialBins, vector<double>(peakChoices, 0.0));

        for (int r = 0; r < radialBins; r++) {
            vector<int> order(thetaBins);
            for (int theta = 0; theta < thetaBins; theta++) {
                order[theta] = theta;
            }

            std::sort(
                order.begin(),
                order.end(),
                [&](int a, int b) {
                    return massGrid[r][a] > massGrid[r][b];
                }
            );

            for (int peak = 0; peak < peakChoices; peak++) {
                topTheta[r][peak] = order[peak];
                topMass[r][peak] = massGrid[r][order[peak]];
            }
        }

        for (int start = 0; start <= radialBins - fitLength; start++) {
            for (int combination = 0; combination < combinations; combination++) {
                vector<double> x;
                vector<double> y;
                vector<int> chosenThetaBins;
                x.reserve(fitLength);
                y.reserve(fitLength);
                chosenThetaBins.reserve(fitLength);

                bool usable = true;
                for (int i = 0; i < fitLength; i++) {
                    int peakChoice = (combination >> i) & 1;
                    int rBin = start + i;
                    if (topMass[rBin][peakChoice] <= 0.0) {
                        usable = false;
                        break;
                    }

                    double radiusCenter =
                        maxRadius * (static_cast<double>(rBin) + 0.5) / radialBins;
                    double thetaCenter =
                        fullCircle * (static_cast<double>(topTheta[rBin][peakChoice]) + 0.5)
                        / thetaBins;

                    if (!y.empty()) {
                        thetaCenter = unwrapNear(thetaCenter, y.back(), fullCircle);
                    }

                    x.push_back(std::log(radiusCenter));
                    y.push_back(thetaCenter);
                    chosenThetaBins.push_back(topTheta[rBin][peakChoice]);
                }

                if (!usable) {
                    continue;
                }

                testsEvaluated++;
                LogSpiralFit fit = fitLine(x, y);
                fit.startRadialBin = start;
                fit.thetaBins = chosenThetaBins;

                if (fit.valid && fit.rSquared > best.rSquared) {
                    best = fit;
                }
            }
        }

        if (!best.valid || testsEvaluated <= 0) {
            return LogSpiralFit{};
        }

        best.adjustedPValue = std::min(1.0, best.pValue * testsEvaluated);
        return best;
    }

    static double unwrapNear(double theta, double reference, double fullCircle) {
        while (theta - reference > 0.5 * fullCircle) {
            theta -= fullCircle;
        }
        while (theta - reference < -0.5 * fullCircle) {
            theta += fullCircle;
        }
        return theta;
    }

    static LogSpiralFit fitLine(const vector<double>& x, const vector<double>& y) {
        LogSpiralFit fit;
        int n = static_cast<int>(x.size());
        if (n < 3 || y.size() != x.size()) {
            return fit;
        }

        double meanX = 0.0;
        double meanY = 0.0;
        for (int i = 0; i < n; i++) {
            meanX += x[i];
            meanY += y[i];
        }
        meanX /= n;
        meanY /= n;

        double sxx = 0.0;
        double sxy = 0.0;
        double syy = 0.0;
        for (int i = 0; i < n; i++) {
            double dx = x[i] - meanX;
            double dy = y[i] - meanY;
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }

        if (sxx <= 0.0 || syy <= 0.0) {
            return fit;
        }

        fit.slope = sxy / sxx;
        fit.intercept = meanY - fit.slope * meanX;

        double sse = 0.0;
        for (int i = 0; i < n; i++) {
            double predicted = fit.slope * x[i] + fit.intercept;
            sse += (y[i] - predicted) * (y[i] - predicted);
        }

        fit.rSquared = std::clamp(1.0 - sse / syy, 0.0, 1.0);
        fit.pValue = slopePValue(fit.rSquared, n);
        fit.adjustedPValue = fit.pValue;
        fit.pitchAngleDegrees =
            std::atan2(1.0, std::max(std::abs(fit.slope), 1e-9)) * 180.0 / std::acos(-1.0);
        fit.valid = true;
        return fit;
    }

    static double slopePValue(double rSquared, int n) {
        if (n <= 2) {
            return 1.0;
        }
        if (rSquared >= 0.999999) {
            return 0.0;
        }

        double t = std::sqrt(rSquared * (n - 2) / std::max(1e-12, 1.0 - rSquared));

        if (n == 5) {
            double x = t / std::sqrt(3.0);
            double positiveArea =
                (std::atan(x) + x / (1.0 + x * x)) / std::acos(-1.0);
            return std::clamp(1.0 - 2.0 * positiveArea, 0.0, 1.0);
        }

        return std::clamp(std::exp(-0.5 * t * t), 0.0, 1.0);
    }

    void drawLogSpiralFitOverlay(
        const LogSpiralFit& fit,
        const sf::Vector2f& plotCenter,
        float plotRadius,
        int radialBins,
        int thetaBins,
        double fullCircle
    ) {
        if (!fit.valid || fit.thetaBins.size() != 5) {
            return;
        }

        vector<sf::Vertex> vertices;
        vertices.reserve(fit.thetaBins.size());

        for (size_t i = 0; i < fit.thetaBins.size(); i++) {
            int rBin = fit.startRadialBin + static_cast<int>(i);
            float radius =
                plotRadius * (static_cast<float>(rBin) + 0.5f) / radialBins;
            double theta =
                fullCircle * (static_cast<double>(fit.thetaBins[i]) + 0.5) / thetaBins;
            sf::Vector2f point(
                plotCenter.x + radius * static_cast<float>(std::cos(theta)),
                plotCenter.y - radius * static_cast<float>(std::sin(theta))
            );

            vertices.push_back(sf::Vertex{point, sf::Color(255, 120, 255, 245)});

            sf::CircleShape marker(2.5f);
            marker.setOrigin({2.5f, 2.5f});
            marker.setPosition(point);
            marker.setFillColor(sf::Color(255, 120, 255, 245));
            window_.draw(marker);
        }

        window_.draw(vertices.data(), vertices.size(), sf::PrimitiveType::LineStrip);
    }

    void drawLogSpiralFitText(const LogSpiralFit& fit, const sf::Vector2f& position) {
        if (!fontLoaded_) {
            return;
        }

        ostringstream text;
        text << fixed << setprecision(3);
        if (fit.valid) {
            text << "log spiral fit  R2=" << fit.rSquared
                 << "  p*=" << fit.adjustedPValue << '\n'
                 << "pitch=" << setprecision(1) << fit.pitchAngleDegrees
                 << " deg  bins " << fit.startRadialBin
                 << "-" << fit.startRadialBin + 4;
        } else {
            text << "log spiral fit\nnot enough peaks";
        }

        sf::Text label(font_, text.str(), 11);
        label.setPosition(position);
        label.setFillColor(fit.valid && fit.adjustedPValue < 0.05
            ? sf::Color(255, 245, 190)
            : sf::Color(180, 190, 205));
        window_.draw(label);
    }

    static sf::Color densityRankColor(int rank) {
        if (rank < 12) {
            return sf::Color(18, 24, 34, 230);
        }

        if (rank == 12) {
            return sf::Color(255, 245, 235, 240);
        }

        if (rank == 13) {
            return sf::Color(255, 215, 45, 238);
        }

        return sf::Color(255, 65, 55, 240);
    }

    static sf::ConvexShape makeAnnularSector(
        const sf::Vector2f& center,
        float innerRadius,
        float outerRadius,
        double startAngle,
        double endAngle
    ) {
        const int segments = 5;
        sf::ConvexShape shape;
        shape.setPointCount(2 * (segments + 1));

        for (int i = 0; i <= segments; i++) {
            double t = static_cast<double>(i) / segments;
            double angle = startAngle + (endAngle - startAngle) * t;
            shape.setPoint(
                i,
                sf::Vector2f(
                    center.x + outerRadius * static_cast<float>(std::cos(angle)),
                    center.y - outerRadius * static_cast<float>(std::sin(angle))
                )
            );
        }

        for (int i = 0; i <= segments; i++) {
            double t = static_cast<double>(segments - i) / segments;
            double angle = startAngle + (endAngle - startAngle) * t;
            shape.setPoint(
                segments + 1 + i,
                sf::Vector2f(
                    center.x + innerRadius * static_cast<float>(std::cos(angle)),
                    center.y - innerRadius * static_cast<float>(std::sin(angle))
                )
            );
        }

        return shape;
    }

    void zoom(float factor) {
        view_.zoom(factor);
        window_.setView(view_);
    }

    static sf::Vector2f toScreenPoint(const Vec2& position) {
        return sf::Vector2f(
            static_cast<float>(position.x),
            static_cast<float>(-position.y)
        );
    }

    static float particleRadius(double mass, double maxMass) {
        double normalized = maxMass > 0.0 ? mass / maxMass : 0.0;
        return static_cast<float>(0.06 + 0.22 * std::sqrt(normalized));
    }

    static sf::Color colorByMass(double mass, double maxMass, bool isGas) {
        double normalized = maxMass > 0.0 ? mass / maxMass : 0.0;
        normalized = std::clamp(normalized, 0.0, 1.0);
        double brightness = 0.35 + 0.65 * std::sqrt(normalized);

        if (isGas) {
            return sf::Color(
                static_cast<uint8_t>(190 + 65 * brightness),
                static_cast<uint8_t>(165 + 80 * brightness),
                static_cast<uint8_t>(70 + 55 * brightness)
            );
        }

        return sf::Color(
            static_cast<uint8_t>(120 + 135 * brightness),
            static_cast<uint8_t>(150 + 105 * brightness),
            static_cast<uint8_t>(190 + 65 * brightness)
        );
    }
};

static double unwrapAngleDelta(double delta) {
    const double pi = std::acos(-1.0);
    while (delta > pi) {
        delta -= 2.0 * pi;
    }
    while (delta < -pi) {
        delta += 2.0 * pi;
    }
    return delta;
}

int main() {
    const int N = SimulationDefaults::bodyCount;
    int configurationChoice;
    double tvelMean;
    double tvelStd;
    double rvelMean;
    double rvelStd;

    cout << "Choose configuration (1-4):" << endl;
    cout << "1. R=48, position-independent Gaussian velocity" << endl;
    cout << "2. R=48, circular-reference Gaussian multiplier" << endl;
    cout << "3. R=120, position-independent Gaussian velocity" << endl;
    cout << "4. R=120, circular-reference Gaussian multiplier" << endl;
    cout << "Configuration: ";
    cin >> configurationChoice;

    SimulationConfig config;
    config.bodyCount = N;
    applyConfiguration(config, configurationChoice);

    if (config.velocityMode == VelocityInitializationMode::CircularReferenceGaussian) {
        cout << "Tangential circular multiplier mean: ";
    } else {
        cout << "Tangential velocity mean: ";
    }
    cin >> tvelMean;

    if (config.velocityMode == VelocityInitializationMode::CircularReferenceGaussian) {
        cout << "Tangential circular multiplier std: ";
    } else {
        cout << "Tangential velocity std: ";
    }
    cin >> tvelStd;

    if (config.velocityMode == VelocityInitializationMode::CircularReferenceGaussian) {
        cout << "Radial circular multiplier mean: ";
    } else {
        cout << "Radial velocity mean: ";
    }
    cin >> rvelMean;

    if (config.velocityMode == VelocityInitializationMode::CircularReferenceGaussian) {
        cout << "Radial circular multiplier std: ";
    } else {
        cout << "Radial velocity std: ";
    }
    cin >> rvelStd;

    config.tangentVelocityMean = tvelMean;
    config.tangentVelocityStd = tvelStd;
    config.radialVelocityMean = rvelMean;
    config.radialVelocityStd = rvelStd;

    GravitySimulator simulator(config);
    ViewConfig view;
    view.scale = viewScaleForConfiguration(config.configurationId);
    DistributionGrid grid(view, simulator.minMass(), simulator.maxMass());
    SpiralAnalyzer analyzer(simulator);
    DistributionLogger logger(grid, N, tvelMean, tvelStd, rvelMean, rvelStd);
    SimulationWindow display(tvelMean, tvelStd, rvelMean, rvelStd, config.configurationId);

    int steps = 2000;
    double dt = 0.2;
    int ignoreSteps = 1000;
    int analyzeInterval = 1;
    int snapshotInterval = 200;
    int drawInterval = 5;
    double initialAngularMomentum = simulator.angularMomentumZ();
    double meanA2 = 0.0;
    double stdA2 = 0.0;
    double relativeStdA2 = 0.0;
    double patternSpeed = 0.0;
    double previousPhase2 = 0.0;
    bool hasPreviousPhase2 = false;
    SpiralMetrics latestMetrics;
    vector<double> a2History;
    vector<double> m2LogSpiralHistory;

    logger.prepare(simulator, steps, dt, ignoreSteps, analyzeInterval);
    latestMetrics = analyzer.computeMetrics();
    simulator.computeForces();

    cout << "Number of bodies = " << N << endl;
    cout << "Using dt = " << dt << endl;
    cout << "Using gas gamma = " << config.gasVelocityDampingGamma << endl;
    cout << "Using dark matter alpha = "
         << config.linearAttractionGravityMultiplier << endl;
    cout << "SFML display: +/- zoom, space pause, left-drag pan." << endl;

    int t = 0;
    while (t < steps && display.isOpen()) {
        display.handleEvents();

        if (display.isPaused()) {
            display.draw(
                simulator,
                t,
                steps,
                dt,
                latestMetrics
            );
            this_thread::sleep_for(chrono::milliseconds(16));
            continue;
        }

        simulator.step(dt);
        if (simulator.gasDampingInterval() > 0 &&
            (t + 1) % simulator.gasDampingInterval() == 0) {
            simulator.applyGasDamping();
        }
        simulator.applyMergePolicy();
        display.updateTrail(simulator);

        latestMetrics = analyzer.computeMetrics();
        FourierSnapshot snapshot = latestMetrics.snapshot;

        if (hasPreviousPhase2) {
            patternSpeed = unwrapAngleDelta(latestMetrics.phase2 - previousPhase2) / dt;
        }
        previousPhase2 = latestMetrics.phase2;
        hasPreviousPhase2 = true;

        if (t >= ignoreSteps && t % analyzeInterval == 0) {
            a2History.push_back(snapshot.a2);
            m2LogSpiralHistory.push_back(latestMetrics.m2LogSpiral.score);
            analyzer.summarize(a2History, meanA2, stdA2, relativeStdA2);
        }

        if (t % snapshotInterval == 0) {
            logger.writeSnapshot(simulator, t, snapshot, dt);
        }

        if (t % drawInterval == 0) {
            display.draw(
                simulator,
                t,
                steps,
                dt,
                latestMetrics
            );
        }

        t++;
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

    cout << "Result folder = " << logger.runDirectory() << endl;
    cout << "Summary file = " << logger.path() << endl;

    return 0;
}
