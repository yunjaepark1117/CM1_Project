#ifndef DISTRIBUTION_GRID_H
#define DISTRIBUTION_GRID_H

#include <string>
#include <vector>

#include "Body.h"

struct ViewConfig {
    int width = 80;
    int height = 30;
    Vec2 viewCenter = Vec2(0.0, 0.0);
    double scale = 30.0;
};

class DistributionGrid {
public:
    DistributionGrid(const ViewConfig& view, double minMass, double maxMass)
        : view_(view), maxMass_(maxMass) {
        (void)minMass;
    }

    std::vector<std::vector<double> > massGrid(const std::vector<Body>& bodies) const {
        std::vector<std::vector<double> > grid(
            view_.height,
            std::vector<double>(view_.width, 0.0)
        );

        for (const Body& body : bodies) {
            int x;
            int y;

            if (project(body.position(), x, y)) {
                grid[y][x] += body.mass();
            }
        }

        return grid;
    }

    std::vector<std::string> screen(
        const std::vector<Body>& bodies,
        const Vec2& markerPosition
    ) const {
        std::vector<std::string> output(view_.height, std::string(view_.width, ' '));
        std::vector<std::vector<double> > grid = massGrid(bodies);

        for (int y = 0; y < view_.height; y++) {
            for (int x = 0; x < view_.width; x++) {
                output[y][x] = symbolForMass(grid[y][x]);
            }
        }

        int x;
        int y;
        if (project(markerPosition, x, y)) {
            output[y][x] = 'X';
        }

        return output;
    }

    char symbolForMass(double mass) const {
        if (mass <= 0.0) {
            return ' ';
        }

        if (mass < maxMass_) {
            return '.';
        }

        if (mass < 2.0 * maxMass_) {
            return '^';
        }

        if (mass < 5.0 * maxMass_) {
            return '*';
        }

        if (mass < 15.0 * maxMass_) {
            return 'O';
        }

        return '#';
    }

    const ViewConfig& view() const {
        return view_;
    }

private:
    ViewConfig view_;
    double maxMass_;

    bool project(const Vec2& position, int& x, int& y) const {
        Vec2 relativePosition = position - view_.viewCenter;

        x = static_cast<int>(
            view_.width / 2 + relativePosition.x / view_.scale * view_.width
        );
        y = static_cast<int>(
            view_.height / 2 - relativePosition.y / view_.scale * view_.height
        );

        return x >= 0 && x < view_.width && y >= 0 && y < view_.height;
    }

    Vec2 unproject(int x, int y) const {
        double worldX =
            (static_cast<double>(x) - view_.width / 2.0) * view_.scale / view_.width
            + view_.viewCenter.x;

        double worldY =
            -(static_cast<double>(y) - view_.height / 2.0) * view_.scale / view_.height
            + view_.viewCenter.y;

        return Vec2(worldX, worldY);
    }
};

#endif
