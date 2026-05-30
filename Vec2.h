#ifndef VEC2_H
#define VEC2_H

class Vec2 {
public:
    double x;
    double y;
    double z;

    Vec2() : x(0.0), y(0.0), z(0.0) {}
    Vec2(double x, double y) : x(x), y(y), z(0.0) {}
    Vec2(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec2 operator+(const Vec2& other) const {
        return Vec2(x + other.x, y + other.y, z + other.z);
    }

    Vec2 operator-(const Vec2& other) const {
        return Vec2(x - other.x, y - other.y, z - other.z);
    }

    Vec2 operator*(double scalar) const {
        return Vec2(x * scalar, y * scalar, z * scalar);
    }

    Vec2 operator/(double scalar) const {
        return Vec2(x / scalar, y / scalar, z / scalar);
    }

    Vec2& operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    double normSquared() const {
        return x * x + y * y + z * z;
    }

    double normSquaredXY() const {
        return x * x + y * y;
    }
};

#endif
