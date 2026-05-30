#ifndef BODY_H
#define BODY_H

#include "Vec2.h"

class Body {
public:
    Body(double mass, const Vec2& position, const Vec2& velocity)
        : mass_(mass),
          position_(position),
          velocity_(velocity),
          force_(0.0, 0.0, 0.0),
          isGas_(false) {}

    Body(double mass, const Vec2& position, const Vec2& velocity, bool isGas)
        : mass_(mass),
          position_(position),
          velocity_(velocity),
          force_(0.0, 0.0, 0.0),
          isGas_(isGas) {}

    double mass() const {
        return mass_;
    }

    const Vec2& position() const {
        return position_;
    }

    const Vec2& velocity() const {
        return velocity_;
    }

    const Vec2& force() const {
        return force_;
    }

    bool isGas() const {
        return isGas_;
    }

    void setMass(double mass) {
        mass_ = mass;
    }

    void setPosition(const Vec2& position) {
        position_ = position;
    }

    void setVelocity(const Vec2& velocity) {
        velocity_ = velocity;
    }

    void setGas(bool isGas) {
        isGas_ = isGas;
    }

    void resetForce() {
        force_ = Vec2(0.0, 0.0, 0.0);
    }

    void addForce(const Vec2& force) {
        force_ += force;
    }

    void resetAfterMerge(const Vec2& position, const Vec2& velocity, double mass) {
        mass_ = mass;
        position_ = position;
        velocity_ = velocity;
        resetForce();
    }

    void dampVelocity(double factor) {
        velocity_ = velocity_ * factor;
    }

    void applyHalfKick(double dt) {
        Vec2 acceleration = force_ / mass_;
        velocity_ += acceleration * (0.5 * dt);
    }

    void drift(double dt) {
        position_ += velocity_ * dt;
    }

private:
    double mass_;
    Vec2 position_;
    Vec2 velocity_;
    Vec2 force_;
    bool isGas_;
};

#endif
