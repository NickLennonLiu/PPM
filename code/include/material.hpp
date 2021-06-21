#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>
#include <iostream>
#include "ray.hpp"
#include "hit.hpp"

class Material {
public:

    explicit Material(const Vector3f &color, int t = 0) :
            Color(color), type(t) {

    }

    virtual ~Material() = default;

    virtual Vector3f getColor() const {
        return Color;
    }

    int getType() const {
        return type;
    }

    Vector3f Color;
    int type; // 0: Matte 1: Mirror, 2: Glass
};


#endif // MATERIAL_H
