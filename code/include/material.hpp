#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>
#include <iostream>
#include "ray.hpp"
#include "hit.hpp"

class Material {
public:

    explicit Material(const Vector3f &d_color, const Vector3f &s_color = Vector3f::ZERO, float s = 0, int t = 0) :
            diffuseColor(d_color), specularColor(s_color), shininess(s), type(t) {

    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    float getShininess() const {
        return shininess;
    }


    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;
        Vector3f Ri = 2.0 * Vector3f::dot(dirToLight, hit.getNormal()) * hit.getNormal() - dirToLight;
        
        shaded += lightColor * 
                (diffuseColor * clamp(Vector3f::dot(dirToLight, hit.getNormal())) + 
                 specularColor * pow(clamp(Vector3f::dot(-(ray.getDirection().normalized()), Ri.normalized())), shininess));
        return shaded;
    }

    int type; // 0: Matte 1: Mirror, 2: Glass

protected:
    Vector3f diffuseColor;
    Vector3f specularColor;
    float shininess;
    

    inline float clamp(float x)
    {
        return (x > 0) ? x : 0; 
    }
};


#endif // MATERIAL_H
