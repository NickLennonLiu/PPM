#ifndef LIGHT_H
#define LIGHT_H

#include <Vector3f.h>
#include "object3d.hpp"
#include "halton.hpp"
#include "ray.hpp"

class Light {
public:
    Light() = default;

    virtual ~Light() = default;

    virtual void genp(Ray &pr, Vector3f *f, int i) = 0;
    
    Vector3f color;
    Vector3f flux;

};


class DirectionalLight : public Light {
    Vector3f direction;
    Vector3f position;
    Vector3f axis1, axis2;
public:
    DirectionalLight() = delete;

    DirectionalLight(const Vector3f &d, const Vector3f &c) {
        direction = d.normalized();
        color = c;
    }

    ~DirectionalLight() override = default;

    virtual void genp(Ray &pr, Vector3f *f, int i) override
    {
        (*f) = flux * (D_PI * 4.0); // flux
        auto p = 2.0 * D_PI * halton(0, i);
        auto t = 2.0 * acos(sqrt(1. - halton(1, i)));
        auto st = sin(t);
        pr = Ray(position, Vector3f(cos(p) * st, cos(t), sin(p) * st));
    }
};

class PointLight : public Light {
    Vector3f position;
public:
    PointLight() = delete;

    PointLight(const Vector3f &p, const Vector3f &c, const Vector3f& f = Vector3f({100, 100, 100})) {
        position = p;
        color = c;
        flux = f;
    }

    ~PointLight() override = default;

    virtual void genp(Ray &pr, Vector3f *f, int i) override
    {
        (*f) = flux * (D_PI * 4.0); // flux
        auto p = 2.0 * D_PI * halton(0, i);
        auto t = 2.0 * acos(sqrt(1. - halton(1, i)));
        auto st = sin(t);
        pr = Ray(position, Vector3f(cos(p) * st, cos(t), sin(p) * st));
    }
};

#endif // LIGHT_H
