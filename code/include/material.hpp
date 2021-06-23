#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>
#include <iostream>
#include "ray.hpp"
#include "hit.hpp"
#include <vector>
using namespace std;

class Material {
public:

    explicit Material(const Vector3f &color, int t = 0, int tt = 0, vector<Vector3f>* texture = nullptr) :
            Color(color), type(t), texture_type(tt), uv(texture) {
        
    }

    virtual ~Material()
    {
        //delete uv;
    }

    virtual Vector3f getColor() const {
        return Color;
    }

    Vector3f getUV(int u, int v) const {
        return uv[u][v];
    }

    int getType() const {
        return type;
    }

    Vector3f Color;
    int type; // 0: Matte 1: Mirror, 2: Glass
    int texture_type;
    vector<Vector3f>* uv;
};

class TriangleMaterial : public Material{
    
};

class SphereMaterial : public Material {

};


#endif // MATERIAL_H
