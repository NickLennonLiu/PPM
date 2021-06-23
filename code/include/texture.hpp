#include "FractalNoise.hpp"
#include "Perlin.hpp"
#include "material.hpp"
#include <cmath>

float cut(float x){
    return x - (int)x;
}

//https://www.researchgate.net/publication/340769728_REALISTIC_RENDERING_OF_WOOD
Vector3f wood(const Vector3f& x, Perlin* gen)
{
    float g = fabs(sin(x.x() + gen->noise(x.x(), x.y(), x.z()))) * 10;
    return cut(g) * Vector3f(0.345, 0.188, 0.074);
}

Vector3f marble(const Vector3f& x, Perlin* gen)
{
    return (cos(x.x() + 30 *gen->noise(x.x(), x.y(), x.z())) + 1)/2 * Vector3f(0.902, 0.894, 0.847);
}

Vector3f rainbow(const Vector3f& x, Perlin* gen)
{
    float z = abs(x.z()*5);
    float zz = z - floor(z);
    switch((int)z % 6)
    {
        case 0: return Vector3f(1,0, 1 - zz);
        case 1: return Vector3f(1,zz, 0);
        case 2: return Vector3f(1-zz, 1, 0);
        case 3: return Vector3f(0, 1, zz);
        case 4: return Vector3f(0, 1 - zz, 1);
        case 5: return Vector3f(zz, 0, 1);
        default: return Vector3f(0,0,0);
    }
}