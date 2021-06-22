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