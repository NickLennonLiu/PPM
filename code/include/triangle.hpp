#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

using namespace std;

class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m) 
	: Object3D(m) 
	{
		vertices[0] = a;
		vertices[1] = b;
		vertices[2] = c;
		normal = (Vector3f::cross(a - b, a - c)).normalized();
	}

	bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
		Vector3f solve = Matrix3f(ray.getDirection(), 
									vertices[0] - vertices[1], 
									vertices[0] - vertices[2], true).inverse()
						* (vertices[0] - ray.getOrigin()) ;
		float t = solve[0], beta = solve[1], gamma = solve[2];
		if (0 <= beta && 0 <= gamma && (beta + gamma) <= 1)
			if(t < hit.getT() && t > tmin)
			{
				Vector3f n = Vector3f::dot(ray.getDirection(), normal) > 0 ? -normal : normal;
				hit.set(t, material, n);
				return true;
			}
        return false;
	}
	Vector3f normal;
	Vector3f vertices[3];
protected:
	
};

#endif //TRIANGLE_H
