#ifndef AABB_HPP
#define AABB_HPP

#include "object3d.hpp"
#include <vecmath.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <iostream>
#include "plane.hpp"
using namespace std;


class AABB : public Object3D
{
    static float dir[3][3];
    // 3个平面确定一个包围盒
    float axis_planes[3][2]; // 0:x, 1:y, 2:z // 0: min, 1: max
public:
    AABB()
    {

    }
    AABB(std::vector<std::pair<float, float>> axis, Material *material = nullptr)
    : Object3D(material)
    {
        axis_planes[0][0] = axis[0].first, axis_planes[0][1] = axis[0].second;
        axis_planes[1][0] = axis[1].first, axis_planes[1][1] = axis[1].second;
        axis_planes[2][0] = axis[2].first, axis_planes[2][1] = axis[2].second;
    }
    bool intersect(const Ray &r, Hit &h, float tmin) override
    {
        //cout << "Start intersecting" << endl;
        Vector3f o = r.getOrigin();
        float tm = 1e38;
        int idx = 0;
        bool re = false;
        Hit ans;
        for(int i = 0; i < 3; ++i)
        {
            int close = (fabsf(axis_planes[i][0] - o[i]) > fabsf(axis_planes[i][1] - o[i])) ? 1 : 0;
            Hit ht;
            float temp = get_axis_plane_t(i, close, r, ht);
            Vector3f inte = r.pointAtParameter(ht.getT());
            bool inside = true;
            for (int j = 0; j < 3; ++j)
            {
                if (j == i)
                    continue;
                if ((inte[j] > axis_planes[j][1]) || (inte[j] < axis_planes[j][0]))
                {
                        
                    inside = false;
                    break;
                }
            }

            if((temp < tm) && inside)
            {
                tm = temp;
                idx = i;
                ans = ht;
                re = true;
            }
        }
        
        if(ans.getT() > h.getT() || ans.getT() < tmin)
            return false;
        h = ans;
        //h.set(h.getT(), material, h.getNormal());
        return re;
    }

    void debug()
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 2; ++j)
                cout << axis_planes[i][j] << " ";
            cout << endl;
        }
    }

protected:
    float get_axis_plane_t(int axis, int direction, const Ray& r, Hit &h)
    {
        Vector3f normal(dir[axis][0], dir[axis][1], dir[axis][2]);
        float D = axis_planes[axis][direction];
        Plane p(normal, D, material);
        p.intersect(r, h, 0);
        return h.getT();
    }
};

float AABB::dir[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

#endif // AABB_HPP