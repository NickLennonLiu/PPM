#ifndef KDNODE_H
#define KDNODE_H

#include "object3d.hpp"
#include "hit.hpp"
#include <vector>

class KdNode
{

public:
    KdNode *lc, *rc;
    AABB bbox;
    KdNode(AABB b, KdNode*lc = nullptr, KdNode* rc = nullptr)
    : lc(lc), rc(rc), bbox(b) {}

    bool intersect(const Ray &r, Hit &h, float tmin)
    {
        Hit tmp;
        Hit th;
        if(!bbox.intersect(r, th, tmin))
            return false;
        // 根节点
        if(lc == nullptr && rc == nullptr)
        {
            return bbox.content->intersect(r,h,tmin);
        }
        bool re = false;
        th = tmp;
        if(lc->bbox.intersect(r, th, tmin))
            re |= lc->intersect(r, h, tmin);
        th = tmp;
        if(rc->bbox.intersect(r, th, tmin))
            re |= rc->intersect(r, h, tmin);
        return re;
    }

    static bool cmpX(const AABB &a, const AABB &b)
    {
        if (a.axis_planes[0][1] < b.axis_planes[0][1])
            return true;
        if (a.axis_planes[0][1] > b.axis_planes[0][1])
            return false;
        return a.axis_planes[0][0] < b.axis_planes[0][0];
    }

    static bool cmpY(const AABB &a, const AABB &b)
    {
        if (a.axis_planes[1][1] < b.axis_planes[1][1])
            return true;
        if (a.axis_planes[1][1] > b.axis_planes[1][1])
            return false;
        return a.axis_planes[1][0] < b.axis_planes[1][0];
    }

    static bool cmpZ(const AABB &a, const AABB &b)
    {
        if (a.axis_planes[2][1] < b.axis_planes[2][1])
            return true;
        if (a.axis_planes[2][1] > b.axis_planes[2][1])
            return false;
        return a.axis_planes[2][0] < b.axis_planes[2][0];
    }

    static KdNode *split(int start, int end, int dir, std::vector<AABB> &bboxs)
    {
        if (end == start)
            return nullptr;
        if (end == start + 1)
            return new KdNode(bboxs[start]);
        if (dir == 0)
            sort(bboxs.begin() + start, bboxs.begin() + end, cmpX);
        else if (dir == 1)
            sort(bboxs.begin() + start, bboxs.begin() + end, cmpY);
        else
            sort(bboxs.begin() + start, bboxs.begin() + end, cmpZ);
        KdNode *lc = split(start, start + (end - start) / 2, (dir + 1) % 3, bboxs),
               *rc = split(start + (end - start) / 2, end, (dir + 1) % 3, bboxs);
        AABB box = AABB::merge(lc->bbox, rc->bbox);
        std::cout << box;
        return new KdNode(box, lc, rc);
    }
};

#endif