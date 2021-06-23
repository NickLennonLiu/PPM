#include "curve.hpp"
#include <vector>
#include "object3d.hpp"
#include "kdnode.hpp"
#include "triangle.hpp"

struct SurfacePoint {
    Vector3f V;
    Vector3f uT, vT;
    Vector3f N;
};

class Surface : public Object3D
{
protected:
    int m, n; // u - m, v - n
    std::vector<std::vector<Vector3f>> controls;
    std::vector<std::vector<SurfacePoint>> samplePoints;
    KdNode* root;
    std::vector<Triangle> triangles;
    AABB box;
    typedef std::tuple<unsigned, unsigned, unsigned> Tup3u;
    // 表面所存储的顶点、顶点对应的法向量、
    struct Face {
        std::vector<Vector3f> VV;
        std::vector<Vector3f> VN;
        std::vector<Tup3u> VF;
    } faces;


public:
    explicit Surface(std::vector<std::vector<Vector3f>> points, Material* material) 
    : Object3D(material)
    , controls(std::move(points)) {
        m = controls.size();
        n = controls[0].size();
    }

    virtual bool intersect(const Ray&r, Hit &h, float tmin) override {
        return false;
    }

    virtual void discretize(int resolution, std::vector<std::vector<SurfacePoint>>& data) = 0;


};

class BezierSurface : public Surface {
public:
    explicit BezierSurface(const std::vector<std::vector<Vector3f>> &points, Material* material) : Surface(points, material) {
        if ((n < 4 || n % 3 != 1) || (m < 4 || m % 3 != 1)) {
            printf("Number of control points of axis of BezierSurface must be 3n+1!\n");
            exit(0);
        }
        discretize(30, samplePoints); // 在初始化的时候即计算离散点
        Meshify();
    }

    void Meshify()
    {
        int row = samplePoints.size();
        int col = samplePoints[0].size();
        //cout << row << " " << col << endl;
        triangles.clear();
        for(int i = 0; i < row-1; ++i)
        {
            for(int j = 0; j < col-1; ++j)
            {
                Triangle face1(samplePoints[i][j].V, 
                                samplePoints[i+1][j].V,
                                samplePoints[i][j+1].V,
                                material, 
                                samplePoints[i][j].N,
                                samplePoints[i+1][j].N,
                                samplePoints[i][j+1].N);
                Triangle face2(samplePoints[i+1][j+1].V, 
                                samplePoints[i+1][j].V,
                                samplePoints[i][j+1].V,
                                material, 
                                samplePoints[i+1][j+1].N,
                                samplePoints[i+1][j].N,
                                samplePoints[i][j+1].N);
                triangles.push_back(face1);
                triangles.push_back(face2);
            }
        }
        std::vector<AABB> bboxs;
        for(auto &i: triangles)
        {
            AABB tbox = i.bbox();
            box = AABB::merge(box, tbox);
            bboxs.push_back(tbox);
        }
        root = KdNode::split(0, bboxs.size(), 0, bboxs);
        box.content = this;
    }

    void discretize(int resolution, std::vector<std::vector<SurfacePoint>> & data) override {
        data.clear();
        std::vector<double> uknots = Bernstein::bezier_knot(m-1),
                            vknots = Bernstein::bezier_knot(n-1);
        Bernstein ubezier = Bernstein(m-1, m-1, uknots),
                  vbezier = Bernstein(n-1, n-1, vknots);
        double ustart = ubezier.vstart(), uend = ubezier.vend(),
            vstart = vbezier.vstart(), vend = vbezier.vend();
        double ustep = (uend - ustart) / resolution,
                vstep = (vend - vstart) / resolution;
        std::vector<double> us, uds, vs, vds;
        us.clear();
        uds.clear();
        vs.clear();
        vds.clear();
        //while(ustart < uend)
        for(float u = ustart; u < uend; u += ustep)
        {
            std::vector<SurfacePoint> row;
            int ulsk = ubezier.evalute(u, us, uds);
            for(float v = vstart; v < vend; v+= vstep)
            //while(vstart < vend)
            {
                int vlsk = vbezier.evalute(v, vs, vds);                   
                SurfacePoint p;
                for(int i = 0; i < us.size(); ++i)
                    for(int j = 0; j < vs.size(); ++j)
                    {
                        Vector3f point = controls[ulsk+i][vlsk+j];
                        p.V += point * us[i] * vs[j];
                        p.uT += point * uds[i] * vs[j];
                        p.vT += point * us[i] * vds[j]; 
                        p.N = Vector3f::cross(p.uT, p.vT).normalized();
                        //cout << point <<" "<< us[i]<<" " << uds[i]<<" " << vs[i]<<" " << vds[i] << endl;
                    }
                row.push_back(p);   
            }
            data.push_back(row);
        }
    }

    AABB bbox() override
    {
        return box;
    }

    ~BezierSurface() override {

    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        return root->intersect(r, h, tmin);
    }
};