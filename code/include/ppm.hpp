//-------------------------------------------------------------------------------------------
// File : main.cpp
// Desc : expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)
//        Original Code by T.Hachisuka (https://cs.uwaterloo.ca/~thachisu/smallppm_exp.cpp)
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
// Includes
//-------------------------------------------------------------------------------------------
#include <cmath>
#include <cstdlib>
#include <cstdio>
//#include <bitmap.h>
//#include <math_util.h>
#include <list>
#include <vector>
#include <chrono>
//#include <sphere.h>
//#include <hitrecord.h>
#include "hitrecord.hpp"
//#include "aabb.hpp"
#include "sphere.hpp"
#include "bbox.hpp"
#include "material.hpp"
#include "halton.hpp"
#include "group.hpp"
#include "image.hpp"

#define PHOTON_COUNT_MUTIPLIER 1000

namespace /* anonymous */
{

    //-------------------------------------------------------------------------------------------
    // Constant Values
    //-------------------------------------------------------------------------------------------
    static const double ALPHA = 0.7; // the alpha parameter of PPM

    //-------------------------------------------------------------------------------------------
    // Global Variables
    //-------------------------------------------------------------------------------------------
    std::list<HitRecord *> hitpoints;
    std::vector<std::list<HitRecord *>> hash_grid;
    double hash_s;
    BoundingBox hpbbox;
    Group* group;
    SceneParser* parser;
    Image* _img;
    int _w;
    int unreach = 0;
    /*
    Sphere sph[] =
        {
            // Scene: radius, position, color, material
            Sphere(Vector3f(1e5 + 1, 40.8, 81.6), 1e5, Material::Matte),   //Right
            Sphere(Vector3f(-1e5 + 99, 40.8, 81.6), 1e5, Material::Matte), //Left
            Sphere(Vector3f(50, 40.8, 1e5), 1e5, Material::Matte),         //Back
            Sphere(Vector3f(50, 40.8, -1e5 + 170), 1e5, Material::Matte),     //Front
            Sphere(Vector3f(50, 1e5, 81.6), 1e5, Material::Matte),         //Bottomm
            Sphere(Vector3f(50, -1e5 + 81.6, 81.6), 1e5, Material::Matte), //Top
            Sphere(Vector3f(27, 16.5, 47), 16.5, Material::Mirror),        //Mirror
            Sphere(Vector3f(73, 16.5, 88), 16.5, Material::Glass),         //Glass
            Sphere(Vector3f(50, 8.5, 60), 8.5, Material::Matte),           //Middle
    };
    */

} // namespace

//-------------------------------------------------------------------------------------------
//      ハッシュキーを取得します.
//-------------------------------------------------------------------------------------------
inline unsigned int get_hash(
    const int ix, const int iy, const int iz)
{
    return (unsigned int)((ix * 73856093) ^
                          (iy * 19349663) ^
                          (iz * 83492791)) %
           hash_grid.size();
}

//-------------------------------------------------------------------------------------------
//      ハッシュグリッドを構築します.
//-------------------------------------------------------------------------------------------
void build_hash_grid(
    const int w, const int h)
{
    // heuristic for initial radius
    auto size = hpbbox.maxi - hpbbox.mini;
    auto irad = ((size.x() + size.y() + size.z()) / 3.0) / ((w + h) / 2.0) * 4.0;

    // determine hash table size
    // we now find the bounding box of all the measurement points inflated by the initial radius
    hpbbox.reset();
    auto photon_count = 0;
    for (auto itr = hitpoints.begin(); itr != hitpoints.end(); ++itr)
    {
        auto hp = (*itr);
        hp->r2 = irad * irad;
        hp->n = 0;
        hp->flux = Vector3f();

        photon_count++;
        hpbbox.merge(hp->pos - irad);
        hpbbox.merge(hp->pos + irad);
    }

    // make each grid cell two times larger than the initial radius
    hash_s = 1.0 / (irad * 2.0);

    // build the hash table
    hash_grid.resize(photon_count * 3.5);
    cout << hash_grid.size() << " --- size" << endl;
    hash_grid.shrink_to_fit();
    for (auto itr = hitpoints.begin(); itr != hitpoints.end(); ++itr)
    {
        auto hp = (*itr);
        auto min = ((hp->pos - irad) - hpbbox.mini) * hash_s;
        auto max = ((hp->pos + irad) - hpbbox.mini) * hash_s;

        auto minX = abs(int(min.x()));
        auto maxX = abs(int(max.x()));

        auto minY = abs(int(min.y()));
        auto maxY = abs(int(max.y()));

        auto minZ = abs(int(min.z()));
        auto maxZ = abs(int(max.z()));

        for (int iz = minZ; iz <= maxZ; iz++)
        {
            for (int iy = minY; iy <= maxY; iy++)
            {
                for (int ix = minX; ix <= maxX; ix++)
                {
                    int hv = get_hash(ix, iy, iz);
                    hash_grid[hv].push_back(hp);
                    
                }
            }
        }
    }
}

//-------------------------------------------------------------------------------------------
//      交差判定を行います.
//-------------------------------------------------------------------------------------------
/*   使用框架中的相交
inline bool intersect(const Ray &r, double &t, int &id)
{

    int n = sizeof(sph) / sizeof(sph[0]);
    auto d = D_INF;
    t = D_INF;
    for (int i = 0; i < n; i++)
    {
        d = sph[i].intersect(r, Hit(), 0);
        if (d < t)
        {
            t = d;
            id = i;
        }
    }

    return (t < D_INF);
}
*/

//-------------------------------------------------------------------------------------------
//     生成Photon Ray
//-------------------------------------------------------------------------------------------
/*
void genp(Ray &pr, Vector3f *f, int i, Light* light)
{
    // generate a photon ray from the point light source with QMC
    (*f) = Vector3f(100, 100, 100) * (D_PI * 4.0); // flux
    auto p = 2.0 * D_PI * halton(0, i);
    auto t = 2.0 * acos(sqrt(1. - halton(1, i)));
    auto st = sin(t);
    pr = Ray(light->position, Vector3f(cos(p) * st, cos(t), sin(p) * st));
}
*/

//-------------------------------------------------------------------------------------------
//      光线跟踪
//-------------------------------------------------------------------------------------------
void trace(const Ray &r, int dpt, bool m, const Vector3f &fl, const Vector3f &adj, int i, int depth)
{
    double t;
    int id;
    Hit h;
    dpt++;
    if (!group->intersect(r, h, 1e-1) || (dpt >= 20))
        return;

    auto d3 = dpt * 3;

    //auto x = r.getOrigin() + r.getDirection() * t, n = (x - obj.pos).normal;    
    auto x = r.pointAtParameter(h.getT()), n = h.getNormal();
    auto material = h.getMaterial();
    auto f = material->getColor();
    auto nl = ((Vector3f::dot(r.getDirection(), n)) < 0) ? n : n * -1;
    auto p = (f.x() > f.y() && f.x() > f.z()) ? f.x() : (f.y() > f.z()) ? f.y() : f.z();

    if (material->getType() == 0) // Matte
    {
        if (m)
        {
            // eye ray
            // store the measurment point
            auto hp = new HitRecord;
            hp->f = f * adj;
            hp->pos = x;
            hp->nrm = n;
            hp->idx = i;
            hitpoints.push_back(hp);

            // find the bounding box of all the measurement points
            hpbbox.merge(x);
        }
        else
        {
            // photon ray
            // find neighboring measurement points and accumulate flux via progressive density estimation
            auto hh = (x - hpbbox.mini) * hash_s;
            auto ix = abs(int(hh.x()));
            auto iy = abs(int(hh.y()));
            auto iz = abs(int(hh.z()));
            // strictly speaking, we should use #pragma omp critical here.
            // it usually works without an artifact due to the fact that photons are
            // rarely accumulated to the same measurement points at the same time (especially with QMC).
            // it is also significantly faster.
            {
                auto list = hash_grid[get_hash(ix, iy, iz)];
                for (auto itr = list.begin(); itr != list.end(); itr++)
                {
                    auto hp = (*itr);
                    auto v = hp->pos - x;
                    // check normals to be closer than 90 degree (avoids some edge brightning)
                    if ((Vector3f::dot(hp->nrm, n) > 1e-3) && (Vector3f::dot(v, v) <= (hp->r2)))
                    {
                        // unlike N in the paper, hp->n stores "N / ALPHA" to make it an integer value
                        auto g = (hp->n * ALPHA + ALPHA) / (hp->n * ALPHA + 1.0);
                        hp->r2 = hp->r2 * g;
                        hp->n++;
                        hp->flux = (hp->flux + (hp->f *  fl) / D_PI) * g;
                    }
                }
            }

            // use QMC to sample the next direction
            auto r1 = 2.0 * D_PI * halton(d3 - 1, i);
            auto r2 = halton(d3 + 0, i);
            auto r2s = sqrt(r2);
            auto w = nl;
            auto u = (Vector3f::cross((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w)).normalized();
            auto v = Vector3f::cross(w, u);
            auto d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalized();

            if (halton(d3 + 1, i) < p)
                trace(Ray(x, d), dpt, m, (f * fl) * (1. / p), (f * adj), i, depth+1);
        }
    }
    else if (material->getType() == 1) // Mirror
    {
        trace(Ray(x, Vector3f::reflect(r.getDirection(), n)), dpt, m, (f * fl), (f * adj), i, depth+1);
    }
    else    // Glass
    {
        Ray lr(x, Vector3f::reflect(r.getDirection(), n));
        auto into = Vector3f::dot(n, nl) > 0.0;
        auto nc = 1.0;
        auto nt = 1.5;
        auto nnt = (into) ? nc / nt : nt / nc;
        auto ddn = Vector3f::dot(r.getDirection(), nl);
        auto cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

        // total internal reflection
        if (cos2t < 0)
            return trace(lr, dpt, m, (f * fl), (f * adj), i, depth+1);

        auto td = (r.getDirection() * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
        auto a = nt - nc;
        auto b = nt + nc;
        auto R0 = a * a / (b * b);
        auto c = 1 - (into ? -ddn : Vector3f::dot(td, n));
        auto Re = R0 + (1 - R0) * c * c * c * c * c;
        auto P = Re;
        Ray rr(x, td);
        auto fa = f * adj;
        auto ffl = f * fl;

        if (m)
        {
            // eye ray (trace both rays)
            trace(lr, dpt, m, ffl, fa * Re, i, depth+1);
            trace(rr, dpt, m, ffl, fa * (1.0 - Re), i, depth+1);
        }
        else
        {
            // photon ray (pick one via Russian roulette)
            (halton(d3 - 1, i) < P)
                ? trace(lr, dpt, m, ffl, fa * Re, i, depth+1)
                : trace(rr, dpt, m, ffl, fa * (1.0 - Re), i, depth+1);
        }
    }
}

//-------------------------------------------------------------------------------------------
//      eye rayを追跡します.
//-------------------------------------------------------------------------------------------
void trace_ray(int w, int h, Camera* camera)
{
    auto start = std::chrono::system_clock::now();
    hitpoints.clear();
    // trace eye rays and store measurement points
    cout << w * h << endl;
    for (int y = 0; y < h; y++)
    {
        std::fprintf(stdout, "\rHitPointPass %5.2f%%", 100.0 * y / (h - 1));
        for (int x = 0; x < w; x++)
        {
            auto idx = x + y * w;   
            trace(camera->generateRay({x, y}), 0, true, Vector3f(), Vector3f(1, 1, 1), idx, 0);
        }
    }
    std::fprintf(stdout, "\n");

    auto end = std::chrono::system_clock::now();
    auto dif = end - start;
    std::fprintf(stdout, "Ray Tracing Pass : %lld(msec)\n", std::chrono::duration_cast<std::chrono::milliseconds>(dif).count());

    start = std::chrono::system_clock::now();

    // build the hash table over the measurement points
    build_hash_grid(w, h);

    end = std::chrono::system_clock::now();
    dif = end - start;
    std::fprintf(stdout, "Build Hash Grid: %lld(msec)\n", std::chrono::duration_cast<std::chrono::milliseconds>(dif).count());
}

//-------------------------------------------------------------------------------------------
//      photon rayを追跡します.
//-------------------------------------------------------------------------------------------
void trace_photon(int s)
{
    auto start = std::chrono::system_clock::now();

    // trace photon rays with multi-threading
    auto vw = Vector3f(1, 1, 1);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < s; i++)
    {
        auto p = 100.0 * (i + 1) / s;
        std::fprintf(stdout, "\rPhotonPass %5.2f%%", p);
        int m = PHOTON_COUNT_MUTIPLIER * i;
        Ray r({0,0,0}, {0,0,0});
        Vector3f f;

        for (int j = 0; j < PHOTON_COUNT_MUTIPLIER; j++)
        {
            for (int li = 0; li < parser->getNumLights(); ++li)
            {
                //genp(r, &f, m + j, parser->getLight(li));
                parser->getLight(li)->genp(r, &f, m+j);
                trace(r, 0, false, f, vw, m + j, 0);
            }
            
        }
    }

    std::fprintf(stdout, "\n");

    auto end = std::chrono::system_clock::now();
    auto dif = end - start;
    std::fprintf(stdout, "Photon Tracing Pass : %lld(sec)\n", std::chrono::duration_cast<std::chrono::seconds>(dif).count());
}

//-------------------------------------------------------------------------------------------
//      密度推定を行います.
//-------------------------------------------------------------------------------------------
void density_estimation(Vector3f *color, int num_photon)
{
    // density estimation
    for (auto itr = hitpoints.begin(); itr != hitpoints.end(); ++itr)
    {
        auto hp = (*itr);
        auto i = hp->idx;
        color[i] = color[i] + hp->flux * (1.0 / (D_PI * hp->r2 * num_photon * PHOTON_COUNT_MUTIPLIER));
    }
}

//-------------------------------------------------------------------------------------------
//      メインエントリーポイントです.
//-------------------------------------------------------------------------------------------
int ppm(int w, int h, int s, Image* img, SceneParser* _parser)
{
    //auto w = 1280;  // 画像の横幅.
    //auto h = 1080;  // 画像の縦幅.
    //auto s = 10000; // s * 1000 photon paths will be traced (s * PHOTON_COUNT_MULTIPLIER).
    auto c = new Vector3f[w * h];
    group = _parser->getGroup();
    group->getKdTree();
    hpbbox.reset();
    parser = _parser;
    _img = img;
    _w = w;
    trace_ray(w, h, _parser->getCamera());
    trace_photon(s);
    density_estimation(c, s);

    const auto kGamma = 2.2;
    // TODO(edit Image.SetPixel() to have gamma argument)
    //save_to_bmp("image.bmp", w, h, &c[0].x, kGamma);

    for(int i = 0; i < (w*h); ++i)    // idx = x + y * w
        img->SetPixel(i % w, i / w, c[i]);

    cout << unreach << endl;
    delete[] c;
    c = nullptr;

    return 0;
}
