#ifndef GROUP_H
#define GROUP_H


#include "object3d.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include <vector>


class Group : public Object3D {

public:

    Group() {

    }

    explicit Group (int num_objects) {
        this->num_objects = num_objects;
    }

    ~Group() override {    
    }

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        bool result = false;
        for(int i = 0; i < num_objects; ++i)
        {
            result |= objs[i]->intersect(r, h, tmin);
        }
        return result;
    }

    void addObject(int index, Object3D *obj) {
        objs.insert(objs.begin() + index, obj);
    }

    int getGroupSize() {
        return num_objects;
    }

private:

    std::vector<Object3D*> objs;    
    int num_objects;
};

#endif
	
