/*
 * File:   Box.h
 * Author: savins
 *
 * Created on 4 de enero de 2016, 23:30
 */

#ifndef BOX_H
#define BOX_H
#include "Point3DDouble.h"
#include "Point3DInteger.h"

class Box {
public:
    Box();
    Point3DDouble lowestPoint;
    Point3DInteger extent;
    double res;
    float* points;
    virtual ~Box();
private:

};

#endif /* BOX_H */

