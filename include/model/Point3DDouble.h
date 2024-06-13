/*
 * File:   Point3D.h
 * Author: savins
 *
 * Created on 4 de enero de 2016, 23:31
 */

#ifndef POINT3DDOUBLE_H
#define POINT3DDOUBLE_H

#include <string>
#include <sstream>
using namespace std;
class Point3DDouble {
public:
    Point3DDouble();
    Point3DDouble(double x, double y, double z);
    virtual ~Point3DDouble();
    double x;
    double y; 
    double z;
    string toString();

private:

};

Point3DDouble operator-(Point3DDouble& a, Point3DDouble& b);
    

#endif /* POINT3D_H */

