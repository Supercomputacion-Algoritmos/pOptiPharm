
/*
 * File:   Point3D.cpp
 * Author: savins
 *
 * Created on 4 de enero de 2016, 23:31
 */

#include "model/Point3DDouble.h"

Point3DDouble::Point3DDouble() {
  this->x = 0;
  this->y = 0;
  this->z = 0;
}

Point3DDouble::Point3DDouble(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

Point3DDouble::~Point3DDouble() {}
string Point3DDouble::toString() {
  std::stringstream ss;
  ss << "(" << this->x << ", " << this->y << ", " << this->z << ")";
  std::string s = ss.str();
  return s;
}

Point3DDouble operator-(Point3DDouble &a, Point3DDouble &b) {
  Point3DDouble point;
  point.x = a.x - b.x;
  point.y = a.y - b.y;
  point.z = a.z - b.z;
  return point;
}
