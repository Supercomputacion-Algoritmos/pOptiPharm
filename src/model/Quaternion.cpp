/*
 * File:   Quaternion.cpp
 * Author: savins
 *
 * Created on 9 de enero de 2016, 22:00
 */

#include "model/Quaternion.h"

Quaternion::Quaternion() {
  this->w = 0;
  this->x = 0;
  this->y = 0;
  this->z = 0;
}

Quaternion::Quaternion(double w, double x, double y, double z) {
  this->w = w;
  this->x = x;
  this->y = y;
  this->z = z;
}

Quaternion::~Quaternion() {}

Quaternion Quaternion::conjugate(Quaternion &a) {
  Quaternion q;
  q.w = a.w;
  q.x = -a.x;
  q.y = -a.y;
  q.z = -a.z;
  return q;
}

Quaternion Quaternion::point3DToQuaternion(Point3DDouble a) {
  Quaternion q;
  q.w = 0;
  q.x = a.x;
  q.y = a.y;
  q.z = a.z;
  return q;
}

Quaternion Quaternion::XYZToQuaternion(double x, double y, double z) {
  Quaternion q;
  q.w = 0;
  q.x = x;
  q.y = y;
  q.z = z;
  return q;
}

Quaternion operator+(Quaternion &a, Quaternion &b) {
  Quaternion q;
  q.w = a.w + b.w;
  q.x = a.x + b.x;
  q.y = a.y + b.y;
  q.z = a.z + b.z;
  return q;
}

Quaternion operator-(Quaternion &a, Quaternion &b) {
  Quaternion q;
  q.w = a.w - b.w;
  q.x = a.x - b.x;
  q.y = a.y - b.y;
  q.z = a.z - b.z;
  return q;
}

Quaternion operator*(Quaternion &a, Quaternion &b) {
  Quaternion q;
  q.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
  q.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
  q.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
  q.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
  return q;
}

Quaternion Quaternion::result(Quaternion q, Quaternion p) {
  Quaternion result;
  Quaternion punto;
  punto.w = 0;
  punto.x = p.x;
  punto.y = p.y;
  punto.z = p.z;
  result.w = punto.w * q.w - punto.x * q.x - punto.y * q.y - punto.z * q.z;
  result.x = punto.w * q.x + punto.x * q.w - punto.y * q.z + punto.z * q.y;
  result.y = punto.w * q.y + punto.x * q.z + punto.y * q.w - punto.z * q.x;
  result.z = punto.w * q.z - punto.x * q.y + punto.y * q.x + punto.z * q.w;
  return result;
}
