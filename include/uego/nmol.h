#include <cstdint>
#ifndef NREAL_H
#define NREAL_H

#include "functions/VolumeOverlap.h"
#include "searchsp.h"
////////////////////////////////////////////////////////////
// $Id: nreal.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// nreal.h
// declaration for n dimensional real search spaces
////////////////////////////////////////////////////////////
// modification history:
//	Jelasity 98 02 15 inner repr. is (0,1)^dim
////////////////////////////////////////////////////////////

//--------------------------------------------------------------------------
//---- base class for n dim real search spaces
//--------------------------------------------------------------------------
// the inner representation works on a (0,1)^n cube, but Value() sees
// the real user-given values in x. This is implemented through UpdateValue().
// Value() is never called directly. The inner representation does not affect
// the behaviour of public functions.
//--------------------------------------------------------------------------

class NDimMolElement : public SearchSpElement {

protected:
  virtual double Value();
  virtual double PreciseValue(); // No optimization in objective function
  double Gauss(double, double);
  void Add(double *, signed char = 1);

  static long dim;
  double *x; // the space element

public:
  NDimMolElement(long); // creates an uninitialized 'empty' element
  virtual ~NDimMolElement() {
    if (x != NULL)
      delete[] x;
  };

  virtual double UpdateValue();
  virtual double PreciseUpdateValue();

  virtual double Distance(SearchSpElement *, SearchSpElement *s = NULL);

  // ----- implementations of SearchSpElement virtual functions
  virtual SearchSpElement *RandNew();

  virtual SearchSpElement *generateMutationNew(short,
                                               SearchSpElement *s = NULL);

  virtual SearchSpElement *RandNew(short, SearchSpElement *s = NULL);

  virtual SearchSpElement *MutateNew(short, SearchSpElement *s = NULL);
  virtual SearchSpElement *BetweenNew(SearchSpElement *,
                                      SearchSpElement *s = NULL);

  virtual SearchSpElement *MutateBetweenNew(SearchSpElement *,
                                            SearchSpElement *s = NULL);

  virtual long Optimize(short, long); // returns funct. evals

  virtual void UpdateFrom(SearchSpElement *);

  virtual double Diameter(Ini *);
  virtual double v(double);

  virtual void Save(FILE *);
  virtual void Save(FILE *, bool);
  virtual void GetX(double *y) {
    for (int k = 0; k < dim; k++)
      y[k] = x[k];
  };
  virtual void SetX(double *y) {
    for (int k = 0; k < dim; k++)
      x[k] = (y[k] - INI.Lowb(k)) / (INI.Upb(k) - INI.Lowb(k));
  }
};

#endif
