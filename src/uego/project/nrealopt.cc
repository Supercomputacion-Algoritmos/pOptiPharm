#include "uego/uego.h"

////////////////////////////////////////////////////////////
// $Id: nrealopt.cc,v 2.5 1998/03/17 23:14:53 jelasity Exp $
// nrealopt.cc
// NDimRealElement::Optimize() : a stochastic hillclimber
// for real spaces
////////////////////////////////////////////////////////////
// modification history:
//	Jelasity 98 03 10 < -> <= to handle platos in optimize
////////////////////////////////////////////////////////////

long NDimMolElement::Optimize(short radind, long maxeval) {

  const long Scnt = 5, Fcnt = 3;
  const double ct = .5, ex = 2.0,
               rad = INI.R(radind) / INI.R(0) * sqrt(dim), // normalized rad.
      sigmaub = rad, sigmalb = (rad / 1000.0 > 1e-5 ? 1e-5 : rad / 1000.0);

  NDimMolElement *tmp;
  long j, i, scnt, fcnt;
  double sigma, *b, *epsi;

  FailFlag = 1 == 1;

  // --- memory allocations ----------------------------
  tmp = new NDimMolElement(dim);
  b = new double[2 * dim];
  if (tmp == NULL || tmp->Fail() || b == NULL) {
    message("No memory in NDimRealElement::Optimize()", MSG_ERROR);
    if (tmp != NULL)
      delete tmp;
    return 0;
  };
  epsi = b + dim;

  scnt = 0;
  fcnt = 0;
  sigma = sigmaub;

  for (i = 0; i < dim; i++)
    b[i] = 0.0;
  j = 0;

  // while ((j < maxeval) && (fcnt <= 20)) //  && (sigma < sigmalb )
  while ((j < maxeval)) {

    if (scnt > Scnt)
      sigma = ex * sigma;
    if (fcnt > Fcnt)
      sigma = ct * sigma;
    if (sigma < sigmalb || fcnt > 20) // restart search
    {
      // scnt = 0;
      // fcnt = 0;
      sigma = sigmaub;
      for (i = 0; i < dim; i++)
        b[i] = 0.0;
    }
    if (sigma > sigmaub)
      sigma = sigmaub;
    //	printf("\t rad=%f\t sigma=%f, j<maxeval=%d\n" ,rad ,sigma,j );

    // CAMBIO SAVINS ROTACION CIRCULAR
    // POR ESO EL BUCLE COMIENZA EN 2
    epsi[0] = Gauss(b[0], sigma);
    if (epsi[0] > rad || epsi[0] < -rad)
      epsi[0] = fmod(((fmod(epsi[0], rad)) + rad), rad);

    for (i = 1; i < dim; i++) {
      epsi[i] = Gauss(b[i], sigma);
      if (epsi[i] > rad)
        epsi[i] = rad;
      else if (epsi[i] < -rad)
        epsi[i] = -rad;
    };

    // --- try tmp = center + epsi --------------------------
    // printf("Suma\n");
    tmp->UpdateFrom(this);
    tmp->Add(epsi);
    tmp->UpdateValue();
    j++;

    // printf("CVB=%lf  CVA=%lf",  tmp->CurrValue(), CurrValue() );
    if (tmp->CurrValue() >= CurrValue()) {
      UpdateFrom(tmp);
      for (i = 0; i < dim; ++i)
        b[i] = 0.4 * epsi[i] + 0.2 * b[i];
      ++scnt;
      fcnt = 0;
      if (j >= maxeval)
        break;
    } else {
      // printf("Resta\n");
      //  --- try tmp = center - epsi --------------------
      tmp->UpdateFrom(this);
      tmp->Add(epsi, -1);
      tmp->UpdateValue();
      j++;
      // printf("CVB=%lf  CVA=%lf",  tmp->CurrValue(), CurrValue() );
      if (tmp->CurrValue() >= CurrValue()) {
        UpdateFrom(tmp);
        for (i = 0; i < dim; ++i)
          b[i] = b[i] - 0.4 * epsi[i];
        ++scnt;
        fcnt = 0;
      } else {
        // --- center is still better than tmp -------------
        for (i = 0; i < dim; ++i)
          b[i] = 0.5 * b[i];
        fcnt = fcnt + 1;
        scnt = 0;
      }
    }
  };

  delete tmp;
  delete[] b;
  FailFlag = 1 == 0;

  // printf("rad=%f, spec. level=%d, master level=%d, INI.R(0)=%f, maxeval=%d, "
  //        "j=%d\n",
  //        rad, radind, INI.getLevelMaster(), INI.R(0), maxeval, j);
  return j;
  // return maxeval;
};
