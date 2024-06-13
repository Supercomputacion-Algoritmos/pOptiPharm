#include "uego/uego.h"

////////////////////////////////////////////////////////////
// $Id: sspopt.cc,v 2.5 1998/03/17 23:14:54 jelasity Exp $
// sspopt.cc
// definition of SearchSpElement::Optimize()
// abstract implementation of the stochastic hillclimber
////////////////////////////////////////////////////////////
// modification history:
//	Jelasity 98 01 10: RandNew changed to MutateNew
////////////////////////////////////////////////////////////

long SearchSpElement::Optimize(short radind, long maxeval) {

  SearchSpElement *spe;

  FailFlag = 1 == 1;

  for (long i = 0; i < maxeval; ++i) {
    spe = MutateNew(radind);
    if (Fail())
      return -1;
    else
      FailFlag = 1 == 1;

    if (spe->CurrValue() >= CurrValue())
      UpdateFrom(spe);
    delete spe;
  };

  FailFlag = 1 == 0;

  return maxeval;
};
