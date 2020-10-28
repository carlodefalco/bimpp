// Copyright (C) 2020  Carlo de Falco
//
// This file is part of bimpp
//
//  bimpp is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  bimpp is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bimpp; If not, see <http://www.gnu.org/licenses/>.
//
//  author: Carlo de Falco     <carlo.defalco _AT_ polimi.it>

/*! \file bimutil.cpp
  \brief Utility Functions for bimpp.
*/

#include <cmath>

void
bimu_bernoulli (double x, double &bp, double &bn)
{
  const double xlim = 1.0e-2;
  double ax  = fabs (x);

  bp  = 0.0;
  bn  = 0.0;

  //  X=0
  if (x == 0.0)
    {
      bp = 1.0;
      bn = 1.0;
      return;
    }

  // ASYMPTOTICS
  if (ax > 80.0)
    {
      if (x > 0.0)
        {
          bp = 0.0;
          bn = x;
        }
      else
        {
          bp = -x;
          bn = 0.0;
        }
      return;
    }

  // INTERMEDIATE VALUES
  if (ax <= 80 &&  ax > xlim)
    {
      bp = x / (exp (x) - 1.0);
      bn = x + bp;
      return;
    }

  // SMALL VALUES
  if (ax <= xlim &&  ax != 0.0)
    {
      double jj = 1.0;
      double fp = 1.0;
      double fn = 1.0;
      double df = 1.0;
      double segno = 1.0;
      while (fabs (df) > 1.0e-16)
        {
          jj += 1.0;
          segno = -segno;
          df = df * x / jj;
          fp = fp + df;
          fn = fn + segno * df;
        }
      bp = 1 / fp;
      bn = 1 / fn;
      return;
    }

};

void
bimu_bernoulli_derivative (double x, double &bpp, double &bnp)
{
  const double xlim = 1.0e-5;
  double ax  = fabs (x);

  double bp  = 0.0;
  double bn  = 0.0;
  bimu_bernoulli (x, bp, bn);

  bpp  = 0.0;
  bnp  = 0.0;

  //  X=0
  if (x == 0.0)
    {
      bpp = -.5;
      bnp =  .5;
      return;
    }

  // INTERMEDIATE VALUES
  if (ax > xlim)
    {
      bpp = (bp / x) * (1 - bn);
      bnp = - (bn / x) * (1 - bp);
      return;
    }

  // SMALL VALUES
  if (ax <= xlim &&  ax != 0.0)
    {
      bpp = -.5 + x / 6.0 - pow (x, 3) / 180.0;
      bnp =  .5 + x / 6.0 - pow (x, 3) / 180.0;
      return;
    }

};
