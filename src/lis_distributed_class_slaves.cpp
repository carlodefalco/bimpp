/*
  Copyright (C) 2011 Carlo de Falco
  This software is distributed under the terms
  the terms of the GNU/GPL licence v3
*/
/*! \file lis_distributed_class.cpp
  \brief interface for linear solver built with for lis library.
*/

#include "lis_distributed_class.h"
#include <stdlib.h>
#include <sstream>
#include <string>
#include <cstring>

void
lis_distributed::cleanup_slaves ()
{
  delete [] data;
  delete [] rhs;
 
  if (have_initial_guess)
    delete [] initial_guess;
}
