// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once



#define USE_PRECOMPILED_HEADER

#ifdef USE_PRECOMPILED_HEADER
#include "targetver.h"

#include "mpi.h"


#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include "PBLAS.h"
#include "BLACS.h"




#include "boost/lexical_cast.hpp"
#include <tclap/CmdLine.h>

#define USE_PROFILER
#define USE_BOOST_CHRONO
#include "tic-toc-profiler.hpp"

#include "H5inDMD.h"

#include "mkl.h"
//#include <stdio.h>
//#include <tchar.h>

#include "Vandermonde.h"

#include "ColumnSquaredNorm.h"

#endif

// TODO: reference additional headers your program requires here
