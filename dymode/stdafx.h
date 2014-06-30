// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once



#define USE_PRECOMPILED_HEADER

#ifdef USE_PRECOMPILED_HEADER
#include "targetver.h"

// External dependencies
#include "mpi.h"
#include "mkl.h"

#include <iostream>
#include <fstream>

// Internal dependencies
#include "peigen.h"

#define USE_PROFILER
#define USE_BOOST_CHRONO
#include "tic-toc-profiler.hpp"

using namespace std;
using namespace Eigen;
using namespace peigen;

// Dymode helper functions
#include "ColumnSquaredNorm.h"
#include "H5inDMD.h"
#include "GEOinDMD.h"
#include "options.h"

using namespace phdfp;

#endif

// TODO: reference additional headers your program requires here
