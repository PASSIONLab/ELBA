//==================================================================
// Title:  Cuda x-drop seed-and-extend alignment algorithm
// Author: G. Guidi, A. Zeni
// Date:   6 March 2019
//==================================================================

#ifndef __LOGAN_CUH__
#define __LOGAN_CUH__

#include<vector>
#include<iostream>
#include<chrono>
#include<numeric>
#include<functional>
#include<iterator>

#include"seed.hpp"
#include"score.hpp"

#ifdef ADAPTABLE
#include"functions.hpp"
#else
#include"functions_static.hpp"
#endif

#endif
