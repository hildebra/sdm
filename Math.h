#pragma once
#include <stdio.h>
//#include <tchar.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
//#include <string.h>
#include <string>
#include <map>
#include <stdlib.h>
#include <algorithm>
#include <time.h>
#include <random>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#define _gziprea//d
#else
#define _gzipread
#endif



//#ifdef _gzipread
//#include "gzstream.h"
//#endif

const bool verbose=1;


int getRand(int until);

void swap(int &x,int &y);
