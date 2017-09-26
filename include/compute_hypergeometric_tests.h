#ifndef COMPUTE_HYPERGEO_TEST_H
#define COMPUTE_HYPERGEO_TEST_H

#include<iostream>
#include<fstream>
#include<climits>
#include<vector>
#include<string>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<cmath>
#include<cctype>
#include <random>
#include <ctime>

#include "hygestat.h"

extern bool DEBUG;
extern pthread_mutex_t mutex_output;
extern pthread_mutex_t mutex_cerr;
extern pthread_mutex_t mutex_read1;
extern pthread_mutex_t mutex_read_map;
extern unsigned int WINDOWSIZE;
extern unsigned int WINDOWADVANCE;
extern unsigned int USER_WINDOWADVANCE;

#endif
