#ifndef HEADER_H
#define HEADER_H

#include<iostream>
#include<fstream>
#include<climits>
#include<vector>
#include<string>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<cctype>
#include <random>
#include <ctime>
#include <unistd.h>
#include <getopt.h>
#include <libgen.h>
#include <sys/stat.h>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_cdf.h>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/type_traits/is_empty.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <parallel/algorithm>
#define WAITUSER std::cerr<<std::endl<<"Press Return";std::cin.get();
#define X_CHROMOSOME 1000
#define Y_CHROMOSOME 1001
#define WINDOWSIZEDEFAULT 10250
#define WINDOWADVANCEDEFAULT WINDOWSIZEDEFAULT+1
#define WAIT cin.get();
#define PHRED_CUTOFF 0
#define PHRED_QUALITY 30
#define LENGTH_CUTOFF 34
#define MAXTHREADS 100

#define  DEF_TELOMERE       false         // Are counting the number of barcoded telomeres?
#define  DEF_TELOMERE_FILE  false   // create telomere count files ?
#define  DEF_DEBUG          false          // output debugging information?
#define  DEF_TIME_POINT     false       // Are you running time point?
#define  DEF_BOWTIE         false          // have we ran bowtie previously? (saving computing time)
#define  DEF_RAREFACTION    false     // Are you doing rarefaction
#define  DEF_PRE_QUALITY    false  // output information for the prequality check?
#define  DEF_POS_QUALITY    false  // output the information about the mapping
#define  DEF_rand_maxi      20  // number of random files (odd number for rand_maxi is best);



#endif
