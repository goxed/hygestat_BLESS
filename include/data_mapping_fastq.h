#ifndef DATA_MAPPING_FASTQ_H
#define DATA_MAPPING_FASTQ_H

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

extern std::string close_bar;
extern std::string pre_qual;
extern std::string pos_qual;
extern std::string data_name;
bool   hits_done(std::string file);
extern bool TELOMERE_FILE;
extern bool DEBUG;
extern bool PRE_QUALITY;
extern bool POS_QUALITY ;
#endif
