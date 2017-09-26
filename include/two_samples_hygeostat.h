#ifndef TWO_SAMPLE_HYGEOSTAT_H
#define TWO_SAMPLE_HYGEOSTAT_H

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


extern std::string chrReadsFileName;
extern std::string outputFileName;
extern std::string org_name;
extern std::string data_dir;
extern std::string Dir1;
extern std::string Dir2;
extern std::string Dir3;
extern std::string Dir4;
extern std::string Dir5;
extern std::string Dir11;
extern std::string Dir22;
extern std::string Dir33;
extern std::string Dir44;
extern std::string Dir55;
extern bool no_bedgraph_file;
extern pthread_mutex_t mutex_output;
extern pthread_mutex_t mutex_cerr;
extern pthread_mutex_t mutex_read1;
extern pthread_mutex_t mutex_read_map;
extern unsigned int WINDOWSIZE;
extern unsigned int WINDOWADVANCE;
extern unsigned int USER_WINDOWADVANCE;
extern bool DEBUG;          // output debugging information?
std::string create_wig_file_with_log_pvalue(std::string file, std::string fastq_t, std::string fastq_c, int resolution);
#endif

