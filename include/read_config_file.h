#ifndef READ_CONFIG_FILE_H
#define READ_CONFIG_FILE_H

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

extern bool TIME_POINT;
extern pthread_mutex_t mutex_output;
extern pthread_mutex_t mutex_cerr;
extern pthread_mutex_t mutex_read1;
extern pthread_mutex_t mutex_read_map;
extern unsigned int WINDOWSIZE;
extern unsigned int resolution;
extern unsigned int WINDOWADVANCE;
extern unsigned int USER_WINDOWADVANCE;
extern bool TELOMERE ;
extern bool DEBUG;
extern bool BOWTIE;
extern bool RAREFACTION;
extern bool PRE_QUALITY;
extern bool POS_QUALITY;
extern int rand_maxi,N_P;
extern std::string genome_mouse;
extern std::string genome_human ;
extern std::string genome_yeast ;
extern std::string output_dir, genome_dir,data_dir;
extern std::string nature_of_analysis,org_name;
extern std::string nature_of_data1, nature_of_data2, nature_of_data3, nature_of_data4, nature_of_data5, data_name;
extern std::string fbandsFileName;
extern std::string cbandsFileName;
extern std::string time_course_file;
extern std::string fastqFile1;
extern std::string fastqFile2;
extern std::string fastqFile3;
extern std::string fastqFile4;
extern std::string fastqFile5;
extern std::string outputFileName;
extern std::string conf_file;
extern std::string s_windowSize;
extern std::string s_windowAdvance;
extern std::vector <int> wind_run;
extern std::string wind_type;
extern std::vector <std::string> treat_fastq_tim;
extern std::vector <std::string> cont_fastq_tim;
extern std::vector <std::string> nat_treat_tim;
extern std::vector <std::string> nat_cont_tim;
extern std::string mapDir;
extern bool no_mappability;
std::string mappability_dir="./";

#endif

