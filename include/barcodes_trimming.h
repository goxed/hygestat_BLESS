
#ifndef BARCODES_TRIMMING_H
#define BARCODES_TRIMMING_H

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
extern bool TELOMERE ;
extern std::string close_bar;
extern std::string close_link;
extern std::string distant_bar;
extern std::string distant_link ;

unsigned int bowtie_function(std::string fastafile);
unsigned int bowtie_function_genomic(std::string fastqfile);
void         extract_barcode(std::string fastFile);
bool         contain_telomere(std::string line);

#endif


