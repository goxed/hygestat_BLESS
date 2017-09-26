#ifndef EVALUATE_MAPPABILITY_H
#define EVALUATE_MAPPABILITY_H

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



extern std::string mapDir;
extern std::string org_name;
extern bool no_mappability;
extern bool DEBUG;          // output debugging information?

struct ChrMapp{
std::vector< std::vector <unsigned int>> chrSize;
std::vector< std::vector <std::string>> chrName;
};


std::vector <std::string> chrN_human{
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY"
};
std::vector <unsigned int> chrZ_human{
        248956422,
        243199373,
        198022430,
        191154276,
        180915260,
        171115067,
        159138663,
        146364022,
        141213431,
        135534747,
        135006516,
        133851895,
        115169878,
        107349540,
        102531392,
        90354753,
        81195210,
        78077248,
        59128983,
        63025520,
        48129895,
        51304566,
        155270560,
        59373566
};
std::vector <std::string> chrN_mouse{
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chrX",
        "chrY"
};
std::vector <unsigned int> chrZ_mouse{
        195471971,
        182113224,
        160039680,
        156508116,
        151834684,
        149736546,
        145441459,
        129401213,
        124595110,
        130694993,
        122082543,
        120129022,
        120421639,
        124902244,
        104043685,
        98207768,
        94987271,
        90702639,
        61431566,
        171031299,
        91744698
};
std::vector <std::string> chrN_yeast{
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16"
};
std::vector <unsigned int> chrZ_yeast{
        230218,
        813184,
        316620,
        1531933,
        576874,
        270161,
        1090940,
        562643,
        439888,
        745751,
        666816,
        1078177,
        924431,
        784333,
        1091291,
        948066
};

#endif

