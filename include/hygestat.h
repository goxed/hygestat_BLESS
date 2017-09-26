#ifndef HYGESTAT_H
#define HYGESTAT_H
#include "header.h"


void *threadFunc(void *arg);
void *threadFunc1(void *arg);
void unzipFile(std::string InFile);
void extract_barcode(std::string fastFile);
void read_config_file(std::string infile);
void sampling_btt_file(std::string reads_file,long min_map,int rand_cut);
std::vector <std::string> return_hit_dir(std::string fastq, std::string nature_of_data);
long get_min_mappable (std::vector <std::string> all_fastq);
std::vector<std::string> pick_n_random_lines( const std::string& path, size_t n , int rand_cut);
//inline bool is_empty(std::ifstream& pFile){ return pFile.peek() == std::ifstream::traits_type::eof();}
void set_ouput_files_and_stat(void);
void stat_rarefaction(std::string file);
double median(std::vector<double> vec);
std::string extract_telomere(std::string line);
bool contain_telomere(std::string line);
void extract_telomere_file(std::string fastFile);
size_t count_line_in_file(std::string file);

unsigned int bowtie_function(std::string fastafile, std::string nature_of_data);
unsigned int bowtie_function_genomic(std::string fastqfile);
unsigned int extract_chrnum(std::string s_chrNum, bool XY);
std::string convert_chrnum_to_string2(unsigned int chrNum);
std::string convert_chrnum_to_string(unsigned int chrNum);
void process_fast(std::string fastqFileName);

void function_hits_bless(std::string file_from_bowtie,std::string nature);
void function_hits_genomic(std::string file_from_bowtie);
void cytoband_function(std::string read_directory);
void cytoband_function1(std::string reads_directory) ;
void   print_usage();
int    mkpath(const char *s, mode_t mode);
int parse_options(int argc, char** argv);
void clean_up_unwanted_files();
void ListAndRemoveFilesRecursively(const char *dir, const char* ext);

unsigned int getMappableWindowEnd(unsigned int startChrPos,
                                  unsigned int mappableWindowSize, std::vector<bool>&mappabilityMap);
double hgcalc_gte(unsigned int ps, unsigned int sa,unsigned int sp,unsigned int ss);

int getPhredOffset(std::string fastqFileName);



inline bool fexists(const char *filename);
template <class T>
inline std::string to_string(const T& t){
  std::stringstream ss;
  ss<<t;
  return ss.str();
}

enum fastqFileSequenceType{fastaHeader, fastaLine, qualityHeader, qualityLine};
struct CYTOBAND{
  unsigned int chrNum;
  unsigned int chrStart;
  unsigned int chrEnd;
  std::string name;
};
struct CHR_READS{
  unsigned int chrNum;
  unsigned int reads;
};

struct SEQUENCE_READS{
  unsigned int chrPos;
  unsigned int hits;
};

struct MAPPABILITY_DATA{
  unsigned int readSize;
  unsigned int chrSize;
  unsigned int chrNum;
  unsigned int chrMapableBases;
  std::vector<bool> mappabilityMap;
};

struct WINDOW_DATA{
  unsigned int chrNum;
  unsigned int windowReads1;
  unsigned int windowReads2;
  unsigned int windowStart;
  unsigned int windowEnd;
  double pval;
  double qval;
  double bval;
  double pvalCNVIns;
  double pvalCNVDel;
  double correctionCoefficientSequencing;
  double correctionCoefficientCNV;
  double pvalCNV2Ins;
  double pvalCNV2Del;
  double correctionCoefficientCNV2;
  double pvalCorrected;
  double qvalCorrected;
};


struct SortByChrPos {
  bool operator () (const SEQUENCE_READS & lhs , const SEQUENCE_READS & rhs) const
  {
    return lhs.chrPos < rhs.chrPos;
  }
};

struct chrData{
  //std::vector<unsigned int> chrNum;
  //std::vector<unsigned int> chrPos;
  unsigned int chrPos;
  unsigned int chrNum;
  std::string strand;
  //std::vector<int> strand;
};



struct threadData{
  unsigned int chrNum;
  unsigned int chrReads;
  std::string chrReadsFileName;
  std::string Dir1;
  std::string Dir2;
  std::string Dir3;
  std::string Dir4;
  std::string Dir5;
  std::vector<MAPPABILITY_DATA> *mappabilityMapList;
  std::vector<WINDOW_DATA> *outData;
  unsigned int mapableBasesTotal;
  unsigned int mappableBases;
  double averageWindowReads;
  double averageWindowReadsCNV;
};

inline void terminate(int retval)
{
  std::cerr << std::endl << "Terminating! because of Error " << retval << std::endl;
  exit(retval);
}

inline std::string convert(std::string sentence)
{
  for(int i(0); i<sentence.size(); ++i)

    {
      sentence[i] = toupper(sentence[i]);
    }
  return sentence;
}

inline bool ChrSortFun(unsigned int c1, unsigned int c2){
  return (c1 < c2);
}

inline bool ChrSortFun1(struct chrData c1, struct chrData c2){
  return (c1.chrNum < c2.chrNum);
}
inline int checkFile(std::string fileName) {
  std::ifstream fileOpenTest;
  int isFileOpen=-1;
  fileOpenTest.open(fileName.c_str());
  if (fileOpenTest.is_open())
    {
      isFileOpen=1;
      fileOpenTest.close();
    }
  return isFileOpen;
}

//**********************check if file exists ************************************************//
inline bool fexists(const char *filename)
{
  std::ifstream ifile(filename);
  return ifile;
}
//************************************************************************* //
inline bool sortSeqReadsFun (struct SEQUENCE_READS s1, struct SEQUENCE_READS s2) {
//**************************************************************************//
  return (s1.chrPos < s2.chrPos);
}
//**************************************************************************//
inline bool SortOutDataPvalFun (struct WINDOW_DATA s1, struct WINDOW_DATA s2){
//***************************************************************************
  return (s1.pval<s2.pval);
}
//******************************************************************************//
inline bool SortOutDataPvalCorrectedFun (struct WINDOW_DATA s1, struct WINDOW_DATA s2){
//*******************************************************************************//
  return (s1.pvalCorrected<s2.pvalCorrected);
}
//******************************************************************************//
inline bool SortOutDataWindowStartFun (struct WINDOW_DATA s1, struct WINDOW_DATA s2){
  return (s1.windowStart<s2.windowStart);
}
//******************************************************************************//
inline bool SortOutDataChrNumFun (struct WINDOW_DATA s1, struct WINDOW_DATA s2){
  return (s1.chrNum<s2.chrNum);
}
//------------------- Functions using local structures -------------------------------
inline std::string group_integers(long long int num){
  std::stringstream ss_num;
  ss_num<<num;
  std::string s_num="";
  ss_num>>s_num;
  std::string::size_type s=s_num.size();
  if (s<4) return s_num;
  std::string g_num="";
  for (unsigned int p=0; p<s; ++p){
    g_num=g_num+s_num[p];
    if ((s-p-1)&&(s-p-1)%3==0){
      g_num=g_num+",";
    }
  }
  return g_num;
}
void write_lab_info();
int getChrReads(std::string FileName, std::vector<CHR_READS> &chrList );
int getMappabilityMap(std::string chrNumFileName, std::vector <MAPPABILITY_DATA> &mappabilityMapList);
double calcAverageWindowReads(unsigned chrNum, std::string chrReadsFileName, std::string Dir3,
                               std::vector<MAPPABILITY_DATA>&mappabilityMapList);
double getSequenceReads(std::string Dir, std::vector<SEQUENCE_READS> &sequenceReads,unsigned int chrNum);
int calcHypergeometricPvalWith(unsigned int chrNum,unsigned int chrReads, double averageWindowReads,
                                double averageWindowReadsCNV, unsigned mapableBasesTotal,
                                std::string chrReadsFileName, std::string Dir1, std::string Dir2, std::string Dir3,
                                std::string Dir4, std::string Dir5, std::vector<MAPPABILITY_DATA>&mappabilityMapList,
                                 std::vector<WINDOW_DATA> &outData);
int calcHypergeometricPvalWithout(unsigned int chrNum,unsigned int chrReads,std::string chrReadsFileName,
                                  std::string Dir1,std::string Dir2,std::vector<MAPPABILITY_DATA>&mappabilityMapList,
                                   std::vector<WINDOW_DATA> &outData);
int calcHypergeometricPvalOneSample(unsigned int chrNum,unsigned int chrReads,unsigned int mappableBases,
                                  std::string chrReadsFileName,std::string Dir1,std::vector<MAPPABILITY_DATA>&mappabilityMapList,
                                  std::vector<WINDOW_DATA> &outData);
int hgStatsTwoSampleWith(std::string Dir1, std::string Dir2, std::string Dir3, std::string Dir4, std::string Dir5,
                                 std::string chrReadsFileName, std::string outputFileName, unsigned   WINDOWSIZE,
                                 unsigned WINDOWADVANCE);
int hgStatsTwoSampleWithout(std::string Dir1, std::string Dir2, std::string chrReadsFileName, std::string outputFileName,
                                 unsigned   WINDOWSIZE, unsigned WINDOWADVANCE);
unsigned int bh(std::vector<double>&pvalList,std::vector<double>&qValList);
unsigned int bonf(std::vector<double>&pvalList,std::vector<double>&qValList);
int hgStatsOneSample(std::string Dir1, std::string chrReadsFileName, std::string outputFileName, unsigned WINDOWSIZE,
                       unsigned WINDOWADVANCE);
//void two_samples_hygeostat(std::string fast_file_t, std::string fast_file_c, std::string fast_file_t_cel, std::string fast_file_cel, std::string fast_file_fib, std::string nature_of_analysis);
void data_mapping_fastq(std::string fast_file, std::string nature_of_data);
void two_samples_hygeostat(std::string fast_file_t, std::string nat_dat_t, std::string fast_file_c,  std::string nat_dat_c,
			   std::string fast_file_t_cel, std::string nat_dat_t_cel, std::string fast_file_cel, std::string nat_dat_cel,
			   std::string fast_file_fib, std::string nat_dat_fib, std::string nature_of_analysis,int resolution);

int return_wind_type(std::vector <int> wind_t);
void stat_funct();

#define CORRECTION_SAMPLE_1       202
#define CORRECTION_SAMPLE_2       203
#define CORRECTION_SAMPLE_3       204
#define TIME_SERIE_FILE           205
#define TIME_SERIE_RAND           206
#define TIME_SERIE_RARE           207
#define NO_UPDATE_CHECK           208
#define CORRECTION_SAMPLE_1_NAT   209
#define CORRECTION_SAMPLE_2_NAT   210
#define CORRECTION_SAMPLE_3_NAT   211


#endif
