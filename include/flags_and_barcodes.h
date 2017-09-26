#ifndef FLAGS_AND_BARCODES_H
#define FLAGS_AND_BARCODES_H

bool TELOMERE        = DEF_TELOMERE;
bool TELOMERE_FILE   = DEF_TELOMERE_FILE;
bool DEBUG           = DEF_DEBUG;
bool TIME_POINT      = DEF_TIME_POINT;
bool BOWTIE          = DEF_BOWTIE;
bool RAREFACTION     = DEF_RAREFACTION;
bool PRE_QUALITY     = DEF_PRE_QUALITY;
bool POS_QUALITY     = DEF_POS_QUALITY;
bool no_mappability  = true;
int  rand_mini=1,rand_maxi=DEF_rand_maxi,rand_step=1;
size_t min_mapp_reads;

std::string close_bar        = "TCGAGGTAGTA";
std::string close_link       = "CCCTAGCGTAACTCTCGAGGTAGTA";
std::string distant_bar      = "TCGAGACGACG";
//std::string distant_link ="CCCTAGCGTAACTCTCGAGACGACG"; // keep the same linker as close..........
std::string distant_link     = "CTAGCGTAACTCTCGAGACGACG";
std::string telomere_ver     = "TTAGGG";
std::string telomere_yea     = "ACTGGTGT"; //Candida guillermondii
std::string XhoI             = "TCGAG";
int N_P = 16; // number of processors
std::string chrReadsFileName = "chrReads.txt";
std::string btt              = "btt ";
std::string bowt_arg         = " -l61 -m1 -n0 -r -p"+to_string(N_P)+" ";
std::string bowt_arg_q       = " -l61 -m1 -n0 -q -p"+to_string(N_P)+" ";
std::string bowt             = "bowtie ";
std::string bowt_q           = "bowtie ";
std::string genome_mouse     = "";
std::string genome_human     = "";
std::string genome_yeast     = "";
std::string link1            = "hit.sh ";
std::string link2            = "hit1.sh ";
std::string reads_dist       = "chr_reads_dist.sh";
std::string pre_qual         = "fastqc ";                       // fastqc for pre-quality check
std::string pos_qual         = "samstat ";            // samstat for pos-quality mapping
std::string mapDir           = "";
std::string mappaDir         = "";
pthread_mutex_t mutex_output       = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_cerr         = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_read1        = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_read_map     = PTHREAD_MUTEX_INITIALIZER;
unsigned int    WINDOWSIZE         = WINDOWSIZEDEFAULT;//make them local in next revision
unsigned int    WINDOWADVANCE      = WINDOWADVANCEDEFAULT;
unsigned int    USER_WINDOWADVANCE = WINDOWADVANCE-1;
std::string resolution;
bool        resol_flag        = false;
bool        no_bedgraph_file  = false;
bool        dont_keep_files   = false;
int         run_type=2,time_serie_rand=20;
std::string org_name="",nature_of_analysis="without_correction",nature_of_data1="", nature_of_data2="", data_name="";
std::string nature_of_data3="", nature_of_data4="", nature_of_data5="";
std::string cbandsFileName="",fbandsFileName="";
std::string human, HUMAN, mouse, MOUSE, yeast, YEAST,time_course_file;
std::string with_correction,bless, genomic, without_correction,outputFileName="output.txt";
std::string fastqFile1,fastqFile2,fastqFile3,fastqFile4,fastqFile5,conf_file;
std::string s_windowSize,s_windowAdvance;
std::vector <std::string> treat_fastq_tim,cont_fastq_tim;
std::vector <std::string> nat_treat_tim,nat_cont_tim;
std::vector <int>    wind_run;
std::string wind_type;
std::string Dir1="./";
std::string Dir2="./";
std::string Dir3="./";
std::string Dir4="./";
std::string Dir5="./";
std::string Dir11="./";
std::string Dir22="./";
std::string Dir33="./";
std::string Dir44="./";
std::string Dir55="./";
std::string btt_file= "", bowt_file="";
std::vector<chrData> data;
std::string output_dir="./",genome_dir="./",genome_type,fastq_control,fastq_treatment, fasta_or_fastq="fastq";
std::string nature_control,nature_treatment,data_dir="./",dont_run_bowtie="false",telomere="false",time_serie="false";
std::string pre_quality_check="false",post_mapping_check="false",cytoband,fragile_band;
std::string input_file,corSample1,corSample2,corSample3,output_file="";
std::string no_update_check,correction="",verbose,quiet_mode, time_serie_rare;
std::string corSample1Nat,corSample2Nat,corSample3Nat;
std::string time_serie_file;
unsigned int total_R1=0, total_R2   =0, total_R3   =0, total_R4   =0, total_R5   =0;
unsigned int barc_R1 =0, barc_R2    =0, barc_R3    =0, barc_R4    =0, barc_R5    =0;
unsigned int map_R1  =0, map_R2     =0, map_R3     =0, map_R4     =0, map_R5     =0;
double perc_bar1     =0.0, perc_bar2=0.0, perc_bar3=0.0, perc_bar4=0.0, perc_bar5=0.0;
double perc_map1     =0.0, perc_map2=0.0, perc_map3=0.0, perc_map4=0.0, perc_map5=0.0;
#endif
