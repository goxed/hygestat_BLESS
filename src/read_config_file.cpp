#include "../include/read_config_file.h"

//*********************************************FUNCTION********************************************

void read_config_file(std::string infile){
  //**************************************************************************************************

  std::string    tel,bow,s_rare,rar,posCheck;
  std::string    s_np,deb,preCheck;
  size_t found;
  std::ifstream  in(infile.c_str());
  std::string    temp;
  getline(in,temp);
  getline(in,org_name);
  getline(in,temp);
  getline(in,nature_of_analysis);
  getline(in,temp);
  getline(in,fastqFile1);
  found = fastqFile1.find(".gz");
  if(found!=std::string::npos){
   unzipFile(fastqFile1);
   fastqFile1.erase(fastqFile1.length()-3);
  }
  getline(in,temp);
  getline(in,nature_of_data1);
  getline(in,temp);
  getline(in,fastqFile2);
  found = fastqFile2.find(".gz");
  if(found!=std::string::npos){
   unzipFile(fastqFile2);
   fastqFile2.erase(fastqFile2.length()-3);
  }
  getline(in,temp);
  getline(in,nature_of_data2);
  if(nature_of_analysis=="with_correction"){
    getline(in,temp);
    getline(in,fastqFile3);
    found = fastqFile3.find(".gz");
    if(found!=std::string::npos){
     unzipFile(fastqFile3);
     fastqFile3.erase(fastqFile3.length()-3);
    }
    getline(in,temp);
    getline(in,nature_of_data3);
    getline(in,fastqFile4);
    found = fastqFile4.find(".gz");
    if(found!=std::string::npos){
     unzipFile(fastqFile4);
     fastqFile4.erase(fastqFile4.length()-3);
    }
    getline(in,temp);
    getline(in,nature_of_data4);
    getline(in,fastqFile5);
    found = fastqFile5.find(".gz");
    if(found!=std::string::npos){
     unzipFile(fastqFile5);
     fastqFile5.erase(fastqFile5.length()-3);
    }
    getline(in,temp);
    getline(in,nature_of_data5);
  }else{ for(int i =0;i<10;i++){ getline(in,temp);} }
  if(org_name=="yeast"|| org_name=="YEAST"){
    for(int ii=0;ii<4;ii++) {
      getline(in,temp);
    }
  }else{
    getline(in,temp);
    getline(in,fbandsFileName);
    getline(in,temp);
    getline(in,cbandsFileName);
  }
  getline(in,temp);
  getline(in,outputFileName);
  getline(in,temp);
  getline(in,s_windowSize);
  getline(in,temp);
  getline(in,s_windowAdvance);
  std::stringstream ss_windowAdvance(s_windowAdvance);
  ss_windowAdvance>>USER_WINDOWADVANCE;
  getline(in,temp);
  std::string temp_file;
  getline(in,temp_file);
  if((temp_file!="NO")&&(temp_file!="no")) {
    TIME_POINT = true;
    if (checkFile(temp_file)!=1) terminate(1);
    std::ifstream in2(temp_file.c_str());
    std::string line;
    std::vector <std::string> temp;
    while(getline(in2,line)) {
      boost::split(temp,line,boost::is_space());
      std::string tr   = temp[0];
      std::string ct   = temp[1];
      std::string na_t = temp[2];
      std::string na_c = temp[3];
      treat_fastq_tim.push_back(tr);
      cont_fastq_tim.push_back(ct);
      nat_treat_tim.push_back(na_t);
      nat_cont_tim.push_back(na_c);
    }
  }
  getline(in,temp);
  getline(in,rar);
  RAREFACTION = rar=="false"? false:true;
  getline(in,temp);
  getline(in,s_rare);
  std::stringstream ss_rare(s_rare);
  ss_rare>>rand_maxi;
  getline(in,temp);
  getline(in,data_name);
  getline(in,temp);
  getline(in,genome_dir);
  getline(in,temp);
  getline(in,mappability_dir);
  getline(in,temp);
  getline(in,data_dir);
  getline(in,temp);
  getline(in,bow);
  BOWTIE = bow == "false"? false:true;
  getline(in,temp);
  getline(in,tel);
  TELOMERE = tel == "false"? false:true;
  getline(in,temp);
  getline(in,preCheck);
  PRE_QUALITY = preCheck == "false"? false:true;
  getline(in,temp);
  getline(in,posCheck);
  POS_QUALITY = posCheck == "false"? false:true;
  getline(in,temp);
  getline(in,s_np);
  std::stringstream ss_np(s_np);
  ss_np>>N_P;
  getline(in,temp);
  getline(in,deb);
  DEBUG = deb == "false"? true:false;
  in.close();

  boost::filesystem::path data_dir(mappability_dir);
  if(boost::filesystem::is_directory(mappability_dir)){
     no_mappability = false;
     mapDir         = mappability_dir;
  }else{std::cerr<<"No mappability data provided, Will run without \n";}


  if      (org_name=="human") {genome_human  = genome_dir;}
  else if (org_name=="mouse") {genome_mouse  = genome_dir;}
  else if (org_name=="yeast") {genome_yeast  = genome_dir;}
  else { std::cerr<< "Error: please provide a correct genome type (human,mouse, or yeast)\n"; exit(1);}

  std::vector <std::string> wind;
  split(wind,s_windowSize,boost::is_space());

  for(unsigned int i = 0; i < wind.size(); i++){
  if(!wind[i].empty()) wind_run.push_back(stoi(wind[i]));
  }
}

void unzipFile(std::string InFile){
    char tmp[255];
    std::string OutFile, outf;
    sprintf(tmp, "%s %s %s", "gunzip", "-k", InFile.c_str());
    OutFile = system(tmp);
  }

void write_lab_info(){

    std::cerr << "                           *      *  *     *  ******   ******   ******  *******     *           *******    \n";
    std::cerr << "                          *      *    *   *  *    *   *        *          *        * *            *        \n";
    std::cerr << "                         *      *      * *  *        *        *          *        *   *          *         \n";
    std::cerr << "                        ********        *  *  ***   ******   ******     *        *******        *          \n";
    std::cerr << "                       *      *        *  *    *   *             *     *        *       *      *           \n";
    std::cerr << "                      *      *        *  *    *   *             *     *        *         *    *            \n";
    std::cerr << "                     *      *        *  ******   ******   ******     *        *           *  *             \n";
    std::cerr << "                     ######################################################################################\n";
    std::cerr << "                     #                                                                                    #\n";
    std::cerr << "                     # Running hygestat_windows for BLESS  data (version hgw.1.2.3, September 2016)       #\n";
    std::cerr << "                     # This software comes without any warranty and .....                                 #\n";
    std::cerr << "                     # Copyright Rowicka Lab @ UTMB.                                                      #\n";
    std::cerr << "                     # Please contact Jules Nde jbndeken@utmb.edu for technical help                      #\n";
    std::cerr << "                     #                                                                                    #\n";
    std::cerr << "                     ######################################################################################\n";

}
