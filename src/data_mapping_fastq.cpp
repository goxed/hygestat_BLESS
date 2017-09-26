#include "../include/data_mapping_fastq.h"

//*********************************************FUNCTION********************************************
std::vector <std::string> return_hit_dir(std::string fq_file, std::string nature_of_data) {
//**************************************************************************************************
  std::vector <std::string> dir;
  std::string D_1,D_2;
  if(nature_of_data=="bless"){
    if(!hits_done(fq_file)){ function_hits_bless(fq_file,"both");}
    D_1=fq_file+"_close_barcode/";
    D_2=fq_file+"_distant_barcode/";
  }else if(nature_of_data=="genomic"){
    if(!hits_done(fq_file)){function_hits_genomic(fq_file);}
    D_1=fq_file+"_no_barcode/";
    D_2=fq_file+"_no_barcode/";
  }else{ std::cerr<<"MAJOR ERROR in return hits.....\n";}
  dir.push_back(D_1);
  dir.push_back(D_2);
  return dir;
}



//*********************************************FUNCTION********************************************
void data_mapping_fastq(std::string fast_file, std::string nature_of_data){
//**************************************************************************************************
// 1st step : preprocess all fastq files -- output are fasta files
boost::filesystem::path s_qc("./"+fast_file+"c.html");
char tmp[257];
std::string fastfile = fast_file+".dump.fa";
if((PRE_QUALITY)&&(!exists(s_qc))){ // if we decide to do the pre-quality check
      char pre[255];
      std::sprintf(pre,"%s %s",pre_qual.c_str(),fast_file.c_str());
      int t_cmd1 = std::system(pre); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
  }
  if(data_name=="fastq"){
  process_fast(fast_file);
  } else if(data_name=="fasta" && nature_of_data=="bless"){
  std::sprintf(tmp, "%s %s %s", "cp ", fast_file.c_str(), fastfile.c_str());
  std::system(tmp);
  }
  if(nature_of_data=="bless"){
    extract_barcode(fast_file);
    if(TELOMERE_FILE) extract_telomere_file(fast_file);
    bowtie_function(fast_file,nature_of_data);
      boost::filesystem::path m_qc("./"+fast_file+".samstat.html");
    if((POS_QUALITY)&&(!exists(m_qc))) { // if we decide to do the pos-quality check
      char pos[255], conv[255],cle[255];
      std::string bt = fast_file+".dump.fa."+close_bar+".fasta.bt.btt";
      std::string bt_bam = fast_file+".dump.fa."+close_bar+".fasta.bam";
      std::sprintf(conv,"%s %s %s ","cp ",bt.c_str(), bt_bam.c_str());
      std::sprintf(pos,"%s %s %s",pos_qual.c_str(), bt_bam.c_str(), fast_file.c_str());
      int c_cmd1 = std::system(conv); if(c_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
      int p_cmd1 = std::system(pos);  if(p_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
  }

  }else if(nature_of_data=="genomic")  {
    if(data_name=="fasta"){
    bowtie_function(fast_file,nature_of_data);
    } else if(data_name=="fastq"){
    bowtie_function_genomic(fast_file);
    }
  }else{
    std::cerr<<"This case is still in progress ....."<<std::endl;
  }
}

bool hits_done(std::string file){
    bool done = false;
    std::string output1=file+"_no_barcode";
    std::string output2=file+ "_close_barcode";
    std::string output3=file+"_distant_barcode";
    boost::filesystem::path   path_1(output1);
    boost::filesystem::path   path_2(output2);
    boost::filesystem::path   path_3(output3);
    if (exists(path_1) || exists(path_2) || exists(path_3)){
        done = true;
    }
  return done;
}
