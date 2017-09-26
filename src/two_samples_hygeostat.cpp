
#include "../include/two_samples_hygeostat.h"

//*********************************************FUNCTION********************************************
void two_samples_hygeostat(std::string fast_file_t, std::string nat_dat_t, std::string fast_file_c,  std::string nat_dat_c,
      std::string fast_file_t_cel, std::string nat_dat_t_cel, std::string fast_file_cel, std::string nat_dat_cel,
      std::string fast_file_fib, std::string nat_dat_fib, std::string nature_of_analysis, int resolution) {
//**************************************************************************************************
  WINDOWSIZE    = resolution;
  WINDOWADVANCE = WINDOWSIZE+1;
  //WINDOWADVANCE = USER_WINDOWADVANCE > WINDOWSIZE ? USER_WINDOWADVANCE: WINDOWSIZE+1;
  std::vector <std::string> D_t, D_c,D_t_cel,D_cel,D_fib;
  D_t     = return_hit_dir(fast_file_t,nat_dat_t);
  D_c     = return_hit_dir(fast_file_c,nat_dat_c);
  if(nat_dat_t=="bless") cytoband_function(fast_file_t);
  else if (nat_dat_t=="genomic") cytoband_function1(fast_file_t);
  else{std::cerr<<"No cytoband (either genome is yeast or bad file provided.Stop the program if the latter)\n";}

  if(nat_dat_c=="bless" ) cytoband_function(fast_file_c);
  else if (nat_dat_c=="genomic") cytoband_function1(fast_file_c);
  else{std::cerr<<"No cytoband (either genome is yeast or bad file provided. Stop the program if the latter)\n";}

  Dir1    = D_t[0];
  Dir11   = D_t[1];
  Dir2    = D_c[0];
  Dir22   = D_c[1];
  if(nature_of_analysis=="with_correction"){
    D_t_cel = return_hit_dir(fast_file_t_cel,nat_dat_t_cel);
    D_cel   = return_hit_dir(fast_file_cel,nat_dat_cel);
    D_fib   = return_hit_dir(fast_file_fib,nat_dat_fib);
    Dir3    = D_t_cel[0];
    Dir33   = D_t_cel[1];
    Dir4    = D_cel[0];
    Dir44   = D_cel[1];
    Dir5    = D_fib[0];
    Dir55   = D_fib[1];
  }else if (nature_of_analysis=="without_correction"){ if(DEBUG) std::cerr <<"CASE_WITHOUT_CORRECTION :::: No need to correct the data \n"; }
  else{std::cerr<<"This case is still in progress ....."<<std::endl; }

   std::string f_t = fast_file_t.substr(fast_file_t.find_last_of("/")+1);
   std::string f_c = fast_file_c.substr(fast_file_c.find_last_of("/")+1);
   std::vector <std::string> temp1, temp2;
   boost::split(temp1,f_t,boost::is_any_of("."));
   boost::split(temp2,f_c,boost::is_any_of("."));
   std::string out1=data_dir+"/"+temp1[0]+"_vs_"+temp2[0]+"_close_barcode";
   std::string out2=data_dir+"/"+temp1[0]+"_vs_"+temp2[0]+"_distant_barcode";
   std::string out3=data_dir+"/"+temp1[0]+"_vs_"+temp2[0]+"_no_barcode";
   boost::filesystem::path   path_1(out1);
   boost::filesystem::path   path_2(out2);
   boost::filesystem::path   path_3(out3);
   char tmp[255],tmp1[255],tmp2[255],tmp3[255], tmp4[255];
  if (nat_dat_c=="bless" || nat_dat_t=="bless"){
  if (!boost::filesystem::exists(path_1)){
        std::sprintf(tmp2,"%s %s","mkdir ",out1.c_str());
        std::sprintf(tmp3,"%s %s","mkdir ",out2.c_str());
        int t_cmd2 = std::system(tmp2); if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
        int t_cmd3 = std::system(tmp3); if(t_cmd3 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    }
    } else{
        if (!boost::filesystem::exists(path_3)){
        std::sprintf(tmp4,"%s %s","mkdir ",out3.c_str());
        int t_cmd4 = std::system(tmp4); if(t_cmd4 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
     }
    }

  if(nat_dat_t=="bless" || nat_dat_c=="bless"){
  if(fast_file_t!=fast_file_c){
    //"Two Sample test";
    if(nature_of_analysis=="without_correction"){
        std::string new_out_file = to_string(resolution) + "_" + outputFileName;
        std::cerr<< "Computing hygestat for window : "<<WINDOWSIZE<<std::endl;
        hgStatsTwoSampleWithout(Dir1,Dir2,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp, "%s %s %s %s", "mv ", wig_file.c_str(), new_out_file.c_str(),out1.c_str());
        }else{
                std::sprintf(tmp, "%s %s %s", "mv ", new_out_file.c_str(),out1.c_str());
        }
        int t_cmd = std::system(tmp); if(t_cmd == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";

        hgStatsTwoSampleWithout(Dir11,Dir22,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp1,"%s %s %s %s", "mv ",  wig_file.c_str(), new_out_file.c_str(),out2.c_str());
        }else{
                 std::sprintf(tmp1,"%s %s %s", "mv ", new_out_file.c_str(),out2.c_str());
        }
        int t_cmd1 = std::system(tmp1); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    } else if(nature_of_analysis=="with_correction"){
        std::string new_out_file = to_string(resolution) + "_" + outputFileName;
        std::cerr<< "Computing hygestat for window : "<<WINDOWSIZE<<std::endl;
        hgStatsTwoSampleWith(Dir1,Dir2,Dir3,Dir4,Dir5,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp, "%s %s %s %s", "mv ", wig_file.c_str(), new_out_file.c_str(),out1.c_str());
        }else{
                std::sprintf(tmp, "%s %s %s", "mv ", new_out_file.c_str(),out1.c_str());
        }
       int t_cmd = std::system(tmp); if(t_cmd == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
       hgStatsTwoSampleWith(Dir11,Dir22,Dir33,Dir44,Dir55,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
       if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp1,"%s %s %s %s", "mv ",  wig_file.c_str(), new_out_file.c_str(),out2.c_str());
        }else{
                std::sprintf(tmp1,"%s %s %s", "mv ", new_out_file.c_str(),out2.c_str());
        }
       int t_cmd1 = std::system(tmp1); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    }


  } else if(fast_file_t==fast_file_c){
    std::string Y,y,N,n;
    std::cerr<<std::endl<<"One Sample test provided, enrichement will be computed versus average. ";
    std::cerr<<std::endl<<"Do you want to continue (Y/N)? ";
    char test;
    std::string test_cond;
    std::stringstream ss;
    std::cin.get(test);
    ss<<test;
    ss>>test_cond;
    if(test_cond=="Y" || test_cond=="y"){
        std::string new_out_file = to_string(resolution) + "_" + outputFileName;
        std::cerr<< "Computing hygestat one sample for window : "<<WINDOWSIZE<<std::endl;
        hgStatsOneSample(Dir1,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp, "%s %s %s %s", "mv ", wig_file.c_str(), new_out_file.c_str(),out1.c_str());
        }else{
                std::sprintf(tmp, "%s %s %s", "mv ", new_out_file.c_str(),out1.c_str());
        }
        int c_cmd = std::system(tmp); if(c_cmd == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
        hgStatsOneSample(Dir1,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
       if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp1,"%s %s %s %s", "mv ",  wig_file.c_str(), new_out_file.c_str(),out2.c_str());
        }else{
                std::sprintf(tmp1,"%s %s %s", "mv ", new_out_file.c_str(),out2.c_str());
        }
      int t_cmd1 = std::system(tmp1); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    } else if(test_cond=="N" || test_cond=="n"){
      std::cerr<<"Thanks for using InstSeq from Rowicka lab (UTMB Health). "<<std::endl;
    }
  }
 } else{
       if(fast_file_t!=fast_file_c){
    //"Two Sample test";
    if(nature_of_analysis=="without_correction"){
        std::string new_out_file = to_string(resolution) + "_" + outputFileName;
        std::cerr<< "Computing hygestat for window : "<<WINDOWSIZE<<std::endl;
        hgStatsTwoSampleWithout(Dir1,Dir2,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp, "%s %s %s %s", "mv ", wig_file.c_str(), new_out_file.c_str(),out3.c_str());
        }else{
                std::sprintf(tmp, "%s %s %s", "mv ", new_out_file.c_str(),out3.c_str());
        }
        int t_cmd = std::system(tmp); if(t_cmd == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    } else if(nature_of_analysis=="with_correction"){
        std::string new_out_file = to_string(resolution) + "_" + outputFileName;
        std::cerr<< "Computing hygestat for window : "<<WINDOWSIZE<<std::endl;
        hgStatsTwoSampleWith(Dir1,Dir2,Dir3,Dir4,Dir5,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp, "%s %s %s %s", "mv ", wig_file.c_str(), new_out_file.c_str(),out3.c_str());
        }else{
                std::sprintf(tmp, "%s %s %s", "mv ", new_out_file.c_str(),out3.c_str());
        }
       int t_cmd = std::system(tmp); if(t_cmd == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    }

  } else if(fast_file_t==fast_file_c){
    std::string Y,y,N,n;
    std::cerr<<std::endl<<"One Sample test provided, enrichement will be computed versus average. ";
    std::cerr<<std::endl<<"Do you want to continue (Y/N)? ";
    char test;
    std::string test_cond;
    std::stringstream ss;
    std::cin.get(test);
    ss<<test;
    ss>>test_cond;
    if(test_cond=="Y" || test_cond=="y"){
        std::string new_out_file = to_string(resolution) + "_" + outputFileName;
        std::cerr<< "Computing hygestat one sample for window : "<<WINDOWSIZE<<std::endl;
        hgStatsOneSample(Dir1,chrReadsFileName,new_out_file,WINDOWSIZE,WINDOWADVANCE);
        if(!no_bedgraph_file) {
                std::string wig_file = create_wig_file_with_log_pvalue(new_out_file,fast_file_t,fast_file_c,resolution);
                std::sprintf(tmp, "%s %s %s %s", "mv ", wig_file.c_str(), new_out_file.c_str(),out3.c_str());
        }else{
                std::sprintf(tmp, "%s %s %s", "mv ", new_out_file.c_str(),out3.c_str());
        }
        int c_cmd = std::system(tmp); if(c_cmd == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    } else if(test_cond=="N" || test_cond=="n"){
      std::cerr<<"Thanks for using InstSeq from Rowicka lab (UTMB Health). "<<std::endl;
    }
  }
  }
}



std::string create_wig_file_with_log_pvalue(std::string file, std::string fastq_t, std::string fastq_c, int resolution) {
    fastq_t.substr(fastq_t.find_last_of("/")+1);
    fastq_c.substr(fastq_c.find_last_of("/")+1);
    std::vector <std::string> temp1, temp2;
    boost::split(temp1,fastq_t,boost::is_any_of("."));
    boost::split(temp2,fastq_c,boost::is_any_of("."));
    std::string                 outfile = temp1[0]+"_VS_"+temp2[0]+"_at_"+to_string(resolution)+".bedgraph";
    if(fastq_t == fastq_c) outfile = temp1[0]+"_at_"+to_string(resolution)+".bedgraph";
    std::ifstream in(file.c_str());
    std::ofstream out(outfile.c_str());
    out<<"browser pack refGene encodeRegions\n";
    out<<"browser hide all\n";
    out<<"browser full altGraph\n";
    std::string name = temp1[0]+"_VS_"+temp2[0]+"_at_"+to_string(resolution);
    out<<"track type=bedGraph name=\""<<name<<"\" description=\"DSB detected at "<< resolution<<"bp resolution\" visibility=full color=0,64,255 altColor=255,0,128 autoScale=on gridDefault=on graphType=line smoothingWindow=16 priority=20\n";
    std::string line;
    std::vector <std::string> temp;
    while(getline(in,line)) {
      boost::split(temp,line,boost::is_space());
      double p    = std::atof(temp[5].c_str());
      double pval = p == 2 ? 0 : -log10(p);
      out<<temp[0]<<"\t"<<temp[1]<<"\t"<<temp[2]<<"\t"<<pval<<std::endl;
    }
return  outfile;
}






