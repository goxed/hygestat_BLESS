

#include "../include/bowtie_alignments.h"
//***************** bowtie function ***************************
unsigned int bowtie_function(std::string fastafile, std::string nature_of_data) {
//--------------------------------------------------------------

   std::string file11=fastafile;
   std::string file12=fastafile;
   std::string file13=fastafile;
   std::string genome="";

   file11.insert(file11.size(),".dump.fa."+close_bar+".fasta");

   file12.insert(file12.size(),".dump.fa."+distant_bar+".fasta");

   if(org_name=="mouse" || org_name=="MOUSE"){
   genome = genome_mouse;
   }
   else if(org_name=="human" || org_name=="HUMAN"){
   genome = genome_human;
   }
   else if(org_name=="yeast" || org_name=="YEAST"){
   genome = genome_yeast;
   }
   else{
    if(DEBUG) std::cout<<"organism non identify: human by default "<<std::endl;
    genome = genome_human;
    }

   if(TELOMERE_FILE) {
      std::string telom;
      std::string file_tel=fastafile;
      std::string file_out_tel = file_tel+".bt";
      if((org_name=="mouse")||(org_name=="human")) telom = telomere_ver;
      else if (org_name=="yeast")  telom =telomere_yea;
      else std::cerr<< "Please give the telomere sequence\n";
      file_tel.insert(file_tel.size(),".dump.fa."+telom+".fasta");
      char cmd1_tel[250],cmd2_tel[250];
      std::sprintf(cmd1_tel,"%s" "%s" "%s" "%s" "%s" "%s",bowt.c_str(),genome.c_str(),bowt_arg.c_str(),file_tel.c_str(),">",file_out_tel.c_str());
      int t_cmd1 = std::system(cmd1_tel); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the mapping file...\n";
      std::sprintf(cmd2_tel, "%s" "%s", btt.c_str(), file_out_tel.c_str());
      int t_cmd2 = std::system(cmd2_tel); if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
   }

   if(DEBUG) {
            std::cerr<<"running bowtie using "<<N_P<<" CPU"<<std::endl;
            std::cerr<<"checking if processor is available..."<<std::endl;
            if(std::system(NULL)) puts ("ok");
            else exit (EXIT_FAILURE);
            }
   char cmd1[250];
   char cmd2[250];
   char cmd11[250];
   char cmd22[250];

   std::string file_out1 = file11+".bt";
   std::string file_out2 = file12+".bt";
   std::string file_out3 = file13+".bt";
   if(nature_of_data=="bless"){
   std::sprintf(cmd1,"%s" "%s" "%s" "%s" "%s" "%s",bowt.c_str(),genome.c_str(),bowt_arg.c_str(),file11.c_str(),">",file_out1.c_str());

   int r_cmd1 = std::system(cmd1); if(r_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";

   std::sprintf(cmd11, "%s" "%s", btt.c_str(), file_out1.c_str());

   int r_cmd11 = std::system(cmd11);if(r_cmd11 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";

   std::sprintf(cmd2,"%s" "%s"  "%s" "%s" "%s" "%s",bowt.c_str(),genome.c_str(),bowt_arg.c_str(),file12.c_str(),">",file_out2.c_str());

   if (DEBUG) std::cout<<"Executing command DIR..."<<std::endl;

   int r_cmd2 =  std::system(cmd2);if(r_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";

    std::sprintf(cmd22, "%s" "%s", btt.c_str(), file_out2.c_str());

   int r_cmd22 =  std::system(cmd22); if(r_cmd22 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
   } else if(nature_of_data=="genomic"){
   std::sprintf(cmd1,"%s" "%s" "%s" "%s" "%s" "%s",bowt.c_str(),genome.c_str(),bowt_arg.c_str(),file13.c_str(),">",file_out3.c_str());

   int r_cmd1 = std::system(cmd1); if(r_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";

   std::sprintf(cmd11, "%s" "%s", btt.c_str(), file_out3.c_str());

   int r_cmd11 = std::system(cmd11);if(r_cmd11 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
   }
    return 0;
}

//*************** bowtie genomic ******************************//
unsigned int bowtie_function_genomic(std::string fastqfile) {
//--------------------------------------------------------------
    if(DEBUG) {
            std::cerr<<"running bowtie using "<<N_P<<" CPU"<<std::endl;
            std::cerr<<"checking if processor is available..."<<std::endl;
            if(std::system(NULL)) puts ("ok");
            else exit (EXIT_FAILURE);
            }
   char cmd1[250];
   char cmd11[250];
   std::string genome="";
   if(org_name=="mouse" || org_name=="MOUSE"){
   genome = genome_mouse;
   }
   else if(org_name=="human" || org_name=="HUMAN"){
   genome = genome_human;
   }
   else if(org_name=="yeast" || org_name=="YEAST"){
   genome = genome_yeast;
   }
   else{
    if(DEBUG) std::cout<<"organism non identify: human by default "<<std::endl;
    genome = genome_human;
    }
   std::string file_out1 = fastqfile+".bt";
   std::sprintf(cmd1,"%s" "%s" "%s" "%s" "%s" "%s",bowt.c_str(),genome.c_str(),bowt_arg_q.c_str(),fastqfile.c_str(),">",file_out1.c_str());
   int r_cmd1 =  std::system(cmd1); if(r_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
   std::sprintf(cmd11, "%s" "%s", btt.c_str(), file_out1.c_str());
   int r_cmd11 =  std::system(cmd11);if(r_cmd11 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    return 0;
}

