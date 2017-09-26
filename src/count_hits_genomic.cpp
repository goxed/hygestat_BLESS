
#include "../include/count_hits_genomic.h"
//********************running hits for genomic data*************************************//
void function_hits_genomic(std::string file_from_bowtie)   {
//**************************************************************************************//
    std::string input1 =file_from_bowtie+".bt.btt";
    std::string output1=file_from_bowtie+"_no_barcode";
    //std::string output1=file_from_bowtie+"_"+to_string(resolution)+"_no_barcode";
    char cm1[255];
    char cm2[255];
    char cm3[255];
    boost::filesystem::path   path_c(output1);
    std::sprintf(cm1,"%s ""%s ""%s",link2.c_str(),input1.c_str(),org_name.c_str());
    int t_cmd1 = std::system(cm1); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
     if (!exists(path_c)){
        std::sprintf(cm2,"%s %s","mkdir ",output1.c_str());
        int t_cmd2 = std::system(cm2); if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
        std::sprintf(cm3,"%s %s %s","mv ","chr* ",output1.c_str());
        int t_cmd3 = std::system(cm3); if(t_cmd3 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
     }
}
