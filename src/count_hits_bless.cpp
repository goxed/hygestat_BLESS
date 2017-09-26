#include "../include/count_hits_bless.h"

//********************************Running hits for bless data ****************************//
void function_hits_bless(std::string file_from_bowtie,std::string nature)     {
//****************************************************************************************//

    std::string input1=file_from_bowtie+".dump.fa.TCGAGGTAGTA.fasta.bt.btt";
    std::string input2=file_from_bowtie+".dump.fa.TCGAGACGACG.fasta.bt.btt";

    std::string output1=file_from_bowtie + "_close_barcode";
    std::string output2=file_from_bowtie +"_distant_barcode";

//    std::string output1=file_from_bowtie+"_"+to_string(resolution) + "_close_barcode";
//    std::string output2=file_from_bowtie+"_"+to_string(resolution) +"_distant_barcode";

    boost::filesystem::path   path_c(output1);
    boost::filesystem::path   path_d(output2);
    char cm1[255],cm2[255],cm3[255],cm4[255],cm5[255],cm6[255];
    char cm1_tel[255],cm2_tel[255],cm3_tel[255];
    std::sprintf(cm1,"%s %s %s ",link1.c_str(),input1.c_str(),org_name.c_str());
    std::sprintf(cm2,"%s %s"    ,"mkdir ",output1.c_str());
    std::sprintf(cm3,"%s %s %s" ,"mv","chr* ",output1.c_str());

    std::sprintf(cm4,"%s %s %s ",link1.c_str(),input2.c_str(),org_name.c_str());
    std::sprintf(cm5,"%s %s"    ,"mkdir ",output2.c_str());
    std::sprintf(cm6,"%s %s %s" ,"mv ","chr* ",output2.c_str());

    if(TELOMERE_FILE) {
        std::string tel;
        if((org_name=="mouse")||(org_name=="human")) tel = telomere_ver;
        else if (org_name=="yeast")  tel =telomere_yea;
        else std::cerr<< "Please give the telomere sequence\n";

        std::string input_tel  = file_from_bowtie+".dump.fa."+tel+".fasta.bt.btt";
        std::string output_tel = file_from_bowtie+"_telomere_barcode";

        std::sprintf(cm1_tel,"%s %s %s ",link1.c_str(),input_tel.c_str(),org_name.c_str());
        std::sprintf(cm2_tel,"%s %s"    ,"mkdir ",output_tel.c_str());
        std::sprintf(cm3_tel,"%s %s %s" ,"mv","chr* ",output_tel.c_str());
        int t_cmd1 = std::system(cm1_tel); if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
        int t_cmd2 = std::system(cm2_tel); if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
        int t_cmd3 = std::system(cm3_tel); if(t_cmd3 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
    }

    if(nature=="close"){
     //   if (!exists(path_c)){
            int t_cmd1 = std::system(cm1);  if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd2 = std::system(cm2);  if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd3 = std::system(cm3);  if(t_cmd3 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
       // }
    }else if(nature=="distant"){
       // if (!exists(path_d)){
            int t_cmd4 = std::system(cm4);  if(t_cmd4 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd5 = std::system(cm5);  if(t_cmd5 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd6 = std::system(cm6);  if(t_cmd6 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
      //  }
    }else{
     //  if((!exists(path_c))&&(!exists(path_d))){
            int t_cmd1 = std::system(cm1);  if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd2 = std::system(cm2);  if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd3 = std::system(cm3);  if(t_cmd3 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd4 = std::system(cm4);  if(t_cmd4 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd5 = std::system(cm5);  if(t_cmd5 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            int t_cmd6 = std::system(cm6);  if(t_cmd6 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
     //  }
    }
}
