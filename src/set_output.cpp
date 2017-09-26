//#include "../include/header.hpp"
#include "../include/set_output.h"
//***************************************************************************************************
void set_ouput_files_and_stat(void) {
//---------------------------------------------------------------------------------------------------
    if(fexists("bless_and_sequency_quality.txt")) std::system("rm bless_and_sequency_quality.txt");
    if(fexists("telomere_stat.txt")) int r_cm1 = std::system("rm telomere_stat.txt");
    std::ofstream bl;
    bl.open ("bless_and_sequency_quality.txt");
    bl<<"Sample_name \t"<<"Total_number_of_reads\t"<<"Number_of_mappable_reads"<< "Number_of_barcoded_reads \t"<<"Not_barcoded_reads\t"<<"No_call_reads\n";
    if(TELOMERE){
        std::ofstream fs;
        fs.open ("telomere_stat.txt");
        fs<<"Fasta_file \t"<< "Number_of_barcoded_reads \t"<<"Number_of_telomeres_barcoded\n"; // open output file for headers
      }

}

