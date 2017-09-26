
#include "../include/barcodes_trimming.h"

// ***************** Extract bar code **************************//
void extract_barcode(std::string fastFile){
//----------------------------------------------------------------
     std::string    file1 = fastFile+".dump.fa";
     std::ifstream  fastaf(file1.c_str());
     std::string    outpu1 = fastFile+".dump.fa."+close_bar+".fasta";
     std::ofstream  out_close(outpu1.c_str());
     std::string    outpu2 = fastFile+".dump.fa."+distant_bar+".fasta";
     std::ofstream  out_distant(outpu2.c_str());
     std::string    line;
     std::fstream   fs;
     int            n_barcoded=0, tel_barcoded=0;
     size_t         bef;
     while(getline(fastaf,line)){
       if(line.compare(0,close_bar.size(),close_bar)==0){
            if((TELOMERE)&&(contain_telomere(line))) tel_barcoded++; // count the number of bless barcoded containg telomeres
            out_close<<line.substr(close_bar.size())<<std::endl;
            n_barcoded++;
       } else if (line.compare(0,close_link.size(),close_link)==0) {
            if((TELOMERE)&&(contain_telomere(line))) tel_barcoded++;
            out_close<<line.substr(close_link.size())<<std::endl;
            n_barcoded++;
       }
       if(line.compare(0,distant_bar.size(),distant_bar)==0){
            if((TELOMERE)&&(contain_telomere(line))) tel_barcoded++;
            out_distant<<line.substr(distant_bar.size())<<std::endl;
            n_barcoded++;
       } else if (line.compare(0,distant_link.size(),distant_link)==0) {
            if((TELOMERE)&&(contain_telomere(line))) tel_barcoded++;
            out_distant<<line.substr(distant_link.size())<<std::endl;
            n_barcoded++;
       }
     }
     fastaf.close();
     if(TELOMERE) {
       fs.open ("telomere_stat.txt", std::fstream::in | std::fstream::out | std::fstream::app);
       fs<<fastFile.erase(fastFile.find(".fastq"))<<" \t "<<n_barcoded<<" \t "<<tel_barcoded<<"\n";
     }

}
