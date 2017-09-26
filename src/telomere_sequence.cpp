
#include "../include/telomere_sequence.h"

//***************** Extract bar code **************************
std::string extract_telomere(std::string line){
  //----------------------------------------------------------------
  std::string telom;

  if((org_name=="mouse")||(org_name=="human")) telom = telomere_ver;
  else if (org_name=="yeast")  telom =telomere_yea;
  else std::cerr<< "Please give the telomere sequence\n";

  while(line.find(telom)!=std::string::npos)  line.erase(line.find(telom),telom.size());

  return line;

}


bool contain_telomere(std::string line){
    std::string telom;
    if((org_name=="mouse")||(org_name=="human")) telom = telomere_ver;
    else if (org_name=="yeast")  telom =telomere_yea;
    else std::cerr<< "Please give the telomere sequence";
    bool tel = false;
    if(line.find(telom)!=std::string::npos) tel = true;
    return tel;
}
// ***************** Extract telomere bar code **************************//
void extract_telomere_file(std::string fastFile){
  //----------------------------------------------------------------
  std::string telom;
  std::string    file1=fastFile+".dump.fa";
  std::ifstream  fastaf(file1.c_str());
  std::string    outpu1 = fastFile+".dump.fa."+telom+".fasta";
  std::string line;

  if((org_name=="mouse")||(org_name=="human")) telom = telomere_ver;
  else if (org_name=="yeast")  telom =telomere_yea;
  else std::cerr<< "Please give the telomere sequence";

  
  std::ofstream  telom_close(outpu1.c_str());

  
  while(getline(fastaf,line)){
    line = extract_telomere(line);
    telom_close<<line<<std::endl;
  }
  fastaf.close();
}
