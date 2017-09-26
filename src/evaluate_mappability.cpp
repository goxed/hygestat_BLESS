#include "../include/evaluate_mappability.h"


//***********************************function ******************************************************//
int getMappabilityMap(std::string chrNumFileName, std::vector <MAPPABILITY_DATA> &mappabilityMapList){
//**************************************************************************************************//
  #define READSIZE "45"
  std::ifstream controlChrFile(chrNumFileName.c_str());
  unsigned int mappableBases=0;
  std::vector<unsigned>chrNumList;
  while (!controlChrFile.eof()){
    std::string line="";
    getline (controlChrFile,line);
    std::string::size_type tab_location;
    tab_location=line.find_first_of('\t');
    unsigned int chrNum = extract_chrnum(line.substr(0,tab_location),true);
    if(chrNum!=0)chrNumList.push_back(chrNum);
  }
/*  omp_set_num_threads(20);
 #pragma omp parallel
  {
  #pragma omp for*/
  for(unsigned c=0;c<chrNumList.size();++c)
  {
    unsigned chrNum=chrNumList[c];

    if (chrNum!=0){
      unsigned int chrSize=0;
      unsigned int chrMapableBases=0;

      if(no_mappability){
            MAPPABILITY_DATA mapDataTmp;
            std::vector <std::string> chrNames;
            std::vector <unsigned int> chrSizes;
            if(org_name=="human" || org_name=="HUMAN"){
                    chrNames = chrN_human;
                    chrSizes = chrZ_human;
                    mappableBases= accumulate(chrSizes.begin(), chrSizes.end(), 0);
            }
            else if(org_name=="mouse" || org_name=="MOUSE"){
                    chrNames = chrN_mouse;
                    chrSizes = chrZ_mouse;
                    mappableBases= accumulate(chrSizes.begin(), chrSizes.end(), 0);
            }
            else if(org_name=="yeast" || org_name=="YEAST"){
                    chrNames = chrN_yeast;
                    chrSizes = chrZ_yeast;
                    mappableBases= accumulate(chrSizes.begin(), chrSizes.end(), 0);

            }
            else{std::cout<<std::endl<<"organism non identified"; terminate(5);}

       for(unsigned int jj=0;jj<chrNames.size();jj++){
            unsigned int chr_n = extract_chrnum(chrNames[jj],true);
            if(chr_n==chrNum){
                chrSize = chrSizes[jj];
                chrMapableBases = chrSize;
            }
       }
   std::vector<bool> numb(chrMapableBases,1);
   mapDataTmp.mappabilityMap = numb;
   mapDataTmp.chrSize=chrSize;
   mapDataTmp.chrNum=chrNum;
   mapDataTmp.readSize=45;
   mapDataTmp.chrMapableBases=chrMapableBases;
   mappabilityMapList.push_back(mapDataTmp);

  }
   else{

      std::string mapFileName=mapDir+"map_chr"+convert_chrnum_to_string2(chrNum)+"_"+"45"+".txt";
      if (checkFile(mapFileName)!=1){
	       std::cout << std::endl << "Cannot open " << mapFileName; WAITUSER;
	       terminate(chrNum);
      }
      MAPPABILITY_DATA mapDataTmp;
      boost::iostreams::stream_buffer<boost::iostreams::file_source> mapFileSource(mapFileName.c_str());
      std::istream mapFile(&mapFileSource);

      std::stringstream buffer;
      buffer <<  mapFile.rdbuf();
      std::string content(buffer.str());
      for (unsigned int n=0;n!=content.size();n+=2){
            bool mapVal=0;
            if (content[n]=='2'){
                mapVal=0;
            }
            else if (content[n]=='1'){
                mapVal=1;
                mappableBases++;
                chrMapableBases++;
            }
            else{ if(DEBUG) std::cerr << std::endl << "Warning nonstandard mapvalue of " << content[n] << " Assuming non mappable!"; }
	  mapDataTmp.mappabilityMap.push_back(mapVal);
	  ++chrSize;
      }
      mapDataTmp.chrSize=chrSize;
      mapDataTmp.chrNum=chrNum;
      mapDataTmp.readSize=45;
      mapDataTmp.chrMapableBases=chrMapableBases;
      mappabilityMapList.push_back(mapDataTmp);
      if(DEBUG) std::cerr << std::endl << "Loading Mapability Map of Chromosome " << convert_chrnum_to_string(chrNum) << " Size " << chrSize << ", Mapable " << chrMapableBases<<std::endl;
    }
   }
  }
   if(DEBUG) std::cerr<<std::endl<<"Total mappable bases = "<<mappableBases;
  return mappableBases;
}

//******************************************************************************//
unsigned int getMappableWindowEnd(unsigned int startChrPos, unsigned int mappableWindowSize, std::vector<bool>&mappabilityMap){
  unsigned int windowEnd=startChrPos;
  if(startChrPos>mappabilityMap.size()){
    return (mappabilityMap.size()-1);
  }
   for (unsigned int mappableCount=0;windowEnd<=mappabilityMap.size();++windowEnd){
      if (mappabilityMap[windowEnd]==1){
	++mappableCount;
      }
      if (mappableCount==mappableWindowSize) break;
    }
    return windowEnd;
}
