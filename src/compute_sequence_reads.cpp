
#include "../include/compute_sequence_reads.h"
//******************************************************************************//
double getSequenceReads(std::string Dir, std::vector<SEQUENCE_READS> &sequenceReads,unsigned int chrNum){
//----------------------------------------------------------------------------------------------
  if (checkFile(Dir+"/chr"+convert_chrnum_to_string(chrNum)+"+.hits.txt")!=1){ std::cerr << Dir+"/chr"+convert_chrnum_to_string(chrNum)+"+.hits.txt"; terminate(chrNum);}
  if (checkFile(Dir+"/chr"+convert_chrnum_to_string(chrNum)+"-.hits.txt")!=1){ std::cerr << Dir+"/chr"+convert_chrnum_to_string(chrNum)+"-.hits.txt"; terminate(chrNum);}
  std::ifstream chrFileW;
  chrFileW.open((Dir+"/chr"+convert_chrnum_to_string(chrNum)+"+.hits.txt").c_str());
  std::ifstream chrFileC;
  chrFileC.open((Dir+"/chr"+convert_chrnum_to_string(chrNum)+"-.hits.txt").c_str());
  std::vector <SEQUENCE_READS> reads;
  double sequenceReadsCount=0;
  while (!chrFileW.eof()){
    std::string line="";
    getline(chrFileW,line);
    std::stringstream ss(line);
    SEQUENCE_READS tmp;
    ss >> tmp.hits >> tmp.chrPos;
    if(!chrFileW.eof()){
    reads.push_back(tmp);
    }
    sequenceReadsCount+=tmp.hits;
  }
  while (!chrFileC.eof()){
    std::string line="";
    getline(chrFileC,line);
    std::stringstream ss(line);
    SEQUENCE_READS tmp;
    ss >> tmp.hits >> tmp.chrPos;
    if(!chrFileC.eof()){
    reads.push_back(tmp);
    }
    sequenceReadsCount+=tmp.hits;
  }
  sort(reads.begin(),reads.end(),sortSeqReadsFun);
// for (std::vector <SEQUENCE_READS>::iterator i=reads.begin();i<reads.end();++i)
//     std::cout << std::endl << (*i).chrPos;
//   WAITUSER;
  sequenceReads=reads;
  if(nature_of_analysis=="without_correction"){
   return 0;
   }
  else if(nature_of_analysis=="with_correction"){
  return sequenceReadsCount;
  }
}

//********************************************************************************************************
double calcAverageWindowReads(unsigned chrNum, std::string chrReadsFileName, std::string Dir3, std::vector<MAPPABILITY_DATA>&mappabilityMapList){
//--------------------------------------------------------------------------------------------------------
      const unsigned int WINDOWLIMIT=WINDOWSIZE*10;
      unsigned int numMappableWindows=0;
      unsigned int numWindows=0;
      std::vector <CHR_READS> ChrList1;
       if(DEBUG) std::cerr<<std::endl<<"Total Reads in Correction Map = "<<getChrReads(Dir3+chrReadsFileName,ChrList1);
      unsigned int currentChrNum=chrNum;
       if(DEBUG) std::cerr<<std::endl<<"Chr"<<chrNum;
      unsigned int Reads1=0;
      for (std::vector <CHR_READS>::iterator iC1=ChrList1.begin();iC1!=ChrList1.end();++iC1){
	  if ((*iC1).chrNum==currentChrNum)
	  {
	    Reads1=(*iC1).reads;
	  }
      }
      std::vector <SEQUENCE_READS> SequenceReads1;
      getSequenceReads(Dir3,SequenceReads1,currentChrNum);
      std::vector <SEQUENCE_READS>::iterator iCorrectionMapStart=SequenceReads1.begin();
      std::vector<bool> mappabilityMap;
      for (std::vector<MAPPABILITY_DATA>::iterator iM=mappabilityMapList.begin();iM!=mappabilityMapList.end();++iM){
	if ((*iM).chrNum==currentChrNum){
	    mappabilityMap = (*iM).mappabilityMap;
	  }
      }
      if (mappabilityMap.size()==0){std::cerr << "Mappability map not found Chr " <<currentChrNum; terminate(10);}
      unsigned int limit = mappabilityMap.size()-1;
       if(DEBUG) std::cerr << std::endl << convert_chrnum_to_string(currentChrNum) << " lim = " <<  limit << " "; //WAITUSER;
      unsigned totalReads=0;
      for (unsigned int chrPos=0; ; \
	  chrPos=getMappableWindowEnd(chrPos, WINDOWADVANCE, mappabilityMap))
      {
	    numMappableWindows++;
	    unsigned int windowReads1=0;
    	    unsigned int windowStart=chrPos;
    	    unsigned int windowEnd=getMappableWindowEnd(windowStart, WINDOWSIZE, mappabilityMap);
	    if (windowEnd>=limit) break;
	    if ((windowEnd-windowStart) >= WINDOWLIMIT){continue;}
    	    unsigned int chrPos1=0;
	    {
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iCorrectionMapStart; iC!=SequenceReads1.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		      if(firstHitFound==0){
			firstHitFound=1;
			iCorrectionMapStart=iC;
		      }
		      hitsFound=1;
		      windowReads1+=(*iC).hits;
		      totalReads+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//std::cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }

      }
      if(DEBUG) std::cerr<<std::endl<<"Number of Mapable Windows="<<numMappableWindows;
      if(DEBUG) std::cerr<<std::endl<<"Total number of reads="<<totalReads;
       if(DEBUG)std::cerr<<std::endl<<"Average number of reads / window="<<totalReads*1.0/numMappableWindows;
//       WAITUSER
      return totalReads*1.0/numMappableWindows;
  }
      //}

