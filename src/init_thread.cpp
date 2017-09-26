
#include "../include/init_thread.h"
//*****************************************************
void *threadFunc(void *arg)
{
   struct threadData *threadArgs;
   threadArgs = (struct threadData *) arg;
   unsigned int chrNum=threadArgs->chrNum;
   unsigned int chrReads=threadArgs->chrReads;
  // if(nature_of_analysis=="with_correction"){
   double averageWindowReads=threadArgs->averageWindowReads;
   double averageWindowReadsCNV=threadArgs->averageWindowReadsCNV;
   unsigned int mapableBasesTotal=threadArgs->mapableBasesTotal;
   //std::string chrReadsFileName=threadArgs->chrReadsFileName;
   //std::string Dir1=threadArgs->Dir1;
   //std::string Dir2=threadArgs->Dir2;
   std::string Dir3=threadArgs->Dir3;
   std::string Dir4=threadArgs->Dir4;
   std::string Dir5=threadArgs->Dir5;
  // }
   std::string Dir1=threadArgs->Dir1;
   std::string Dir2=threadArgs->Dir2;
   std::string chrReadsFileName=threadArgs->chrReadsFileName;
   std::vector<MAPPABILITY_DATA>*mappabilityMapList=threadArgs->mappabilityMapList;
   std::vector<WINDOW_DATA>*outData=threadArgs->outData;
   if(nature_of_analysis=="with_correction"){
   calcHypergeometricPvalWith(chrNum,chrReads,averageWindowReads,averageWindowReadsCNV,mapableBasesTotal,chrReadsFileName,Dir1,Dir2,Dir3,Dir4,Dir5,*mappabilityMapList,*outData);
   } else if(nature_of_analysis=="without_correction"){
       calcHypergeometricPvalWithout(chrNum,chrReads,chrReadsFileName,Dir1,Dir2,*mappabilityMapList,*outData);
  }
   pthread_exit(NULL);
}
//********************************************************
void *threadFunc1(void *arg) {
//-----------------------------------------------------
   struct threadData *threadArgs;
   threadArgs = (struct threadData *) arg;
   unsigned int chrNum=threadArgs->chrNum;
   unsigned int chrReads=threadArgs->chrReads;
   std::string chrReadsFileName=threadArgs->chrReadsFileName;
   std::string Dir1=threadArgs->Dir1;
   std::vector<MAPPABILITY_DATA>*mappabilityMapList=threadArgs->mappabilityMapList;
   std::vector<WINDOW_DATA>*outData=threadArgs->outData;
   if(nature_of_analysis=="with_correction"){
   unsigned int mapableBasesTotal=threadArgs->mapableBasesTotal;
   calcHypergeometricPvalOneSample(chrNum,chrReads,mapableBasesTotal,chrReadsFileName,Dir1,*mappabilityMapList,*outData);
   }
   else if(nature_of_analysis=="without_correction"){
   unsigned int mappableBases=threadArgs->mappableBases;
   calcHypergeometricPvalOneSample(chrNum,chrReads,mappableBases,chrReadsFileName,Dir1,*mappabilityMapList,*outData);
   }
   pthread_exit(NULL);
}
