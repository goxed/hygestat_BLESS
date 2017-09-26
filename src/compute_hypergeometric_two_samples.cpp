
#include "../include/compute_hypergeometric_two_samples.h"

//*************************************************************************************************************************
int hgStatsTwoSampleWith(std::string Dir1, std::string Dir2, std::string Dir3, std::string Dir4, std::string Dir5, std::string chrReadsFileName, std::string outputFileName, unsigned   WINDOWSIZE, unsigned WINDOWADVANCE) {
//-------------------------------------------------------------------------------------------------------------------------
    std::vector <CHR_READS> ChrList1, ChrList2;
    getChrReads(Dir1+chrReadsFileName,ChrList1);
    getChrReads(Dir2+chrReadsFileName,ChrList2);
    if (ChrList1.size()!=ChrList2.size()){
      if(DEBUG) std::cerr << std::endl << ChrList1.size() << "  " << ChrList2.size() << " Chromosome size mismatch"; WAITUSER;
    }
    std::vector <MAPPABILITY_DATA> mappabilityMapList;
    unsigned int mapableBasesTotal=getMappabilityMap(Dir1+chrReadsFileName, mappabilityMapList);
    if (ChrList2.size()>MAXTHREADS){std::cerr<<std::endl<<"maxthreads exceeded!";}
    std::vector<WINDOW_DATA>outData;
    pthread_t a_thread[MAXTHREADS];
    unsigned int threadNum=0;
    struct threadData threadArgs[MAXTHREADS];
    std::string logFileName=outputFileName+"_log_hgwindows.txt";
    std::ofstream logFile(logFileName.c_str());
    logFile<<"Two sample test"<<std::endl;
    logFile<<"Using "<<WINDOWSIZE<<" bp mappable window"<< " and " <<WINDOWADVANCE<< " bp WINDOW ADVANCE " <<std::endl;
    for (std::vector <CHR_READS>::reverse_iterator i=ChrList2.rbegin();i!=ChrList2.rend();++i,++threadNum)//iterating over chromosomes
    {
      threadArgs[threadNum].chrNum=(*i).chrNum;
      threadArgs[threadNum].chrReads=(*i).reads;
      threadArgs[threadNum].chrReadsFileName=chrReadsFileName;
      threadArgs[threadNum].mapableBasesTotal=mapableBasesTotal;
      threadArgs[threadNum].Dir1=Dir1;
      threadArgs[threadNum].Dir2=Dir2;
      threadArgs[threadNum].Dir3=Dir3;
      threadArgs[threadNum].Dir4=Dir4;
      threadArgs[threadNum].Dir5=Dir5;
      threadArgs[threadNum].mappabilityMapList=&mappabilityMapList;
      threadArgs[threadNum].outData=&outData;
//       threadArgs[threadNum].averageWindowReads=calcAverageWindowReads((*i).chrNum, chrReadsFileName, Dir3, mappabilityMapList);
//       threadArgs[threadNum].averageWindowReadsCNV=calcAverageWindowReads((*i).chrNum, chrReadsFileName, Dir4, mappabilityMapList);
     if(DEBUG) std::cerr<<std::endl<<threadArgs[threadNum].chrNum;
     if(DEBUG) std::cerr<<std::endl<<threadArgs[threadNum].chrReads;
      int p=0;
      if(p=pthread_create(&a_thread[threadNum], NULL, threadFunc, (void *)(&threadArgs[threadNum]))){
	if(DEBUG) std::cerr<<std::endl<<"Thread Creation Failed! " << p ;
      };      /*pthread_join(a_thread[threadNum],NULL);*///uncomment for single threaded execution
     }
     for (unsigned int t=0;t<threadNum;++t){
      pthread_join(a_thread[t],NULL);
     }
    if(DEBUG) std::cerr<<std::endl<<"Hypergeometric done";
    logFile<<"Outdata Size="<<outData.size();
    {//pval
        omp_set_num_threads(12);
      __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataPvalFun);
      std::vector<double> pvalList;
      for(std::vector<WINDOW_DATA>::iterator iw=outData.begin();iw!=outData.end();++iw){
	if((*iw).pval<2.0){
	  pvalList.push_back((*iw).pval);
	}
      }
      std::vector<double>qValList;
      bh(pvalList,qValList);
      if(pvalList.size()!=qValList.size()){
	if(DEBUG) std::cerr<<std::endl<<"Oops BH procedure failed!";
	terminate(0);
      }
      logFile<<"pval List Size="<<pvalList.size();
      std::vector<WINDOW_DATA>::iterator iw=outData.begin();
      std::vector<double>::iterator q=qValList.begin();
      for(;iw!=outData.end(),q!=qValList.end();++iw,++q){
	(*iw).qval=(*q);
      }
      __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataWindowStartFun);
      __gnu_parallel::stable_sort(outData.begin(),outData.end(),SortOutDataChrNumFun);
    }

    {//pvalCorrected
          omp_set_num_threads(12);
      __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataPvalCorrectedFun);
      std::vector<double> pvalList;
      for(std::vector<WINDOW_DATA>::iterator iw=outData.begin();iw!=outData.end();++iw){
	if((*iw).pvalCorrected<2.0){
	  pvalList.push_back((*iw).pvalCorrected);
	}
      }
      std::vector<double>qValList;
      bh(pvalList,qValList);
      if(pvalList.size()!=qValList.size()){
	if(DEBUG) std::cerr<<std::endl<<"Oops BH procedure failed on pval corrected!";
	terminate(01);
      }
      logFile<<"pval Corrected List Size="<<pvalList.size();
      std::vector<WINDOW_DATA>::iterator iw=outData.begin();
      std::vector<double>::iterator q=qValList.begin();
      for(;iw!=outData.end(),q!=qValList.end();++iw,++q){
	(*iw).qvalCorrected=(*q);
      }
      __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataWindowStartFun);
      __gnu_parallel::stable_sort(outData.begin(),outData.end(),SortOutDataChrNumFun);
    }

    std::vector<WINDOW_DATA>::iterator iw=outData.begin();
    boost::iostreams::stream_buffer<boost::iostreams::file_sink> outputFileSink(outputFileName.c_str());
    std::ostream outputFile(&outputFileSink);
    outputFile<<"chrnum"<<"\t"<<"windowStart"<<"\t"<<"windowEnd"<<"\t" \
	      <<"rdBLESStrt"<<"\t"<<"rdCntr"<<"\t"<<"pval"<<"\t"<<"qval"<<"\t"<<"SeqNormlz"<<"\t" \
	      <<"CNVcell"<<"\t"<<"CNVcellNormlz"<<"\t"<<"CNVcellINSpval"<<"\t"<<"CNVcellDELpval"<<"\t" \
	      <<"CNVtrtcell"<<"\t"<<"CNVtrtcellNormlz"<<"\t"<<"CNVtrtcellINSpval"<<"\t"<<"CNVtrtcellDELpval"<<std::endl;
    for(;iw!=outData.end();++iw){
      outputFile<<"chr"<<convert_chrnum_to_string((*iw).chrNum)<<"\t"<<(*iw).windowStart<<"\t"<<(*iw).windowEnd<<"\t" \
      <<(*iw).windowReads1<<"\t"<<(*iw).windowReads2<<"\t"<<(*iw).pvalCorrected<<"\t"<<(*iw).qvalCorrected<<"\t"<<(*iw).correctionCoefficientSequencing<<"\t";
     if((*iw).pvalCNVIns * outData.size()<=0.05){
	outputFile<<"INS_CNV1";
      }
      if((*iw).pvalCNVDel * outData.size()<=0.05){
	outputFile<<"DEL_CNV1";
      }
      if((*iw).pvalCNVDel * outData.size() >0.05 && (*iw).pvalCNVIns * outData.size() >0.05){
	outputFile<<"OK_CNV1";
      }
      outputFile<<"\t"<<(*iw).correctionCoefficientCNV<<"\t"<<(*iw).pvalCNVIns<<"\t"<<(*iw).pvalCNVDel<<"\t";
      if((*iw).pvalCNV2Ins * outData.size()<=0.05){
	outputFile<<"INS_CNV2";
      }
      if((*iw).pvalCNV2Del * outData.size()<=0.05){
	outputFile<<"DEL_CNV2";
      }
      if((*iw).pvalCNV2Del * outData.size() >0.05 && (*iw).pvalCNV2Ins * outData.size() >0.05){
	outputFile<<"OK_CNV2";
      }
      outputFile<<"\t"<<(*iw).correctionCoefficientCNV2<<"\t"<<(*iw).pvalCNV2Ins<<"\t"<<(*iw).pvalCNV2Del<<"\t";
      outputFile<<std::endl;

    }
    return 0;
  }
//**************************************** Hypergeometric two samples ******************************************************************************//
int hgStatsTwoSampleWithout(std::string Dir1, std::string Dir2, std::string chrReadsFileName, std::string outputFileName, unsigned   WINDOWSIZE, unsigned WINDOWADVANCE) {
//--------------------------------------------------------------------------------------------------------------------------------------------------//
    std::vector <CHR_READS> ChrList1, ChrList2;
    getChrReads(Dir1+chrReadsFileName,ChrList1);
    getChrReads(Dir2+chrReadsFileName,ChrList2);
    if (ChrList1.size()!=ChrList2.size()){
      std::cerr << std::endl << ChrList1.size() << "  " << ChrList2.size() << " Chromosome size mismatch"; WAITUSER;
    }
    std::vector <MAPPABILITY_DATA> mappabilityMapList;
    getMappabilityMap(Dir1+chrReadsFileName, mappabilityMapList);
    if (ChrList2.size()>MAXTHREADS){std::cerr<<std::endl<<"maxthreads exceeded!";}
    std::vector<WINDOW_DATA>outData;
    pthread_t a_thread[MAXTHREADS];
    unsigned int threadNum=0;
    struct threadData threadArgs[MAXTHREADS];
    std::string logFileName=outputFileName+"_log_hgwindows.txt";
    std::ofstream logFile(logFileName.c_str());
    logFile<<"Two sample test"<<std::endl;
    logFile<<"Using "<<WINDOWSIZE<<" bp mappable window"<< " and " <<WINDOWADVANCE<< " bp WINDOW ADVANCE " <<std::endl;
    for (std::vector <CHR_READS>::reverse_iterator i=ChrList2.rbegin();i!=ChrList2.rend();++i,++threadNum)//iterating over chromosomes
    {
      threadArgs[threadNum].chrNum=(*i).chrNum;
      threadArgs[threadNum].chrReads=(*i).reads;
      threadArgs[threadNum].chrReadsFileName=chrReadsFileName;
      threadArgs[threadNum].Dir1=Dir1;
      threadArgs[threadNum].Dir2=Dir2;
      threadArgs[threadNum].mappabilityMapList=&mappabilityMapList;
      threadArgs[threadNum].outData=&outData;
      if(DEBUG) std::cerr<<std::endl<<threadArgs[threadNum].chrNum;
      if(DEBUG) std::cerr<<std::endl<<threadArgs[threadNum].chrReads;
      int p=0;
      if(p=pthread_create(&a_thread[threadNum], NULL, threadFunc, (void *)(&threadArgs[threadNum]))){
	             std::cerr<<std::endl<<"Thread Creation Failed! " << p ;
      };      //pthread_join(a_thread[threadNum],NULL);//uncomment for single threaded execution
    }
     for (unsigned int t=0;t<threadNum;++t){
      pthread_join(a_thread[t],NULL);
    }
    if(DEBUG) std::cerr<<std::endl<<"Hypergeometric done";
    omp_set_num_threads(N_P);
    __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataPvalFun);
    logFile<<"Outdata Size="<<outData.size();
    std::vector<double> pvalList;
    for(std::vector<WINDOW_DATA>::iterator iw=outData.begin();iw!=outData.end();++iw){
      if((*iw).pval<2.0){
	pvalList.push_back((*iw).pval);
      }
    }
    std::vector<double> pvalList1;
    for(std::vector<WINDOW_DATA>::iterator iw=outData.begin();iw!=outData.end();++iw){
        pvalList1.push_back((*iw).pval);
    }
    std::vector<double>bValList;
    std::vector<double>qValList;
    bh(pvalList,qValList);
    bonf(pvalList1,bValList);
    if(pvalList.size()!=qValList.size() || pvalList1.size()!=bValList.size()){
      std::cerr<<std::endl<<"Oops BH procedure failed!";
      terminate(0);
    }
    logFile<<"pval List Size="<<pvalList.size();
    std::vector<WINDOW_DATA>::iterator iw=outData.begin();
    std::vector<double>::iterator q=qValList.begin();
    for(;iw!=outData.end(),q!=qValList.end();++iw,++q){
      (*iw).qval=(*q);
    }
    iw=outData.begin();
    std::vector<double>::iterator b=bValList.begin();
    for(;iw!=outData.end(),b!=bValList.end();++iw,++b){
      (*iw).bval=(*b);
    }
    __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataWindowStartFun);
    __gnu_parallel::stable_sort(outData.begin(),outData.end(),SortOutDataChrNumFun);
    iw=outData.begin();
    boost::iostreams::stream_buffer<boost::iostreams::file_sink> outputFileSink(outputFileName.c_str());
    std::ostream outputFile(&outputFileSink);
    for(;iw!=outData.end();++iw){
      outputFile<<"chr"<<convert_chrnum_to_string((*iw).chrNum)<<"\t"<<(*iw).windowStart<<"\t"<<(*iw).windowEnd<<"\t" \
      <<(*iw).windowReads1<<"\t"<<(*iw).windowReads2<<"\t"<<(*iw).pval<<"\t"<<(*iw).qval<<"\t"<<(*iw).bval<<std::endl;
    }
    return 0;
 }
//******************** One sample hypergeometric test **************************************************************
int hgStatsOneSample(std::string Dir1, std::string chrReadsFileName, std::string outputFileName, unsigned WINDOWSIZE, unsigned WINDOWADVANCE) {
  //------------------------------------------------------------------------------------------------------------------
    std::vector <CHR_READS> ChrList1;
    getChrReads(Dir1+chrReadsFileName,ChrList1);
    std::vector <MAPPABILITY_DATA> mappabilityMapList;
    unsigned int mapableBasesTotal=getMappabilityMap(Dir1+chrReadsFileName, mappabilityMapList);
    if (ChrList1.size()>MAXTHREADS){std::cerr<<std::endl<<"maxthreads exceeded!"; terminate(MAXTHREADS);};
    std::vector<WINDOW_DATA>outData;
    pthread_t a_thread[MAXTHREADS];
    unsigned int threadNum=0;
    struct threadData threadArgs[MAXTHREADS];
    std::string logFileName=outputFileName+"_log_hgwindows.txt";
    std::ofstream logFile(logFileName.c_str());
    logFile<<"One sample test"<<std::endl;
    logFile<<"Using "<<WINDOWSIZE<<" bp mappable window"<< " and " <<WINDOWADVANCE<< " bp WINDOW ADVANCE " <<std::endl;
    for (std::vector <CHR_READS>::reverse_iterator i=ChrList1.rbegin();i!=ChrList1.rend();++i,++threadNum)//iterating over chromosomes
    {
      threadArgs[threadNum].chrNum=(*i).chrNum;
      threadArgs[threadNum].chrReads=(*i).reads;
      threadArgs[threadNum].chrReadsFileName=chrReadsFileName;
      threadArgs[threadNum].Dir1=Dir1;
      threadArgs[threadNum].mappabilityMapList=&mappabilityMapList;
      threadArgs[threadNum].outData=&outData;
      threadArgs[threadNum].mapableBasesTotal=mapableBasesTotal;
      std::cerr<<std::endl<<threadArgs[threadNum].chrNum;
      std::cerr<<std::endl<<threadArgs[threadNum].chrReads;
      int p=0;
      if(p=pthread_create(&a_thread[threadNum], NULL, threadFunc1, (void *)(&threadArgs[threadNum]))){
	std::cerr<<std::endl<<"Thread Creation Failed! " << p ;
      };
      //pthread_join(a_thread[threadNum],NULL);//uncomment for single threaded execution
    }
     for (unsigned int t=0;t<threadNum;++t){
      pthread_join(a_thread[t],NULL);
    }
    std::cerr<<std::endl<<"Hypergeometric done";
    omp_set_num_threads(N_P);
    __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataPvalFun);
    logFile<<"Outdata Size="<<outData.size();
    std::vector<double> pvalList;
    for(std::vector<WINDOW_DATA>::iterator iw=outData.begin();iw!=outData.end();++iw){
      if((*iw).pval<2.0){
	pvalList.push_back((*iw).pval);
      }
    }
    std::vector<double> pvalList1;
    std::vector<double>qValList;
    bh(pvalList,qValList);
    if(pvalList.size()!=qValList.size()){
      std::cerr<<std::endl<<"Oops BH procedure failed!";
      terminate(0);
    }
    std::cerr<<std::endl<<"BH done\n";
    std::vector<WINDOW_DATA>::iterator iw=outData.begin();
    std::vector<double>::iterator q=qValList.begin();
    for(;iw!=outData.end(),q!=qValList.end();++iw,++q){
      (*iw).qval=(*q);
    }
     __gnu_parallel::sort(outData.begin(),outData.end(),SortOutDataWindowStartFun);
    __gnu_parallel::stable_sort(outData.begin(),outData.end(),SortOutDataChrNumFun);
    iw=outData.begin();
    boost::iostreams::stream_buffer<boost::iostreams::file_sink> outputFileSink(outputFileName.c_str());
    std::ostream outputFile(&outputFileSink);
    for(;iw!=outData.end();++iw){
      outputFile<<"chr"<<convert_chrnum_to_string((*iw).chrNum)<<"\t"<<(*iw).windowStart<<"\t"<<(*iw).windowEnd<<"\t" \
      <<(*iw).windowReads1<<"\t"<<(*iw).windowReads2<<"\t"<<(*iw).pval<<"\t"<<(*iw).qval<<std::endl;
    }
    return 0;
}
