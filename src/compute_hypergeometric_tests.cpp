
#include "../include/compute_hypergeometric_tests.h"
//********************************* Hypergeometric test ***************************************//
double hgcalc_gte(unsigned int ps, unsigned int sa,unsigned int sp,unsigned int ss){
//-------------------------------------------------------------------------------------------------------------
// ps = /*pop. size*/
// sa = /*sample size*/
// sp = /*Success Pop.*/
// ss = /*Success Sample*/
  double pval=1.0;
  if (ss>=1&&ps>=sa&&ps>=sp&&sa>=ss){
      pval=gsl_cdf_hypergeometric_Q(ss-1,sp,ps-sp,sa);
  }
  return pval;
}

//**********************************************************************************
unsigned int bh(std::vector<double>&pvalList,std::vector<double>&qValList)  {
  //*****************************
      unsigned int count=0;
      double previous_bhVal=0.0;
      unsigned listSize=pvalList.size();
      for(std::vector<double>::iterator p=pvalList.begin();p!=pvalList.end();++p,++count){
	#define DEBUGPVAL
	#ifdef DEBUGPVAL
	  if ((*p)>1.0 || (*p)<0.0){
	   if(DEBUG) std::cerr<<std::endl<<(*p)<<" bad pval!";
	  }
	#endif
	double bhVal=(*p)*listSize/(count+1);
	if(bhVal>1.0){
	  bhVal=1.0;
	}
	if(bhVal<previous_bhVal){
	  bhVal=previous_bhVal;
	}
	previous_bhVal=bhVal;
	qValList.push_back(bhVal);
      }
      return(qValList.size());
 }

unsigned int bonf(std::vector<double>&pvalList,std::vector<double>&bValList){
      double previous_bfVal=0.0;
      unsigned listSize=pvalList.size();
      for(std::vector<double>::iterator p=pvalList.begin();p!=pvalList.end();++p){
	/*#define DEBUGPVAL
  	#ifdef DEBUGPVAL
  		  if ((*p)>1.0 || (*p)<0.0){
  		  	    cerr<<endl<<(*p)<<" bad pval!";
  		  	    	  }
  		  	    	  	#endif */
	double bfVal=(*p)*listSize;
	if(bfVal>1.0){
	  bfVal=1.0;
	}
	if(bfVal<previous_bfVal){
	  bfVal=previous_bfVal;
	}
	previous_bfVal=bfVal;
	bValList.push_back(bfVal);
      }
      return(bValList.size());
 }

//********************************************** hyper with corr *************************************************//
 int calcHypergeometricPvalWith(unsigned int chrNum,unsigned int chrReads, double averageWindowReads, double averageWindowReadsCNV, unsigned mapableBasesTotal, \
  std::string chrReadsFileName, std::string Dir1, std::string Dir2, std::string Dir3, std::string Dir4, std::string Dir5, std::vector<MAPPABILITY_DATA>&mappabilityMapList, std::vector<WINDOW_DATA> &outData) {
 //---------------------------------------------------------------------------------------------------------------

      const unsigned int WINDOWLIMIT=WINDOWSIZE*10;
      unsigned int numMappableWindows=0;
      unsigned int numWindows=0;
      std::vector <CHR_READS> ChrList1;
      std::vector <CHR_READS> ChrList2;
      std::vector <CHR_READS> ChrList3;
      std::vector <CHR_READS> ChrList4;
      std::vector <CHR_READS> ChrList5;
      double numReads1=getChrReads(Dir1+chrReadsFileName,ChrList1);
      double numReads2=getChrReads(Dir2+chrReadsFileName,ChrList2);
      double numReads3=getChrReads(Dir3+chrReadsFileName,ChrList3);
      double numReads4=getChrReads(Dir4+chrReadsFileName,ChrList4);
      double numReads5=getChrReads(Dir5+chrReadsFileName,ChrList5);
      unsigned int currentChrNum=chrNum;
      if(DEBUG) std::cerr<<std::endl<<"Chr"<<chrNum;
      unsigned int Reads1=0;
      for (std::vector <CHR_READS>::iterator iC1=ChrList1.begin();iC1!=ChrList1.end();++iC1){
	  if ((*iC1).chrNum==currentChrNum)
	  {
	    Reads1=(*iC1).reads;
	  }
      }
      unsigned int Reads2=0;
      for (std::vector <CHR_READS>::iterator iC2=ChrList2.begin();iC2!=ChrList2.end();++iC2){
	  if ((*iC2).chrNum==currentChrNum)
	  {
	    Reads2=(*iC2).reads;
	  }
      }
      unsigned int Reads3=0;
      for (std::vector <CHR_READS>::iterator iC3=ChrList3.begin();iC3!=ChrList3.end();++iC3){
	  if ((*iC3).chrNum==currentChrNum)
	  {
	    Reads3=(*iC3).reads;
	  }
      }
      unsigned int Reads4=0;
      for (std::vector <CHR_READS>::iterator iC4=ChrList4.begin();iC4!=ChrList4.end();++iC4){
	  if ((*iC4).chrNum==currentChrNum)
	  {
	    Reads4=(*iC4).reads;
	  }
      }
      unsigned int Reads5=0;
      for (std::vector <CHR_READS>::iterator iC5=ChrList5.begin();iC5!=ChrList5.end();++iC5){
	  if ((*iC5).chrNum==currentChrNum)
	  {
	    Reads5=(*iC5).reads;
	  }
      }
	//cout <<  "\t" << Reads1 << "\t" << Reads2+Reads1;  WAITUSER
	std::vector<bool> mappabilityMap;
	unsigned mapableBases=0;
/*	pthread_mutex_lock(&mutex_read_map);  */
	for (std::vector<MAPPABILITY_DATA>::iterator iM=mappabilityMapList.begin();iM!=mappabilityMapList.end();++iM){
	  if ((*iM).chrNum==currentChrNum){
	    mappabilityMap = (*iM).mappabilityMap;
	    mapableBases=(*iM).chrMapableBases;
	  }
	}
	std::vector <SEQUENCE_READS> SequenceReads1;
	std::vector <SEQUENCE_READS> SequenceReads2;
	std::vector <SEQUENCE_READS> SequenceReads3;
	std::vector <SEQUENCE_READS> SequenceReads4;
	std::vector <SEQUENCE_READS> SequenceReads5;
	double numChrReads1=getSequenceReads(Dir1,SequenceReads1,currentChrNum);
	double numChrReads2=getSequenceReads(Dir2,SequenceReads2,currentChrNum);
	double numChrReads3=getSequenceReads(Dir3,SequenceReads3,currentChrNum);
	double numChrReads4=getSequenceReads(Dir4,SequenceReads4,currentChrNum);
	double numChrReads5=getSequenceReads(Dir5,SequenceReads5,currentChrNum);
	double coverageChrSequence1=numChrReads1/mapableBases;
	double coverageChrSequence2=numChrReads2/mapableBases;
	double coverageChrSequence3=numChrReads3/mapableBases;
	double coverageChrSequence4=numChrReads4/mapableBases;
	double coverageChrSequence5=numChrReads5/mapableBases;
	 if(DEBUG){
	std::cerr<<std::endl<<numChrReads1<<" Reads in "<<Dir1<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence1 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence1*WINDOWSIZE ;
	std::cerr<<std::endl<<numChrReads2<<" Reads in "<<Dir2<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence2 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence2*WINDOWSIZE ;
	std::cerr<<std::endl<<numChrReads3<<" Reads in "<<Dir3<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence3 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence3*WINDOWSIZE ;
	std::cerr<<std::endl<<numChrReads4<<" Reads in "<<Dir4<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence4 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence4*WINDOWSIZE ;
	std::cerr<<std::endl<<numChrReads4<<" Reads in "<<Dir5<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence5 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence5*WINDOWSIZE ;
    }
// 	WAITUSER
// 	unsigned int limit = (SequenceReads2[SequenceReads2.size()-1].chrPos > \
	SequenceReads1[SequenceReads1.size()-1].chrPos ? SequenceReads2[SequenceReads2.size()-1].chrPos: \
	SequenceReads1[SequenceReads1.size()-1].chrPos);
	std::vector <SEQUENCE_READS>::iterator iControlStart=SequenceReads1.begin();
	std::vector <SEQUENCE_READS>::iterator iExpStart=SequenceReads2.begin();
	std::vector <SEQUENCE_READS>::iterator iCorrectionMapStart=SequenceReads3.begin();
	std::vector <SEQUENCE_READS>::iterator iCNVMapStart=SequenceReads4.begin();
	std::vector <SEQUENCE_READS>::iterator iCNV2MapStart=SequenceReads5.begin();
// 	WAITUSER;
	//cout << std::endl << SequenceReads2[0].chrPos; WAITUSER;
	//cout << std::endl << SequenceReads1[0].chrPos; WAITUSER;

	if (mappabilityMap.size()==0){std::cerr << "Mappability map not found Chr " <<currentChrNum; terminate(10);}
	double coverageSequence1=numReads1/mapableBasesTotal;
	double coverageSequence2=numReads2/mapableBasesTotal;
	double coverageSequence3=numReads3/mapableBasesTotal;
	double coverageSequence4=numReads4/mapableBasesTotal;
	double coverageSequence5=numReads5/mapableBasesTotal;
     if(DEBUG){
	std::cerr<<std::endl<<"Total Reads in data1 = "<< numReads1<< " Coverage = " << coverageSequence1 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageSequence1*WINDOWSIZE ;
	std::cerr<<std::endl<<"Total Reads in data2 = "<< numReads2<< " Coverage = " << coverageSequence2 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageSequence2*WINDOWSIZE ;
	std::cerr<<std::endl<<"Total Reads in data3 = "<< numReads3<< " Coverage = " << coverageSequence3 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageSequence3*WINDOWSIZE ;
	std::cerr<<std::endl<<"Total Reads in data4 = "<< numReads4<< " Coverage = " << coverageSequence4 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageSequence4*WINDOWSIZE ;
	std::cerr<<std::endl<<"Total Reads in data5 = "<< numReads5<< " Coverage = " << coverageSequence5 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageSequence5*WINDOWSIZE ;
    }
	unsigned int limit = mappabilityMap.size()-1;
	pthread_mutex_lock(&mutex_cerr);
	 if(DEBUG) std::cerr << std::endl << "chr" <<convert_chrnum_to_string(currentChrNum) << " lim = " <<  limit << " "; /*WAITUSER*/;
	pthread_mutex_unlock(&mutex_cerr);
	for (unsigned int chrPos=0; ; \
	  chrPos=getMappableWindowEnd(chrPos, WINDOWADVANCE, mappabilityMap))
	{
	    numMappableWindows++;
	    unsigned int windowReads1=0;
	    unsigned int windowReads2=0;
	    unsigned int windowReads3=0;
	    unsigned int windowReads4=0;
        unsigned int windowReads5=0;
	    unsigned int windowStart=chrPos;
	    unsigned int windowEnd=getMappableWindowEnd(windowStart, WINDOWSIZE, mappabilityMap);
	    //cout << std::endl << chrPos << "\t" <<windowEnd << "\t" << getMappableWindowEnd(chrPos, WINDOWADVANCE, mappabilityMap) ;
	    if (windowEnd>=limit) break;//WILL avoid Last Window
	    if ((windowEnd-windowStart) >= WINDOWLIMIT){ /*cout << std::endl << windowEnd-windowStart << "\t" << WINDOWLIMIT;*/ continue;}//avoid large window, esp in

	    {//window1
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iControlStart; iC!=SequenceReads1.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		      if(firstHitFound==0){
			firstHitFound=1;
			iControlStart=iC;
		      }
		      hitsFound=1;
		      windowReads1+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }
	    {
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iExpStart; iC!=SequenceReads2.end(); ++iC)
		{
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))
		    {
      		      if(firstHitFound==0){
			firstHitFound=1;
			iExpStart=iC;
		      }
		      hitsFound=1;
		      windowReads2+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << "\t" << windowStart << "\t"  << windowEnd-1 << "\t" << windowReads2;
		if (hitsFound==0){
		}
	    }
	    {
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iCorrectionMapStart; iC!=SequenceReads3.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		      if(firstHitFound==0){
			firstHitFound=1;
			iCorrectionMapStart=iC;
		      }
		      hitsFound=1;
		      windowReads3+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }
	    {
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iCNVMapStart; iC!=SequenceReads4.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		      if(firstHitFound==0){
			firstHitFound=1;
			iCNVMapStart=iC;
		      }
		      hitsFound=1;
		      windowReads4+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }
	    {
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iCNV2MapStart; iC!=SequenceReads5.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		      if(firstHitFound==0){
			firstHitFound=1;
			iCNV2MapStart=iC;
		      }
		      hitsFound=1;
		      windowReads5+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }
	    double correctionCoefficientSequencing=coverageChrSequence3*WINDOWSIZE/windowReads3;//expectedReads/actualReads
        double correctionCoefficientCNV=(windowReads3*1.0/Reads3) / (windowReads4*1.0/Reads4)  ;//windowReads3*correctionCoefficientSequencing*
        double correctionCoefficientCNV2= (windowReads3*1.0/Reads3) / (windowReads5*1.0/Reads5) ;//windowReads3*correctionCoefficientSequencing*
        if(DEBUG){
            std::cerr<<std::endl<<numChrReads3<<" Reads in "<<Dir3<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence3 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence3*WINDOWSIZE ;
            std::cerr<<std::endl<<numChrReads4<<" Reads in "<<Dir4<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence4 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence4*WINDOWSIZE ;
            std::cerr<<std::endl<<numChrReads5<<" Reads in "<<Dir5<<"chr"<<convert_chrnum_to_string2(chrNum)<<" coverage = "<< coverageChrSequence5 << ": expected reads in "<<WINDOWSIZE << " mapable window = "<< coverageChrSequence5*WINDOWSIZE ;
            std::cerr<<std::endl<<"Windows Reads 3= "<<windowReads3;
            std::cerr<<std::endl<<"Windows Reads 4= "<<windowReads4;
            std::cerr<<std::endl<<"Windows Reads 5= "<<windowReads5;
            std::cerr<<std::endl<<"Correction coefficient Sequencing="<<correctionCoefficientSequencing;
            std::cerr<<std::endl<<"Correction coefficient CNV="<<correctionCoefficientCNV;
            std::cerr<<std::endl<<"Correction coefficient CNV2 ="<<correctionCoefficientCNV2;
            }
	    double pvalCNVIns=1.0;
	    double pvalCNVDel=1.0;
	    if (windowReads4>=1 && windowReads3>=0){
	      unsigned int ps=Reads3+Reads4;/*Pop. Size*/
	      unsigned int sa=windowReads3+windowReads4;/*Sample Size*/
	      unsigned int sp=Reads4;/*Success Population*/
	      unsigned int ss=windowReads4;/*Success Sample*/
// 	      double pval= gsl_cdf_hypergeometric_Q(windowReads1-1,Reads1,Reads2,windowReads2+windowReads1);
	      pvalCNVIns=hgcalc_gte(ps,sa,sp,ss);
// 	      pthread_mutex_lock(&mutex_output);
// 	      cout << std::endl << convert_chrnum_to_string2(currentChrNum) << "\t" << windowReads4 << "\t" << windowReads3+windowReads4 << \
// 	      "\t" << Reads4 << "\t" << Reads3+Reads4<<"\t"<<windowStart<<"\t"<<windowEnd<<"\t"<<pvalCNVIns<<"\t"<<1+windowEnd-windowStart<<"\t";
// 	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,pval});
// 	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }

	    if (windowReads3>=1 && windowReads4>=0){
	      unsigned int ps=Reads3+Reads4;/*Pop. Size*/
	      unsigned int sa=windowReads3+windowReads4;/*Sample Size*/
	      unsigned int sp=Reads3;/*Success Population*/
	      unsigned int ss=windowReads3;/*Success Sample*/
// 	      double pval= gsl_cdf_hypergeometric_Q(windowReads1-1,Reads1,Reads2,windowReads2+windowReads1);
	      pvalCNVDel=hgcalc_gte(ps,sa,sp,ss);
// 	      pthread_mutex_lock(&mutex_output);
// 	      cout << pvalCNVDel;
// 	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,pval});
// 	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }
// 	    WAITUSER
	    if(DEBUG)  std::cerr<<std::endl<<"pvalIns="<<pvalCNVIns<<"\tpvalDel="<<pvalCNVDel;
	    /////////////////////////////////////////////////////////////
    	    double pvalCNV2Ins=1.0;
	    double pvalCNV2Del=1.0;

	    if (windowReads5>=1 && windowReads4>=0){
	      unsigned int ps=Reads4+Reads5;/*Pop. Size*/
	      unsigned int sa=windowReads4+windowReads5;/*Sample Size*/
	      unsigned int sp=Reads5;/*Success Population*/
	      unsigned int ss=windowReads5;/*Success Sample*/
// 	      double pval= gsl_cdf_hypergeometric_Q(windowReads1-1,Reads1,Reads2,windowReads2+windowReads1);
	      pvalCNV2Ins=hgcalc_gte(ps,sa,sp,ss);
// 	      pthread_mutex_lock(&mutex_output);
// 	      cout << std::endl << convert_chrnum_to_string2(currentChrNum) << "\t" << windowReads4 << "\t" << windowReads3+windowReads4 << \
// 	      "\t" << Reads4 << "\t" << Reads3+Reads4<<"\t"<<windowStart<<"\t"<<windowEnd<<"\t"<<pvalCNVIns<<"\t"<<1+windowEnd-windowStart<<"\t";
// 	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,pval});
// 	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }

	    if (windowReads4>=1 && windowReads5>=0){
	      unsigned int ps=Reads4+Reads5;/*Pop. Size*/
	      unsigned int sa=windowReads4+windowReads5;/*Sample Size*/
	      unsigned int sp=Reads4;/*Success Population*/
	      unsigned int ss=windowReads4;/*Success Sample*/
// 	      double pval= gsl_cdf_hypergeometric_Q(windowReads1-1,Reads1,Reads2,windowReads2+windowReads1);
	      pvalCNV2Del=hgcalc_gte(ps,sa,sp,ss);
// 	      pthread_mutex_lock(&mutex_output);
// 	      cout << pvalCNVDel;
// 	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,pval});
// 	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }
    	  if(DEBUG)   std::cerr<<std::endl<<"pvalIns2="<<pvalCNV2Ins<<"\tpvalDel2="<<pvalCNV2Del;
	    /////////////////////////////////////////////////////////////
	    if (windowReads1>=1 && windowReads2>=0){
	      numWindows++;
	      unsigned int ps=Reads1+Reads2;/*Pop. Size*/
	      unsigned int sa=windowReads1+windowReads2;/*Sample Size*/
	      unsigned int sp=Reads1;/*Success Population*/
	      unsigned int ss=windowReads1;/*Success Sample*/
// 	      double pval= gsl_cdf_hypergeometric_Q(windowReads1-1,Reads1,Reads2,windowReads2+windowReads1);
	      double pval=hgcalc_gte(ps,sa,sp,ss);
	      double pvalCorrected=hgcalc_gte(ps,windowReads1*correctionCoefficientCNV2+windowReads2*correctionCoefficientCNV,sp,ss*correctionCoefficientCNV2);
	       if(DEBUG) std::cerr<<std::endl<<pval<<"\t"<<pvalCorrected;
	      pthread_mutex_lock(&mutex_output);
// 	      cout << std::endl << convert_chrnum_to_string2(currentChrNum) << "\t" << windowReads1 << "\t" << windowReads2+windowReads1 << \
// 	      "\t" << Reads1 << "\t" << Reads2+Reads1<<"\t"<<windowStart<<"\t"<<windowEnd<<"\t"<<pval<<"\t"<<1+windowEnd-windowStart<<"\t";
	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,pval,1.0,pvalCNVIns,pvalCNVDel,correctionCoefficientSequencing,correctionCoefficientCNV,pvalCNV2Ins,pvalCNV2Del,correctionCoefficientCNV2,pvalCorrected,1.0});
	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }
	    else{	    //avoid holes in the data list, pval=2.0 for filtering later
	      pthread_mutex_lock(&mutex_output);
	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,2.0,1.0,pvalCNVIns,pvalCNVDel,correctionCoefficientSequencing,correctionCoefficientCNV,pvalCNV2Ins,pvalCNV2Del,correctionCoefficientCNV2,2.0,1.0});
	      pthread_mutex_unlock(&mutex_output);
	    }
	    //cout << std::endl << chrPos1 << "\t" << chrPos2 << "\t" << "chrPos " << chrPos ;
	    //cout << std::endl << windowReads1;
	    //std::vector <SEQUENCE_READS>::iterator iExpSequence=SequenceReads2.begin();

	}
	pthread_mutex_lock(&mutex_cerr);
	 if(DEBUG) std::cerr << std::endl << convert_chrnum_to_string(currentChrNum) << " Done"; //WAITUSER;
	pthread_mutex_unlock(&mutex_cerr);
	return numMappableWindows;
  }

  int calcHypergeometricPvalWithout(unsigned int chrNum,unsigned int chrReads,std::string chrReadsFileName,std::string Dir1,std::string Dir2,std::vector<MAPPABILITY_DATA>&mappabilityMapList,std::vector<WINDOW_DATA> &outData)
{
      const unsigned int WINDOWLIMIT=WINDOWSIZE*10;
      unsigned int numMappableWindows=0;
      unsigned int numWindows=0;
/*      std::cerr<<std::endl<<WINDOWSIZE<<"\t"<<WINDOWADVANCE;WAITUSER*/
      std::vector <CHR_READS> ChrList1, ChrList2;
      int chr_data1 = getChrReads(Dir1+chrReadsFileName,ChrList1);
      int chr_data2 = getChrReads(Dir2+chrReadsFileName,ChrList2);
      if(DEBUG) std::cerr<<"Total Reads in data1 = "<< chr_data1<<std::endl;
      if(DEBUG) std::cerr<<"Total Reads in data2 = "<< chr_data2<<std::endl;

      unsigned int currentChrNum=chrNum;
      if(DEBUG) std::cerr<<std::endl<<"Chr"<<chrNum;
      unsigned int Reads1=0;
      for (std::vector <CHR_READS>::iterator iC1=ChrList1.begin();iC1!=ChrList1.end();++iC1){
	  if ((*iC1).chrNum==currentChrNum)
	  {
	    Reads1=(*iC1).reads;
	  }
      }
      unsigned int Reads2=0;
      for (std::vector <CHR_READS>::iterator iC2=ChrList2.begin();iC2!=ChrList2.end();++iC2){
	  if ((*iC2).chrNum==currentChrNum)
	  {
	    Reads2=(*iC2).reads;
	  }
      }
	//cout <<  "\t" << Reads1 << "\t" << Reads2+Reads1;  WAITUSER
	std::vector <SEQUENCE_READS> SequenceReads1, SequenceReads2;
	getSequenceReads(Dir1,SequenceReads1,currentChrNum);
	getSequenceReads(Dir2,SequenceReads2,currentChrNum);
// 	unsigned int limit = (SequenceReads2[SequenceReads2.size()-1].chrPos > \
	SequenceReads1[SequenceReads1.size()-1].chrPos ? SequenceReads2[SequenceReads2.size()-1].chrPos: \
	SequenceReads1[SequenceReads1.size()-1].chrPos);
	std::vector <SEQUENCE_READS>::iterator iControlStart=SequenceReads1.begin();
	std::vector <SEQUENCE_READS>::iterator iExpStart=SequenceReads2.begin();
	//cout << std::endl << SequenceReads2[0].chrPos; WAITUSER;
	//cout << std::endl << SequenceReads1[0].chrPos; WAITUSER;
	std::vector<bool> mappabilityMap;
/*	pthread_mutex_lock(&mutex_read_map);  */
	for (std::vector<MAPPABILITY_DATA>::iterator iM=mappabilityMapList.begin();iM!=mappabilityMapList.end();++iM){
	  if ((*iM).chrNum==currentChrNum){
	    mappabilityMap = (*iM).mappabilityMap;
	  }
	}
	if (mappabilityMap.size()==0){std::cerr << "Mappability map not found Chr " <<currentChrNum; terminate(1);}
	unsigned int limit = mappabilityMap.size()-1;
	pthread_mutex_lock(&mutex_cerr);
	 if(DEBUG) std::cerr << std::endl << convert_chrnum_to_string(currentChrNum) << " lim = " <<  limit << " "; //WAITUSER;
	pthread_mutex_unlock(&mutex_cerr);
	for (unsigned int chrPos=0; ; \
	  chrPos=getMappableWindowEnd(chrPos, WINDOWADVANCE, mappabilityMap))
	{
	    numMappableWindows++;
	    unsigned int windowReads1=0;
	    unsigned int windowReads2=0;
	    unsigned int windowStart=chrPos;
	    unsigned int windowEnd=getMappableWindowEnd(windowStart, WINDOWSIZE, mappabilityMap);
	    //cout << std::endl << chrPos << "\t" <<windowEnd << "\t" << getMappableWindowEnd(chrPos, WINDOWADVANCE, mappabilityMap) ;
	    if (windowEnd>=limit-WINDOWSIZE) break;//WILL avoid Last Window
	    if ((windowEnd-windowStart) >= WINDOWLIMIT){ /*cout << std::endl << windowEnd-windowStart << "\t" << WINDOWLIMIT;*/ continue;}//avoid large window, esp in
	    //centromeric regions
	    unsigned int chrPos1=0;
	    unsigned int chrPos2=0;
	    //cout << std::endl << chrPos << "\t" <<getMappableWindowEnd(chrPos,WINDOWADVANCE ,mappabilityMap) << "\t" <<getMappableWindowEnd(chrPos,1000,mappabilityMap);WAITUSER;
	    {//window1
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iControlStart; iC!=SequenceReads1.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		      if(firstHitFound==0){
			firstHitFound=1;
			iControlStart=iC;
		      }
		      hitsFound=1;
		      windowReads1+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }
	    {
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iExpStart; iC!=SequenceReads2.end(); ++iC)
		{
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))
		    {
      		      if(firstHitFound==0){
			firstHitFound=1;
			iExpStart=iC;
		      }
		      hitsFound=1;
		      windowReads2+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << "\t" << windowStart << "\t"  << windowEnd-1 << "\t" << windowReads2;
		if (hitsFound==0){
		}
	    }
	    if (windowReads1>=1 && windowReads2>=0){
	      numWindows++;
	      unsigned int ps=Reads1+Reads2;/*Pop. Size*/
	      unsigned int sa=windowReads1+windowReads2;/*Sample Size*/
	      unsigned int sp=Reads1;/*Success Population*/
	      unsigned int ss=windowReads1;/*Success Sample*/
// 	      double pval= gsl_cdf_hypergeometric_Q(windowReads1-1,Reads1,Reads2,windowReads2+windowReads1);
	      double pval=hgcalc_gte(ps,sa,sp,ss);
	      pthread_mutex_lock(&mutex_output);
// 	      cout << std::endl << convert_chrnum_to_string2(currentChrNum) << "\t" << windowReads1 << "\t" << windowReads2+windowReads1 << \
// 	      "\t" << Reads1 << "\t" << Reads2+Reads1<<"\t"<<windowStart<<"\t"<<windowEnd<<"\t"<<pval<<"\t"<<1+windowEnd-windowStart<<"\t";
	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,pval});
	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }
	    else{	    //avoid holes in the data list, pval=2.0 for filtering later
	      pthread_mutex_lock(&mutex_output);
	      outData.push_back({currentChrNum,windowReads1,windowReads2,windowStart,windowEnd,2.0,1.0});
	      pthread_mutex_unlock(&mutex_output);
	    }
	    //cout << std::endl << chrPos1 << "\t" << chrPos2 << "\t" << "chrPos " << chrPos ;
	    //cout << std::endl << windowReads1;
	    //std::vector <SEQUENCE_READS>::iterator iExpSequence=SequenceReads2.begin();

	}
	pthread_mutex_lock(&mutex_cerr);
if(DEBUG) std::cerr << std::endl << convert_chrnum_to_string(currentChrNum) << " Done"; //WAITUSER;
	pthread_mutex_unlock(&mutex_cerr);
	return numMappableWindows;
  }
int calcHypergeometricPvalOneSample(unsigned int chrNum,unsigned int chrReads,unsigned int mappableBases,std::string chrReadsFileName,std::string Dir1,std::vector<MAPPABILITY_DATA>&mappabilityMapList,std::vector<WINDOW_DATA> &outData)
{
      const unsigned int WINDOWLIMIT=WINDOWSIZE*10;
/*      std::cerr<<std::endl<<WINDOWSIZE<<"\t"<<WINDOWADVANCE;WAITUSER*/
      std::vector <CHR_READS> ChrList1;
      unsigned totalReads=getChrReads(Dir1+chrReadsFileName,ChrList1);
    if(DEBUG)  std::cerr<<std::endl<<"Total Reads in data1 = "<<totalReads;
      unsigned int currentChrNum=chrNum;
    if(DEBUG)  std::cerr<<std::endl<<"Chr"<<chrNum;
      unsigned int Reads1=0;
      for (std::vector <CHR_READS>::iterator iC1=ChrList1.begin();iC1!=ChrList1.end();++iC1){
	  if ((*iC1).chrNum==currentChrNum)
	  {
	    Reads1=(*iC1).reads;
	  }
      }
	std::vector <SEQUENCE_READS> SequenceReads1;
	getSequenceReads(Dir1,SequenceReads1,currentChrNum);
/*	unsigned int limit = SequenceReads1[SequenceReads1.size()-1].chrPos;*/
	std::vector <SEQUENCE_READS>::iterator iControlStart=SequenceReads1.begin();
	std::vector<bool> mappabilityMap;
	for (std::vector<MAPPABILITY_DATA>::iterator iM=mappabilityMapList.begin();iM!=mappabilityMapList.end();++iM){
	  if ((*iM).chrNum==currentChrNum){
	    mappabilityMap = (*iM).mappabilityMap;
	  }
	}
	if (mappabilityMap.size()==0){std::cerr << "Mappability map not found Chr " <<currentChrNum; terminate(1); }
	unsigned int limit = mappabilityMap.size()-1;
	pthread_mutex_lock(&mutex_cerr);
	if(DEBUG)std::cerr << std::endl << convert_chrnum_to_string(currentChrNum) << " lim = " <<  limit << " "; //WAITUSER;
	pthread_mutex_unlock(&mutex_cerr);
	for (unsigned int chrPos=0; ; \
	  chrPos=getMappableWindowEnd(chrPos, WINDOWADVANCE, mappabilityMap))
	{
	    unsigned int windowReads1=0;
	    unsigned int windowStart=chrPos;
	    unsigned int windowEnd=getMappableWindowEnd(windowStart, WINDOWSIZE, mappabilityMap);
	    if (windowEnd>=limit-WINDOWSIZE) break;//WILL avoid Last Window
	   if ((windowEnd-windowStart) >= WINDOWLIMIT){ /*cout << std::endl << windowEnd-windowStart << "\t" << WINDOWLIMIT;*/ continue;}//avoid large window, esp in
	    //centromeric regions
	    unsigned int chrPos1=0;
	    {//window Control reads
		bool hitsFound=0;
		bool firstHitFound=0;
		for (std::vector <SEQUENCE_READS>::iterator iC=iControlStart; iC!=SequenceReads1.end(); ++iC)
		{//find all reads in the current window
		    if (((*iC).chrPos>=windowStart) && ((*iC).chrPos<=windowEnd))//
		    {
		  /*    if(firstHitFound==0){
			firstHitFound=1;
			iControlStart=iC;
		      }
		      hitsFound=1;    */
		      windowReads1+=(*iC).hits;
		    }
		    if ((*iC).chrPos > windowEnd) break;//
		}
		//cout << std::endl << windowStart << "\t"  << windowEnd-1 << "\t"<<windowReads1;
		if (hitsFound==0){
		}
	    }
	    if (windowReads1>=1){
	      unsigned int ps=mappableBases;/*Pop. Size*/
	      unsigned int sa=WINDOWSIZE;/*Sample Size*/
	      unsigned int sp=totalReads;/*Success Population*/
	      unsigned int ss=windowReads1;/*Success Sample*/
	      double pval=hgcalc_gte(ps,sa,sp,ss);
	      pthread_mutex_lock(&mutex_output);
	      outData.push_back({currentChrNum,windowReads1,windowReads1,windowStart,windowEnd,pval});
	      pthread_mutex_unlock(&mutex_output);
	      //WAITUSER
	    }
	    //avoid holes in the data list, pval=2.0 for filtering later
	    else{
      	      pthread_mutex_lock(&mutex_output);
      	      outData.push_back({currentChrNum,windowReads1,windowReads1,windowStart,windowEnd,2.0,1.0});
	      pthread_mutex_unlock(&mutex_output);
	    }
	}
	pthread_mutex_lock(&mutex_cerr);
	if(DEBUG) std::cerr << std::endl << convert_chrnum_to_string(currentChrNum) << " Done"; //WAITUSER;
	pthread_mutex_unlock(&mutex_cerr);
	return 0;
}
