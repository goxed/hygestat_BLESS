
#include "../include/process_fastq.h"
//**********************Running process fastq for trimming*****************************************//
void process_fast(std::string fastqFileName){
//--------------------------------------------------------------------------------------------------
  boost::iostreams::stream_buffer<boost::iostreams::file_source> fastqFileSource(fastqFileName.c_str());
  std::istream fastqFile (&fastqFileSource);
  int phredOffset=getPhredOffset(fastqFileName);
  if(DEBUG) std::cerr<< "Phred Offset:"<<phredOffset;
  int phredCutoff=PHRED_CUTOFF;
  if(DEBUG) std::cerr <<std::endl<< "Phred cutoff for trimming reads:" << phredCutoff;
  unsigned int lengthCutoff=LENGTH_CUTOFF;
  if(DEBUG) std::cerr <<std::endl<< "Trimmed reads will be >=" << lengthCutoff<< " contiguous bases";
  std::string outputFastaFileName;
  if(phredCutoff>0){
    outputFastaFileName=fastqFileName+".qct.fa";
  }
  else if (phredCutoff==0){
    outputFastaFileName=fastqFileName+".dump.fa";
  }
  boost::iostreams::stream_buffer<boost::iostreams::file_sink> fastaFileSink(outputFastaFileName.c_str());
  std::ostream outputFastaFile (&fastaFileSink);

  if (phredOffset==0){
    if(DEBUG) std::cerr<<std::endl<<"Sorry Phred Offset cannot be figured out";
    terminate(2);
  }
  unsigned int lineCnt=0;
  long long int qualVals=0;
  long long int numbp=0;
  long long int qualValsCutoff=0;
  unsigned int numGoodReads=0;
  std::vector<std::vector<unsigned int> >qV;
  std::stringstream outFileStream;
  std::stringstream outFastaFileStream;
  const std::string Hiseq2500CENTSignature2="P6KDKXP1";
  const std::string Hiseq2500CENTSignature="M01523";
  const std::string Hiseq4000CENTSignature="K00248";
  const std::string genomeAnalyzerHiseqSignature="HWI";
  const std::string solexaSignature="slxa";
  const std::string gaiiSignature="USI";
  const std::string hiseqStringUT="UT";
  const std::string hiseqStringCENT="HISEQ";
  const std::string hiseqStringCENT2="SN870";
  while (!fastqFile.eof()){
    static std::string fasta="";
    static std::string quality="";
    std::string line="";
    getline(fastqFile,line);
    static std::string f="",q="";
    ++lineCnt;
    char testChar=line[0];
    static fastqFileSequenceType fastqFileSequence;
    if ((testChar=='@')||(line.find(genomeAnalyzerHiseqSignature)!=std::string::npos||line.find(Hiseq4000CENTSignature)!=std::string::npos||line.find(gaiiSignature)!=std::string::npos|| \
    line.find(solexaSignature)!=std::string::npos||line.find(hiseqStringUT)!=std::string::npos||line.find(hiseqStringCENT)!=std::string::npos||line.find(hiseqStringCENT2)!=std::string::npos|| \
    line.find(Hiseq2500CENTSignature)!=std::string::npos||line.find(Hiseq2500CENTSignature2)!=std::string::npos)){
#ifdef DEBUG
      std::cerr << std::endl <<"Header\t"<< "line [0]:"<<testChar;
#endif
      if(fastqFileSequence==fastaHeader)
	fastqFileSequence=fastaLine;
    }
    else if ((testChar=='+')||((line.size()==1)||line.find(gaiiSignature)!=std::string::npos||line.find(genomeAnalyzerHiseqSignature)!=std::string::npos || \
    line.find(Hiseq4000CENTSignature)!=std::string::npos||line.find(solexaSignature)!=std::string::npos||line.find(hiseqStringUT)!=std::string::npos||line.find(hiseqStringCENT)!=std::string::npos|| \
    line.find(hiseqStringCENT2)!=std::string::npos||line.find(Hiseq2500CENTSignature)!=std::string::npos||line.find(Hiseq2500CENTSignature2)!=std::string::npos)){
#ifdef DEBUG
      std::cerr << std::endl << "Quality Header\t"<<"line [0]:"<<testChar;
#endif
      if(fastqFileSequence==qualityHeader)
	fastqFileSequence=qualityLine;
    }
    else {
      if(fastqFileSequence==fastaLine){
	fastqFileSequence=qualityHeader;
#ifdef DEBUG
	std::cerr<<std::endl<<"Data";
	outFileStream << std::endl << "F: " << line;
#endif
	f=line;
	numbp+=line.size();
      }
      else if(fastqFileSequence==qualityLine){
	fastqFileSequence=fastaHeader;
#ifdef DEBUG
	std::cerr<<std::endl<<"Quality Data";
	outFileStream << std::endl << "Q: "<< line;
#endif
	q=line;
	std::string phred="";
	std::string::size_type cutoff=line.size();
	bool cutoffdecided=0;
	if(qV.size()<line.size())qV.resize(line.size());
	for (std::string::size_type l=0;l<q.size();++l){
	  char c=q[l];
	  int p=c-phredOffset;
  	  qualVals+=p;
	  qV[l].push_back(p);
	  if(p>PHRED_QUALITY)qualValsCutoff++;
	  if (p<-5){
	    std::cerr << "Quality score not compatible with this program. Terminating!";
	    terminate(3);
	  }
	  else if (p<PHRED_CUTOFF&&cutoffdecided==0){
	    cutoff=l;
	    cutoffdecided=1;
	  }
	}
#ifdef DEBUG
	outFileStream<<"\tReadLen cut=" <<cutoff;
#endif
	if (cutoff>=LENGTH_CUTOFF) {
	  numGoodReads++;
#ifdef DEBUG
	  outFileStream << "\tgood";
#endif
	  outFastaFileStream<<f.substr(0,cutoff)<<std::endl;
	  outputFastaFile<<outFastaFileStream.str();
	  outFastaFileStream.str("");
	}
      }
      else {
	if (line.size()>1)
	//  std::cerr << std::endl << "Possible Error in fastq file line or sequencer name is not in database. Please modify lines 61, 68, add sequencer name and re-compile:" <<lineCnt << "\t"<<(int)line[0]<<"\t"<<line; /*WAITUSER*/
	fastqFileSequence=fastaHeader;
      }
    }
  }
#ifdef DEBUG
  std::string outputFileName=fastqFileName+".qc";
  boost::iostreams::stream_buffer<boost::iostreams::file_sink> fastqFileSink(outputFileName.c_str());
  std::ostream outputFastqFile (&fastqFileSink);
#endif
#ifdef DEBUG
  outputFastqFile<<outFileStream.str();
#endif
  unsigned int numReads=lineCnt/4;
  float averagePhred=qualVals*1.0/numbp;
  unsigned int readLen=numbp/numReads;
  unsigned int phredQuality=PHRED_QUALITY;
   if(DEBUG) std::cerr<<std::endl<<"Total yield:"<<group_integers(numbp) <<" (bases) and "<<group_integers(lineCnt/4) <<" (reads),"<<" Read Length: "<<readLen<<std::endl \
      <<"Bases with quality >"<<phredQuality<<": "<<group_integers(qualValsCutoff)<<",\t"<<(qualValsCutoff*100.0/numbp)<<"%";
  std::string  IlluminaVerdict="No";
  if ((qualValsCutoff*100.0/numbp)>=85.0){
    if(DEBUG) std::cerr << "\t(Consistent with Illumina's spec of >85% Phred 30 bases)";
    IlluminaVerdict="Yes";
  }
  else{
     if(DEBUG) std::cerr << "\t(! Not consistent with Illumina's spec of >85% Phred 30 bases) !";
    IlluminaVerdict="No";
  }
   if(DEBUG)std::cerr << std::endl <<"Good Quality reads "<<group_integers(numGoodReads);
  //  <<averagePhred <<":average Phred,\t"<<qualValsCutoff<<":quality bp (> Phred " << phredQuality<< "),\t"<<(qualValsCutoff*100.0/numbp)<<":% quality bp,\t" <<std::endl;
  float percentageGoodReads=100.0*numGoodReads/numReads;
  std::string verdict="NA";
  if (percentageGoodReads>=80.0)verdict="Good Quality";
  else if (percentageGoodReads>=50.0) verdict="Fair Quality";
  else if (percentageGoodReads<50.0) verdict="Poor Quality";
   if(DEBUG)std::cerr << std::endl << verdict;
  std::string resultFileName=fastqFileName+".csv";
  fastqFileName=fastqFileName.substr(fastqFileName.find_last_of("/")+1);
  std::ofstream logFile(resultFileName.c_str());
  logFile<<"Sample name,Total reads,Reads with Phred >=20; Len. >=34,%Good,Sample Verdict,% Bases with quality > "<<phredQuality<<", Consistent with Illumina's Spec. >85%";
  logFile<<std::endl;
  logFile<<"\""<<fastqFileName<<"\",";
  logFile<<"\""<<group_integers(lineCnt/4)<<"\",";
  logFile<<"\""<<group_integers(numGoodReads)<<"\",";
  logFile<<percentageGoodReads<<"%,";
  logFile<<verdict<<",";
  logFile<<(qualValsCutoff*100.0/numbp)<<"%,";
  logFile<<IlluminaVerdict<<",";
  logFile<<std::endl;
  logFile<<"Median Quality Scores";
  //omp_set_num_threads(12);
  for(std::vector<std::vector<unsigned int> >::iterator q=qV.begin();q!=qV.end();++q){
    //__gnu_parallel::
    sort((*q).begin(),(*q).end());
    std::vector<unsigned int>::iterator qb=(*q).begin()+(*q).size()/2;
    logFile<<","<<(*qb);
  }
  logFile<<std::endl;
  logFile<<"Num Scores for each Base Position";
  for(std::vector<std::vector<unsigned int> >::iterator q=qV.begin();q!=qV.end();++q){
    logFile<<","<<(*q).size();
  }
  logFile<<std::endl;
  /* for (unsigned n=0;n<10;++n)*/{
    logFile<<",";
  }
  logFile<<"Results produced by Instant-Seq quality module; iseq_qm v0.4b; www.iseq.utmb.edu" ;

  logFile.close();
  /*WAITUSER*/;
  //return 0;
}

