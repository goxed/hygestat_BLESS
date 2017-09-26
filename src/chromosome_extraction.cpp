
#include "../include/chromosome_extraction.h"

//************  Chr extraction ********************************//
unsigned int extract_chrnum(std::string s_chrNum, bool XY){
//--------------------------------------------------------------
  unsigned int len=s_chrNum.size();
  unsigned int count =0;
  unsigned int chrNum=0;
  for (unsigned int n=len; n>0; --n) {
     char c = s_chrNum[n-1];
     if (c >= '0' && c<='9'){
       ++count;
       if     (count==1) {
	 chrNum+= (c-'0');
       }else if (count==2) {
	 chrNum+= (c-'0')*10;
       }else if (count==3) {
	 chrNum+= (c-'0')*100;
       }else if (count > 3)
	 chrNum=0;
     }
  }

  if (chrNum==0 && XY==0) {
    for (unsigned int n=len; n>0; --n) {
      char c = s_chrNum[n-1];
      if (c == 'X' || c== 'V' || c=='I'){
	++count;
	if (count==1) {
	  if      (c=='X') chrNum+=10;
	  else if (c=='V') chrNum+=5;
	  else if (c=='I') chrNum+=1;
	}else if (count>1 && count <=4) {
	  if (c < char(s_chrNum[n])){
	    if (c=='X')      chrNum-=10;
	    else if (c=='V') chrNum-=5;
	    else if (c=='I') chrNum-=1;
	  }
	  else if (c >= char(s_chrNum[n])){
	    if (c=='X')      chrNum+=10;
	    else if (c=='V') chrNum+=5;
	    else if (c=='I') chrNum+=1;
	  }
	}else if (count >4) chrNum=0;
      }
    }
  } else if (chrNum==0 && XY==1)//XY chromosome
    for (unsigned int n=len; n>0; --n){
      char c=s_chrNum[n-1];
      if (c=='X')      chrNum=X_CHROMOSOME;
      else if (c=='Y') chrNum=Y_CHROMOSOME;
    }
  return chrNum;
}

//***********************************************************
std::string convert_chrnum_to_string2(unsigned int chrNum){
//--------------------------------------------------------

  std::string s_chrNum="";
  std::stringstream ss_chrNum;
  if (chrNum <10){
    std::string s_chrNumTemp="";
    ss_chrNum << chrNum;
    ss_chrNum >> s_chrNumTemp;
    s_chrNum = "0" + s_chrNumTemp;
  }
  else if (chrNum <X_CHROMOSOME){
    ss_chrNum << chrNum;
    ss_chrNum >> s_chrNum;
  }
  else {
    char c_chrNum=chrNum-X_CHROMOSOME+'X';
    s_chrNum=c_chrNum;
  }

  return s_chrNum;
}
//*******************************************************
std::string convert_chrnum_to_string(unsigned int chrNum){
//------------------------------------------------------
  std::string s_chrNum="";
  std::stringstream ss_chrNum;
  if (chrNum <1000){
    ss_chrNum << chrNum;
    ss_chrNum >> s_chrNum;
  }
  else {
    char c_chrNum=chrNum-1000+'X';
    s_chrNum=c_chrNum;
  }

  return s_chrNum;
}


//***********************************************************************************
int getChrReads(std::string FileName, std::vector<CHR_READS> &chrList ){
//-----------------------------------------------------------------------------------
  std::ifstream controlChrFile(FileName.c_str());
  if(!controlChrFile.is_open()){
        std::cerr<<"Can't find chromosome file,"<< FileName<<" exiting the program...\n";
         exit(0);
  }
  unsigned int totalReads=0;
  while(!controlChrFile.eof()){
    std::string line;
    getline (controlChrFile,line);
    std::string::size_type tab_location;
    tab_location=line.find_first_of('\t');
    CHR_READS crtemp;
    unsigned int chrNum = extract_chrnum(line.substr(0,tab_location),true);
    if (chrNum!=0){
      crtemp.chrNum=chrNum;
      std::stringstream ss_reads(line.substr(tab_location+1));
      unsigned int reads;
      ss_reads >> reads;
      crtemp.reads=reads;
      totalReads+=reads;
      chrList.push_back(crtemp);
    }
  }
  return totalReads;
}


//******************************************************
int getPhredOffset(std::string fastqFileName) {
  //------------------------------------------------------
  int phredOffset=64;
  std::ifstream fastqFileTest(fastqFileName.c_str());
  unsigned int lineMax=4*1000;
  unsigned int lineCount=0;
  while (!fastqFileTest.eof()&&lineCount<=lineMax){
    static std::string fasta="";
    static std::string quality="";
    static std::string f="",q="";
    std::string line="";
    fastqFileTest>>line;
    char testChar=line[0];
    static fastqFileSequenceType fastqFileSequence;
    if (testChar=='@') {
      if(fastqFileSequence==fastaHeader)
	fastqFileSequence=fastaLine;
    }
    else if (testChar=='+'){
      if(fastqFileSequence==qualityHeader)
	fastqFileSequence=qualityLine;
    }
    else {
      if(fastqFileSequence==fastaLine){
	fastqFileSequence=qualityHeader;
      }
      else if(fastqFileSequence==qualityLine){
	fastqFileSequence=fastaHeader;
	q=line;
	std::string phred="";
	std::string::size_type cutoff=line.size();
	for (std::string::size_type l=0;l<q.size();++l){
	  char c=q[l];
	  int p64=c-64;
	  int p33=c-33;
	  if (p64<-5){
	    phredOffset=33;
	  }
	  if (p33<0){
	  }
	}
      }
    }
    lineCount++;
  }
  return phredOffset;
}
