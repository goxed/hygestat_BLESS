
#include "../include/compute_cytoband.h"

//******************function cytoband for human and mouse data ******************************//
void cytoband_function(std::string read_directory)  {
//*******************************************************************************************//
  std::string sequenceDataDir1="./";
  std::string sequenceDataDir2="./";
  std::string sequenceDataDir3="./";
//  sequenceDataDir1=read_directory+"_"+to_string(WINDOWSIZE)+"_close_barcode";
//  sequenceDataDir2=read_directory+"_"+to_string(WINDOWSIZE)+"_distant_barcode";
//  sequenceDataDir3=read_directory+"_"+to_string(WINDOWSIZE)+"_no_barcode";

  sequenceDataDir1=read_directory+"_close_barcode";
  sequenceDataDir2=read_directory+"_distant_barcode";
  sequenceDataDir3=read_directory+"_no_barcode";
  std::ifstream cbandsFile;
  std::ifstream fragileCytoBandsFile;
  std::vector <CYTOBAND> cytoBands;
  std::vector <CYTOBAND> fragileBands;

  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE"){
    //std::ifstream cbandsFile;
    cbandsFile.open(cbandsFileName.c_str());
    //std::vector <CYTOBAND> cytoBands;

    while (!cbandsFile.eof()){
      std::string line;
      getline (cbandsFile, line);
      unsigned int columns=0;
      size_t tab;
      std::string s_chrNum, s_chrStart, s_chrEnd, s_name;
      do{
	tab=line.find_first_of('\t');
	columns++;
	switch(columns){
	case 1: s_chrNum=line.substr(0,tab);
	  break;
	case 2: s_chrStart=line.substr(0,tab);
	  break;
	case 3: s_chrEnd=line.substr(0,tab);
	  break;
	case 4: s_name=line.substr(0,tab);
	  break;
	default:
	  break;
	}
        line = line.substr(tab+1, line.size());
      }while(tab!=ULONG_MAX);
      //std::cout << std::endl << extract_chrnum(s_chrNum) << "\t" << s_chrStart << "\t" << \
      s_chrEnd << "\t" << s_name; WAITUSER;
      if (extract_chrnum(s_chrNum, true)>0){
	CYTOBAND cytoBandTmp;
	std::stringstream ss_chrStart(s_chrStart);
	std::stringstream ss_chrEnd(s_chrEnd);
	unsigned int chrStart;
	unsigned int chrEnd;
	ss_chrStart >> chrStart;
	ss_chrEnd >> chrEnd;
	cytoBandTmp.chrNum=extract_chrnum(s_chrNum, true);
	cytoBandTmp.chrStart=chrStart;
	cytoBandTmp.chrEnd=chrEnd;
	cytoBandTmp.name=s_name;
	cytoBands.push_back(cytoBandTmp);
      }
    }
    if(DEBUG) std::cout << std::endl << cytoBands.size() << " Cytobands parsed"; //WAITUSER;
    ///////Processing Fragile Cytobands///////////////////
    //std::ifstream fragileCytoBandsFile;

    fragileCytoBandsFile.open(fbandsFileName.c_str());
    //std::vector <CYTOBAND> fragileBands;
    while(!fragileCytoBandsFile.eof()){
      std::string line;
      getline(fragileCytoBandsFile,line);
      unsigned int pqPos=UINT_MAX;
      for(size_t n=0;n<=line.size();++n)
	if (line[n]=='p'||line[n]=='q'){
	  pqPos=n;
	  break;
	}
      if (pqPos!=UINT_MAX){
	//std::cout << std::endl << pqPos << line[pqPos] ; WAITUSER;
	std::string s_chrNum=line.substr(0,pqPos);
	unsigned int chrNum=extract_chrnum(s_chrNum, true);
	std::string name=line.substr(pqPos) ;
	//std::cout << std::endl << chrNum << " " << name ;//WAITUSER
	std::vector <CYTOBAND>::iterator i;
	unsigned int numbands=0;
	bool foundband=false;
	std::string::size_type fcbDot=name.find_first_of('.');
	for (i=cytoBands.begin();i<cytoBands.end();++i){
	  if ((*i).chrNum==chrNum){
	    if(fcbDot!=std::string::npos){//there's a dot in fcb
	      if ((*i).name.substr(0,name.size())==name){
		//std::cout <<  "\tSubstring match " << (*i).name.substr(0,name.size()) << " " << name <<;WAITUSER;
		fragileBands.push_back(*i);
		++numbands;
		foundband=true;
		//break;
	      }
	    }
	    else if (fcbDot==std::string::npos){//did not find a dot in fcb
	      std::string cbName=(*i).name;
	      std::string::size_type cbDot=cbName.find_first_of('.');
	      if (cbDot==std::string::npos){//did not find a dot in cb
		if (cbName.substr(0,cbDot)==name){
		  fragileBands.push_back(*i);
		  ++numbands;
		  foundband=true;
		  break;
		}
	      }else {
		if (cbName.substr(0,cbDot)==name){
		  fragileBands.push_back(*i);
		  ++numbands;
		  foundband=true;
		}
	      }
	    }
	  }
	}
        if (foundband==false){ std::cout << std::endl  << name << " not found in chr " << chrNum ; WAITUSER; }
        else {//std::cout << std::endl << name << " found " << numbands << " times";
        }
      }
    }
    if (fragileCytoBandsFile.is_open()) fragileCytoBandsFile.close();
    if (cbandsFile.is_open()) cbandsFile.close();
    {
      std::vector <CYTOBAND>::iterator i;
      std::vector <std::string> fragileCytoNames;
      for (i=fragileBands.begin();i<fragileBands.end();++i){ //dump fragile cyto names to a vector
	//if(chrPos>=(*i).chrStart && chrPos<=(*i).chrEnd)
	//std::cout << std::endl << (*i).chrNum << "\t" << (*i).name << "\t" << (*i).chrStart << "\t" << (*i).chrEnd;// WAITUSER;
	std::stringstream s;
	s<<convert_chrnum_to_string((*i).chrNum)<<(*i).name;
	fragileCytoNames.push_back(s.str());
      }
      sort(fragileCytoNames.begin(),fragileCytoNames.end());
      std::ofstream fcbFile ("fragilecytobands.txt");
      for (std::vector<std::string>::iterator i=fragileCytoNames.begin();i<fragileCytoNames.end();++i){
	fcbFile << (*i) << std::endl;
      }
      if (fcbFile.is_open())fcbFile.close();
    }
  }

  /////////Processing bowtie hits files//////////////////////////
  unsigned int totalNumFragileSites=0;
  unsigned int totalNumFragileSites1=0;
  unsigned int totalReads=0;
  unsigned int totalReads1=0;
  std::vector <int> chrNumList;
  std::vector <int> cytobandReads(cytoBands.size(),0);
  if(org_name=="human" || org_name=="HUMAN"){
    for (unsigned int chrNum=1; chrNum<=22; ++chrNum)
      chrNumList.push_back(chrNum);
    chrNumList.push_back(X_CHROMOSOME); chrNumList.push_back(Y_CHROMOSOME);
    //std::vector <int> cytobandReads(cytoBands.size(),0);
  }
  else if(org_name=="mouse" || org_name=="MOUSE"){
    for (unsigned int chrNum=1; chrNum<=19; ++chrNum)
      chrNumList.push_back(chrNum);
    chrNumList.push_back(X_CHROMOSOME); chrNumList.push_back(Y_CHROMOSOME);
  }
  else if(org_name=="yeast" || org_name=="YEAST"){
    for (unsigned int chrNum=1; chrNum<=16; ++chrNum)
      chrNumList.push_back(chrNum);
    //chrNumList.push_back(M_CHROMOSOME);
  }
  std::string chrReadFileName=sequenceDataDir1+"/chrReads.txt";
  std::ofstream chrReadsFile(chrReadFileName.c_str());
  for (std::vector <int>::iterator i=chrNumList.begin(); i<chrNumList.end();++i){
    unsigned int chrNum=*i;
    std::stringstream ss_chrNum;
    ss_chrNum << chrNum;
    std::string s_chrNum;
    ss_chrNum >> s_chrNum;
    if (chrNum>=X_CHROMOSOME) s_chrNum=chrNum-X_CHROMOSOME+'X';
    unsigned int chrReadsTwin=0;
    for (unsigned int s=0; s<2; ++s){
      char strand='+';
      if (s==0) strand='+';
      else if(s==1) strand='-';
       if(DEBUG) std::cout << std::endl << "chr"+s_chrNum+strand;//WAITUSER;
      std::string readsFileName=sequenceDataDir1+"/chr"+s_chrNum+strand+".hits.txt";
      if (checkFile(readsFileName)==-1) {
	   if(DEBUG) std::cout << std::endl << "File " << readsFileName;
	terminate (2);
      }
      std::ifstream readsFile;
      readsFile.open(readsFileName.c_str());
      unsigned int chrNumFragileSites=0;
      unsigned int chrReads=0;
      while (!readsFile.eof()){
	std::string line;
	getline(readsFile,line);
	unsigned int chrPos, hits;
	size_t c=line.size();
	for (c=line.size(); c!=ULONG_MAX; --c)
	  {
	    if(line[c]==' ')
	      break;
	  }
	//std::cout << std::endl << c << " " << ULONG_MAX;WAITUSER;
	if (c!=ULONG_MAX){
	  // std::cout << std::endl << line.substr(0,c) << " " << line.substr(c);// WAITUSER
	  //std::cout << std::endl << c;
	  std::stringstream ss_chrPos(line.substr(c));
	  std::stringstream ss_hits(line.substr(0,c));
	  ss_chrPos>>chrPos;
	  ss_hits>>hits;
	  chrReads+=hits;
	  std::vector <CYTOBAND>::iterator i;
	  unsigned int bandNum=0;
	  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE"){
	    for (i=cytoBands.begin();i<cytoBands.end();++i,++bandNum){
	      if ((*i).chrNum==chrNum && chrPos>=(*i).chrStart && chrPos<=(*i).chrEnd){
		//std::cout << std::endl <<  chrPos << "\t" << (*i).name;// WAITUSER;
		//chrNumFragileSites+=hits;
		cytobandReads[bandNum]+=hits;

	      }
	    }
	  }
	}
      }
       if(DEBUG) std::cout << std::endl << "Total Reads = " << chrReads ;//WAITUSER;
      chrReadsTwin+=chrReads;
       if(DEBUG) std::cout << std::endl << chrNumFragileSites*1.0 / chrReads;
      if (readsFile.is_open()) readsFile.close();
      totalNumFragileSites+=chrNumFragileSites;
      totalReads+=chrReads;
    }
     if(DEBUG) std::cout << std::endl << "Total chrReadstwin = " << chrReadsTwin;
    chrReadsFile << convert_chrnum_to_string(chrNum) << "\t" << chrReadsTwin << std::endl;
  }
   if(DEBUG) std::cout << std::endl << "Total Reads close = " << totalReads;

  std::string chrReadFileName1=sequenceDataDir2+"/chrReads.txt";
  std::ofstream chrReadsFile1(chrReadFileName1.c_str());
  for (std::vector <int>::iterator i1=chrNumList.begin(); i1<chrNumList.end();++i1){
    unsigned int chrNum=*i1;
    std::stringstream ss_chrNum;
    ss_chrNum << chrNum;
    std::string s_chrNum;
    ss_chrNum >> s_chrNum;
    if (chrNum>=X_CHROMOSOME) s_chrNum=chrNum-X_CHROMOSOME+'X';
    unsigned int chrReadsTwin=0;
    for (unsigned int s=0; s<2; ++s){
      char strand='+';
      if (s==0) strand='+';
      else if(s==1) strand='-';
       if(DEBUG) std::cout << std::endl << "chr"+s_chrNum+strand;//WAITUSER;
      std::string readsFileName=sequenceDataDir2+"/chr"+s_chrNum+strand+".hits.txt";
      if (checkFile(readsFileName)==-1) {
	 if(DEBUG) std::cout << std::endl << "File " << readsFileName;
	terminate (2);
      }
      std::ifstream readsFile;
      readsFile.open(readsFileName.c_str());
      unsigned int chrNumFragileSites=0;
      unsigned int chrReads1=0;
      while (!readsFile.eof()){
	std::string line;
	getline(readsFile,line);
	unsigned int chrPos, hits;
	size_t c=line.size();
	for (c=line.size(); c!=ULONG_MAX; --c)
	  {
	    if(line[c]==' ')
	      break;
	  }
	//std::cout << std::endl << c << " " << ULONG_MAX;WAITUSER;
	if (c!=ULONG_MAX){
	  // std::cout << std::endl << line.substr(0,c) << " " << line.substr(c);// WAITUSER
	  //std::cout << std::endl << c;
	  std::stringstream ss_chrPos(line.substr(c));
	  std::stringstream ss_hits(line.substr(0,c));
	  ss_chrPos>>chrPos;
	  ss_hits>>hits;
	  chrReads1+=hits;
	  std::vector <CYTOBAND>::iterator i;
	  unsigned int bandNum=0;
	  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE"){
	    for (i=cytoBands.begin();i<cytoBands.end();++i,++bandNum){
	      if ((*i).chrNum==chrNum && chrPos>=(*i).chrStart && chrPos<=(*i).chrEnd){
		//std::cout << std::endl <<  chrPos << "\t" << (*i).name;// WAITUSER;
		//chrNumFragileSites+=hits;
		cytobandReads[bandNum]+=hits;

	      }
	    }
	  }
	}
      }
       if(DEBUG) std::cout << std::endl << "Total Reads distant = " << chrReads1 ;//WAITUSER;
      chrReadsTwin+=chrReads1;
       if(DEBUG) std::cout << std::endl << chrNumFragileSites*1.0 / chrReads1;
      if (readsFile.is_open()) readsFile.close();
      totalNumFragileSites1+=chrNumFragileSites;
      totalReads1+=chrReads1;
    }
     if(DEBUG) std::cout << std::endl << "Total chrReadstwin = " << chrReadsTwin;
    chrReadsFile1 << convert_chrnum_to_string(chrNum) << "\t" << chrReadsTwin << std::endl;
  }
   if(DEBUG) std::cout << std::endl << "Total Reads = " << totalReads1;

  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE") {
     if(DEBUG) std::cout << std::endl << "cytobands";
    std::ofstream cytoBandsReadsFile;
    std::string cbrFileName=sequenceDataDir1+"/cbrFile.txt";
    cytoBandsReadsFile.open(cbrFileName.c_str());
    std::ofstream cytoBandsReadsFile1;
    std::string cbrFileName1=sequenceDataDir2+"/cbrFile.txt";
    cytoBandsReadsFile1.open(cbrFileName1.c_str());

    for(unsigned int n=0; n<cytoBands.size(); ++n){
      if (cytoBands[n].chrNum==X_CHROMOSOME || cytoBands[n].chrNum==Y_CHROMOSOME)

         if(DEBUG) std::cout << std::endl << char(cytoBands[n].chrNum-X_CHROMOSOME+'X') << cytoBands[n].name << "\t" <<  cytobandReads[n];
        else if(DEBUG) std::cout << std::endl << cytoBands[n].chrNum << cytoBands[n].name << "\t" <<  cytobandReads[n];
      cytoBandsReadsFile << convert_chrnum_to_string(cytoBands[n].chrNum) << "\t" << cytoBands[n].name << "\t" <<  cytobandReads[n] << std::endl ;
      cytoBandsReadsFile1 << convert_chrnum_to_string(cytoBands[n].chrNum) << "\t" << cytoBands[n].name << "\t" <<  cytobandReads[n] << std::endl ;
    }
     if(DEBUG) std::cout << std::endl << totalNumFragileSites*1.0 / totalReads << std::endl;
    if (chrReadsFile.is_open()) chrReadsFile.close();
    if (cytoBandsReadsFile.is_open()) cytoBandsReadsFile.close();
     if(DEBUG) std::cout << std::endl << totalNumFragileSites1*1.0 / totalReads1 << std::endl;
    if (chrReadsFile1.is_open()) chrReadsFile1.close();
    if (cytoBandsReadsFile1.is_open()) cytoBandsReadsFile1.close();
  }
}

//*******************************************************************************
void cytoband_function1(std::string reads_directory)     {
//----------------------------------------------------------------------------
  std::string sequenceDataDir3="./";
  sequenceDataDir3=reads_directory+"_no_barcode";
  std::ifstream cbandsFile;
  std::ifstream fragileCytoBandsFile;
  std::vector <CYTOBAND> cytoBands;
  std::vector <CYTOBAND> fragileBands;

  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE"){
    //std::ifstream cbandsFile;
    cbandsFile.open(cbandsFileName.c_str());
    //std::vector <CYTOBAND> cytoBands;
    while (!cbandsFile.eof()){
      std::string line;
      getline (cbandsFile, line);
      unsigned int columns=0;
      size_t tab;
      std::string s_chrNum, s_chrStart, s_chrEnd, s_name;
      do{
	tab=line.find_first_of('\t');
	columns++;
	switch(columns){
	case 1: s_chrNum=line.substr(0,tab);
	  break;
	case 2: s_chrStart=line.substr(0,tab);
	  break;
	case 3: s_chrEnd=line.substr(0,tab);
	  break;
	case 4: s_name=line.substr(0,tab);
	  break;
	default:
	  break;
	}
	line = line.substr(tab+1, line.size());
      }while(tab!=ULONG_MAX);
      //std::cout << std::endl << extract_chrnum(s_chrNum) << "\t" << s_chrStart << "\t" << \
      s_chrEnd << "\t" << s_name; WAITUSER;
      if (extract_chrnum(s_chrNum, true)>0){
	CYTOBAND cytoBandTmp;
	std::stringstream ss_chrStart(s_chrStart);
	std::stringstream ss_chrEnd(s_chrEnd);
	unsigned int chrStart;
	unsigned int chrEnd;
	ss_chrStart >> chrStart;
	ss_chrEnd >> chrEnd;
	cytoBandTmp.chrNum=extract_chrnum(s_chrNum, true);
	cytoBandTmp.chrStart=chrStart;
	cytoBandTmp.chrEnd=chrEnd;
	cytoBandTmp.name=s_name;
	cytoBands.push_back(cytoBandTmp);
      }
    }
     if(DEBUG) std::cout << std::endl << cytoBands.size() << " Cytobands parsed"; //WAITUSER;
    ///////Processing Fragile Cytobands///////////////////
    //std::ifstream fragileCytoBandsFile;
    fragileCytoBandsFile.open(fbandsFileName.c_str());
    //std::vector <CYTOBAND> fragileBands;
    while(!fragileCytoBandsFile.eof()){
      std::string line;
      getline(fragileCytoBandsFile,line);
      unsigned int pqPos=UINT_MAX;
      for(size_t n=0;n<=line.size();++n)
	if (line[n]=='p'||line[n]=='q'){
	  pqPos=n;
	  break;
	}
      if (pqPos!=UINT_MAX){
	//std::cout << std::endl << pqPos << line[pqPos] ; WAITUSER;
	std::string s_chrNum=line.substr(0,pqPos);
	unsigned int chrNum=extract_chrnum(s_chrNum, true);
	std::string name=line.substr(pqPos) ;
	//std::cout << std::endl << chrNum << " " << name ;//WAITUSER
	std::vector <CYTOBAND>::iterator i;
	unsigned int numbands=0;
	bool foundband=false;
	std::string::size_type fcbDot=name.find_first_of('.');
	for (i=cytoBands.begin();i<cytoBands.end();++i){
	  if ((*i).chrNum==chrNum){
	    if(fcbDot!=std::string::npos){//there's a dot in fcb
	      if ((*i).name.substr(0,name.size())==name){
		//std::cout <<  "\tSubstring match " << (*i).name.substr(0,name.size()) << " " << name <<;WAITUSER;
		fragileBands.push_back(*i);
		++numbands;
		foundband=true;
		//break;
	      }
	    }
	    else if (fcbDot==std::string::npos){//did not find a dot in fcb
	      std::string cbName=(*i).name;
	      std::string::size_type cbDot=cbName.find_first_of('.');
	      if (cbDot==std::string::npos){//did not find a dot in cb
		if (cbName.substr(0,cbDot)==name){
		  fragileBands.push_back(*i);
		  ++numbands;
		  foundband=true;
		  break;
		}
	      }
	      else {
		if (cbName.substr(0,cbDot)==name){
		  fragileBands.push_back(*i);
		  ++numbands;
		  foundband=true;
		}
	      }
	    }
	  }
	}
	if (foundband==false){ std::cout << std::endl  << name << " not found in chr " << chrNum ; WAITUSER;
	}
	else {//std::cout << std::endl << name << " found " << numbands << " times";
	}
      }

    }
    if (fragileCytoBandsFile.is_open()) fragileCytoBandsFile.close();
    if (cbandsFile.is_open()) cbandsFile.close();

    {
      std::vector <CYTOBAND>::iterator i;
      std::vector <std::string> fragileCytoNames;
      for (i=fragileBands.begin();i<fragileBands.end();++i){ //dump fragile cyto names to a vector
	//if(chrPos>=(*i).chrStart && chrPos<=(*i).chrEnd)
	//std::cout << std::endl << (*i).chrNum << "\t" << (*i).name << "\t" << (*i).chrStart << "\t" << (*i).chrEnd;// WAITUSER;
	std::stringstream s;
	s<<convert_chrnum_to_string((*i).chrNum)<<(*i).name;
	fragileCytoNames.push_back(s.str());
      }
      sort(fragileCytoNames.begin(),fragileCytoNames.end());
      std::ofstream fcbFile ("fragilecytobands.txt");
      for (std::vector<std::string>::iterator i=fragileCytoNames.begin();i<fragileCytoNames.end();++i){
	fcbFile << (*i) << std::endl;
      }
      if (fcbFile.is_open())fcbFile.close();
    }
  }

  /////////Processing bowtie hits files//////////////////////////
  unsigned int totalNumFragileSites=0;
  unsigned int totalReads=0;
  std::vector <int> chrNumList;
  std::vector <int> cytobandReads(cytoBands.size(),0);
  if(org_name=="human" || org_name=="HUMAN"){
    for (unsigned int chrNum=1; chrNum<=22; ++chrNum)
      chrNumList.push_back(chrNum);
    chrNumList.push_back(X_CHROMOSOME); chrNumList.push_back(Y_CHROMOSOME);
    //std::vector <int> cytobandReads(cytoBands.size(),0);
  }
  else if(org_name=="mouse" || org_name=="MOUSE"){
    for (unsigned int chrNum=1; chrNum<=19; ++chrNum)
      chrNumList.push_back(chrNum);
    chrNumList.push_back(X_CHROMOSOME); chrNumList.push_back(Y_CHROMOSOME);
  }
  else if(org_name=="yeast" || org_name=="YEAST"){
    for (unsigned int chrNum=1; chrNum<=16; ++chrNum)
      chrNumList.push_back(chrNum);
    //chrNumList.push_back(M_CHROMOSOME);
  }
  std::string chrReadFileName=sequenceDataDir3+"/chrReads.txt";
  std::ofstream chrReadsFile(chrReadFileName.c_str());
  for (std::vector <int>::iterator i=chrNumList.begin(); i<chrNumList.end();++i){
    unsigned int chrNum=*i;
    std::stringstream ss_chrNum;
    ss_chrNum << chrNum;
    std::string s_chrNum;
    ss_chrNum >> s_chrNum;
    if (chrNum>=X_CHROMOSOME) s_chrNum=chrNum-X_CHROMOSOME+'X';
    unsigned int chrReadsTwin=0;
    for (unsigned int s=0; s<2; ++s){
      char strand='+';
      if (s==0) strand='+';
      else if(s==1) strand='-';
      if(DEBUG) std::cout << std::endl << "chr"+s_chrNum+strand;//WAITUSER;
      std::string readsFileName=sequenceDataDir3+"/chr"+s_chrNum+strand+".hits.txt";
      if (checkFile(readsFileName)==-1) {
	 if(DEBUG) std::cout << std::endl << "File " << readsFileName;
	terminate (2);
      }
      std::ifstream readsFile;
      readsFile.open(readsFileName.c_str());
      unsigned int chrNumFragileSites=0;
      unsigned int chrReads=0;
      while (!readsFile.eof()){
	std::string line;
	getline(readsFile,line);
	unsigned int chrPos, hits;
	size_t c=line.size();
	for (c=line.size(); c!=ULONG_MAX; --c)
	  {
	    if(line[c]==' ')
	      break;
	  }
	//std::cout << std::endl << c << " " << ULONG_MAX;WAITUSER;
	if (c!=ULONG_MAX){
	  // std::cout << std::endl << line.substr(0,c) << " " << line.substr(c);// WAITUSER
	  //std::cout << std::endl << c;
	  std::stringstream ss_chrPos(line.substr(c));
	  std::stringstream ss_hits(line.substr(0,c));
	  ss_chrPos>>chrPos;
	  ss_hits>>hits;
	  chrReads+=hits;
	  std::vector <CYTOBAND>::iterator i;
	  unsigned int bandNum=0;
	  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE"){
	    for (i=cytoBands.begin();i<cytoBands.end();++i,++bandNum){
	      if ((*i).chrNum==chrNum && chrPos>=(*i).chrStart && chrPos<=(*i).chrEnd){
		//std::cout << std::endl <<  chrPos << "\t" << (*i).name;// WAITUSER;
		//chrNumFragileSites+=hits;
		cytobandReads[bandNum]+=hits;

	      }
	    }
	  }
	}
      }
     // std::cout << std::endl << "Total Reads = " << chrReads ;//WAITUSER;
      chrReadsTwin+=chrReads;
    //  std::cout << std::endl << chrNumFragileSites*1.0 / chrReads;
      if (readsFile.is_open()) readsFile.close();
      totalNumFragileSites+=chrNumFragileSites;
      totalReads+=chrReads;
    }
   // std::cout << std::endl << "Total chrReadstwin = " << chrReadsTwin;
    chrReadsFile << convert_chrnum_to_string(chrNum) << "\t" << chrReadsTwin << std::endl;
  }
//  std::cout << std::endl << "Total Reads close = " << totalReads;

  if(org_name=="human" || org_name=="HUMAN" || org_name=="mouse" || org_name=="MOUSE") {
   // std::cout << std::endl << "cytobands";
    std::ofstream cytoBandsReadsFile;
    std::string cbrFileName=sequenceDataDir3+"/cbrFile.txt";
    cytoBandsReadsFile.open(cbrFileName.c_str());

    for(unsigned int n=0; n<cytoBands.size(); ++n){
      if (cytoBands[n].chrNum==X_CHROMOSOME || cytoBands[n].chrNum==Y_CHROMOSOME)

        std::cout << std::endl << char(cytoBands[n].chrNum-X_CHROMOSOME+'X') << cytoBands[n].name << "\t" <<  cytobandReads[n];
      else
	std::cout << std::endl << cytoBands[n].chrNum << cytoBands[n].name << "\t" <<  cytobandReads[n];
      cytoBandsReadsFile << convert_chrnum_to_string(cytoBands[n].chrNum) << "\t" << cytoBands[n].name << "\t" <<  cytobandReads[n] << std::endl ;
    }
    std::cout << std::endl << totalNumFragileSites*1.0 / totalReads << std::endl;
    if (chrReadsFile.is_open()) chrReadsFile.close();
    if (cytoBandsReadsFile.is_open()) cytoBandsReadsFile.close();
  }
}
