
#include "../include/hygestat.h"
#include "../include/flags_and_barcodes.h"

const char *short_options = "o:O:G:g:F:f:N:n:d:b:t:s:a:A:w:c:y:Y:vqp:i:hr:m:k";

//************************** MAIN FUNCTION ****************************************//
int main (int argc, char **argv) {
//*********************************************************************************//
    std::string cmdline;
    std::vector <std::string> all_fastq;
    std::vector <std::string> nat_all_fastq;
    std::fstream fs;
    int  windows_type;

    write_lab_info();
    if (argc < 3) { print_usage(); return 1; }
    for (int i = 0; i < argc; i++) {cmdline+=argv[i]; cmdline+=" "; }
    int parse_ret = parse_options(argc,argv);
    if (parse_ret) return parse_ret;
    fprintf(stderr, "Command line:\n%s\n", cmdline.c_str());
    set_ouput_files_and_stat();
    if(!TIME_POINT){   // Run without time points
      all_fastq.push_back(fastqFile1);
      all_fastq.push_back(fastqFile2);
      nat_all_fastq.push_back(nature_of_data1);
      nat_all_fastq.push_back(nature_of_data2);
      if(nature_of_analysis=="with_correction"){
        all_fastq.push_back(fastqFile3);
        all_fastq.push_back(fastqFile4);
        all_fastq.push_back(fastqFile5);
        nat_all_fastq.push_back(nature_of_data3);
        nat_all_fastq.push_back(nature_of_data4);
        nat_all_fastq.push_back(nature_of_data5);
      }

      fs.open ("bless_and_sequencing_quality.txt", std::fstream::in | std::fstream::out | std::fstream::app);
      for (unsigned j=0;j<all_fastq.size();j++)  {
       if(!BOWTIE) data_mapping_fastq(all_fastq[j], nat_all_fastq[j]);
            if( nat_all_fastq[j]=="bless"){
                std::string bttfile_close     = all_fastq[j]+".dump.fa."+close_bar+".fasta.bt.btt";
                std::string bttfile_dist      = all_fastq[j]+".dump.fa."+distant_bar+".fasta.bt.btt";
                size_t num_reads_close   = count_line_in_file(bttfile_close);
                size_t num_reads_distant = count_line_in_file(bttfile_dist);
                size_t total_reads       = count_line_in_file(all_fastq[j]+".dump.fa");
               // fs<< all_fastq[j].erase(all_fastq[j].find(".fastq"))<<" \t "<<total_reads<<" \t "<<num_reads_close+num_reads_distant<<"\n"; // uncompleted!
            }
       }
       windows_type = return_wind_type(wind_run);
       switch ( windows_type ) {
       case 10:
            two_samples_hygeostat(fastqFile1, nature_of_data1,
                                  fastqFile2, nature_of_data2,
			                      fastqFile3, nature_of_data3,
                                  fastqFile4, nature_of_data4,
                                  fastqFile5, nature_of_data5,
                                  nature_of_analysis,
			                      wind_run[0]
                          );
            break;
       case 20:
            for(unsigned int i = 0; i < wind_run.size();i++){
                    two_samples_hygeostat(fastqFile1, nature_of_data1,
                                          fastqFile2, nature_of_data2,
                                          fastqFile3, nature_of_data3,
                                          fastqFile4, nature_of_data4,
			                              fastqFile5, nature_of_data5,
                                          nature_of_analysis,
                                          wind_run[i]
                          );
            }
            break;
       case 30:
            for( int w = wind_run[0]; w <= wind_run[1];w+=wind_run[2]){
                    two_samples_hygeostat(fastqFile1, nature_of_data1,
                                          fastqFile2, nature_of_data2,
                                          fastqFile3, nature_of_data3,
                                          fastqFile4, nature_of_data4,
			                              fastqFile5, nature_of_data5,
                                          nature_of_analysis,
                                          w
                          );
            }
            break;
       default:            // Note the colon, not a semicolon
           std::cout<<"Error, bad input resolution, quitting\n";
           break;
       }
    }else{  // if running time points
      for (unsigned i=0;i<treat_fastq_tim.size();i++){
        all_fastq.push_back(treat_fastq_tim[i]);
        all_fastq.push_back(cont_fastq_tim[i]);
        nat_all_fastq.push_back(nat_treat_tim[i]);
        nat_all_fastq.push_back(nat_cont_tim[i]);
      }

      min_mapp_reads=1e15;
      size_t num_reads;
      for (unsigned j=0;j<all_fastq.size();j++){  // do the mapping and find the sample with minimum reads
            std::string bttfile_close,bttfile_dist;
            if(DEBUG)    std::cerr<<"File = "<<all_fastq[j]<<" nature ="<<nat_all_fastq[j]<<std::endl;
            if(!BOWTIE)  data_mapping_fastq(all_fastq[j],nat_all_fastq[j]); // run and create the .btt files
            std::fstream fs;
            fs.open ("bless_and_sequency_quality.txt", std::fstream::in | std::fstream::out | std::fstream::app);
            if(nat_all_fastq[j]=="bless")  {
                bttfile_close  = all_fastq[j]+".dump.fa."+close_bar+".fasta.bt.btt";
                bttfile_dist   = all_fastq[j]+".dump.fa."+distant_bar+".fasta.bt.btt";
                size_t num_reads_close   = count_line_in_file(bttfile_close);
                size_t num_reads_distant = count_line_in_file(bttfile_dist);
                size_t total_reads       = count_line_in_file(all_fastq[j]+".dump.fa");
                num_reads = num_reads_close < num_reads_distant ? num_reads_close : num_reads_distant;
                fs<< all_fastq[j].erase(all_fastq[j].find(".fastq"))<<" \t "<<total_reads<<" \t "<<num_reads_close+num_reads_distant<<"\n"; // uncompleted!
            }else if(nat_all_fastq[j]=="genomic") {
                std::string bttfile = all_fastq[j]+".bt.btt";
                num_reads      = count_line_in_file(bttfile);
               // No need to compute
            }else { if(DEBUG) std::cerr<<"TYPE NOT DEFINED...\n";}

            if((RAREFACTION)&&(num_reads <= min_mapp_reads)) min_mapp_reads = num_reads;
          }
            if((DEBUG)&&(RAREFACTION)) std::cerr<<"RAREFICATION MINIMUM : "<< min_mapp_reads <<std::endl;

         for (unsigned h=0;h<treat_fastq_tim.size();h++){
              two_samples_hygeostat(treat_fastq_tim[h],nat_treat_tim[h],cont_fastq_tim[h] ,nat_cont_tim[h],
                        "na","na","na","na","na","na","without_correction",15000);
              if(RAREFACTION){
                 for (int kk=rand_mini;kk<=rand_maxi;kk+=rand_step){
                      std::string rare = to_string(kk);
                      pick_n_random_lines(treat_fastq_tim[h]+".dump.fa."+close_bar+".fasta.bt.btt", min_mapp_reads,kk); //close
                      std::string new_treat = "sample_"+rare+"_"+treat_fastq_tim[h];
                      std::string new_cont  =  cont_fastq_tim[h];
                      std::string close="sample_"+rare+"_"+treat_fastq_tim[h]+".dump.fa."+close_bar+".fasta.bt.btt"; // avoid sampling the distant will save time
                      std::string dist ="sample_"+rare+"_"+treat_fastq_tim[h]+".dump.fa."+distant_bar+".fasta.bt.btt";
                      char tmp[255];
                      std::sprintf(tmp,"%s %s %s","cp ",close.c_str(),dist.c_str());
                      int t_cmd1 = std::system(tmp);  if(t_cmd1 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
                      two_samples_hygeostat(new_treat,nat_treat_tim[h],new_cont ,nat_cont_tim[h],
                                            "na","na","na","na","na","na","without_correction",15000);
                 }
                  std::string cle="rm -r sample_* *_no_barcode";
                  char stat[255],clean[255],rare_stat[255];
                  std::sprintf(stat,"%s %s","../tools/summarize_rarefaction.sh",treat_fastq_tim[h].c_str());
                  int t_cmd2 = std::system(stat); if(t_cmd2 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
                  std::string file_sum = "/tmp/"+treat_fastq_tim[h]+"_summary.txt";
                  stat_rarefaction(file_sum);
                  int t_cmd3 = std::system(cle.c_str()); if(t_cmd3 == -1) std::cerr<<" ERROR:: Can't create the btt file...\n";
            }
        }
      }
      stat_funct();
   if(dont_keep_files) clean_up_unwanted_files();
} // end main
//--------------------------------------------------------------------------
static struct option long_options[] = {
//--------------------------------------------------------------------------
{"genomeDir",	            required_argument,		 0,			 'G'},
{"genomeType",		    required_argument,		 0,			 'g'},
{"inputFile",		    required_argument,		 0,			 'i'},
{"fastqControl",            required_argument,		 0,			 'F'},
{"natureControl",	    required_argument,		 0,			 'N'},
{"fastqTreatment",          required_argument,		 0,			 'f'},
{"natureTreatment",         required_argument,		 0,			 'n'},
{"data_typ",	            required_argument,		 0,			 'o'},
{"output-file",	            required_argument,		 0,			 'O'},
{"mapDir",                  required_argument,		 0,			 'm'},
{"verbose",	            no_argument,		 0,			 'v'},
{"wigOutput",	            no_argument,		 0,			 'k'},
{"quiet",		    no_argument,		 0,			 'q'},
{"help",		    no_argument,		 0,			 'h'},
//{"no-update-check",         no_argument,             0,          OPT_NO_UPDATE_CHECK},
{"num-threads",		    required_argument,           0,          'p'},
{"dataDir",                 required_argument,       0,          'd'},
{"noBowtie ",               required_argument,       0,          'b'},
{"telomere",		    required_argument,       0,          't'},
{"time-serie",		    required_argument,       0,          's'},
{"time-serie-file",         required_argument,       0,          TIME_SERIE_FILE},
{"time-serie-rare",         required_argument,	 	 0,	         TIME_SERIE_RARE},
{"time-serie-rand",         required_argument,		 0,			 TIME_SERIE_RAND},
{"preCheck",                required_argument,		 0,			 'a'},
{"postCheck",		    required_argument,		 0,			 'A'},
{"resolution",              required_argument,		 0,			 'r'},
{"windowsAdvance",          required_argument,       0,			 'w'},
{"correction",              required_argument,       0,			 'c'},
{"corSample1",              required_argument,       0,			 CORRECTION_SAMPLE_1},
{"corSample1Nat",           required_argument,	     0,	         CORRECTION_SAMPLE_1_NAT},
{"corSample2",              required_argument,	     0,	         CORRECTION_SAMPLE_2},
{"corSample2Nat",           required_argument,       0,          CORRECTION_SAMPLE_2_NAT},
{"corSample3",	            required_argument,		 0,			 CORRECTION_SAMPLE_3},
{"corSample3Nat",	    required_argument,		 0,			 CORRECTION_SAMPLE_3},
{"cytoPath",	            required_argument,		 0,			 'y'},
{"fragilePath",             required_argument,		 0,			 'Y'},
{0, 0, 0, 0} // terminator
};
//---------------------------------------------------------------------
void print_usage()  {
//----------------------------------------------------------------------
    //NOTE: SPACES ONLY, bozo
  //  fprintf(stderr, "hygestat v%s\n", PACKAGE_VERSION);
  //  fprintf(stderr, "linked against Boost version %d\n", BOOST_VERSION);
    fprintf(stderr, "------------------------------------------*********  HYGESTAT MENU  ***********---------------------------------------------- \n");
    fprintf(stderr, "Usage:  hygestat currently support 2 types of usage (1): hygestat -i config_file.txt [default]                                \n");
    fprintf(stderr, "or (2): hygestat [options] <control.fastq> <treatment.fastq>                                                                  \n");
    fprintf(stderr, "Only option for (1) :                                                                                                         \n");
    fprintf(stderr, "  -i/--inputFile         <config_file.txt>           configuration file for method 1                                          \n");
    fprintf(stderr, "General Options for (2) :                                                                                                     \n");
    fprintf(stderr, "  -o/--data-typ        <hygest-results>            write all output files to this directory              [ default:    fastq ]\n");
    fprintf(stderr, "  -O/--output-file       <output.txt>                main output file name                         [ default:     output.txt ]\n");
    fprintf(stderr, "  -p/--num-threads       <1>                         number of threads used during analysis                [ default:      1 ]\n");
    fprintf(stderr, "  -G/--genomeDir         </path/to/genome>           absolute directory path to the reference genome       [ default:     ./ ]\n");
    fprintf(stderr, "  -g/--genomeType        <genome>                    genome type (currently support mouse, human, and yeast genome)           \n");
    fprintf(stderr, "  -F/--fastqControl      <control.fastq>             control fastq file                                    [ no default:     ]\n");
    fprintf(stderr, "  -N/--natureControl     <bless/genomic>             nature of the control data                            [ default: genomic]\n");
    fprintf(stderr, "  -f/--fastqTreatment    <treatment.fastq>           treatment fastq file                                  [ no default:     ]\n");
    fprintf(stderr, "  -n/--natureTreatment   <bless/genomic>             nature of the treatment data                            [ default: bless]\n");
    fprintf(stderr, "  -d/--dataDir           </path/to/fastq/files>      absolute directory path to the data                   [ default:     ./ ]\n");
    fprintf(stderr, "  -m/--mapDir           </path/to/mappability/files> absolute directory path to the data             [ default:no mappability]\n");
    fprintf(stderr, "  -b/--noBowtie          <FALSE/TRUE>                skip alignement with bowtie                           [ default:  FALSE ]\n");
    fprintf(stderr, "  -t/--telomere          <FALSE/TRUE>                estimate telomeres contribution to the analysis       [ default:  FALSE ]\n");
    //Running time serie bless analysis
    fprintf(stderr, "  -s/--time-serie        <FALSE/TRUE>                time serie analysis [provide a configuration file]    [ default:  FALSE ]\n");
    fprintf(stderr, "  --time-serie-file      <config_time_serie.txt>     time serie configuration file                                            \n");
    fprintf(stderr, "  --time-serie-rare      <FALSE/TRUE>                rarefaction correction of time serie data             [ default:  FALSE ]\n");
    fprintf(stderr, "  --time-serie-rand      <20>                        number of rarefaction data [will use number +1]            [ default: 20]\n");
    // pre and post alignment analysis
    fprintf(stderr, "  -a/--preCheck          <FALSE/TRUE>                use FastQC for pre-quality check of all fastq files     [ default: FALSE]\n");
    fprintf(stderr, "  -k/--wigOutput                                     output wig files for visualization                      [ default: no  ]\n");
    fprintf(stderr, "  -A/--postCheck         <FALSE/TRUE>                use samstat for post alignment analysis                 [ default: FALSE]\n");
    // Bless resolution
    fprintf(stderr, "  -r/--resolution        <10250>                       resolution of the dsb detection                       [ default: 10250]\n");
    fprintf(stderr, "  -w/--windowsAdvance    <1>                         windows advance                                             [ default: 1]\n");
    //Running bless with/without correction
    fprintf(stderr, "  -c/--correction          <false/true>                run bless with correction                           [ default:  false ]\n");
    fprintf(stderr, "  --corSample1           <correction_sample_1.fastq> correction sample 1                                   [ no default:     ]\n");
    fprintf(stderr, "  --corSample1Nat        <bless/genomic>             nature of sample 1 (bless or genomic)                 [ no default:     ]\n");
    fprintf(stderr, "  --corSample2           <correction_sample_1.fastq> correction sample 2                                   [ no default:     ]\n");
    fprintf(stderr, "  --corSample2Nat        <bless/genomic>      nature of sample 2 (bless or genomic)                        [ no default:     ]\n");
    fprintf(stderr, "  --corSample3           <correction_sample_1.fastq> correction sample 3                                   [ no default:     ]\n");
    fprintf(stderr, "  --corSample3Nat        <bless/genomic>             nature of sample 3 (bless or genomic)                 [ no default:     ]\n");
    // path to fragile band and cyto bands
    fprintf(stderr, "  -y/--cytoPath          <cytoband-human.txt>        absolute path to the cytoband file                    [ no default:     ]\n");
    fprintf(stderr, "  -Y/--fragilePath       <fragile-band-human.txt>    absolute path to the fragile band file                [ no default:     ]\n");

    fprintf(stderr, "\nAdvanced Program Behavior Options:                                                                                          \n");
    fprintf(stderr, "  -v/--verbose                                       log-friendly verbose processing (no progress bar)     [ default:  FALSE ]\n");
    fprintf(stderr, "  -q/--quiet                                         log-friendly quiet processing (no progress bar)       [ default:  FALSE ]\n");
    fprintf(stderr, "  -h/--help                                          print help menu                                       [ no default:     ]\n");
    fprintf(stderr, "  --no-update-check                                  do not contact server to check for update availability[ default:  FALSE ]\n");
} // end function

//----------------------------------------------------------------------------------------------------------------------------
int parse_options(int argc, char** argv)  {
    int option_index = 0;
    int next_option;
	bool F_set = false;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
			case 'o':
				fasta_or_fastq  = optarg;
				data_name       = fasta_or_fastq;
				break;
			case 'O':
				output_file     = optarg;
				outputFileName  = output_file;
				break;
			case 'G':
				genome_dir      = optarg;
				break;
			case 'g':
			        genome_type     = optarg;
			        org_name        = genome_type;
				break;
			case 'F':
				fastq_control   = optarg;
				fastqFile2      = fastq_control;
				break;
			case 'f':
				fastq_treatment = optarg;
				fastqFile1      = fastq_treatment;
				break;
			case 'N':
				nature_control  = optarg;
				nature_of_data2 = nature_control;
				break;
			case 'n':
				nature_treatment = optarg;
				nature_of_data1  = nature_treatment;
				break;
			case 'm':
				mapDir           = optarg;
				no_mappability   = false;
				break;
			case 'k':
				no_bedgraph_file = false;
				break;
			case 'd':
				data_dir = optarg;
				break;
			case 'b':
				dont_run_bowtie = optarg;
				if(dont_run_bowtie=="true") BOWTIE = true;
				break;
			case 't':
				telomere = optarg;
				if(telomere=="true") {
				    TELOMERE      = true;
				    TELOMERE_FILE = true;
				}
				break;
			case 's':
				time_serie = optarg;
				if(time_serie!="false") TIME_POINT = true;
			//	if(time_serie){
				  //define some values
				//}
				break;
			case 'a':
				pre_quality_check = optarg;
				if(pre_quality_check!="false")   PRE_QUALITY = true;
				break;
			case 'A':
				post_mapping_check = optarg;
				if(post_mapping_check!="false")  POS_QUALITY = true;
				break;
			case 'r':
			    resolution    = optarg;
			    resol_flag    = true;
				break;
			case 'w':
				WINDOWADVANCE = atoi(optarg);
				WINDOWADVANCE += WINDOWSIZE;
				break;
			case 'c':
				correction = optarg;
				//nature_of_analysis = "without_correction";
				if(correction == "true" || correction == "TRUE") nature_of_analysis = "with_correction";
				break;
			case 'y':
				cytoband       = optarg;
				cbandsFileName = cytoband;
				break;
			case 'Y':
				fragile_band   = optarg;
				fbandsFileName = fragile_band;
				break;
			case 'v':
				verbose = true;
				DEBUG   = true;
				break;
			case 'q':
				quiet_mode = true;
				DEBUG      = false;
				break;
			case 'p':
				N_P       = atoi(optarg);
				bowt_arg  = " -l61 -m1 -n0 -r -p"   +to_string(N_P)+" ";
				bowt_arg_q= " -l61 -m1 -n0 -r -q -p"+to_string(N_P)+" ";
				break;
			case 'h':
				print_usage();
				exit(1);
				break;
			case 'i':
				input_file = optarg;
				run_type  = 1;
                if (checkFile(input_file)!=1) exit(1);
				read_config_file(input_file);
				break;
           // case OPT_NO_UPDATE_CHECK:
           //     no_update_check = true;
			//	break;
			case TIME_SERIE_FILE:
			    time_serie_file = optarg;
				break;
			case TIME_SERIE_RARE:
				time_serie_rare = optarg;
				break;
			case TIME_SERIE_RAND:
				time_serie_rand = atoi(optarg);
				rand_maxi = time_serie_rand;
				break;
			case CORRECTION_SAMPLE_1:
                 corSample1= optarg;
                 fastqFile3=corSample1;
				break;
			case CORRECTION_SAMPLE_1_NAT:
				corSample1Nat= optarg;
				nature_of_data3=corSample1Nat;
				break;
			case CORRECTION_SAMPLE_2:
				corSample2= optarg;
				fastqFile4=corSample2;
				break;
			case CORRECTION_SAMPLE_2_NAT:
				corSample2Nat= optarg;
				nature_of_data4=corSample2Nat;
				break;
			case CORRECTION_SAMPLE_3:
				corSample3= optarg;
				fastqFile5=corSample3;
				break;
			case CORRECTION_SAMPLE_3_NAT:
				corSample3Nat= optarg;
				nature_of_data5=corSample3Nat;
				break;
		    default:
				print_usage();
				return 1;
        }

    } while(next_option != -1);
   //check the run type to be 1 (with conf file) or 2 (with inline arguments)
   if (argc > 3 && (run_type == 1 )){
        std::cerr<<"*ERROR  Hygestat with configuration file can't take more than one argument *********\n";
        std::cerr<<"*ERROR  Please use the configuration file only or check usage below        *********\n";
        print_usage();
        exit(1);
    }
    //make sure the configuration file exists if run type =1
    if((run_type == 1) && (checkFile(input_file) != 1)) {
        fprintf (stderr, "Error: cannot found the configuration file %s\n", input_file.c_str());
        exit(1);
    }
    // run type 2
    if(run_type==2){

    if(resol_flag) {
        // split the resolution to allocate windows
        std::vector <std::string> wind;
        split(wind,resolution,boost::is_space());
        for(unsigned int i = 0; i < wind.size(); i++){
            wind_run.push_back(stoi(wind[i]));
        }
     }
    //check if the out dir can be created
    if (output_dir != "") {
        int retcode = mkpath(output_dir.c_str(), 0777);
        if (retcode == -1) {
            if (errno != EEXIST) {
                fprintf (stderr, "Error: cannot create directory %s\n", output_dir.c_str());
                std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
                exit(1);
            }
        }
    }
    //make sure the data directory is not empty
     if((data_dir != "") && (data_dir!="./")) {
       boost::filesystem::path p(data_dir);
       if( (!boost::filesystem::exists(p)) || (!boost::filesystem::is_directory(p)) || boost::filesystem::is_empty(p)){
            fprintf (stderr, "Error: the data directory doesn't exist, is not a directory or is empty %s\n", data_dir.c_str());
            std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
        exit(1);
       }else{
           fastqFile1 = data_dir+"/"+fastqFile1;
           fastqFile2 = data_dir+"/"+fastqFile2;
           if(nature_of_analysis=="with_correction"){
              fastqFile3 = data_dir+"/"+fastqFile3;
              fastqFile4 = data_dir+"/"+fastqFile4;
              fastqFile5 = data_dir+"/"+fastqFile5;
           }
       }
     }
 //make sur the genome dir exists and is not empty
     if(genome_dir != "") {

         if(genome_type == ""){
            std::cerr<< "Error: please provide the genome type (human/mouse/yeast)\n";
            std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
            exit(1);
          }else{
               if      (genome_type == "human") genome_human  = genome_dir;
               else if (genome_type == "mouse") genome_mouse  = genome_dir;
               else if (genome_type == "yeast") genome_yeast  = genome_dir;
               else {
               std::cerr<< "Error: please provide the correct genome type (example: human,mouse, or yeast)\n";
               std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
               exit(1);
               }
             }
        boost::filesystem::path p(genome_dir);
     // if( (!boost::filesystem::exists(p)) || (boost::filesystem::is_empty(p))){
     if(0){
            fprintf (stderr, "Error: please a correct genome file: the file %s is not readable. \n", genome_dir.c_str());
            std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
        exit(1);
       }
     }
     if(fastq_control == ""||fastq_treatment == "" || nature_control == ""|| nature_control == ""){
         std::cerr<< "Error: please a correct fastq files and their nature (bless/genomic)\n";
         std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
        exit(1);
     }else{
        boost::filesystem::path p_c(data_dir +"/"+ fastq_control);
        boost::filesystem::path p_t(data_dir +"/"+ fastq_treatment);
        if( (!boost::filesystem::exists(p_c)) || (!boost::filesystem::exists(p_t))){
          std::cerr<< "Error: the fastq files don't exist in the directory you provided \n";
          std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
           exit(1);
        }
     }
     if((nature_of_analysis == "with_correction") && ((fastqFile3 == "")||(fastqFile3 == "")||(fastqFile3 == ""))){
         std::cerr<< "Error: The run type is : with_correction, therefore you should provide correction data in fastq format \n";
         std::cerr<< "_____________Error: Run hygestat -h for help________________________________ \n";
           exit(1);
     }
   }
    return 0;
}


//--------------------------------------------------------
int mkpath(const char *s, mode_t mode) {
//---------------------------------------------------------
    char *q, *r = NULL, *path = NULL, *up = NULL;
    int rv = -1;

    if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0) return (0);

    if ((path = strdup(s)) == NULL) exit(1);

    if ((q = strdup(s)) == NULL) exit(1);

    if ((r = dirname(q)) == NULL) goto out;

    if ((up = strdup(r)) == NULL) exit(1);

    if ((mkpath(up, mode) == -1) && (errno != EEXIST)) goto out;

    if ((mkdir(path, mode) == -1) && (errno != EEXIST))
        rv = -1;
    else
        rv = 0;

out:
    if (up != NULL)
        free(up);
    free(q);
    free(path);
    return (rv);
}

//------------------------------------------------------------
int return_wind_type(std::vector <int> wind_t){
//--------------------------------------------
    int w_type;
    if(wind_t.size() < 1) {
          std::cerr<<"Error: please provide the window for hygestat\n";
          exit(0);
    }
    else if(wind_t.size() == 1) {
          w_type = 10;
    }
    else if(wind_t.size() == 2) {
          w_type = 20;
    }
    else if(wind_t.size() == 3) {
          if(wind_t[2] <= wind_t[1]){
                 w_type = 30;
          }else {w_type = 20;}
    }else {      w_type = 20;}
 return w_type;
}

//--------------------------------------------------------------------------
void clean_up_unwanted_files() {
//-------------------------------------------------------------------------
    ListAndRemoveFilesRecursively(".", ".fasta");
    ListAndRemoveFilesRecursively(".", ".bt");
    ListAndRemoveFilesRecursively(".", ".btt");
    ListAndRemoveFilesRecursively(".", ".csv");
    ListAndRemoveFilesRecursively(".", ".fa");
    ListAndRemoveFilesRecursively(".", "fastq_close_barcode");
    ListAndRemoveFilesRecursively(".", "fastq_distant_barcode");
    ListAndRemoveFilesRecursively(".", "_quality.txt");
    ListAndRemoveFilesRecursively(".", "_hgwindows.txt");

}

void stat_funct()
{
std::string line;
std::string name1 = fastqFile1.substr(fastqFile1.find_last_of("/")+1);
std::string name2 = fastqFile2.substr(fastqFile2.find_last_of("/")+1);
std::string name3 = fastqFile3.substr(fastqFile3.find_last_of("/")+1);
std::string name4 = fastqFile4.substr(fastqFile4.find_last_of("/")+1);
std::string name5 = fastqFile5.substr(fastqFile5.find_last_of("/")+1);
std::string out   = name1+"_vs_"+name2+".tab";
std::ofstream tab(out.c_str());
std::string treat_all, treat1_all, treat2_all, bar_all1, bar_all2, bar_all3, bar_all4, bar_all5;
std::string con_all, treat3_all, map1, map2, map3, map4, map5;
if(data_name=="fastq"){
treat_all = fastqFile1+".dump.fa"; treat1_all = fastqFile3+".dump.fa"; treat2_all = fastqFile4+".dump.fa";
con_all   = fastqFile2+".dump.fa"; treat3_all = fastqFile5+".dump.fa";
} else{
treat_all = fastqFile1; treat1_all = fastqFile3; treat2_all = fastqFile4;
con_all   = fastqFile2; treat3_all = fastqFile5;
}
 if(nature_of_data1=="bless"){
 if(data_name=="fastq") bar_all1  = treat_all+".TCGAGGTAGTA.fasta";
 else bar_all1  = treat_all+".dump.fa.TCGAGGTAGTA.fasta";
 map1      = bar_all1+".bt";
 } else{map1 = fastqFile1+".bt";}

 if(nature_of_data2=="bless"){
 if(data_name=="fastq") bar_all2  = con_all+".TCGAGGTAGTA.fasta";
 else bar_all2  = con_all+".dump.fa.TCGAGGTAGTA.fasta";
 map2      = bar_all2+".bt";
 } else{map2 = fastqFile2+".bt";}
 if(nature_of_data3=="bless"){
 if(data_name=="fastq") bar_all3  = treat1_all+".TCGAGGTAGTA.fasta";
 else bar_all3  = treat1_all+".dump.fa.TCGAGGTAGTA.fasta";
 map3      = bar_all3+".bt";
 } else{map3 = treat1_all+".bt";}
 if(nature_of_data4=="bless"){
 if(data_name=="fastq") bar_all4  = treat2_all+".TCGAGGTAGTA.fasta";
 else bar_all4  = treat2_all+".dump.fa.TCGAGGTAGTA.fasta";
 map4      = bar_all4+".bt";
 } else{map4 = treat2_all+".bt";}
 if(nature_of_data5=="bless"){
 if(data_name=="fastq") bar_all5  = treat3_all+".TCGAGGTAGTA.fasta";
 else bar_all5  = treat3_all+".dump.fa.TCGAGGTAGTA.fasta";
 map5      = bar_all5+".bt";
 } else{map5 = treat3_all+".bt";}
std::ifstream file1(treat_all.c_str());
std::ifstream file2(con_all.c_str());
std::ifstream file11(bar_all1.c_str());
std::ifstream file22(bar_all2.c_str());
std::ifstream file111(map1.c_str());
std::ifstream file222(map2.c_str());
while (getline(file1,line)) ++total_R1;
while (getline(file11,line)) ++barc_R1;
while (getline(file2,line)) ++total_R2;
while (getline(file22,line)) ++barc_R2;
while (getline(file111,line)) ++map_R1;
while (getline(file222,line)) ++map_R2;
perc_bar1 = 100.0*barc_R1/total_R1;
perc_bar2 = 100.0*barc_R2/total_R2;
perc_map1 = 100.0*map_R1/total_R1;
perc_map2 = 100.0*map_R2/total_R2;
//cout<<"//**********/////////////// " <<total_R1<<"\t"<<map_R1<<"\t"<<perc_bar1<<"\t"<<perc_map1<<"\t"<<total_R2<<"\t"<<map_R2<<"\t"<<perc_bar2<<"\t"<<perc_map2<<" //******************////////"<<endl;
if(nature_of_analysis=="with_corection"){
std::ifstream file3(treat1_all.c_str());
std::ifstream file4(treat2_all.c_str());
std::ifstream file5(treat3_all.c_str());
std::ifstream file33(bar_all3.c_str());
std::ifstream file44(bar_all4.c_str());
std::ifstream file55(bar_all5.c_str());
std::ifstream file333(map3.c_str());
std::ifstream file444(map4.c_str());
std::ifstream file555(map5.c_str());
while (getline(file3,line)) ++total_R3;
while (getline(file33,line)) ++barc_R3;
while (getline(file4,line)) ++total_R4;
while (getline(file44,line)) ++barc_R4;
while (getline(file5,line)) ++total_R5;
while (getline(file55,line)) ++barc_R5;
while (getline(file333,line)) ++map_R3;
while (getline(file444,line)) ++map_R4;
while (getline(file555,line)) ++map_R5;
perc_bar3 = 100.0*barc_R3/total_R3;
perc_bar4 = 100.0*barc_R4/total_R4;
perc_bar5 = 100.0*barc_R5/total_R5;
perc_map3 = 100.0*map_R3/total_R3;
perc_map4 = 100.0*map_R4/total_R4;
perc_map5 = 100.0*map_R5/total_R5;
}
tab <<std::setw(11)<<"Sample_name"<<std::setw(16)<<" #_total_of_reads"<<std::setw(16)<<" #_barcoded_reads"<<std::setw(14)<<" #_mapped_reads"<<std::setw(19)<<" %_of_barcoded_reads"<<std::setw(17)<<" %_of_mapped_reads"<<std::endl<<std::endl;
tab <<std::setw(name1.size())<<name1<<std::setw(15)<<total_R1<<std::setw(15)<<barc_R1<<std::setw(15)<<map_R1<<std::setw(15)<<perc_bar1<<std::setw(20)<<perc_map1<<std::endl<<std::endl;
tab <<std::setw(name2.size())<<name2<<std::setw(15)<<total_R2<<std::setw(15)<<barc_R2<<std::setw(15)<<map_R2<<std::setw(15)<<perc_bar2<<std::setw(20)<<perc_map2<<std::endl<<std::endl;
if(nature_of_analysis=="with_correction"){
tab <<std::setw(name3.size())<<name3<<std::setw(15)<<total_R3<<std::setw(15)<<barc_R3<<std::setw(15)<<map_R3<<std::setw(15)<<perc_bar3<<std::setw(20)<<perc_map3<<std::endl<<std::endl;
tab <<std::setw(name4.size())<<name4<<std::setw(15)<<total_R4<<std::setw(15)<<barc_R4<<std::setw(15)<<map_R4<<std::setw(15)<<perc_bar4<<std::setw(20)<<perc_map4<<std::endl<<std::endl;
tab <<std::setw(name5.size())<<name5<<std::setw(15)<<total_R5<<std::setw(15)<<barc_R5<<std::setw(15)<<map_R5<<std::setw(15)<<perc_bar5<<std::setw(20)<<perc_map5<<std::endl;
}
}

//------------------------------------------------------------------------
void ListAndRemoveFilesRecursively(const char *dir, const char* ext) {
//------------------------------------------------------------------------
    boost::filesystem::recursive_directory_iterator rdi(dir);
    boost::filesystem::recursive_directory_iterator end_rdi;
    std::string ext_str0(ext);
    for (; rdi != end_rdi; rdi++)
    {
         std::string tata((*rdi).path().string());
         if(boost::algorithm::ends_with(tata,ext_str0))
        {
            try
            {
                {
                    remove(rdi->path());
                }
            }
            catch (const std::exception& ex )
            {
                ex;
            }
        }
    }
}
