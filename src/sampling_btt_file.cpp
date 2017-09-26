#include "../include/sampling_btt_file.h"

//*********************************************FUNCTION********************************************
void sampling_btt_file(std::string reads_file,long min_map,int rand_cut){
//**************************************************************************************************
    std::string     rare = to_string(rand_cut);
    std::string  outfile = "sample_"+rare+"_"+reads_file;
    std::string   quest;
    int   line_counter = 0;
    srand(time(NULL));
    int nb_rand = rand_cut;
    std::ifstream       myfile (reads_file.c_str());
    std::ofstream       out(outfile.c_str());
    std::vector<std::string> text_file, text_file_cp;
    std::vector <std::vector<std::string> > sample_data;
    //reduce to chr19
    if(DEBUG) std::cerr<<"File = "<<reads_file<<std::endl;
    if (myfile.is_open()){
        while( getline( myfile, quest ) )
            {
           // if((!quest.empty())&&(quest.find("chr19")!=std::string::npos))
            if((!quest.empty()))
                {
                   text_file.push_back( quest );
                   line_counter++;
                }
            }
        if(DEBUG) std::cerr<<"line counter == "<<line_counter<<std::endl;
        for(int i = 0; i<nb_rand; i++) {
                std::vector<std::string> text_out;
                text_out.clear();
                std::vector<std::string> text_file_cp;
                 text_file_cp.clear();
                 text_file_cp= text_file;
                int init=0;
                while(init<min_map)
                    {
                        int line = rand() % line_counter;
                        text_out.push_back(text_file_cp[line]);
                        text_file_cp.erase(text_file_cp.begin() + line);
                        init++;
                    }
                sample_data.push_back(text_out);
             }
         myfile.close();
    } else {exit(0);}

  int select_file = rand()%nb_rand;
  for (unsigned l=0;l<sample_data[select_file].size();l++){
     out<<sample_data[select_file][l]<<std::endl;
  }

}
