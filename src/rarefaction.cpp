#include "../include/rarefaction.h"

//************************randomly pick n lines from a file ******************************************//
 std::vector<std::string> pick_n_random_lines( const std::string& path, size_t n, int rand_cut ) {
//****************************************************************************************************//
   static std::mt19937 rng( time(0) ) ;
   static const auto random = []( size_t min, size_t max ) {
     return std::uniform_int_distribution<size_t>(min,max)(rng) ;
   } ;
   std::string     rare = to_string(rand_cut);
   std::vector<std::string> selected ;
   omp_set_num_threads(N_P);
   std::ifstream file(path) ;
   std::string line ;
   size_t nlines = 0 ;
   while( getline( file, line ) ) {
     if( !line.empty() ) {
       if( selected.size() < n ) selected.push_back(line) ; // select the first n
       else // replace a random selected line with probability n / #lines
	 if( random(0,nlines) < n ) selected[ random(0,n-1) ] = line ;
       ++nlines ;
     }
   }
   shuffle( selected.begin(), selected.end(), rng ) ;
   // __gnu_parallel::random_shuffle( selected.begin(), selected.end(), rng ) ;
   std::string outfile = "sample_"+rare+"_"+path;
   std::ofstream out(outfile.c_str());
   for (unsigned i=0;i<selected.size();i++) out<<selected[i]<<std::endl;
   return selected ;
 }

//**************************************************
double median(std::vector<double> vec) {
//**************************************************
        typedef std::vector<int>::size_type vec_sz;

        vec_sz size = vec.size();
        if (size == 0)
                throw std::domain_error("median of an empty vector");

        sort(vec.begin(), vec.end());

        vec_sz mid = size/2;

        return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];
}
//**************************************************
void stat_rarefaction(std::string file){
//***************************************************
std::vector<std::string> new_file;
boost::split(new_file,file,boost::is_any_of("/"));
std::string out = "stat_file_"+new_file[2];
std::ifstream in(file.c_str());
std::ofstream ou(out.c_str());
std::vector <std::vector <std::string>> rare;
std::string line;

int i=0,j;
while(getline(in,line) && !line.empty()) {
   std::vector <std::string> temp;
   split(temp,line,boost::is_space());
   std::vector <std::string> temp2;
      for (size_t tt=0;tt<temp.size();tt++){
        temp2.push_back(temp[tt]);
      }
     rare.push_back(temp2);
    }
int nb_rare = ((rare[0].size()-7)/3)+1;
for (unsigned h=0;h<rare.size();h++){
        std::vector < double> treat;
        std::vector <double> pval;
        std::vector <double> qval;
        std::string ts1 = rare[h][3];
        treat.push_back(atof(ts1.c_str()));
        std::string tp1 = rare[h][5];
        pval.push_back(atof(tp1.c_str()));
        std::string tq1 = rare[h][6];
        qval.push_back(atof(tq1.c_str()));
        for (unsigned k=7;k<rare[h].size();k+=3){
            std::string ts = rare[h][k];
            treat.push_back(atof(ts.c_str()));
            std::string tp = rare[h][k+1];
            pval.push_back(atof(tp.c_str()));
            std::string tq = rare[h][k+2];
            qval.push_back(atof(tq.c_str()));
        }
   double min_t =*min_element(treat.begin(),treat.end());
   double max_t =*max_element(treat.begin(),treat.end());
   double med_t = median(treat);
   double min_p =*min_element(pval.begin(),pval.end());
   double max_p =*max_element(pval.begin(),pval.end());
   double med_p = median(pval);
   double min_q =*min_element(qval.begin(),qval.end());
   double max_q =*max_element(qval.begin(),qval.end());
   double med_q = median(qval);
   ou<<rare[h][0]<<"\t"<<rare[h][1]<<"\t"<<rare[h][2]<<"\t"<<rare[h][4]<<"\t";
   ou<<min_t<<"\t"<<max_t<<"\t"<< med_t<<"\t";
   ou<<min_p<<"\t"<<max_p<<"\t"<< med_p<<"\t";
   ou<<min_q<<"\t"<<max_q<<"\t"<< med_q<<"\n";
  }
 }

//******** count line in file and return size_t *******************
size_t count_line_in_file(std::string file){
//-----------------------------------------------------------------
  std::ifstream in(file.c_str());
  size_t line_count=0;
  std::string line;
  while(getline(in,line)) if((!line.empty())) ++line_count;
  return line_count;
}

