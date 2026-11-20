#include <zlib.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
struct Agg { long long rows=0, positive_rows=0, depth_sum=0; };
int main(int argc,char**argv){
  if(argc!=6){std::cerr<<"usage: scan_coverage_positions coverage.gz barcode_map.tsv positive_positions.tsv matched_barcodes.tsv human_suffix_match(0/1)\n";return 2;}
  const bool suffix_mode=std::atoi(argv[5])!=0;
  std::unordered_map<std::string,int> types,suffix_types; std::string line;
  std::ifstream bm(argv[2]); if(!bm)throw std::runtime_error("cannot open barcode map");
  while(std::getline(bm,line)){
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    auto t=line.find('\t'); if(t==std::string::npos)continue;
    std::string bc=line.substr(0,t); int idx=std::stoi(line.substr(t+1));
    types[bc]=idx;
    if(suffix_mode){
      auto u=bc.rfind('_'); std::string suffix=(u==std::string::npos?bc:bc.substr(u+1));
      auto ins=suffix_types.emplace(suffix,idx);
      if(!ins.second&&ins.first->second!=idx)throw std::runtime_error("ambiguous metadata barcode suffix");
    }
  }
  std::unordered_set<std::string> matched;
  std::unordered_map<int,Agg> pos;
  gzFile gz=gzopen(argv[1],"rb"); if(!gz)throw std::runtime_error("cannot open coverage file");
  char buf[1024]; long long lines=0,selected=0; int err=0;
  while(gzgets(gz,buf,sizeof(buf))){
    ++lines;
    if(!std::strchr(buf,'\n')&&!gzeof(gz))throw std::runtime_error("coverage row exceeds buffer");
    char*a=std::strchr(buf,','); if(!a)throw std::runtime_error("missing first comma"); *a=0;
    char*b=std::strchr(a+1,','); if(!b)throw std::runtime_error("missing second comma"); *b=0;
    std::string bc(a+1); auto it=types.find(bc); int idx=-1;
    if(it!=types.end())idx=it->second;
    else if(suffix_mode){auto u=bc.rfind('_');std::string suffix=(u==std::string::npos?bc:bc.substr(u+1));auto si=suffix_types.find(suffix);if(si!=suffix_types.end())idx=si->second;}
    if(idx<0)continue;
    int p=std::atoi(buf), d=std::atoi(b+1); if(d<0)throw std::runtime_error("negative depth");
    matched.insert(bc); ++selected;
    auto &x=pos[p]; ++x.rows; x.depth_sum+=d; if(d>0)++x.positive_rows;
  }
  const char*msg=gzerror(gz,&err); if(err!=Z_OK&&err!=Z_STREAM_END)throw std::runtime_error(msg); gzclose(gz);
  std::ofstream out(argv[3]); if(!out)throw std::runtime_error("cannot write positions");
  out<<"position\tn_coverage_records\tn_positive_coverage_records\ttotal_depth\n";
  long long npos=0; for(const auto&kv:pos)if(kv.second.positive_rows>0){++npos;out<<kv.first<<'\t'<<kv.second.rows<<'\t'<<kv.second.positive_rows<<'\t'<<kv.second.depth_sum<<'\n';}
  std::ofstream mb(argv[4]); if(!mb)throw std::runtime_error("cannot write matched barcodes");
  mb<<"coverage_barcode\tcelltype_index\n";
  for(const auto&bc:matched){auto it=types.find(bc);int idx=-1;if(it!=types.end())idx=it->second;else if(suffix_mode){auto u=bc.rfind('_');std::string suffix=(u==std::string::npos?bc:bc.substr(u+1));idx=suffix_types.at(suffix);}mb<<bc<<'\t'<<idx<<'\n';}
  std::cerr<<"lines="<<lines<<" selected_coverage_records="<<selected<<" matched_barcodes="<<matched.size()<<" positions_with_positive_depth="<<npos<<"\n";
}
