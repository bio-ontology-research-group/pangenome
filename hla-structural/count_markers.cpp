// Exact canonical 31-mer counting. Output is little-endian uint32 in marker order.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;
const uint64_t MASK=(1ULL<<62)-1;
int base(char c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}}
template<class F> void scan(const string&s,const string*q,F hit){
 uint64_t f=0,r=0;int valid=0;
 for(size_t i=0;i<s.size();i++){
  int b=base(s[i]);if(b<0||(q&&((*q)[i]-33<20))){f=r=0;valid=0;continue;}
  f=((f<<2)|b)&MASK;r=(r>>2)|(uint64_t(3-b)<<60);valid++;
  if(valid>=31)hit(min(f,r));
 }
}
int main(int argc,char**argv){
 if(argc!=5){cerr<<"mode(fasta|fastq) markers input(- permitted) output\n";return 2;}
 string mode=argv[1];if(mode!="fasta"&&mode!="fastq")return 2;
 ifstream mf(argv[2]);if(!mf){cerr<<"markers unavailable\n";return 2;}
 unordered_map<uint64_t,size_t> idx;uint64_t key;while(mf>>key)idx.emplace(key,idx.size());
 vector<uint32_t> counts(idx.size(),0);ifstream inf;istream*in=&cin;
 if(string(argv[3])!="-"){inf.open(argv[3]);if(!inf)return 2;in=&inf;}
 string line,seq,name,plus,qual,prev;unordered_set<size_t> fragment;uint64_t records=0,fragments=0;
 auto flush=[&](){if(!prev.empty()){for(auto k:fragment)counts[k]++;fragments++;fragment.clear();}};
 if(mode=="fasta"){
  auto process=[&](){scan(seq,nullptr,[&](uint64_t k){auto it=idx.find(k);if(it!=idx.end())counts[it->second]++;});seq.clear();};
  while(getline(*in,line)){if(!line.empty()&&line[0]=='>'){process();records++;}else seq+=line;}process();
 }else{
  while(getline(*in,name)){
   if(!getline(*in,seq)||!getline(*in,plus)||!getline(*in,qual)||seq.size()!=qual.size()||name.empty()||name[0]!='@'){cerr<<"invalid FASTQ\n";return 3;}
   name=name.substr(1,name.find_first_of(" \t")-1);
   if(name.size()>2&&name[name.size()-2]=='/'&&(name.back()=='1'||name.back()=='2'))name.resize(name.size()-2);
   if(name!=prev){flush();prev=name;}
   scan(seq,&qual,[&](uint64_t k){auto it=idx.find(k);if(it!=idx.end())fragment.insert(it->second);});records++;
  }flush();
 }
 ofstream out(argv[4],ios::binary);if(!out)return 2;
 for(uint32_t c:counts){for(int j=0;j<4;j++)out.put(char((c>>(8*j))&255));}
 cerr<<"records="<<records<<" fragments="<<fragments<<" markers="<<counts.size()<<"\n";
 return out?0:2;
}
