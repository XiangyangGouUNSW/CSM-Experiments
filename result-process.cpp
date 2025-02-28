#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
using namespace std;

int main(int argc, char* argv[]) 
{
	//string file_path = "./record.txt";
	string file_path = argv[1];
	ifstream fin(file_path.c_str());
	string query;
	unsigned int qid;
	double index_time, processing_time, inter_num;
	double total_it=0, total_pt=0, total_inter=0, total_r = 0, total_cv = 0, total_ce = 0, total_m=0;
	unsigned long long flag, total_result, cand_v, cand_e, memory;
	unsigned int failed = 0, success =0;
	while(fin>>query)
	{
		fin>>qid>>index_time>>flag>>total_result>>processing_time>>inter_num;
		fin>>cand_v>>cand_e>>memory;
		
		cout<<"processing query "<<qid<<' '<<flag<<' '<<processing_time<<endl;
		total_it += index_time;
		total_cv += cand_v;
		total_ce += cand_e;
		total_m += memory;
		
		if(flag==0){
			total_pt += 3600;
			failed++;
		}
		else
		{
			total_pt += processing_time;
			total_inter += inter_num;
			total_r += total_result;
			success++;
		}
	}
	
	cout<<"failed query "<<failed<<endl;
	cout<<"success query "<<success<<endl;
	unsigned int total_q = success + failed;
	cout<<"average indexing time "<<total_it /total_q<<endl;
	cout<<"average processing time "<<total_pt /total_q<<endl;
	cout<<"average result number "<<total_r / success<<endl;
	cout<<"averge intermediate result size "<<total_inter / success<<endl;
	cout<<"average candidate vertices number "<<total_cv / total_q<<endl;
	cout<<"average candidate edges number "<<total_ce / total_q<<endl;
	cout<<"average memory usage "<<total_m/total_q /1024/1024<<" MB"<<endl;

	return 0;
}
