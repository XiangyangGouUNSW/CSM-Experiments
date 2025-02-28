#include<fstream>
#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
using namespace std;

int main(int argc, char* argv[]) 
{
	string file_path1 = argv[1];
	string file_path2 = argv[2];
//	string file_path2 = "./record.txt";
//	string file_path1 = "./NSP-record.txt";
	ifstream fin1(file_path1.c_str());
	ifstream fin2(file_path2.c_str());
	string query;
	unsigned int qid;
	double index_time, processing_time1, processing_time2, inter_num;
	double impr = 0;
	unsigned long long flag, total_result, cand_v, cand_e, memory;
	unsigned int qnum = 0;
	while(fin1>>query)
	{
		fin1>>qid>>index_time>>flag>>total_result>>processing_time1>>inter_num;
		fin1>>cand_v>>cand_e>>memory;
		fin2>>query>>qid>>index_time>>flag>>total_result>>processing_time2>>inter_num;
		fin2>>cand_v>>cand_e>>memory;

		if(processing_time1>=3600&&processing_time2>=3600)
		    continue;
		qnum++;
		cout<<"processing query "<<qid<<' '<<processing_time1<<' '<<processing_time2<<endl;
		impr += 1-(processing_time1/processing_time2);
	}
	cout<<qnum<<" average improvement "<<impr/qnum<<endl;
	return 0;
}
