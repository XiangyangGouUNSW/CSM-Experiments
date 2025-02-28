#include<iostream>
#include<fstream>
#include<vector>
#include<unordered_map>
#include<time.h>
#include"setting.h"
using namespace std;


bool read_query_graph(query_graph* qg, ifstream& fquery)
{
	unsigned int vertice_num;
	if (fquery >> vertice_num) {
		qg->clear();
		qg->resize(vertice_num);
		for (int i = 0; i < vertice_num; i++)
		{
			unsigned int src, src_label, label_num, edge_num, dst, dst_label, edge_label;
			fquery >> src >> src_label >> label_num;
			for (int j = 0; j < label_num; j++)
			{
				fquery >> dst_label >> edge_label >> edge_num;
				for (int k = 0; k < edge_num; k++)
				{
					fquery >> dst;
					if (dst < src) // to avoid inserting the same edge twice.
						qg->insert_edge(src, src_label, dst, dst_label, edge_label);
				}
			}
		}
		return true;
	}
	else
		return false;
}


int main(int argc, char* argv[])
{
	string initial_path = argv[1];
	string update_path = argv[2];
	string query_path = argv[3];
	string output_path = argv[4];


/*	string initial_path =  "D://data//initial-lkml.txt";
	string update_path = "D://data//updated-lkml.txt";
	string query_path = "D://data//sparse-query-lkml.txt";
	string output_path = "./NSP-";*/

	bool record_match = false; // this parameter records if we need to report all match results find in each update. If it is false, we will only report the result count. 
	//Otherwise we will store the these results in a table. Note that report the result will slow down the algortihm and greatly narrow the gap between different algorithms, as the cost of reporting results may
	//become bottleneck and this cost is the same for all methods. We set this parameter to false in our experiments,  which is the same the experimental setting of previous works (we have cheked their open source code).


	ifstream fdata(initial_path.c_str());
	ifstream fupdate(update_path.c_str());
	ifstream fquery(query_path.c_str());

	// here we assume that the number of vertices are known and ids are consecutive inegers start from 0. If they are not, class data_graph
	// will transform them into so.
	
	string outpath1 = output_path + "result.txt";
	string outpath2 = output_path + "record.txt";

	
	ofstream fout(outpath1.c_str());
	ofstream frec(outpath2.c_str());
	
	query_graph* qg = new query_graph;
	unsigned int query_cursor = 0;
	unsigned int target_query = 2;
	while (read_query_graph(qg, fquery))
	{
		if (query_cursor > target_query)
			break;
		cout << "processing query NO" << query_cursor << endl; 
		fout << "query " << query_cursor << endl;
		frec << "query " << query_cursor << endl;
		query_cursor++;

		unsigned int vertice_number;
		fdata >> vertice_number;

		unsigned int src, src_label, dst, dst_label, edge_label;
		data_graph* g = new data_graph(vertice_number);
		unsigned int cnt = 0;
		while (fdata >> src >> src_label >> dst >> dst_label >> edge_label) {
			cnt++;
			if (cnt % 100000 == 0)
				cout << cnt << " edges has been inserted" << endl;
			g->insert_edge(src, src_label, dst, dst_label, edge_label);
		}

		cout << "load initial data finished " << endl;
		solution* su = new solution(g, qg, enable_global_index, enable_edge_view);
		double initialize_time = su->initialization();
		cout << "initialization time " << initialize_time << endl;
		if (enable_global_index) {
			pair<unsigned int, unsigned int> p = su->idx->get_size();
			cout <<"The candidate cnt, candidate edge cnt, and memory usage of the initial index are "<< p.first << ' ' << p.second << ' ' << su->idx->memory_compute() << endl;
		}
		frec << initialize_time << ' ';
		

		bool solved = true;
		unsigned int total_match = 0;

		cnt = 0;

		start_clock = clock();
		inter_result =0;
		filtered_edge = 0;
		label_filtered = 0;
		neighbor_sum = 0;
		expanding_num = 0;
		filtered_expanding = 0; 
		while (fupdate >> src >> src_label >> dst >> dst_label >> edge_label)
		{
			vector<unordered_map<unsigned int, match_result>> match;
			unsigned int match_size = 0;
			su->insert_edge(src, src_label, dst, dst_label, edge_label, match, match_size, record_match);
			if(match_size!=0)
				fout <<cnt<<' '<< match_size << endl;
			total_match += match_size;
			clock_t current = clock();
			double time_passd = (double)(current - start_clock) / CLOCKS_PER_SEC;
			if (time_passd >= time_limit)
			{
				solved = false;
				break;
			}
			cnt++;
			if (cnt % 10000 == 0)
				cout << cnt << " edges processed " << endl;
		}
		clock_t current = clock();
		double time_passd = (double)(current - start_clock) / CLOCKS_PER_SEC;
		frec << solved << ' ' << ' '<<total_match<<' '<<time_passd <<' '<<((double)inter_result/(cnt- label_filtered))<<endl;
		
		if(enable_global_index){
		pair<unsigned int, unsigned int> p = su->idx->get_size();
		frec <<p.first<<' '<<p.second<<' '<<su->idx->memory_compute()<<endl;
		}
		else
			frec<<0<<' '<<0<<' '<<0<<endl;
		if (solved){
			cout << "solved query, eclipsed time " << time_passd <<endl;
		}
		else
			cout << "unsolved query " << endl;
		fout << -1 << endl;
		fout << endl;

		fdata.clear();
		fdata.seekg(0, ios::beg);
		fupdate.clear();
		fupdate.seekg(0, ios::beg);
		delete g;
		delete su;
	}

}
