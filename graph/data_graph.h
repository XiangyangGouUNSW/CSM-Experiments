#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include<map>
#include "../utility/struct.h"
#include<assert.h>
#define merge_long_long(s, d) (((long long)s<<32)|d)
using namespace std;

class data_graph
{
private:
	unordered_map<unsigned long long, vector<int>> edge_label_map;

public:
	vector<node> g; // this vector record the information and neighborhood of vertices in the graph. Note that we use the indexes of vertices in this vector as its
	// incoded-id in the algorithms, so that we can transform original vertice ids to continuous integer which supports fast retrieve.   
	unordered_map<unsigned int, unsigned int> vertice_id_map; // reverse mapping from the original vertice id to the encoded id/
	vector<unsigned int> original_ids;
	bool original_continuous_id;    // Note that if the node ids in the original dataset are continuous integers, the trasform is not needed.This variable indicates this case.
	unordered_map<unsigned int, vector<unsigned int>> label_vertice_map; // map each vertice label to the set of vertices with this label. 
	//unordered_map<unsigned long long, vector<int>> edge_label_map; // this map enables a fast check for the existence and labels of an edge. It maps a vertex-ID pair (merged into a long long type, the first 32 bits represent the source and the last 32bit for the destination) into
	// a set of labels attached to this edge//
	unsigned int edge_number; // the number of edges, increased/decreased as edges inserted/deleted
	unsigned int vertice_number;//the number of vertices. If the vertice ids in the original dataset are continuous integers, the total number of vertices should be given when initiate the data graph, 
	vector<map<unsigned long long, int>> NLF;
	data_graph(unsigned int vertice_number_=0) {
		edge_number = 0;
		vertice_number = vertice_number_;
		g.resize(vertice_number);
		NLF.resize(vertice_number);
		if (vertice_number != 0) 
			original_continuous_id = true;
		else
			original_continuous_id = false;
	}
	~data_graph()
	{
		for (int i = 0; i < g.size(); i++)
			g[i].clear();
		g.clear();
		vertice_id_map.clear();
		for (auto iter = label_vertice_map.begin(); iter != label_vertice_map.end(); iter++)
			iter->second.clear();
		label_vertice_map.clear();
		for (auto iter = edge_label_map.begin(); iter != edge_label_map.end(); iter++)
			iter->second.clear();
		edge_label_map.clear();
		for (int i = 0; i < NLF.size(); i++)
			NLF[i].clear();
		NLF.clear();
	}
	unsigned int memory_compute()
	{
		unsigned int mem1 = sizeof(g) + g.capacity() * sizeof(node);
		for (int i = 0; i < g.size(); i++)
		{
			for (auto iter = g[i].label_edge_map.begin(); iter != g[i].label_edge_map.end(); iter++)
				mem1 += sizeof(unsigned long long) + sizeof(vector<int>) + sizeof(void*) *2 + iter->second.capacity() * sizeof(int);
		}
		cout << "memory of graph " << ((double)mem1) / 1024 / 1024 << " MB" << endl;

		unsigned int mem2 = sizeof(label_vertice_map);
		mem2 += label_vertice_map.bucket_count() * sizeof(void*);
		mem2 += label_vertice_map.size() * (sizeof(unsigned int) + sizeof(vector<int>) + sizeof(void*));
		for (auto iter = label_vertice_map.begin(); iter != label_vertice_map.end(); iter++)
			mem2 += iter->second.capacity() * sizeof(int);

		cout << "memory of label vertice map " << ((double)mem2) / 1024 / 1024 << " MB" << endl;

		unsigned int mem3 = sizeof(edge_label_map);
		mem3 += edge_label_map.bucket_count() * sizeof(void*);
		mem3 += edge_label_map.size() * (sizeof(unsigned long long) + sizeof(vector<int>) + sizeof(void*));
		for (auto iter = edge_label_map.begin(); iter != edge_label_map.end(); iter++)
			mem3 += iter->second.capacity() * sizeof(int);

		cout << "memory of edge label map " << ((double)mem3) / 1024 / 1024 << " MB" << endl;

		unsigned int mem4 = sizeof(NLF);
		mem4 += NLF.capacity() * sizeof(map<unsigned long long, int>);
		for (int i = 0; i < NLF.size(); i++)
			mem4 += NLF[i].size() * (sizeof(unsigned long long) + sizeof(unsigned int) + sizeof(void*));

		cout << "memory of NLF " << ((double)mem4) / 1024 / 1024 << " MB" << endl;
		return mem1 + mem2 + mem3 + mem4;
	}
	void set_max_vertice_number(unsigned int vertice_number_)
	{
		vertice_number = vertice_number_;
		original_continuous_id = true;
		g.resize(vertice_number);
	}

	bool insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{

		if (!original_continuous_id)
		{
			if (vertice_id_map.find(src) == vertice_id_map.end())
			{
				vertice_id_map[src] = g.size();
				original_ids.push_back(src);
				vertice_number++;
				src = g.size();
				g.push_back(node(src, src_label));
				NLF.push_back(map<unsigned long long, int>());
			}
			else
				src = vertice_id_map[src];

			if (vertice_id_map.find(dst) == vertice_id_map.end())
			{
				vertice_id_map[dst] = g.size();
				original_ids.push_back(dst);
				vertice_number++;
				dst = g.size();
				g.push_back(node(dst, dst_label));
				NLF.push_back(map<unsigned long long, int>());
			}
			else
				dst = vertice_id_map[dst];
		}
		else
		{
			g[src].id = src;
			g[src].label = src_label;
			g[dst].id = dst;
			g[dst].label = dst_label;
		}

		if (label_vertice_map.find(src_label) == label_vertice_map.end())
			label_vertice_map[src_label] = vector<unsigned int>();
		if (label_vertice_map.find(src_label) == label_vertice_map.end())
			label_vertice_map[src_label] = vector<unsigned int>();

		unsigned long long label_pair = merge_long_long(dst_label, edge_label);
		unsigned long long reverse_label_pair = merge_long_long(src_label, edge_label);


		bool exist = false;
		vector<unsigned int>::iterator iter = lower_bound(g[src].label_edge_map[label_pair].begin(), g[src].label_edge_map[label_pair].end(), dst);
		if (iter == g[src].label_edge_map[label_pair].end() || *iter != dst)
			g[src].label_edge_map[label_pair].emplace(iter, dst);
		else
			exist = true;

		iter = lower_bound(g[dst].label_edge_map[reverse_label_pair].begin(), g[dst].label_edge_map[reverse_label_pair].end(), src);
		if (iter == g[dst].label_edge_map[reverse_label_pair].end() || *iter != src)
			g[dst].label_edge_map[reverse_label_pair].emplace(iter, src);
		else
			exist = true;
		if (!exist) {
			edge_number++;
			if (NLF[src].find(label_pair) != NLF[src].end())
				NLF[src][label_pair]++;
			else
				NLF[src][label_pair] = 1;
			if (NLF[dst].find(reverse_label_pair) != NLF[dst].end())
				NLF[dst][reverse_label_pair]++;
			else
				NLF[dst][reverse_label_pair] = 1;

			iter = lower_bound(label_vertice_map[src_label].begin(), label_vertice_map[src_label].end(), src);
			if (iter == label_vertice_map[src_label].end() || *iter != src)
				label_vertice_map[src_label].emplace(iter, src);

			iter = lower_bound(label_vertice_map[dst_label].begin(), label_vertice_map[dst_label].end(), dst);
			if (iter == label_vertice_map[dst_label].end() || *iter != dst)
				label_vertice_map[dst_label].emplace(iter, dst);
			unsigned int v1, v2;
			if (src >= dst)
			{
				v1 = src;
				v2 = dst;
			}
			else
			{
				v2 = src;
				v1 = dst;
			}
			unsigned long long edge_pair = merge_long_long(v1, v2);
			vector<int>::iterator iter2 = lower_bound(edge_label_map[edge_pair].begin(), edge_label_map[edge_pair].end(), edge_label);
			edge_label_map[edge_pair].emplace(iter2, edge_label);

		/*	edge_pair = merge_long_long(dst, src);
			iter2 = lower_bound(edge_label_map[edge_pair].begin(), edge_label_map[edge_pair].end(), edge_label);
			edge_label_map[edge_pair].emplace(iter2, edge_label);*/
		}
		return !exist;
	}	

	void delete_edge(unsigned int src, unsigned int dst, int edge_label)
	{
		if (!original_continuous_id)
		{
			src = vertice_id_map[src];
			dst = vertice_id_map[dst];
		}
		unsigned int src_label = g[src].label;
		unsigned int dst_label = g[dst].label;
		unsigned long long label_pair = merge_long_long(dst_label, edge_label);
		unsigned long long reverse_label_pair = merge_long_long(src_label, edge_label);

		vector<unsigned int>::iterator iter = lower_bound(g[src].label_edge_map[label_pair].begin(), g[src].label_edge_map[label_pair].end(), dst);
		if (iter != g[src].label_edge_map[label_pair].end() && *iter == dst) {
			g[src].label_edge_map[label_pair].erase(iter);
			if (g[src].label_edge_map[label_pair].empty())
			{
				g[src].label_edge_map.erase(label_pair);
				if (g[src].label_edge_map.empty()) // if there are no edges associated with this vertex, it is an isolated vertex and can be viewed as unvalid, but we will leave its position in g and NLF, to avoid re-code other vertices
				{
					vector<unsigned int>::iterator search_iter = lower_bound(label_vertice_map[src_label].begin(), label_vertice_map[src_label].end(), src);
					assert(search_iter != label_vertice_map[src_label].end() && *search_iter == src);
					label_vertice_map[src_label].erase(search_iter);
					if (label_vertice_map[src_label].empty())
						label_vertice_map.erase(src_label);
				}
			}

			NLF[src][label_pair]--;
			if (NLF[src][label_pair] == 0)
				NLF[src].erase(label_pair);

		}


		iter = lower_bound(g[dst].label_edge_map[reverse_label_pair].begin(), g[dst].label_edge_map[reverse_label_pair].end(), src);
		if (iter != g[dst].label_edge_map[reverse_label_pair].end() && *iter == src) {
			g[dst].label_edge_map[reverse_label_pair].erase(iter);
			if (g[dst].label_edge_map[reverse_label_pair].empty())
			{
				g[dst].label_edge_map.erase(reverse_label_pair);
				if (g[dst].label_edge_map.empty())
				{
					vector<unsigned int>::iterator search_iter = lower_bound(label_vertice_map[dst_label].begin(), label_vertice_map[dst_label].end(), dst);
					assert(search_iter != label_vertice_map[dst_label].end() && *search_iter == dst);
					label_vertice_map[dst_label].erase(search_iter);
					if (label_vertice_map[dst_label].empty())
						label_vertice_map.erase(dst_label);
				}
			}

			NLF[dst][reverse_label_pair]--;
			if (NLF[dst][reverse_label_pair] == 0)
				NLF[dst].erase(reverse_label_pair);
		}

		unsigned int v1, v2;
		if (src >= dst)
		{
			v1 = src;
			v2 = dst;
		}
		else
		{
			v2 = src;
			v1 = dst;
		}

		unsigned long long edge_pair = merge_long_long(v1, v2);
		if (edge_label_map.find(edge_pair) != edge_label_map.end())
		{
			vector<int>::iterator iter2 = lower_bound(edge_label_map[edge_pair].begin(), edge_label_map[edge_pair].end(), edge_label);
			if (iter2 != edge_label_map[edge_pair].end() && *iter2 == edge_label)
				edge_label_map[edge_pair].erase(iter2);
			if (edge_label_map[edge_pair].empty())
				edge_label_map.erase(edge_pair);

		}
	/*	edge_pair = merge_long_long(dst, src);
		if (edge_label_map.find(edge_pair) != edge_label_map.end())
		{
			vector<int>::iterator iter2 = lower_bound(edge_label_map[edge_pair].begin(), edge_label_map[edge_pair].end(), edge_label);
			if (iter2 != edge_label_map[edge_pair].end() && *iter2 == edge_label)
				edge_label_map[edge_pair].erase(iter2);
			if (edge_label_map[edge_pair].empty())
				edge_label_map.erase(edge_pair);

		}*/
	}

	void count_label_frequency(unsigned int id, map<unsigned long long, int>& label_freq)
	{
		label_freq = NLF[id];
	}

	void get_labeled_node(int label, vector<unsigned int>& vec)
	{
		if (label_vertice_map.find(label) != label_vertice_map.end())
			vec = label_vertice_map[label];
	}

	bool check_edge_label(unsigned int src, unsigned int dst, int edge_label)
	{
		unsigned int v1, v2;
		if (src >= dst)
		{
			v1 = src;
			v2 = dst;
		}
		else
		{
			v2 = src;
			v1 = dst;
		}
		unsigned long long edge_pair = merge_long_long(v1, v2);
		if (edge_label_map.find(edge_pair) != edge_label_map.end())
		{
			vector<int>::iterator iter = lower_bound(edge_label_map[edge_pair].begin(), edge_label_map[edge_pair].end(), edge_label);
			if (iter != edge_label_map[edge_pair].end() && *iter == edge_label)
				return true;
		}
		return false;
	}

	void get_neighbor(unsigned int id, vector<neighbor_info>& neighbor_list)
	{
		for (map<unsigned long long, vector<unsigned int> >::iterator iter = g[id].label_edge_map.begin(); iter != g[id].label_edge_map.end(); iter++)
		{
			unsigned int dst_label = (iter->first >> 32);
			unsigned int edge_label = (iter->first & 0xFFFFFFFF);
			for (int i = 0; i < iter->second.size(); i++)
				neighbor_list.push_back(neighbor_info(iter->second[i], dst_label, edge_label));
		}
	}

	void get_labeled_neighbor(unsigned int id, int target_neighbor_label, int target_edge_label, vector<unsigned int>& neighbor_list)
	{
		unsigned long long label_pair = merge_long_long(target_neighbor_label, target_edge_label);
		if (g[id].label_edge_map.find(label_pair) != g[id].label_edge_map.end())
			neighbor_list = g[id].label_edge_map[label_pair];
	}

	bool check_labeled_neighbor(unsigned int id, int target_neighbor_label, int target_edge_label)
	{
		unsigned long long label_pair = merge_long_long(target_neighbor_label, target_edge_label);
		if (g[id].label_edge_map.find(label_pair) != g[id].label_edge_map.end())
			return true;
		else
			return false;
	}

	unsigned int get_labeled_degree(unsigned int v, unsigned long long label_pair)
	{
		if (NLF[v].find(label_pair) == NLF[v].end())
			return 0;
		else
			return NLF[v][label_pair];
	}
	unsigned int get_labeled_degree(unsigned int v, unsigned int neighbor_label, unsigned int edge_label)
	{
		unsigned long long label_pair = merge_long_long(neighbor_label, edge_label);
		if (NLF[v].find(label_pair) == NLF[v].end())
			return 0;
		else
			return NLF[v][label_pair];
	}
};
