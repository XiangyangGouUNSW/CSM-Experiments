#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include<map>
#include "../graph/data_graph.h"
#include "../graph/query_graph.h"
#include"../indexes/index-list.h"
#include"../setting.h"
#include<assert.h>
#include<algorithm>
#include"../utility/automorphism-check.h"
#include "../utility/OrderGeneration.h"
#include"../indexes/local-index.h"
using namespace std;

class RapidFlow
{
public:
	data_graph* g;
	query_graph* q;
	my_index* idx;
	unordered_map<unsigned long long, matching_order*> mo_map;
	unordered_map<unsigned long long, indexing_order*> io_map;
	map<unsigned long long, vector<pair<unsigned long long, automatch*>>> auto_set;
	vector<automatch> automatches;
	bool use_global_index;
	RapidFlow(data_graph* dg, query_graph* qg, bool use_gd = false, bool use_edge_view = false)
	{
		g = dg;
		q = qg;
		if(use_gd)
			idx = new my_index(g, q, use_edge_view);
		use_global_index = use_gd;
	}
	~RapidFlow()
	{
		if(use_global_index)
			delete idx;
		for (auto iter = mo_map.begin(); iter != mo_map.end(); iter++)
			delete iter->second;
		mo_map.clear();
		automatches.clear();
		for (auto iter = auto_set.begin(); iter != auto_set.end(); iter++)
			iter->second.clear();
		auto_set.clear();
	}

	void print_auto_match(automatch* aut)
	{
		for (auto iter = aut->mapping.begin(); iter != aut->mapping.end(); iter++)
		{
			cout << "match from " << iter->first << " to " << iter->second << endl;
		}
	}
	void print_autoset()
	{
		for (auto iter = auto_set.begin(); iter != auto_set.end(); iter++)
		{
			unsigned int src = (iter->first >> 32);
			unsigned int dst = (iter->first & 0xFFFFFFFF);
			cout << "matching list of " << src << ' ' << dst << endl;
			for (int i = 0; i < iter->second.size(); i++)
			{
				unsigned int ms = (iter->second[i].first >> 32);
				unsigned int md = (iter->second[i].first & 0xFFFFFFFF);
				cout << "edge " << ms << ' ' << md <<" match detail"<< endl;
				print_auto_match(iter->second[i].second);
				cout << endl;
			}
		}
	}


	double initialization()
	{
		for (int i = 0; i < q->g.size(); i++)
		{
			unsigned int src = q->g[i].id;
			vector<unsigned int> neighbor;
			q->get_unlabeled_neighbor(src, neighbor);
			for (int j = 0; j < neighbor.size(); j++)
			{
				src = q->g[i].id;
				unsigned int dst = neighbor[j];
				if (src > dst)
				{
					unsigned int tmp = src;
					src = dst;
					dst = tmp;
				}
				unsigned long long edge_pair = merge_long_long(src, dst);
				if (mo_map.find(edge_pair) == mo_map.end()) {
					query_graph* tmp_q = new query_graph(q->node_number);
					q->copy_with_node_delete(src, dst, tmp_q);
					matching_order* mo = new matching_order;
					RI_order_generation(tmp_q, mo);
					mo_map[edge_pair] = mo;
					delete tmp_q;
				}
				if (io_map.find(edge_pair) == io_map.end())
				{
					indexing_order* io = new indexing_order;
					RF_indexing_order_generation(q, io, src, dst);
					io_map[edge_pair] = io;
				}
			}
		}
		unordered_map<unsigned int, vector<unsigned int>> match;
		automorphism(q, match);
		generate_autoset(q, match, auto_set, automatches);
		for (auto iter = match.begin(); iter != match.end(); iter++)
			iter->second.clear();
		match.clear();
		//print_autoset();

		double time_used = 0;
		if (use_global_index) {
			clock_t start = clock();
			idx->construct_index();
			clock_t end = clock();
			time_used = double(end - start) / CLOCKS_PER_SEC;
		}
		return time_used;
		
	}

	void result_compare(unordered_map<unsigned int, unsigned int>& tmp_map)
	{
		for (auto iter = tmp_map.begin(); iter != tmp_map.end(); iter++)
		{
			unsigned int u = iter->first;
			unsigned int v = iter->second;
			vector<neighbor_info> vec;
			q->get_neighbor(u, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int u1 = vec[i].id;
				unsigned int label = vec[i].e_label;
				unsigned int v1 = tmp_map[u1];
				if (g->g[v1].label != q->g[u1].label)
					cout << "error map" << endl;
				if (!g->check_edge_label(v, v1, label)) {
					cout << "error map " << endl;
				}
			}
		}
	}
	void result_check(unordered_map<unsigned int, match_result>& tmp_match)
	{
		for (int i = 0; i < tmp_match.begin()->second.auto_header.size(); i++)
		{
			for (int j = 0; j < tmp_match.begin()->second.matched_data.size(); j++)
			{
				unordered_map<unsigned int, unsigned int> tmp_mapping;
				for (auto iter = tmp_match.begin(); iter != tmp_match.end(); iter++)
					tmp_mapping[iter->second.auto_header[i]] = iter->second.matched_data[j];
				result_compare(tmp_mapping);
			}
		}
	}
	void query_dfs(local_index* li, query_graph* qg, matching_order* mo, unordered_set<unsigned int>& used, unsigned int index, vector<vector<unsigned int>>& candidates,
		unordered_map<unsigned int, unsigned int>& tmp_maping, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int qid = mo->order[index];
		int parent = mo->parent[index];
		if (parent != -1) {
			unsigned int pcand = tmp_maping[parent];
			li->get_neighbor(parent, qid, pcand, candidates[qid]);
		}
		else
			li->get_candidates(qid, candidates[qid]);

		if (candidates[qid].empty())
			return;

		expanding_num++;
		vector<vector<unsigned int>> tmp_vecs;
		tmp_vecs.push_back(candidates[qid]);
		for (int i = 0; i < mo->left_neighbors[qid].size(); i++)
		{
			unsigned int u = mo->left_neighbors[qid][i].first;
			if (u == parent) {
				continue;
			}
			unsigned int v = tmp_maping[u];
			vector<unsigned int> vec;
			li->get_neighbor(v, merge_long_long(u, qid), vec);
			if (vec.empty())
			{
				for (int j = 0; j < tmp_vecs.size(); j++)
					tmp_vecs[j].clear();
				candidates[qid].clear();
				return;
			}
			tmp_vecs.push_back(vec);
		}
		//candidates[qid].clear();
		vector<unsigned int> joint;
		leap_frog_join(tmp_vecs, joint);
		for (int j = 0; j < tmp_vecs.size(); j++) {
			neighbor_sum += tmp_vecs[j].size();
			tmp_vecs[j].clear();
		}
		candidates[qid].swap(joint);
		joint.clear();
		filtered_expanding += candidates[qid].size();
		for (int i = 0; i < candidates[qid].size(); i++)
		{
			unsigned int candidate_id = candidates[qid][i];
			if (used.find(candidate_id) != used.end())
				continue;
			if (index + 1 == mo->order.size())
			{
				if (record_match) {
					for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
						match[iter->first].matched_data.push_back(iter->second);
					match[qid].matched_data.push_back(candidate_id);
				}
				match_size++;
			}
			else {
				tmp_maping[qid] = candidate_id;
				used.insert(candidate_id);
				inter_result++;
				query_dfs(li, qg, mo, used, index + 1, candidates, tmp_maping, match, match_size, record_match);
				used.erase(candidate_id);
				tmp_maping.erase(qid);
			}
		}
		candidates[qid].clear();
	}

	void query_processing(local_index* li, matching_order* mo, unordered_map<unsigned int, match_result>& match,
		unsigned int& match_size, bool record_match)
	{
		vector<vector<unsigned int>> candidates;
		candidates.resize(q->g.size());
		li->get_candidates(mo->order[0], candidates[mo->order[0]]);
		for (int i = 0; i < candidates[mo->order[0]].size(); i++)
		{
			unordered_map<unsigned int, unsigned int> tmp_map;
			tmp_map[mo->order[0]] = candidates[mo->order[0]][i];
			unordered_set<unsigned int> used;
			used.insert(candidates[mo->order[0]][i]);
			inter_result++;
			query_dfs(li, li->qg, mo, used, 1, candidates, tmp_map, match, match_size, record_match);
			used.clear();
			tmp_map.clear();
		}
		candidates[mo->order[0]].clear();
		candidates.clear();
	}

	bool NLF_check(unsigned int v, unsigned qv)
	{
		if (use_NLF) {
			for (auto iter = q->NLF[qv].begin(); iter != q->NLF[qv].end(); iter++)
			{
				if (g->NLF[v].find(iter->first) == g->NLF[v].end() || g->NLF[v][iter->first] < iter->second)
					return false;
			}
		}
		return true;
	}

	void query(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match = true)
	{
	//	cout << "begin query" << endl;
		vector<pair<unsigned int, unsigned int> > vec;
		q->get_edge_by_label(src_label, dst_label, edge_label, vec);
		if (vec.empty()) {
			label_filtered++;
			return;
		}
		bool enumerated = false;
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int qsrc = vec[i].first;
			unsigned int qdst = vec[i].second;
		//	cout << "exploring " << src << ' ' << dst << ' ' << qsrc << ' ' << qdst << endl;
			if (use_global_index)
			{
				if (!idx->check_candidate(qsrc, src) || !idx->check_candidate(qdst, dst))
					continue;
			}
			else if (!NLF_check(src, qsrc) || !NLF_check(dst, qdst))
				continue;
			enumerated = true;
			unsigned long long edge = merge_long_long(qsrc, qdst);
			if (auto_set.find(edge) == auto_set.end())
				continue;
			matching_order* mo;
			indexing_order* io;
			if (qsrc > qdst) {
				mo = mo_map[merge_long_long(qdst, qsrc)];
				io = io_map[merge_long_long(qdst, qsrc)];
			}
			else {
				mo = mo_map[edge];
				io = io_map[edge];
			}
			local_index* li = new local_index(src, dst, qsrc, qdst, q, idx, g, io, use_global_index);
			bool construct_result = true;
			if (use_global_index) {
				construct_result = li->construct_index();
			}
			else
				construct_result = li->graph_construct_index();
		//	cout << "local index constructed " << endl;
			if (construct_result) {
				unordered_map<unsigned int, match_result> tmp_match;
				unsigned int tmp_size = 0;
				query_processing(li, mo, tmp_match, tmp_size, record_match);
				if (tmp_size > 0)
				{
					if (record_match) {
						tmp_match[qsrc].matched_data.insert(tmp_match[qsrc].matched_data.end(), tmp_size, src);
						tmp_match[qdst].matched_data.insert(tmp_match[qdst].matched_data.end(), tmp_size, dst);
						for (int j = 0; j < auto_set[edge].size(); j++)
						{
							automatch* am = auto_set[edge][j].second;
							for (auto iter = tmp_match.begin(); iter != tmp_match.end(); iter++)
								iter->second.auto_header.push_back(am->mapping[iter->first]);
						}
						result_check(tmp_match);
						match.push_back(tmp_match);
					}
					match_size = match_size + tmp_size * auto_set[edge].size();
				}
			}
			delete li;
			clock_t currect_clock = clock();
			if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
				break;
		}
		if (!enumerated)
			filtered_edge++;
	}
	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>> &match,
		unsigned int &match_size, bool record_match = true)
	{
	//	cout << "inserting " <<src<<' '<<dst<< endl;
		if (src == dst) // in isormorphism we do not accept circles.
			return;
		bool newly_inserted = g->insert_edge(src, src_label, dst, dst_label, edge_label);
		if (!g->original_continuous_id)
		{
			src = g->vertice_id_map[src];
			dst = g->vertice_id_map[dst];
		}
		if (use_global_index&&newly_inserted)
			idx->insert_edge(src, src_label, dst, dst_label, edge_label);

		if (!indexing_test)
			query(src, src_label, dst, dst_label, edge_label, match, match_size, record_match);
	}

	void delete_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match = true)
	{
		if (src == dst) // in isormorphism we do not accept circles.
			return;
		if (!g->original_continuous_id)
		{
			src = g->vertice_id_map[src];
			dst = g->vertice_id_map[dst];
		}
		if (!indexing_test)
			query(src, src_label, dst, dst_label, edge_label, match, match_size, record_match);
		g->delete_edge(src, dst, edge_label);
		if (use_global_index)
			idx->delete_edge(src, src_label, dst, dst_label, edge_label);
	}
};
