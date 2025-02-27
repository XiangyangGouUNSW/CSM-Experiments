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
#include<assert.h>
#include<algorithm>
#include<chrono>
#include<time.h>
#include"../utility/automorphism-check.h"
#include "../utility/OrderGeneration.h"
using namespace std;

class GraphFlow
{
public:
	data_graph* dg;
	query_graph* qg;
	my_index* idx;
	unordered_map<unsigned long long, matching_order*> mo_map;
	bool use_global_index;
	double index_update_time;
	GraphFlow(data_graph* dg_, query_graph* qg_, bool use_gd = false, bool use_edge_view = false)
	{
		dg = dg_;
		qg = qg_;
		use_global_index = use_gd;
		if (use_gd)
			idx = new my_index(dg, qg, use_edge_view);
		else
			idx = NULL;
		index_update_time = 0;
	}
	~GraphFlow()
	{
		if (use_global_index)
			delete idx;
	}

	double initialization() // return the time of building index;
	{
		double time_used = 0;
		if (use_global_index) {
			auto start = std::chrono::system_clock::now();
			idx->construct_index();
			auto end = std::chrono::system_clock::now();
			auto dur = std::chrono::duration<double>(end - start);
			time_used = dur.count();
		}
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int src = qg->g[i].id;
			vector<unsigned int> neighbor;
			qg->get_unlabeled_neighbor(src, neighbor);
			for (int j = 0; j < neighbor.size(); j++)
			{
				src = qg->g[i].id;
				unsigned int dst = neighbor[j];
				if (src > dst)
				{
					unsigned int tmp = src;
					src = dst;
					dst = tmp;
				}
				unsigned long long edge_pair = merge_long_long(src, dst);
				if (mo_map.find(edge_pair) != mo_map.end())
					continue;
				matching_order* mo = new matching_order;
				GraphFlow_order_generation(qg, mo, src, dst);
				mo_map[edge_pair] = mo;
			}
		}
		return time_used;
	}

	bool NLF_check(unsigned int v, unsigned qv)
	{
		if (use_NLF) {
			for (auto iter = qg->NLF[qv].begin(); iter != qg->NLF[qv].end(); iter++)
			{
				if (dg->NLF[v].find(iter->first) == dg->NLF[v].end() || dg->NLF[v][iter->first] < iter->second)
					return false;
			}
		}
		return true;
	}


/*	void query_dfs(matching_order* mo, unordered_set<unsigned int>& used, unsigned int index, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>> &candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int qid = mo->order[index];
		int parent = mo->parent[index];
		int qlabel = qg->g[qid].label;
		if (parent != -1) {
			unsigned int pcand = tmp_maping[parent];
			if (!use_global_index) {
				int edge_label = qg->check_edge_label(parent, qid);
				dg->get_labeled_neighbor(pcand, qlabel, edge_label, candidates[qid]);
				vector<unsigned int> filtered;
				for (int i = 0; i < candidates[qid].size(); i++)
				{
					if (NLF_check(candidates[qid][i], qid))
						filtered.push_back(candidates[qid][i]);
				}
				candidates[qid].swap(filtered);
				filtered.clear();
			}
			else {
				idx->get_neighbor(parent, qid, pcand, candidates[qid]);
			}

			for (int i = 0; i < candidates[qid].size(); i++)
			{
				unsigned int candidate_id = candidates[qid][i];
				if (used.find(candidate_id) != used.end())
					continue;
				bool filtered = false;
				for (int j = 0; j < mo->left_neighbors[qid].size(); j++)
				{
					if (mo->left_neighbors[qid][j].first == parent)
						continue;
					if (!dg->check_edge_label(candidate_id, tmp_maping[mo->left_neighbors[qid][j].first], mo->left_neighbors[qid][j].second))
					{
						filtered = true;
						break;
					}
				}
				if (!filtered)
				{
					if (index + 1 == mo->order.size())
					{
						if (record_match) {
							for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
								match[iter->first].matched_data.push_back(iter->second);
							match[qid].matched_data.push_back(candidate_id);
						}
						//	cout << "find match" << endl;
						match_size++;
					}
					else {
						tmp_maping[qid] = candidate_id;
						used.insert(candidate_id);
						query_dfs(mo, used, index + 1, tmp_maping, candidates, match, match_size, record_match);
						used.erase(candidate_id);
						tmp_maping.erase(qid);
					}
				}
			}
			candidates[qid].clear();
		}
		else
		{
			cout << "error!!disconnected graph" << endl;
		}
	}*/

	void lf_query_dfs(matching_order* mo, unordered_set<unsigned int>& used, unsigned int index, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int qid = mo->order[index];
		int qlabel = qg->g[qid].label;
		vector<vector<unsigned int>> tmp_candidates;
		tmp_candidates.resize(mo->left_neighbors[qid].size());
		for (int i = 0; i < mo->left_neighbors[qid].size(); i++)
		{
			unsigned int neighbor_id = mo->left_neighbors[qid][i].first;
			int edge_label = mo->left_neighbors[qid][i].second;
			if (use_global_index)
				idx->get_neighbor(neighbor_id, qid, tmp_maping[neighbor_id], tmp_candidates[i]);
			else
				dg->get_labeled_neighbor(tmp_maping[neighbor_id], qlabel, edge_label, tmp_candidates[i]);
			if (tmp_candidates[i].empty())
			{
				for (int j = 0; j < tmp_candidates.size(); j++) {
					tmp_candidates[j].clear();
				}
				tmp_candidates.clear();
				return;
			}
		}
		expanding_num++;
		leap_frog_join(tmp_candidates, candidates[qid]);
		for (int i = 0; i < tmp_candidates.size(); i++) {
		//	check_order(tmp_candidates[i]);
			neighbor_sum += tmp_candidates[i].size();
			tmp_candidates[i].clear();
		}
		filtered_expanding += candidates[qid].size();
		tmp_candidates.clear();

		for (int i = 0; i < candidates[qid].size(); i++)
		{
			unsigned int candidate_id = candidates[qid][i];
			if (used.find(candidate_id) != used.end())
				continue;
			if (!use_global_index && !NLF_check(candidate_id, qid))
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
				inter_result++;
				tmp_maping[qid] = candidate_id;
				used.insert(candidate_id);
				lf_query_dfs(mo, used, index + 1, tmp_maping, candidates, match, match_size, record_match);
				used.erase(candidate_id);
				tmp_maping.erase(qid);
			}
		}
		candidates[qid].clear();
	}
	void query_processing(unsigned int src, unsigned int dst, unsigned int qsrc, unsigned int qdst, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match)
	{
		unsigned long long edge_pair = (qsrc > qdst ? merge_long_long(qdst, qsrc) : merge_long_long(qsrc, qdst));
		matching_order* mo = mo_map[edge_pair];
		unordered_map<unsigned int, match_result> tmp_match;
		for (int i = 0; i < qg->g.size(); i++)
			tmp_match[i].auto_header.push_back(i);
		unordered_map<unsigned int, unsigned int> tmp_mapping;
		tmp_mapping[qsrc] = src;
		tmp_mapping[qdst] = dst;
		unordered_set<unsigned int> used;
		used.insert(src);
		used.insert(dst);
		vector < vector<unsigned int> > candidates;
		candidates.resize(qg->g.size());
		lf_query_dfs(mo, used, 2, tmp_mapping, candidates, tmp_match, match_size, record_match);
		if (record_match) {
			if(tmp_match.begin()->second.matched_data.size()>0)
				match.push_back(tmp_match);
		}
		used.clear();
		tmp_mapping.clear();
		for (int i = 0; i < candidates.size(); i++)
			candidates[i].clear();
		candidates.clear();
	}

	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match, unsigned int& match_size, bool record_match = true)
	{
		if (src == dst)
			return;
		bool newly_inserted = dg->insert_edge(src, src_label, dst, dst_label, edge_label);
		if (!dg->original_continuous_id)
		{
			src = dg->vertice_id_map[src];
			dst = dg->vertice_id_map[dst];
		}
		if (use_global_index && newly_inserted)
			idx->insert_edge(src, src_label, dst, dst_label, edge_label);
		if (indexing_test)
			return;
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		if (vec.empty()) {
			label_filtered++;
			return;
		}
		bool enumerated = false;
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int qsrc = vec[i].first;
			unsigned int qdst = vec[i].second;
			if (use_global_index)
			{
				if (!idx->check_candidate(qsrc, src) || !idx->check_candidate(qdst, dst))
					continue;
			}
			else if (!NLF_check(src, qsrc) || !NLF_check(dst, qdst))
				continue;
			unsigned int tmp = match_size;
			enumerated = true;
			query_processing(src, dst, qsrc, qdst, match, match_size, record_match);
			clock_t currect_clock = clock();
			if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
				break;
		}
		if (!enumerated)
			filtered_edge++;
	}

	void delete_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match, unsigned int& match_size, bool record_match = true)
	{
		if (src == dst)
			return;
		if (!dg->original_continuous_id)
		{
			src = dg->vertice_id_map[src];
			dst = dg->vertice_id_map[dst];
		}
		if (!indexing_test) {
			vector<pair<unsigned int, unsigned int> > vec;
			qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int qsrc = vec[i].first;
				unsigned int qdst = vec[i].second;
				if (use_global_index)
				{
					if (!idx->check_candidate(qsrc, src) || !idx->check_candidate(qdst, dst))
						continue;
				}
				else if (!NLF_check(src, qsrc) || !NLF_check(dst, qdst))
					continue;
				query_processing(src, dst, qsrc, qdst, match, match_size, record_match);
				clock_t currect_clock = clock();
				if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
					break;
			}
		}
		dg->delete_edge(src, dst, edge_label);
		if (use_global_index)
			idx->delete_edge(src, src_label, dst, dst_label, edge_label);
	}


};
