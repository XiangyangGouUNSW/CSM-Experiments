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
using namespace std;


class CaliG
{
public:
	data_graph* dg;
	query_graph* qg;
	my_index* idx;
	unordered_map<unsigned long long, ks_order*> mo_map;
	bool use_global_index;
	CaliG(data_graph* dg_, query_graph* qg_, bool use_gd = false, bool use_edge_view = false)
	{
		dg = dg_;
		qg = qg_;
		use_global_index = use_gd;
		if (use_gd) 
			idx = new my_index(dg, qg, use_edge_view);
		else
			idx = NULL;
	}
	~CaliG()
	{
		for (auto iter = mo_map.begin(); iter != mo_map.end(); iter++)
			delete iter->second;
		if (use_global_index)
			delete idx;
	}

	double initialization()
	{
		for (int i = 0; i < qg->g.size(); i++)
		{
			vector<unsigned int> neighbor;
			qg->get_unlabeled_neighbor(qg->g[i].id, neighbor);
			for (int j = 0; j < neighbor.size(); j++)
			{
				unsigned int src = qg->g[i].id;
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
				ks_order* ko = new ks_order;
				CailiG_order_generation(qg, ko, src, dst);
				mo_map[edge_pair] = ko;
			}
		}
		double time_used = 0;
		if (use_global_index) {
			clock_t start = clock();
			idx->construct_index();
			clock_t end = clock();
			time_used = double(end - start) / CLOCKS_PER_SEC;
		}
		return time_used;
	}

	bool NLF_check(unsigned int qv, unsigned v)
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

	void candidate_initialize( unsigned int u, unsigned int parent, vector<vector<unsigned int>>& candidates, 
		unordered_map<unsigned int, unsigned int>& tmp_maping)
	{
		unsigned int p = parent;
		unsigned int u_lable = qg->g[u].label;
		if (use_global_index)
			idx->get_neighbor(p, u, tmp_maping[p], candidates[u]);
		else
		{
			int edge_label = qg->check_edge_label(p, u);
			dg->get_labeled_neighbor(tmp_maping[p], u_lable, edge_label, candidates[u]);
		}
	}

	void candidate_filter(unsigned int u, unsigned int parent, vector<pair<int, int>>& neighbor, unordered_set<unsigned int>& used, 
		unordered_map<unsigned int, unsigned int>& tmp_maping, vector<vector<unsigned int>>& candidates)
	{
		if (candidates[u].empty())
			return;
		vector<unsigned int> filtered_cand;
		vector<vector<unsigned int> >tmp_vecs;
		tmp_vecs.push_back(candidates[u]);
		for (int i = 0; i < neighbor.size(); i++)
		{
			unsigned int u2 = neighbor[i].first;
			unsigned int elabel = neighbor[i].second;
			if (u2 == parent)
				continue;
			vector<unsigned int> vec;
			if (use_global_index)
				idx->get_neighbor(u2, u, tmp_maping[u2], vec);
			else
				dg->get_labeled_neighbor(tmp_maping[u2], qg->g[u].label, elabel, vec);
			if (vec.empty())
			{
				for (int j = 0; j < tmp_vecs.size(); j++)
					tmp_vecs[j].clear();
				candidates[u].clear();
				return;
			}
			tmp_vecs.push_back(vec);
		}
		leap_frog_join(tmp_vecs, filtered_cand);
		candidates[u].clear();
		for (int i = 0; i < filtered_cand.size(); i++)
		{
			unsigned int candidate_id = filtered_cand[i];
			if (used.find(candidate_id) != used.end())
				continue;
			if (!use_global_index)
			{
				if (!NLF_check(u, candidate_id))
					continue;
			}
			candidates[u].push_back(candidate_id);
		}
		filtered_cand.clear();

	}

	bool test_shell(unsigned int u, unsigned int parent, vector<pair<int, int>>& neighbor, unordered_set<unsigned int>& used,
		unordered_map<unsigned int, unsigned int>& tmp_maping, vector<vector<unsigned int>>& candidates)
	{
	//	if (print)
	//		cout << "test shell " << u << endl;
		candidate_initialize(u, parent, candidates, tmp_maping);
		bool test_result = false;
		for (int i = 0; i < candidates[u].size(); i++)
		{
			unsigned int candidate_id = candidates[u][i];
			if (used.find(candidate_id) != used.end())
				continue;
			if (!use_global_index)
			{
				if (!NLF_check(u, candidate_id))
					continue;
			}
			if (neighbor.size() == 1)
			{
				assert(neighbor[0].first == parent);
				test_result = true;
				break;
			}
			bool filtered = false;
			for (int j = 0; j < neighbor.size(); j++)
			{
				if (neighbor[j].first == parent)
					continue;
				if (!dg->check_edge_label(candidate_id, tmp_maping[neighbor[j].first], neighbor[j].second))
				{
					filtered = true;
					break;
				}
			}
			if (!filtered)
			{
				test_result = true;
				break;
			}
		}
		candidates[u].clear();
	//	if (print)
	//		cout << "test shell finished" << endl;
		return test_result;
	}
	void shell_enumerate(vector<unsigned int>& shells, unordered_set<unsigned int>& used, unsigned int index, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int u = shells[index];
		if (index + 1 == shells.size())
		{
			if (record_match)
			{
				for (int i = 0; i < candidates[u].size(); i++) {
					if (used.find(candidates[u][i]) != used.end())
						continue;
					for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
						match[iter->first].matched_data.push_back(iter->second);
					match[u].matched_data.push_back(candidates[u][i]);
					match_size++;
				}
			}
			else
			{
				//sort(candidates[u].begin(), candidates[u].end()); // here we donot need to sort, as the candidates we retrive are always sorted. 
				unsigned int del = 0;
				for (auto iter = used.begin(); iter != used.end(); iter++)
				{
					auto iter2 = lower_bound(candidates[u].begin(), candidates[u].end(), *iter);
					if (iter2 != candidates[u].end() && *iter2 == *iter)
						del++;
				}
				match_size += candidates[u].size() - del;
			}
		}
		else
		{
			for (int i = 0; i < candidates[u].size(); i++)
			{
				unsigned int v = candidates[u][i];
				if (used.find(v) == used.end())
				{
					used.insert(v);
					tmp_maping[u] == v;
					inter_result++;
					shell_enumerate(shells, used, index + 1, tmp_maping, candidates, match, match_size, record_match);
					tmp_maping.erase(u);
					used.erase(v);
				}
			}
		}
	}
	bool shell_match(ks_order* ko, unordered_set<unsigned int>& used,  unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		vector<unsigned int> shell_vec;
		for (auto iter = ko->shell_vertices.begin(); iter != ko->shell_vertices.end(); iter++)
			shell_vec.push_back(*iter);
		bool filtered = false;
		for (int j = 0; j < shell_vec.size(); j++) {
			unsigned int shell = shell_vec[j];
			candidate_initialize(shell, ko->shell_neighbors[shell][0].first, candidates, tmp_maping);
			candidate_filter(shell, ko->shell_neighbors[shell][0].first, ko->shell_neighbors[shell], used, tmp_maping, candidates);
			if (candidates[shell].empty())
			{
				filtered = true;
				break;
			}
		}
		if (filtered)
		{
			for (int k = 0; k < shell_vec.size(); k++)
				candidates[shell_vec[k]].clear();
			return false;
		}
		shell_enumerate(shell_vec, used, 0, tmp_maping, candidates, match, match_size, record_match);
		for (int j = 0; j < shell_vec.size(); j++)
			candidates[shell_vec[j]].clear();
		return true;
	}
	void kernel_dfs(ks_order* ko, unordered_set<unsigned int>& used, unsigned int index, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int qid = ko->mo.order[index];
	//	cout << "kernel dfs for " << qid << endl;
		int parent = ko->mo.parent[index];
		if (parent != -1) {
			candidate_initialize(qid, parent, candidates, tmp_maping);
			candidate_filter(qid, parent, ko->mo.left_neighbors[qid], used, tmp_maping, candidates);
		//	if (print)
		//		cout << "candidate generate for " << qid << " finished" << endl;
			for (int i = 0; i < candidates[qid].size(); i++)
			{
				unsigned int candidate_id = candidates[qid][i];
				/*if (used.find(candidate_id) != used.end())
					continue;
				if (!use_global_index)
				{
					if (!NLF_check(qid, candidate_id))
						continue;
				}*/
				tmp_maping[qid] = candidate_id;
				used.insert(candidate_id);
				if (!ko->shell_to_test[index].empty())
				{
					bool filtered = false;
					for (int j = 0; j < ko->shell_to_test[index].size(); j++)
					{
						unsigned int shell = ko->shell_to_test[index][j];
						if (!test_shell(shell, ko->shell_neighbors[shell][0].first, ko->shell_neighbors[shell], used, tmp_maping, candidates))
						{
							filtered = true;
							break;
						}
					}
					if (filtered) {
						used.erase(candidate_id);
						tmp_maping.erase(qid);
						continue;
					}
				}

				inter_result++;
				if (index + 1 == ko->mo.order.size()) // kernel enummerate phase finished.
					shell_match(ko, used, tmp_maping, candidates, match, match_size, record_match);
				else 
					kernel_dfs(ko, used, index + 1, tmp_maping, candidates, match, match_size, record_match);
		
				used.erase(candidate_id);
				tmp_maping.erase(qid);
			}
			candidates[qid].clear();
		}
		else
		{
			cout << "error!!disconnected graph" << endl;
		}
	}
	void query_processing(unsigned int src, unsigned int dst, unsigned int qsrc, unsigned int qdst, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match)
	{
		unsigned long long edge_pair = (qsrc > qdst ? merge_long_long(qdst, qsrc) : merge_long_long(qsrc, qdst));
		ks_order* ko = mo_map[edge_pair];
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
		if (ko->mo.order.size() == 2)
			shell_match(ko, used, tmp_mapping, candidates, tmp_match, match_size, record_match);
		else
			kernel_dfs(ko, used, 2, tmp_mapping, candidates, tmp_match, match_size, record_match);
		if (record_match)
			match.push_back(tmp_match);
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
		if (!indexing_test) {
	//		if(print)
		//		cout << "index update finished " << endl;
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
			//	if(print)
		//			cout << "matching to " << qsrc << ' ' << qdst << endl;
				if (use_global_index)
				{
					if (!idx->check_candidate(qsrc, src) || !idx->check_candidate(qdst, dst))
						continue;
				}
				else if (!NLF_check(qsrc, src) || !NLF_check(qdst, dst))
					continue;
			//	if (print)
		//			cout << "index check finished " << endl;
				enumerated = true;
				query_processing(src, dst, qsrc, qdst, match, match_size, record_match);
				clock_t currect_clock = clock();
				if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
					break;
			}
			if (!enumerated)
				filtered_edge++;
		}
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
				else if (!NLF_check(qsrc, src) || !NLF_check(qdst, dst))
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
