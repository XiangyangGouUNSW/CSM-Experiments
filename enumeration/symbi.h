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

class Symbi
{
public:
	data_graph* dg;
	query_graph* qg;
	my_index* idx;
	bool use_global_index;
	unordered_map<unsigned int, vector<unsigned int>> equal_vertices;
	unordered_map<unsigned long long, vector<unsigned long long>> equal_edges;
	unordered_map<unsigned long long, matching_order*> mo_map;
	Symbi(data_graph* dg_, query_graph* qg_, bool use_gd = false, bool use_edge_view = false)
	{
		dg = dg_;
		qg = qg_;
		use_global_index = use_gd;
		if (use_gd)
			idx = new my_index(dg, qg, use_edge_view);
		else
			idx = NULL;
	}
	~Symbi()
	{
		if (use_global_index)
			delete idx;
		for (auto iter = equal_vertices.begin(); iter != equal_vertices.end(); iter++)
			iter->second.clear();
		equal_vertices.clear();
		for (auto iter = equal_edges.begin(); iter != equal_edges.end(); iter++)
			iter->second.clear();
		equal_edges.clear();

	}

	void compute_equal_vertices()
	{
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u = qg->g[i].id;
			vector<neighbor_info> vec;
			qg->get_neighbor(u, vec);
			unordered_map<unsigned int, unsigned int> umap;
			for (int i = 0; i < vec.size(); i++)
				umap[vec[i].id] = vec[i].e_label;
			vec.clear();
			unsigned int label = qg->g[u].label;
			vector<unsigned int> cand;
			qg->get_vertex_by_label(label, cand);
			for (int j = 0; j < cand.size(); j++)
			{
				unsigned int u2 = cand[j];
				if (u2 == u)
					continue;
				vector<neighbor_info> vec;
				qg->get_neighbor(u2, vec);
				if (vec.size() != umap.size())
					continue;
				bool filtered = false;
				for (int k = 0; k < vec.size(); k++)
				{
					if (umap.find(vec[k].id) == umap.end() || umap[vec[k].id] != vec[k].e_label)
					{
						filtered = true;
						break;
					}
				}
				if (!filtered)
					equal_vertices[u].push_back(u2);
			}
		}
	}

	void compute_equal_edges()
	{
		unordered_set<unsigned long long> recorded;
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u1 = qg->g[i].id;
			vector<unsigned int> vec;
			qg->get_unlabeled_neighbor(u1, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int u2 = vec[j];
				unsigned long long edge = merge_long_long(u1, u2);
				if (recorded.find(edge) != recorded.end())
					continue;
				equal_edges[edge].push_back(edge);
				recorded.insert(edge);
				if (equal_vertices.find(u2) != equal_vertices.end())
				{
					for (int k = 0; k < equal_vertices[u2].size(); k++)
					{
						unsigned int u3 = equal_vertices[u2][k];
						unsigned long long cand_edge = merge_long_long(u1, u3);
						if (recorded.find(cand_edge) != recorded.end())
							continue;
						recorded.insert(cand_edge);
						equal_edges[edge].push_back(cand_edge);
					 }
				}
				if (equal_vertices.find(u1) != equal_vertices.end())
				{
					for (int k = 0; k < equal_vertices[u1].size(); k++)
					{
						unsigned int u3 = equal_vertices[u1][k];
						unsigned long long cand_edge = merge_long_long(u3, u2);
						if (recorded.find(cand_edge) != recorded.end())
							continue;
						recorded.insert(cand_edge);
						equal_edges[edge].push_back(cand_edge);
					}
				}
			}
		}
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
		compute_equal_vertices();
		compute_equal_edges();

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


	void isolated_vertices_match(unordered_set<unsigned int>& used, vector<unsigned int>& isolated, unsigned int index, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		unsigned int qid = isolated[index];
		for (int i = 0; i < candidates[qid].size(); i++)
		{
			unsigned int v = candidates[qid][i];
			if (used.find(v) != used.end())
				continue;
			if (index + 1 == isolated.size())
			{
				if (record_match) {
					for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
						match[iter->first].matched_data.push_back(iter->second);
					match[qid].matched_data.push_back(v);
				}
				match_size++;
			}
			else {
				tmp_maping[qid] = v;
				used.insert(v);
				isolated_vertices_match(used, isolated, index + 1, tmp_maping, candidates, match, match_size, record_match);
				used.erase(v);
				tmp_maping.erase(qid);
			}
		}
	}

	unsigned int choose_next(unordered_map<unsigned int, unsigned int>& tmp_maping, vector<vector<unsigned int> > &candidates)
	{
		unsigned int min_size = 0xFFFFFFFF;
		int min_u = -1;
		for (int i= 0; i<qg->g.size();i++)
		{
			unsigned int u = qg->g[i].id;
			if (!candidates[u].empty()||tmp_maping.find(u)!=tmp_maping.end())
				continue;
			unsigned int min_u_size = 0xFFFFFFFF;
			vector<neighbor_info> vec;
			qg->get_neighbor(u, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int u2 = vec[j].id;
				if (tmp_maping.find(u2) != tmp_maping.end())
				{
					vector<unsigned int> tmp_vec;
					if (use_global_index) {
						idx->get_neighbor(u2, u, tmp_maping[u2], tmp_vec);
					}
					else
					{
						dg->get_labeled_neighbor(tmp_maping[u2], qg->g[u].label, vec[j].e_label, tmp_vec);
						// if we use NLF check here, the cost will be heavily increased, the the entire throughput will be decreased. Thus we donot include it here, 
							// this may increase the estimation error
					}
					if (tmp_vec.size() < min_u_size)
						min_u_size = tmp_vec.size();
					tmp_vec.clear();
				}
			}
			if (min_u_size < min_size)
			{
				min_size = min_u_size;
				min_u = u;
			}
		}
		return min_u;
	}

	void compute_canidates(unsigned int qid, unordered_set<unsigned int>& used, vector<unsigned int>& match_order, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates)
	{
		int qlabel = qg->g[qid].label;
		if (equal_vertices.find(qid) != equal_vertices.end())
		{
			unsigned int min_size = 0xFFFFFFFF;
			int min_ancestor = -1;
			for (int i = 0; i < equal_vertices[qid].size(); i++)
			{
				unsigned int u = equal_vertices[qid][i];
				if (!candidates[u].empty())
				{
					if (candidates[u].size() < min_size)
					{
						min_size = candidates[u].size();
						min_ancestor = u;
					}
				}
			}
			if (min_ancestor != -1)
			{
				vector<vector<unsigned int>>tmp_vecs;
				for (int i = match_order.size() - 1; i >= 0; i--)
				{
					if (match_order[i] == min_ancestor)
						break;
					int edge_label = qg->check_edge_label(qid, match_order[i]);
					if (edge_label != -1)
					{
						unsigned int neighbor_id = match_order[i];
						vector<unsigned int> vec;
						if (use_global_index)
							idx->get_neighbor(neighbor_id, qid, tmp_maping[neighbor_id], vec);
						else
							dg->get_labeled_neighbor(tmp_maping[neighbor_id], qlabel, edge_label, vec);
						if (vec.empty())
						{
							for (int j = 0; j < tmp_vecs.size(); j++) {
								tmp_vecs[j].clear();
							}
							tmp_vecs.clear();
							return;
						}
						tmp_vecs.push_back(vec);
					}
				}
				if (tmp_vecs.empty()) {
					candidates[qid] = candidates[min_ancestor];
				}
				else
				{
					tmp_vecs.push_back(candidates[min_ancestor]);
				//	vector<unsigned int> joint;
					leap_frog_join(tmp_vecs, candidates[qid]); // here we donot need to check the index or NLF, as they have been checked in the candidates generation of ancestor;
					for (int j = 0; j < tmp_vecs.size(); j++) {
						tmp_vecs[j].clear();
					}
					tmp_vecs.clear();
				}
				return;
			}
		}
		vector<vector<unsigned int>> tmp_vecs;
		vector<neighbor_info> nei;
		qg->get_neighbor(qid, nei);
		for (int i = 0; i < nei.size(); i++)
		{
			unsigned int neighbor_id = nei[i].id;
			int edge_label = nei[i].e_label;
			if (tmp_maping.find(neighbor_id) == tmp_maping.end())
				continue;
			vector<unsigned int> vec;
			if (use_global_index)
				idx->get_neighbor(neighbor_id, qid, tmp_maping[neighbor_id], vec);
			else
				dg->get_labeled_neighbor(tmp_maping[neighbor_id], qlabel, edge_label, vec);
			if (vec.empty())
			{
				for (int j = 0; j < tmp_vecs.size(); j++) {
					tmp_vecs[j].clear();
				}
				tmp_vecs.clear();
				return;
			}
			tmp_vecs.push_back(vec);
		}
		expanding_num++;
	//	vector<unsigned int> joint;
		leap_frog_join(tmp_vecs, candidates[qid]);
		for (int i = 0; i < tmp_vecs.size(); i++) {
			//	check_order(tmp_candidates[i]);
			tmp_vecs[i].clear();
		}
		tmp_vecs.clear();
		filtered_expanding += candidates[qid].size();
	}
	void lf_query_dfs( unordered_set<unsigned int>& used, vector<unsigned int>& isolated, vector<unsigned int> &match_order, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int qid = choose_next(tmp_maping, candidates);
		compute_canidates(qid, used, match_order, tmp_maping, candidates);
		if (candidates[qid].empty())
			return;
		match_order.push_back(qid);
		vector<unsigned int> vec;
		qg->get_unlabeled_neighbor(qid, vec);
		bool iso = true;
		for (int i = 0; i < vec.size(); i++)
		{
			if (tmp_maping.find(vec[i]) == tmp_maping.end())
			{
				iso = false;
				break;
			}
		}
		if (iso)
		{
			isolated.push_back(qid);
			if (match_order.size() == qg->g.size())
				isolated_vertices_match(used, isolated, 0, tmp_maping, candidates, match, match_size, record_match);
			else
				lf_query_dfs(used, isolated, match_order, tmp_maping, candidates, match, match_size, record_match);
			isolated.erase(isolated.begin() + isolated.size() - 1);
			match_order.erase(match_order.begin() + match_order.size() - 1);
			candidates[qid].clear();
			return;
		}
		for (int i = 0; i < candidates[qid].size(); i++)
		{
			unsigned int candidate_id = candidates[qid][i];
			if (used.find(candidate_id) != used.end())
				continue;
			if (!use_global_index)
			{
				if (!NLF_check(candidate_id, qid))
					continue;
			}
			used.insert(candidate_id);
			tmp_maping[qid] = candidate_id;
			if (match_order.size() == qg->g.size()) {
				if (isolated.empty())
				{
					if (record_match) {
						for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
							match[iter->first].matched_data.push_back(iter->second);
					}
					match_size++;
				}
				else
					isolated_vertices_match(used, isolated, 0, tmp_maping, candidates, match, match_size, record_match);
			}
			else {
				inter_result++;
				lf_query_dfs(used, isolated, match_order, tmp_maping, candidates, match, match_size, record_match);
			}
			used.erase(candidate_id);
			tmp_maping.erase(qid);
		}
		match_order.erase(match_order.begin() + match_order.size() - 1);
		candidates[qid].clear();
		return;
	}
	void query_processing(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label,  vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match)
	{
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
			unsigned long long edge = merge_long_long(qsrc, qdst);
			if (equal_edges.find(edge) == equal_edges.end())
				continue;
			if (use_global_index)
			{
				if (!idx->check_candidate(qsrc, src) || !idx->check_candidate(qdst, dst))
					continue;
			}
			else if (!NLF_check(src, qsrc) || !NLF_check(dst, qdst))
				continue;
			enumerated = true;

			unordered_map<unsigned int, match_result> tmp_match;
			for (int i = 0; i < qg->g.size(); i++)
				tmp_match[i] = match_result();
			unsigned int tmp_size = 0;
			unordered_map<unsigned int, unsigned int> tmp_mapping;
			tmp_mapping[qsrc] = src;
			tmp_mapping[qdst] = dst;
			unordered_set<unsigned int> used;
			used.insert(src);
			used.insert(dst);
			vector < vector<unsigned int> > candidates;
			candidates.resize(qg->g.size());
			vector<unsigned int> match_order;
			match_order.push_back(qsrc);
			match_order.push_back(qdst);
			vector<unsigned int> isolated;
			unsigned long long edge_pair = (qsrc > qdst ? merge_long_long(qdst, qsrc) : merge_long_long(qsrc, qdst));
			matching_order* mo = mo_map[edge_pair];
			lf_query_dfs(used, isolated, match_order, tmp_mapping, candidates, tmp_match, tmp_size, record_match);
			if (record_match) {
				if (tmp_match.begin()->second.matched_data.size() > 0) {
					for (int j = 0; j < equal_edges[edge].size();j++)
					{
						unsigned int eq_src = (equal_edges[edge][j] >> 32);
						unsigned int eq_dst = (equal_edges[edge][j] & 0xFFFFFFFF);
						for (auto iter = tmp_match.begin(); iter != tmp_match.end(); iter++)
						{
							unsigned int u = iter->first;
							if (u == qsrc)
								u = eq_src;
							if (u == eq_src) // exchange the equal edges;
								u = qsrc;

							if (u == qdst)
								u = eq_dst;
							if (u == eq_dst) // exchange the equal edges;
								u = qdst;
							iter->second.auto_header.push_back(u);
						}
					}
					match.push_back(tmp_match);
				}
			}
			match_size += tmp_size * equal_edges[edge].size();
			used.clear();
			tmp_mapping.clear();
			for (int i = 0; i < candidates.size(); i++)
				candidates[i].clear();
			candidates.clear();
			clock_t currect_clock = clock();
			if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
				break;
		}
		if (!enumerated)
			filtered_edge++;
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
		query_processing(src, src_label, dst, dst_label, edge_label, match, match_size, record_match);
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
			query_processing(src, src_label, dst, dst_label, edge_label, match, match_size, record_match);
		}
		dg->delete_edge(src, dst, edge_label);
		if (use_global_index)
			idx->delete_edge(src, src_label, dst, dst_label, edge_label);
	}


};
