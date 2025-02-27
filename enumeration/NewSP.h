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
#include<chrono>
#include"../utility/automorphism-check.h"
#include "../utility/OrderGeneration.h"
using namespace std;


class NewSP
{
public:
	data_graph* g;
	query_graph* q;
	my_index* idx;
	unordered_map<unsigned long long, operation_sequence*> op_map;
	bool use_global_index;
	NewSP(data_graph* dg, query_graph* qg, bool use_gd = false, bool use_edge_view = false)
	{
		g = dg;
		q = qg;
		use_global_index = use_gd;
		if (use_global_index)
			idx = new my_index(g, q, use_edge_view);
		else
			idx = NULL;
	}
	~NewSP()
	{
		if(use_global_index)
			delete idx;
		for (auto iter = op_map.begin(); iter != op_map.end(); iter++)
			delete iter->second;
		op_map.clear();
	}
	void load_initial_graph(vector<unsigned int>& src, vector<unsigned int>& dst, vector<int>& src_label, vector<int>& dst_label, vector<int>& edge_label)
	{
		for (int i = 0; i < src.size(); i++)
			g->insert_edge(src[i], src_label[i], dst[i], dst_label[i], edge_label[i]);
	}



	void load_query_graph(unsigned int nodenum, vector<unsigned int>& src, vector<unsigned int>& dst, vector<int>& src_label, vector<int>& dst_label, vector<int>& edge_label)
	{
		q->resize(nodenum);
		for (int i = 0; i < src.size(); i++)
			q->insert_edge(src[i], src_label[i], dst[i], dst_label[i], edge_label[i]);
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
				if (op_map.find(edge_pair) != op_map.end())
					continue;
	

				matching_order* mo = new matching_order;
				//edge_oriented_RI_order(q, mo, src, dst);
				GraphFlow_order_generation(q, mo, src, dst);
				operation_sequence* op = new operation_sequence;
				operation_sequence_generation(q, mo, op, src, dst, use_global_index);
				op_map[edge_pair] = op;
				delete mo;
			}
		}
		cout << "order generation finished" << endl;

		//print_autoset();
		double time_used = 0;
		if (use_global_index) {
			auto start = std::chrono::system_clock::now();
			idx->construct_index();
			auto end = std::chrono::system_clock::now();
			auto dur = std::chrono::duration<double>(end - start);
			time_used = dur.count();
		}
		return time_used;

	}

	bool NLF_check(query_graph* qg, data_graph* dg, unsigned int qid, unsigned int cand)
	{
		if (use_NLF) {
			for (auto iter = qg->NLF[qid].begin(); iter != qg->NLF[qid].end(); iter++)
			{
				if (dg->NLF[cand].find(iter->first) == dg->NLF[cand].end() || dg->NLF[cand][iter->first] < iter->second)
					return false;
			}
		}
		return true;
	}


	bool pre_compute_CPT(query_graph* qg, data_graph* dg, operation_sequence* op, unsigned int qsrc, unsigned int qdst, unordered_set<unsigned int>& used, unsigned int index, vector<int>& tmp_maping, vector<vector<unsigned int>>& CPT,
		vector<vector<unsigned int>>& pre_CPT, vector<bool>& tested)
	{
		unsigned int next_q = op->operation_id[index];
		unsigned int next_p = op->parent[index];
		int type = op->operation_type[index];
		bool cover_test = true;
		if (((type >> 5) & 1) == 1)
			cover_test = false;
		tested[next_q] = true;
		if (op->CPT_reuse.find(next_q) != op->CPT_reuse.end())
		{
			unsigned int ancestor = op->CPT_reuse[next_q];
			for (int i = 0; i < CPT[ancestor].size(); i++) {
				unsigned int candidate_id = CPT[ancestor][i];
				if (used.find(candidate_id) != used.end())
					continue;
				if (use_global_index)
				{
					if (!idx->check_candidate(next_q, candidate_id))
						continue;
				}
				else if (cover_test) {
					if (!NLF_check(qg, dg, next_q, candidate_id))
						continue;
				}
				pre_CPT[next_q].push_back(CPT[ancestor][i]);
			}
			if (pre_CPT[next_q].empty())
				return false;
			
			if (op->compensate_edges.find(next_q) != op->compensate_edges.end())
			{
				vector<vector<unsigned int>>tmp_vec;
				tmp_vec.push_back(pre_CPT[next_q]);
				for (int j = 0; j < op->compensate_edges[next_q].size(); j++)
				{
					int neighbor_id = op->compensate_edges[next_q][j].first;
					int neighbor_label = op->compensate_edges[next_q][j].second;
					if (tmp_maping[neighbor_id] != -1) {
						vector<unsigned int> vec;
						if (use_global_index)
							idx->get_neighbor(neighbor_id, next_q, tmp_maping[neighbor_id], vec);
						else
							dg->get_labeled_neighbor(tmp_maping[neighbor_id], q->g[next_q].label, neighbor_label, vec);
						if (vec.empty())
						{
							for (int i = 0; i < tmp_vec.size(); i++)
								tmp_vec[i].clear();
							tmp_vec.clear();
							vec.clear();
							return false;
						}
						tmp_vec.push_back(vec);
					}
				}
				
				vector<unsigned int> final_cand;
				leap_frog_join(tmp_vec, final_cand);
				pre_CPT[next_q].swap(final_cand);
				final_cand.clear();
				for (int i = 0; i < tmp_vec.size(); i++)
					tmp_vec[i].clear();
				tmp_vec.clear();
			}
		}
		else {
			vector<vector<unsigned int>>tmp_vec;
			for (int i = 0; i < op->left_neighbors[next_q].size(); i++)
			{
				int neighbor_id = op->left_neighbors[next_q][i].first;
				int neighbor_label = op->left_neighbors[next_q][i].second;
				if (tmp_maping[neighbor_id] != -1)
				{
					vector<unsigned int> vec;
					if (use_global_index)
						idx->get_neighbor(neighbor_id, next_q, tmp_maping[neighbor_id], vec);
					else
						dg->get_labeled_neighbor(tmp_maping[neighbor_id], q->g[next_q].label, neighbor_label, vec);
					if (vec.empty())
					{
						for (int i = 0; i < tmp_vec.size(); i++)
							tmp_vec[i].clear();
						tmp_vec.clear();
						vec.clear();
						return false;
					}
					tmp_vec.push_back(vec);
				}
			}
			vector<unsigned int> cand;
			leap_frog_join(tmp_vec, cand);
			for (int i = 0; i < cand.size(); i++)
			{
				unsigned int candidate_id = cand[i];
				if (used.find(candidate_id) != used.end())
					continue;
				if (!use_global_index && cover_test&& !NLF_check(qg, dg, next_q, candidate_id))
					continue;
				pre_CPT[next_q].push_back(candidate_id);
			}
			cand.clear();
			for (int i = 0; i < tmp_vec.size(); i++)
				tmp_vec[i].clear();
			tmp_vec.clear();
		
		}
		if (pre_CPT[next_q].empty())
			return false;
		else
			return true;
	}
	bool compute_CPT(query_graph* qg, data_graph* dg, operation_sequence* op, unsigned int qsrc, unsigned int qdst, unordered_set<unsigned int>& used, unsigned int index, 
		vector<int>& tmp_maping, vector<vector<unsigned int>>& CPT, vector<vector<unsigned int>>& pre_CPT)
	{
		unsigned int qid = op->operation_id[index];
		unsigned int type = op->operation_type[index];
		int parent = op->parent[index];
		if (!CPT[qid].empty())
		{
			assert(qid == qsrc || qid == qdst);
			return true;
		}
		if (pre_CPT[qid].empty()) {
			if (op->CPT_reuse.find(qid) != op->CPT_reuse.end()) // CPT reuse, we can use the candidate set of another query vertex as the initial set of this CPT computation.
			{
				unsigned int ancestor = op->CPT_reuse[qid];
				for (int i = 0; i < CPT[ancestor].size(); i++) {
					if (use_global_index)
					{
						if (!idx->check_candidate(qid, CPT[ancestor][i]))
							continue;
					}
					if (used.find(CPT[ancestor][i]) != used.end())
						continue;
					CPT[qid].push_back(CPT[ancestor][i]);
				}
				if (CPT[qid].empty())
					return false;
				if (op->compensate_edges.find(qid) != op->compensate_edges.end()) // the neighborhood of the recorded ancestor may be different with this query vertex, there may be additional neighbors
				{
					vector<vector<unsigned int>>tmp_vec;
					tmp_vec.push_back(CPT[qid]);
					for (int j = 0; j < op->compensate_edges[qid].size(); j++)
					{
						int neighbor_id = op->compensate_edges[qid][j].first;
						int neighbor_label = op->compensate_edges[qid][j].second;
						if (tmp_maping[neighbor_id] != -1) {
							vector<unsigned int> vec;
							if (use_global_index)
								idx->get_neighbor(neighbor_id, qid, tmp_maping[neighbor_id], vec);
							else
								dg->get_labeled_neighbor(tmp_maping[neighbor_id], q->g[qid].label, neighbor_label, vec);
							if (vec.empty())
							{
								for (int i = 0; i < tmp_vec.size(); i++)
									tmp_vec[i].clear();
								tmp_vec.clear();
								vec.clear();
								CPT[qid].clear();
								return false;
							}
							tmp_vec.push_back(vec);
						}
					}
					vector<unsigned int> final_cand;
					leap_frog_join(tmp_vec, final_cand);
					CPT[qid].swap(final_cand);
					final_cand.clear();
					for (int i = 0; i < tmp_vec.size(); i++)
						tmp_vec[i].clear();
					tmp_vec.clear();
				}
			}
			else {
				expanding_num++;
				vector<vector<unsigned int>>tmp_vec;
				for (int i = 0; i < op->left_neighbors[qid].size(); i++)
				{
					int neighbor_id = op->left_neighbors[qid][i].first;
					int neighbor_label = op->left_neighbors[qid][i].second;
					if (tmp_maping[neighbor_id] != -1)
					{
						vector<unsigned int> vec;
						if (use_global_index)
							idx->get_neighbor(neighbor_id, qid, tmp_maping[neighbor_id], vec);
						else
							dg->get_labeled_neighbor(tmp_maping[neighbor_id], q->g[qid].label, neighbor_label, vec);
						if (vec.empty())
						{
							for (int i = 0; i < tmp_vec.size(); i++)
								tmp_vec[i].clear();
							tmp_vec.clear();
							vec.clear();
							return false;
						}
						neighbor_sum += vec.size();
						tmp_vec.push_back(vec);
					}
				}
				vector<unsigned int> cand;
				leap_frog_join(tmp_vec, cand);
				for (int i = 0; i < cand.size(); i++)
				{
					unsigned int candidate_id = cand[i];
					if (used.find(candidate_id) != used.end())
						continue;
					CPT[qid].push_back(candidate_id);
				}
				filtered_expanding += CPT[qid].size();
				cand.clear();
				for (int i = 0; i < tmp_vec.size(); i++)
					tmp_vec[i].clear();
				tmp_vec.clear();
			}
		}
		else
		{
			vector<pair<unsigned int, int>> multi_exp;
			int cur = index - 1;
			while (cur >= 0 && ((op->operation_type[cur] & 1) == 1))
			{
				unsigned int suc = op->operation_id[cur];
				int e_label = qg->check_edge_label(suc, qid);
				assert(e_label != -1);
				multi_exp.push_back(make_pair(op->operation_id[cur], e_label));
				cur--;
			}
			if (!multi_exp.empty()) {
				vector<vector<unsigned int>>tmp_vec;
				tmp_vec.push_back(pre_CPT[qid]);
				for (int j = 0; j < multi_exp.size(); j++)
				{
					int neighbor_id = multi_exp[j].first;
					int neighbor_label = multi_exp[j].second;
					vector<unsigned int> vec;
					if (use_global_index)
						idx->get_neighbor(neighbor_id, qid, tmp_maping[neighbor_id], vec);
					else
						dg->get_labeled_neighbor(tmp_maping[neighbor_id], q->g[qid].label, neighbor_label, vec);
					if (vec.empty())
					{
						for (int i = 0; i < tmp_vec.size(); i++)
							tmp_vec[i].clear();
						tmp_vec.clear();
						vec.clear();
						return false;
					}
					tmp_vec.push_back(vec);
				}
				vector<unsigned int> final_cand;
				leap_frog_join(tmp_vec, final_cand);
				CPT[qid].clear();
				for (int i = 0; i < final_cand.size(); i++)
				{
					if (used.find(final_cand[i]) == used.end())
						CPT[qid].push_back(final_cand[i]);
				}
				final_cand.clear();
				for (int i = 0; i < tmp_vec.size(); i++)
					tmp_vec[i].clear();
				tmp_vec.clear();
			}
		}
		if (CPT[qid].empty())
			return false;
		else
			return true;
	}

	void nsp_dfs(query_graph* qg, data_graph* dg, operation_sequence* op, unsigned int qsrc, unsigned int qdst, unordered_set<unsigned int>& used, unsigned int index, vector<int>& tmp_maping, 
		 vector<vector<unsigned int>>& CPT, vector<vector<unsigned int>>& pre_CPT, vector<bool>& tested,
		unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match, int threshold)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int qid = op->operation_id[index];
		unsigned int type = op->operation_type[index];
		int parent = op->parent[index];

		if ((type & 1) == 0) // CPT operation  // CPT reuse and cover test omit need to be added 
		{
			bool compute_result = compute_CPT(qg, dg, op, qsrc, qdst, used, index, tmp_maping, CPT, pre_CPT);
			if (!compute_result)
				return;
			int next_q = -1;
			if (((type >> 3) & 1) == 1)
			{
				int cur = index + 1;
				while (cur < op->operation_id.size() && ((op->operation_type[cur] & 1) == 1))
					cur++;
				if (cur < op->operation_id.size()) {
					next_q = op->operation_id[cur];
					compute_result = pre_compute_CPT(qg, dg, op, qsrc, qdst, used, cur, tmp_maping, CPT, pre_CPT, tested);
					if (!compute_result)
					{
						CPT[qid].clear();
						tested[qid] = false;
						return;
					}
				}
			}
			nsp_dfs(qg, dg, op, qsrc, qdst, used, index + 1, tmp_maping, CPT, pre_CPT, tested, match, match_size, record_match, threshold);
			CPT[qid].clear();
			tested[qid] = false;
			if (next_q != -1)
				pre_CPT[next_q].clear();
		}
		else
		{
			for (int i = 0; i < CPT[qid].size(); i++)
			{
				if ((((type >> 2) & 1) != 1) || used.find(CPT[qid][i]) == used.end())
				{
					if (!use_global_index) {
						if (((type >> 5) & 1) != 1 && !tested[qid]) // need cover test
						{
							if (!NLF_check(qg, dg, qid, CPT[qid][i]))
								continue;
						}
					}
					if (index + 1 == op->operation_id.size())
					{
						if (record_match) {
							for (int j = 0; j < tmp_maping.size(); j++) {
								if (tmp_maping[j] != -1)
									match[j].matched_data.push_back(tmp_maping[j]);
							}
							match[qid].matched_data.push_back(CPT[qid][i]);
						}
						match_size++;
					}
					else
					{
							if(qid!=qsrc&&qid!=qdst)
								inter_result++; // record the number of inermediate result, but the match of qsrc and qdst need not to be recorded
							tmp_maping[qid] = CPT[qid][i];
							used.insert(CPT[qid][i]);
							nsp_dfs(qg, dg, op, qsrc, qdst, used, index + 1, tmp_maping, CPT, pre_CPT, tested, match, match_size, record_match, threshold);
							if (threshold != -1 && match_size >= threshold)
								return;
							used.erase(CPT[qid][i]);
							tmp_maping[qid]=-1;
						}
					}
				}
		}
	}

	void query_processing_nsp(unsigned int qsrc, unsigned int qdst, unsigned int src, unsigned int dst,
		query_graph* qg, data_graph* dg, operation_sequence* op, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match, int threshold)
	{
		unordered_set<unsigned int> used;
		vector<int> tmp_map;
		vector<vector<unsigned int>> CPT;
		vector<vector<unsigned int>> pre_CPT;
		vector<bool> tested;
		CPT.resize(qg->g.size());
		pre_CPT.resize(qg->g.size());
		tmp_map.resize(qg->g.size());
		tested.resize(qg->g.size());
		for (int i = 0; i < tmp_map.size(); i++)
		{
			tmp_map[i] = -1;
			tested[i] = false;
		}
		CPT[qsrc].push_back(src);
		CPT[qdst].push_back(dst);
		nsp_dfs(qg, dg, op, qsrc, qdst, used, 0, tmp_map, CPT, pre_CPT, tested, match, match_size, record_match, threshold);
		used.clear();
		CPT.clear();
		tmp_map.clear();
	}

	void query(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match = true, int threshold = -1)
	{
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
			if (use_global_index)
			{
				if (!idx->check_candidate(qsrc, src) || !idx->check_candidate(qdst, dst))
					continue;
			}
			else
			{
				if (!NLF_check(q, g, qsrc, src) || !NLF_check(q, g, qdst, dst))
					continue;
			}
			enumerated = true;
			unsigned long long edge = merge_long_long(qsrc, qdst);
			operation_sequence* op;
			if (qsrc > qdst)
				op = op_map[merge_long_long(qdst, qsrc)];
			else
				op = op_map[edge];
			unordered_map<unsigned int, match_result> tmp_match;
			unsigned int tmp_size = 0;
			int tmp_threshold = -1;
			if (threshold != -1)
				tmp_threshold = tmp_threshold - match_size;
			query_processing_nsp(qsrc, qdst, src, dst, q, g, op, tmp_match, tmp_size, record_match, tmp_threshold);
			if (!tmp_match.empty() || tmp_size > 0)
			{
				if (record_match) {
					for (auto iter = tmp_match.begin(); iter != tmp_match.end(); iter++)
						iter->second.auto_header.push_back(iter->first);
					if (!g->original_continuous_id)
					{
						for (auto iter = tmp_match.begin(); iter != tmp_match.end(); iter++)
						{
							for (int j = 0; j < iter->second.matched_data.size(); j++)
								iter->second.matched_data[j] = g->original_ids[iter->second.matched_data[j]];
						}
					}
					match.push_back(tmp_match);
				}
				match_size = match_size + tmp_size;
				if (threshold != -1 && match_size >= threshold)
					break;
			}
			clock_t currect_clock = clock();
			if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
				break;
		}
		if (!enumerated)
			filtered_edge++;
	}
	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match = true, int threshold = -1)
	{
	//	cout << "inserting edge " << src << ' ' << src_label << ' ' << dst << ' ' << dst_label << endl;
		if (src == dst) // in isormorphism we do not accept circles.
			return;
		bool newly_inserted = g->insert_edge(src, src_label, dst, dst_label, edge_label);
		if (!g->original_continuous_id)
		{
			src = g->vertice_id_map[src];
			dst = g->vertice_id_map[dst];
		}
		if (use_global_index && newly_inserted)
			idx->insert_edge(src, src_label, dst, dst_label, edge_label);
	//	cout << "index update finished" << endl;
		if(!indexing_test)
			query(src, src_label, dst, dst_label, edge_label, match, match_size, record_match, threshold);
	//	cout << "query finished" << endl;
	}

	void delete_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match = true, int threshold = -1)
	{
		if (src == dst) // in isormorphism we do not accept circles.
			return;
		if (!g->original_continuous_id)
		{
			src = g->vertice_id_map[src];
			dst = g->vertice_id_map[dst];
		}
		if (!indexing_test)
			query(src, src_label, dst, dst_label, edge_label, match, match_size, record_match, threshold);

		 g->delete_edge(src, dst, edge_label);
		 idx->delete_edge(src, src_label, dst, dst_label, edge_label);
	}
};
