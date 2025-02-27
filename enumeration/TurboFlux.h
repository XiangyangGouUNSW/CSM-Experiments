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
#include <chrono>
using namespace std;

struct inter_record
{
	map<unsigned int, unsigned long long> child_inter;
	unsigned long long inter;
	inter_record() {
		inter = 0;
	}
	~inter_record()
	{
		child_inter.clear();
	}
	void copy(inter_record& tmp1)
	{
		for (auto iter = tmp1.child_inter.begin(); iter != tmp1.child_inter.end(); iter++)
			child_inter[iter->first] = iter->second;
		inter = tmp1.inter;
	}
};
class TurboFlux
{
public:
	data_graph* dg;
	query_graph* qg;
	my_index* idx;
	matching_order* global_mo;
	query_tree* qt;
	bool use_global_index;
	TurboFlux(data_graph* dg_, query_graph* qg_, bool use_gd = false, bool use_edge_view = false)
	{
		dg = dg_;
		qg = qg_;
		use_global_index = use_gd;
		if (use_gd)
			idx = new my_index(dg, qg, use_edge_view);
		else
			idx = NULL;
		qt = new query_tree();
	}
	~TurboFlux()
	{
		delete global_mo;
		delete qt;
		if (use_global_index)
			delete idx;
	}

	unsigned long long compute_inter(TF_index* tf, vector<vector<inter_record> >& inter, query_tree* qt)
	{
		inter.resize(qg->g.size());
		for (int i = 0; i < qg->g.size(); i++) {
			inter[i].resize(dg->g.size());
			for (int j = 0; j < dg->g.size(); j++)
				inter[i][j].inter = 0;
		}
		vector<unsigned int> visited_degree;
		visited_degree.resize(qg->g.size());
		for (int i = 0; i < visited_degree.size(); i++)
			visited_degree[i] = 0;
		vector<unsigned int> leafs;
		qt->get_leaf_nodes(leafs);
		queue<unsigned int> q;
		for (int i = 0; i < leafs.size(); i++)
		{
			unsigned int u = leafs[i];
			vector<unsigned int> vec;
			tf->get_candidates(u, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int v = vec[j];
				inter[u][v].inter = 1;
			}
			vec.clear();
			q.push(u);
		}
		while (!q.empty())
		{
			unsigned int u = q.front();
			q.pop();
			int parent = qt->tree[u]->parent;
			if (parent != -1)
			{
				visited_degree[parent]++;
				if (visited_degree[parent] == qt->tree[parent]->child.size())
					q.push(parent);
			}
			vector<unsigned int> vec;
			tf->get_candidates(u, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int v = vec[u];
				if (inter[u][v].inter == 0)
				{
					unsigned long long tmp = 1;
					for (auto iter = inter[u][v].child_inter.begin(); iter != inter[u][v].child_inter.end(); iter++)
						tmp = tmp * iter->second;
					inter[u][v].inter = tmp;
				}
				if (parent != -1) {
					vector<unsigned int> pvec;
					tf->get_neighbor(u, parent, v, pvec);
					for (int j = 0; j < pvec.size(); j++)
					{
						unsigned int p = pvec[j];
						if (inter[parent][p].child_inter.find(u) == inter[parent][p].child_inter.end())
							inter[parent][p].child_inter[u] = inter[u][v].inter;
						else
							inter[parent][p].child_inter[u] += inter[u][v].inter;
					}
					pvec.clear();
				}
			}
			vec.clear();
		}
		vector<unsigned int> vec;
		tf->get_candidates(qt->root, vec);
		unsigned long long total_inter = 0;
		for (int i = 0; i < vec.size(); i++)
			total_inter += inter[qt->root][vec[i]].inter;
		return total_inter;
	}

	unsigned long long delete_node(TF_index* tf, vector<vector<inter_record> >& inter, query_tree* qt, unsigned int u)
	{
		queue<unsigned int> q;
		int parent = qt->tree[u]->parent;
		if (parent == -1)
			return 0;
		vector<unsigned int> vec;
		tf->get_candidates(parent, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int p = vec[i];
			inter[parent][p].child_inter[u] = 1;	
		}
		vec.clear();
		q.push(parent);
		while (!q.empty())
		{
			unsigned int u = q.front();
			q.pop();
			int parent = qt->tree[u]->parent;
			if (parent != -1)
			{
				vector<unsigned int> pvec;
				tf->get_candidates(parent, pvec);
				for (int i = 0; i < pvec.size(); i++)
					inter[parent][pvec[i]].child_inter[u]=0;
				pvec.clear();
			}
			vector<unsigned int> vec;
			tf->get_candidates(u, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int v = vec[u];
				unsigned long long tmp = 1;
				for (auto iter = inter[u][v].child_inter.begin(); iter != inter[u][v].child_inter.end(); iter++)
					tmp = tmp * iter->second;
				inter[u][v].inter = tmp;
				if (parent != -1) {
					vector<unsigned int> pvec;
					tf->get_neighbor(u, parent, v, pvec);
					for (int j = 0; j < pvec.size(); j++)
					{
						unsigned int p = pvec[j];
						inter[parent][p].child_inter[u] += inter[u][v].inter;
					}
					pvec.clear();
				}
			}
			vec.clear();
		}
		tf->get_candidates(qt->root, vec);
		unsigned long long total_inter = 0;
		for (int i = 0; i < vec.size(); i++)
			total_inter += inter[qt->root][vec[i]].inter;
		return total_inter;
	}
	unsigned long long try_delete(TF_index* tf, vector<vector<inter_record> >& inter, query_tree* qt, unsigned int u)
	{
		vector<vector<inter_record> > tmp_inter;
		tmp_inter.resize(inter.size());
		for (int i = 0; i < inter.size(); i++)
		{
			tmp_inter[i].resize(inter[i].size());
			for (int j = 0; j < inter[i].size(); j++)
				tmp_inter[i][j].copy(inter[i][j]);
		}
		unsigned long long total = delete_node(tf, tmp_inter, qt, u);
		for (int i = 0; i < tmp_inter.size(); i++)
			tmp_inter[i].clear();
		tmp_inter.clear();
		return total;
	}

	void get_global_order(TF_index* tf, query_tree* qt, vector<unsigned int>& order)
	{
		query_tree* tmp_qt = new query_tree;
		tmp_qt->copy(qt);
		vector<vector<inter_record> > inter;
		vector<unsigned int> reverse_order;
		compute_inter(tf, inter, qt);
		while (reverse_order.size()<qg->g.size())
		{
			vector<unsigned int> leafs;
			tmp_qt->get_leaf_nodes(leafs);
			unsigned long long min_inter = 0xFFFFFFFFFFFFFFFF;
			int min_leaf = -1;
			if (leafs.size() == 1)
				min_leaf = leafs[0];
			else {
				for (int i = 0; i < leafs.size(); i++)
				{
					unsigned int u = leafs[i];
					unsigned long long total = try_delete(tf, inter, qt, u);
					if (total < min_inter)
					{
						min_leaf = u;
						min_inter = total;
					}
				}
			}
			delete_node(tf, inter, qt, min_leaf);
			tmp_qt->delete_tree_node(min_leaf);
			reverse_order.push_back(min_leaf);
		}
		for (int i = reverse_order.size() - 1; i >= 0; i--)
			order.push_back(reverse_order[i]);
		reverse_order.clear();
		for (int i = 0; i < inter.size(); i++)
			inter[i].clear();
		inter.clear();
	}
	double initialization()
	{
		TF_index* tf = new TF_index(dg, qg, true);
		tf->construct_index();
		qt->copy(tf->qt);
		vector<unsigned int> order;
		get_global_order(tf, qt, order);
		global_mo = new matching_order;
		global_mo->pos.resize(qg->g.size());
		for (int i = 0; i < order.size(); i++)
		{
			global_mo->order.push_back(order[i]);
			global_mo->pos[order[i]] = i;
		}
		for (int i = 0; i < order.size(); i++)
		{
			unsigned int u = global_mo->order[i];
			vector<neighbor_info> vec;
			qg->get_neighbor(u, vec);
			int left_pos = global_mo->order.size();
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int u2 = vec[j].id;
				if (global_mo->pos[u2] < i)
				{
					global_mo->left_neighbors[u].push_back(make_pair(u2, vec[j].e_label));
					if (global_mo->pos[u2] < left_pos)
						left_pos = global_mo->pos[u2];
				}
			}
			if (left_pos < global_mo->order.size())
				global_mo->parent.push_back(global_mo->order[left_pos]);
			else
				global_mo->parent.push_back(-1);

		}
		delete tf;
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

	void query_dfs(unordered_set<unsigned int>& used, unsigned int index, unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates,  unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		clock_t currect_clock = clock();
		if ((currect_clock - start_clock) / CLOCKS_PER_SEC > time_limit)
			return;
		unsigned int u = global_mo->order[index];

		vector<neighbor_info> vec;
		qg->get_neighbor(u, vec);
		if (tmp_maping.find(u) != tmp_maping.end()) // this means u has been visited in the backtrack to the tree root, and we need to varify non tree edges. 
		{
			unsigned int v = tmp_maping[u];
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int u1 = vec[i].id;
				if (tmp_maping.find(u1) != tmp_maping.end())
				{
					unsigned int v1 = tmp_maping[u1];
					if(!dg->check_edge_label(v1, v, vec[i].e_label))
						return;
				}
			}
			if (index + 1 == global_mo->order.size())
			{
				if (record_match) {
					for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
						match[iter->first].matched_data.push_back(iter->second);
				}
				match_size++;
			}
			else {
//				inter_result++;
				query_dfs(used, index + 1, tmp_maping, candidates, match, match_size, record_match);
			}
		}
		else
		{
			vector<vector<unsigned int> > tmp_vecs;
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int u1 = vec[i].id;
				if (tmp_maping.find(u1) != tmp_maping.end())
				{
					unsigned int v1 = tmp_maping[u1];
					vector<unsigned int> cand_vec;
					if (use_global_index)
						idx->get_neighbor(u1, u, v1, cand_vec);
					else
						dg->get_labeled_neighbor(v1, qg->g[u].label, vec[i].e_label, cand_vec);
					if (cand_vec.empty())
					{
						for (int j = 0; j < tmp_vecs.size(); j++)
							tmp_vecs[j].clear();
						return;
					}
					tmp_vecs.push_back(cand_vec);
				}
			}
			expanding_num++;
			leap_frog_join(tmp_vecs, candidates[u]);
			for (int i = 0; i < tmp_vecs.size(); i++)
			{
				neighbor_sum += tmp_vecs[i].size();
				tmp_vecs[i].clear();
			}
			tmp_vecs.clear();
			for (int i = 0; i < candidates[u].size(); i++)
			{
				unsigned int v = candidates[u][i];
				if (used.find(v) != used.end())
					continue;
				if (!use_global_index)
				{
					if (!NLF_check(v, u))
						continue;
				}
				filtered_expanding++;
				if (index + 1 == global_mo->order.size())
				{
					if (record_match) {
						for (unordered_map<unsigned int, unsigned int>::iterator iter = tmp_maping.begin(); iter != tmp_maping.end(); iter++)
							match[iter->first].matched_data.push_back(iter->second);
						match[u].matched_data.push_back(v);
					}
					match_size++;
				}
				else {
					inter_result++;
					tmp_maping[u] = v;
					used.insert(v);
					query_dfs(used, index + 1, tmp_maping, candidates, match, match_size, record_match);
					used.erase(v);
					tmp_maping.erase(u);
				}
			}
			candidates[u].clear();
		}
	}

	void tree_backtrack(unordered_set<unsigned int>& used, unsigned int u,  unordered_map<unsigned int, unsigned int>& tmp_maping,
		vector<vector<unsigned int>>& candidates, unordered_map<unsigned int, match_result>& match, unsigned int& match_size, bool record_match)
	{
		int parent = qt->tree[u]->parent;
		if (parent == -1)
			query_dfs(used, 0, tmp_maping, candidates, match, match_size, record_match);
		else if (tmp_maping.find(parent) != tmp_maping.end())
			tree_backtrack(used, parent, tmp_maping, candidates, match, match_size, record_match);
		else
		{
			unsigned int v = tmp_maping[u];
			if (use_global_index)
				idx->get_neighbor(u, parent, v, candidates[parent]);
			else
				dg->get_labeled_neighbor(v, qg->g[parent].label, qt->tree[u]->p_elabel, candidates[parent]);
			for (int j = 0; j < candidates[parent].size(); j++)
			{
				unsigned int vp = candidates[parent][j];
				if (used.find(vp) != used.end())
					continue;
				if (!use_global_index)
				{
					if (!NLF_check(candidates[parent][j], parent))
						continue;
				}
				inter_result++;
				tmp_maping[parent] = vp;
				used.insert(vp);
				tree_backtrack(used, parent, tmp_maping, candidates, match, match_size, record_match);
				used.erase(vp);
				tmp_maping.erase(parent);
			}
			candidates[parent].clear();
		}
	}
	void query_processing(unsigned int src, unsigned int dst, unsigned int qsrc, unsigned int qdst, vector<unordered_map<unsigned int, match_result>>& match,
		unsigned int& match_size, bool record_match)
	{
		unsigned long long edge = merge_long_long(qsrc, qdst);
		unsigned long long reverse_edge = merge_long_long(qdst, qsrc);
		unsigned int u = qsrc;
		if (qt->tree_edges.find(edge) == qt->tree_edges.end() && qt->tree_edges.find(reverse_edge) != qt->tree_edges.end())
			u = qdst;
		vector<vector<unsigned int> >candidates;
		candidates.resize(qg->g.size());
		unordered_map<unsigned int, unsigned int> tmp_mapping;
		unordered_set<unsigned int> used;
		unordered_map<unsigned int, match_result> tmp_match;
		tmp_mapping[qsrc] = src;
		tmp_mapping[qdst] = dst;
		used.insert(src);
		used.insert(dst);
		tree_backtrack(used, u, tmp_mapping, candidates, tmp_match, match_size, record_match);
		if (record_match) {
			if (tmp_match.begin()->second.matched_data.size() > 0)
				match.push_back(tmp_match);
		}
		for (int i = 0; i < candidates.size(); i++)
			candidates[i].clear();
		candidates.clear();
		tmp_mapping.clear();
		used.clear();
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
