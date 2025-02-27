#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include<map>
#include<assert.h>
#include<algorithm>
#include"../utility/struct.h"
#include"../graph/data_graph.h"
#include"../graph/query_graph.h"
#define MAX_INT 0xFFFFFF
using namespace std;


bool my_sort(const pair<unsigned int, unsigned int>& p1, const pair<unsigned int, unsigned int>& p2)
{
	return p1.second < p2.second;
}
void merge_sorted_vec(vector<unsigned int>& vec1, vector<unsigned int>& vec2, vector<unsigned int>& joint)
{
	unsigned int cur1 = 0;
	unsigned int cur2 = 0;
	while (cur1 < vec1.size() && cur2 < vec2.size())
	{
		if (vec1[cur1] < vec2[cur2])
			cur1++;
		else if (vec1[cur1] == vec2[cur2])
			joint.push_back(vec1[cur1]);
		else
			cur2++;
	}
}

class local_index
{
public:
	local_index_template* idx;
	unsigned int target_src;
	unsigned int target_dst;
	unsigned int query_src;
	unsigned int query_dst;
	query_graph* qg;
	my_index* gd;
	data_graph* dg;
	indexing_order* io;
	bool use_gd;


	local_index(unsigned int s, unsigned int d, unsigned int query_s, unsigned int query_d, query_graph* qg_, my_index* gd_, data_graph* dg_, indexing_order* io_, bool use_gd_ = true)
	{
		target_src = s;
		target_dst = d;
		query_src = query_s;
		query_dst = query_d;
		qg = qg_;
		gd = gd_;
		dg = dg_;
		io = io_;
		idx = new local_index_template;
		idx->candidate_map.resize(qg->g.size());
		use_gd = use_gd_;
	}

	~local_index()
	{
		delete idx;
	}

	bool get_initial_candidates()
	{
		vector<unsigned int> src_neighbor_list, dst_neighbor_list;
		qg->get_unlabeled_neighbor(query_src, src_neighbor_list);
		qg->get_unlabeled_neighbor(query_dst, dst_neighbor_list);
		for (int i = 0; i < src_neighbor_list.size(); i++)
		{
			unsigned int neighbor_id = src_neighbor_list[i];
			if (neighbor_id == query_dst)
				continue;
			unsigned long long edge_pair = merge_long_long(query_src, neighbor_id);

			unsigned int set_size = gd->get_neighbor(target_src, edge_pair,  idx->candidate_map[neighbor_id]);
			if (set_size == 0)
				return false;

			vector<unsigned int>::iterator search_iter = lower_bound(idx->candidate_map[neighbor_id].begin(), idx->candidate_map[neighbor_id].end(), target_src);
			if (search_iter != idx->candidate_map[neighbor_id].end() && *search_iter == target_src)
				idx->candidate_map[neighbor_id].erase(search_iter);
			search_iter = lower_bound(idx->candidate_map[neighbor_id].begin(), idx->candidate_map[neighbor_id].end(), target_dst);
			if (search_iter != idx->candidate_map[neighbor_id].end() && *search_iter == target_dst)
				idx->candidate_map[neighbor_id].erase(search_iter);

			if (idx->candidate_map[neighbor_id].empty())
				return false;
		}

		for (int i = 0; i < dst_neighbor_list.size(); i++)
		{
			unsigned int neighbor_id = dst_neighbor_list[i];
			unsigned long long edge_pair = merge_long_long(query_dst, neighbor_id);

			if (neighbor_id == query_src)
				continue;

			if (!idx->candidate_map[neighbor_id].empty())
			{
				vector<unsigned int> vec;
				gd->get_neighbor(target_dst, edge_pair, vec);
				vector<unsigned int> joint;
				vector<unsigned int>::iterator iter = idx->candidate_map[neighbor_id].begin();
				for (int i = 0; i < vec.size(); i++)
				{
					iter = lower_bound(iter, idx->candidate_map[neighbor_id].end(), vec[i]);
					if (iter != idx->candidate_map[neighbor_id].end() && *iter == vec[i])
						joint.push_back(vec[i]); // we can directly push_back as vec is an ordered vector.
					if (iter == idx->candidate_map[neighbor_id].end())
						break;
				}
				if (joint.empty())
					return false;
				idx->candidate_map[neighbor_id].swap(joint);
				joint.clear();
			}
			else
			{
				unsigned int set_size = gd->get_neighbor(target_dst, edge_pair, idx->candidate_map[neighbor_id]);
				if (set_size == 0)
					return false;

				vector<unsigned int>::iterator search_iter = lower_bound(idx->candidate_map[neighbor_id].begin(), idx->candidate_map[neighbor_id].end(), target_src);
				if (search_iter != idx->candidate_map[neighbor_id].end() && *search_iter == target_src)
					idx->candidate_map[neighbor_id].erase(search_iter);
				search_iter = lower_bound(idx->candidate_map[neighbor_id].begin(), idx->candidate_map[neighbor_id].end(), target_dst);
				if (search_iter != idx->candidate_map[neighbor_id].end() && *search_iter == target_dst)
					idx->candidate_map[neighbor_id].erase(search_iter);

				if (idx->candidate_map[neighbor_id].empty())
					return false;
			}
		}
		return true;
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

	bool graph_get_initial_candidates()
	{
		vector<neighbor_info> src_neighbor_list, dst_neighbor_list;
		qg->get_neighbor(query_src, src_neighbor_list);
		qg->get_neighbor(query_dst, dst_neighbor_list);
		for (int i = 0; i < src_neighbor_list.size(); i++)
		{
			unsigned int neighbor_id = src_neighbor_list[i].id;
			if (neighbor_id == query_dst)
				continue;
			unsigned long long edge_pair = merge_long_long(query_src, neighbor_id);

			vector<unsigned int> vec;
			dg->get_labeled_neighbor(target_src, src_neighbor_list[i].v_label, src_neighbor_list[i].e_label, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				if (vec[i] != target_src && vec[i] != target_dst && NLF_check(vec[i], neighbor_id)) {
					idx->candidate_map[neighbor_id].push_back(vec[i]); // here we can directly push back as vec is sorted.
				}
			}
			if (idx->candidate_map[neighbor_id].empty())
				return false;
		}

		for (int i = 0; i < dst_neighbor_list.size(); i++)
		{
			unsigned int neighbor_id = dst_neighbor_list[i].id;
			unsigned long long edge_pair = merge_long_long(query_dst, neighbor_id);

			if (neighbor_id == query_src)
				continue;

			if (!idx->candidate_map[neighbor_id].empty())
			{
				vector<unsigned int> vec;
				dg->get_labeled_neighbor(target_dst, dst_neighbor_list[i].v_label, dst_neighbor_list[i].e_label, vec);
				vector<unsigned int> joint;
				for (int i = 0; i < vec.size(); i++)
				{
				//	if (NLF_check(vec[i], neighbor_id)) { // no need to carry out NLF check again, as candidates in the candidate map has been checked
						vector<unsigned int>::iterator iter = lower_bound(idx->candidate_map[neighbor_id].begin(), idx->candidate_map[neighbor_id].end(), vec[i]);
						if (iter != idx->candidate_map[neighbor_id].end() && *iter == vec[i])
							joint.push_back(vec[i]);
				//	}
				}
				if (joint.empty())
					return false;
				idx->candidate_map[neighbor_id].swap(joint);
			}
			else
			{
				vector<unsigned int> vec;
				dg->get_labeled_neighbor(target_dst, dst_neighbor_list[i].v_label, dst_neighbor_list[i].e_label, vec);
				for (int i = 0; i < vec.size(); i++)
				{
					if (vec[i] != target_src && vec[i] != target_dst && NLF_check(vec[i], neighbor_id)) {
						idx->candidate_map[neighbor_id].push_back(vec[i]);
					}
				}
				if (idx->candidate_map[neighbor_id].empty())
					return false;
			}
		}
		return true;
	}


	bool filter_candidates(unsigned int query_node_id)
	{
		for (int j = 0; j < io->left_neighbors[query_node_id].size(); j++)
		{
			unsigned int query_neighbor_id = io->left_neighbors[query_node_id][j];
			unsigned long long edge_pair = merge_long_long(query_node_id, query_neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(query_neighbor_id, query_node_id);
	//		cout << "exam left neighbor " << query_neighbor_id << endl;
			if (idx->candidate_map[query_node_id].empty())
			{
		//		cout << "query node not initialized " << endl;
				for (int i = 0; i < idx->candidate_map[query_neighbor_id].size(); i++)
				{
					unsigned int neighbor_cand = idx->candidate_map[query_neighbor_id][i];
					vector<unsigned int> neighbor_vec;
					gd->get_neighbor(neighbor_cand, reverse_edge_pair, neighbor_vec);
		//			vector<unsigned int>::iterator search_iter = idx->candidate_map[query_node_id].begin();
					for (int k = 0; k < neighbor_vec.size(); k++)
					{
						if (neighbor_vec[k] != target_src && neighbor_vec[k] != target_dst) {
							vector<unsigned int>::iterator search_iter = lower_bound(idx->candidate_map[query_node_id].begin(), idx->candidate_map[query_node_id].end(), neighbor_vec[k]);
							if (search_iter == idx->candidate_map[query_node_id].end() || *search_iter != neighbor_vec[k])
								idx->candidate_map[query_node_id].emplace(search_iter, neighbor_vec[k]);
						//	idx->edges[edge_pair].insert_edge(neighbor_vec[k], neighbor_cand);
							//idx->edges[reverse_edge_pair].insert_edge(neighbor_cand, neighbor_vec[k]);
						}
					}
				}
			}
			else {
		//		cout << "query node hase been initialized " << endl;
				vector<unsigned int> joint;
				bool* matched = new bool[idx->candidate_map[query_node_id].size()];
				for (int i = 0; i < idx->candidate_map[query_node_id].size(); i++)
					matched[i] = false;
				for (int i = 0; i < idx->candidate_map[query_neighbor_id].size(); i++) // we use the candidate set of its neighbor to filter this node, since its neighbor' canidate set  is smaller in size
				{
					unsigned int neighbor_cand = idx->candidate_map[query_neighbor_id][i];
					vector<unsigned int> neighbor_vec;
					gd->get_neighbor(neighbor_cand, reverse_edge_pair, neighbor_vec);
					if (neighbor_vec.empty())
							continue;
					vector<unsigned int>::iterator last_pos = neighbor_vec.begin();
					for (int k = 0; k < idx->candidate_map[query_node_id].size(); k++)
					{
						if (matched[k])
							continue;
						unsigned int v = idx->candidate_map[query_node_id][k];
						vector<unsigned int>::iterator search_iter = lower_bound(last_pos, neighbor_vec.end(), v);// here the left parameter can be the last search iter, as neighbor vec is sorted
						if (search_iter != neighbor_vec.end() && *search_iter == v) // find a neighbor in the candidates query_neighbor_id
						{
							matched[k] = true;
							joint.push_back(v);
						}
						last_pos == search_iter;
						if (last_pos == neighbor_vec.end())
							break;
					}
				}
				delete[] matched;
				sort(joint.begin(), joint.end());
				idx->candidate_map[query_node_id].swap(joint);
				joint.clear();
			}
			if (idx->candidate_map[query_node_id].empty())
				return false;
		}
		return true;
	}

	bool graph_filter_candidates(unsigned int query_node_id)
	{
		for (int j = 0; j < io->left_neighbors[query_node_id].size(); j++)
		{
			unsigned int query_neighbor_id = io->left_neighbors[query_node_id][j];
			int label = qg->g[query_node_id].label;
			int edge_label = qg->check_edge_label(query_neighbor_id, query_node_id);
			unsigned long long edge_pair = merge_long_long(query_node_id, query_neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(query_neighbor_id, query_node_id);
			if (idx->candidate_map[query_node_id].empty())
			{
				for (int i = 0; i < idx->candidate_map[query_neighbor_id].size(); i++)
				{
					unsigned int neighbor_cand = idx->candidate_map[query_neighbor_id][i];
					vector<unsigned int> neighbor_vec;
					dg->get_labeled_neighbor(neighbor_cand, label, edge_label, neighbor_vec);
					//unsigned int last_cand = 0;
					for (int k = 0; k < neighbor_vec.size(); k++)
					{
						if (NLF_check(neighbor_vec[k], query_node_id) && neighbor_vec[k] != target_src && neighbor_vec[k] != target_dst) {
							vector<unsigned int>::iterator search_iter = lower_bound(idx->candidate_map[query_node_id].begin(), idx->candidate_map[query_node_id].end(), neighbor_vec[k]);
							if (search_iter == idx->candidate_map[query_node_id].end() || *search_iter != neighbor_vec[k])
								idx->candidate_map[query_node_id].emplace(search_iter, neighbor_vec[k]);
						//	idx->candidate_map[query_node_id].push_back(neighbor_vec[k]);
						//	idx->edges[edge_pair].insert_edge(neighbor_vec[k], neighbor_cand);
							//idx->edges[reverse_edge_pair].insert_edge(neighbor_cand, neighbor_vec[k]);
							//last_cand = neighbor_vec[k];
						}
					}
				}
			}
			else {
				vector<unsigned int> joint;
				bool* matched = NULL;
				matched = new bool[idx->candidate_map[query_node_id].size()];
				for (int i = 0; i < idx->candidate_map[query_node_id].size(); i++)
					matched[i] = false;
				for (int i = 0; i < idx->candidate_map[query_neighbor_id].size(); i++) // we use the candidate set of its neighbor to filter this node, since its neighbor' canidate set  is smaller in size 
				{
					unsigned int neighbor_cand = idx->candidate_map[query_neighbor_id][i];
					vector<unsigned int> neighbor_vec;
					dg->get_labeled_neighbor(neighbor_cand, label, edge_label, neighbor_vec);
					vector<unsigned int>::iterator last_pos = neighbor_vec.begin();
					for (int k = 0; k < idx->candidate_map[query_node_id].size(); k++)
					{
						if (matched[k])
							continue;
						unsigned int v = idx->candidate_map[query_node_id][k];
						vector<unsigned int>::iterator search_iter = lower_bound(last_pos, neighbor_vec.end(), v);
						if (search_iter != neighbor_vec.end() && *search_iter == idx->candidate_map[query_node_id][k])
						{
							joint.push_back(idx->candidate_map[query_node_id][k]);
							matched[k] = true;
						}
							last_pos == search_iter;
							if (last_pos == neighbor_vec.end())
								break;
					}
				}
				sort(joint.begin(), joint.end());
				idx->candidate_map[query_node_id].swap(joint);
				delete[]matched;
				joint.clear();
			}
				if (idx->candidate_map[query_node_id].empty())
					return false;
			
		}
		return true;
	}

	void build_edge()
	{
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u = i;
			for (int j = 0; j < io->left_neighbors[u].size(); j++)
			{
				unsigned int u1 = io->left_neighbors[u][j];
				unsigned int edge_label = qg->check_edge_label(u, u1);
				unsigned int label = qg->g[u1].label;
				for (int k = 0; k < idx->candidate_map[u].size(); k++)
				{
					unsigned int v = idx->candidate_map[u][k];
					vector<unsigned int> vec;
					if (use_gd)
						gd->get_neighbor(u, u1, v, vec);
					else
						dg->get_labeled_neighbor(v, label, edge_label, vec);
					for (int m = 0; m < vec.size(); m++)
					{
						unsigned int v1 = vec[m];
						if (idx->check_candidate(u1, v1)) {
							idx->edges[merge_long_long(u, u1)].insert_edge(v, v1);
							idx->edges[merge_long_long(u1, u)].insert_edge(v1, v);
						}
					}

				}
			}
		}
	}

	bool construct_index()
	{
	//	cout << "begin construct" << endl;
		if (get_initial_candidates() == false)
			return false;
	//	cout << "get initial finished" << endl;
		vector<pair<unsigned int, unsigned int> > temporary_cardinality;
		for (int i = 0; i < idx->candidate_map.size(); i++)
		{
			if(!idx->candidate_map[i].empty())
				temporary_cardinality.push_back(make_pair(i, idx->candidate_map[i].size()));
		}
		sort(temporary_cardinality.begin(), temporary_cardinality.end(), my_sort);
		unordered_map<unsigned int, unsigned int> pos;
		for (int i = 0; i < temporary_cardinality.size(); i++)
			pos[temporary_cardinality[i].first] = i;

		for (int i = 0; i < temporary_cardinality.size(); i++)
		{
			unsigned int query_node_id = temporary_cardinality[i].first;
			vector<unsigned int> vec;
			qg->get_unlabeled_neighbor(query_node_id, vec);
			io->left_neighbors[query_node_id].clear();
			for (int j = 0; j < vec.size(); j++)
			{
				if (pos.find(vec[j]) != pos.end() && pos[vec[j]] < i)
					io->left_neighbors[query_node_id].push_back(vec[j]);
			}
		//	cout << "filtering candidates of " << query_node_id << endl;
			if(!filter_candidates(query_node_id))
				return false;
		}
		for (int i = 0; i < io->order_of_remain.size(); i++)
		{
			unsigned int query_node_id = io->order_of_remain[i];
		//	cout << "filtering candidates of " << query_node_id << endl;
			if (!filter_candidates(query_node_id))
				return false;
			pos[query_node_id] = pos.size();
		}
	//	cout << "building edge" << endl;
		build_edge();
		return true;
	}

	bool graph_construct_index()
	{
		if (graph_get_initial_candidates() == false)
			return false;
		vector<pair<unsigned int, unsigned int> > temporary_cardinality;
		for (int i = 0; i < idx->candidate_map.size(); i++)
		{
			if (!idx->candidate_map[i].empty())
				temporary_cardinality.push_back(make_pair(i, idx->candidate_map[i].size()));
		}
		sort(temporary_cardinality.begin(), temporary_cardinality.end(), my_sort);
		unordered_map<unsigned int, unsigned int> pos;
		for (int i = 0; i < temporary_cardinality.size(); i++)
			pos[temporary_cardinality[i].first] = i;

		for (int i = 0; i < temporary_cardinality.size(); i++)
		{
			unsigned int query_node_id = temporary_cardinality[i].first;
			vector<unsigned int> vec;
			qg->get_unlabeled_neighbor(query_node_id, vec);
			io->left_neighbors[query_node_id].clear();
			for (int j = 0; j < vec.size(); j++)
			{
				if (pos.find(vec[j]) != pos.end() && pos[vec[j]] < i)
					io->left_neighbors[query_node_id].push_back(vec[j]);
			}
			if (!graph_filter_candidates(query_node_id))
				return false;
		}
		for (int i = 0; i < io->order_of_remain.size(); i++)
		{
			unsigned int query_node_id = io->order_of_remain[i];
			if (!graph_filter_candidates(query_node_id))
				return false;
			pos[query_node_id] = pos.size();
		}
		build_edge();
		return true;
	}

	bool check_candidate(unsigned int qid, unsigned int cand)
	{
		vector<unsigned int>::iterator vec_iter = lower_bound(idx->candidate_map[qid].begin(), idx->candidate_map[qid].end(), cand);
		if (vec_iter == idx->candidate_map[qid].end() || *vec_iter != cand)
			return false;
		else
			return true;
	}


	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		idx->get_candidates(id, vec);
	}

	void graph_get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
		int edge_label = qg->check_edge_label(src_qid, dst_qid);
		int dst_label = qg->g[dst_qid].label;
		vector<unsigned int> vec;
		dg->get_labeled_neighbor(candidate_id, dst_label, edge_label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			if (idx->check_candidate(dst_qid, vec[i]))
				neighbor_list.push_back(vec[i]);
		}
	}

	void get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		unordered_map<unsigned long long, edge_view>::iterator iter = idx->edges.find(edge_pair);
		if (iter != idx->edges.end())
		{
			vector<unsigned int> vec;
			iter->second.get_neighbor(candidate_id, vec);
			for (int i = 0; i < vec.size(); i++)
			{
					neighbor_list.push_back(vec[i]);
			}
		}
	}

	void get_neighbor( unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
		unsigned int dst_qid = (edge_pair & 0xFFFFFFFF);
		unordered_map<unsigned long long, edge_view>::iterator iter = idx->edges.find(edge_pair);
		if (iter != idx->edges.end())
		{
			vector<unsigned int> vec;
			iter->second.get_neighbor(candidate_id, vec);
			for (int i = 0; i < vec.size(); i++)
			{
					neighbor_list.push_back(vec[i]);
			}
		}
	}
	bool check_edge(unsigned int src_qid, unsigned int dst_qid, unsigned int src, unsigned int dst)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		if (idx->edges.find(edge_pair) != idx->edges.end() && idx->edges[edge_pair].check_edge(src, dst))
			return true;
		else
			return false;
	}
};

