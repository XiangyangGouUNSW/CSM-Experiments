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
#include"../graph/DAG-query.h"
#include "../utility/struct.h"
#include<assert.h>
#include<set>
#include<algorithm>


// we add NLF check in other indexes, but not in CaliG, because the injective matching of CaliG is stronger than NLF.

class CaliG_index
{
public:
	global_index_template* idx;
	unordered_map<unsigned long long, edge_view> in_edges;
	bool use_edge_view;
	vector<vector<bipartite_graph*>> bi_idx;
	unsigned int node_num;
	data_graph* dg;
	query_graph* qg;

	CaliG_index(data_graph* d, query_graph* q, bool use_edge_view_ = false)
	{
		dg = d;
		qg = q;
		node_num = d->vertice_number;
		idx = new global_index_template;
		idx->candidate_map.resize(node_num);
		for (int i = 0; i < node_num; i++) {
			idx->candidate_map[i] = 0;
		}
		bi_idx.resize(qg->g.size());
		for (int i = 0; i < qg->g.size(); i++) {
			bi_idx[i].resize(node_num);
			for(int j=0;j<node_num;j++)
				bi_idx[i][j] = NULL;
		}

		use_edge_view = use_edge_view_;
	}
	~CaliG_index()
	{
		delete idx;
		for (int i = 0; i < bi_idx.size(); i++)
		{
			for (int j = 0; j < bi_idx[i].size(); j++)
			{
				if (bi_idx[i][j])
					delete bi_idx[i][j];
			}
		}
		in_edges.clear();
	}

	unsigned int memory_compute()
	{
		unsigned int mem = 0;
		mem += sizeof(void*) * 3 + sizeof(unsigned int) + sizeof(bool) + sizeof(in_edges) + sizeof(bi_idx);
		mem = byte_align(mem, 8);

		mem += in_edges.bucket_count() * sizeof(void*);
		unsigned int kvsize = sizeof(unsigned long long) + sizeof(edge_view);
		kvsize = byte_align(kvsize, 8);
		mem += in_edges.size() * (kvsize + sizeof(void*));
		for (auto iter = in_edges.begin(); iter != in_edges.end(); iter++)
			mem += iter->second.memory_compute();

		mem += bi_idx.capacity() * sizeof(vector<bipartite_graph*>);
		for (int i = 0; i < bi_idx.size(); i++) {
			mem += bi_idx[i].capacity() * sizeof(void*);
			for (int j = 0; j < bi_idx[i].size(); j++)
			{
				if (bi_idx[i][j])
					mem += bi_idx[i][j]->memory_compute();
			}
		}
		return mem;
	}
	void delete_and_check(unsigned int query_node, unsigned int data_node)
	{
		queue<pair<unsigned int, unsigned int>> q;
		q.push(make_pair(query_node, data_node));
		while (!q.empty()) {
			pair<unsigned int, unsigned int> p = q.front();
			q.pop();
			unsigned int u = p.first;
			unsigned int v = p.second;
		//	cout << "function delete and check " << u << ' ' << v << endl;
			if (use_edge_view)
			{
				vector<unsigned int> vec;
				qg->get_unlabeled_neighbor(u, vec);
				for (int i = 0; i < vec.size(); i++)
				{
					unsigned int uj = vec[i];
					unsigned long long edge = merge_long_long(u, uj);
					unsigned long long reverse_edge = merge_long_long(uj, u);
					vector<unsigned int> data_vec;
					idx->edges[edge].get_neighbor(v, data_vec);
					for (int j = 0; j < data_vec.size(); j++)
					{
						unsigned int vj = data_vec[j];
						if (idx->check_candidate(uj, vj))
						{
							idx->edges[edge].delete_edge(v, vj);
							in_edges[reverse_edge].delete_edge(vj, v);
							if (bi_idx[uj][vj]->delete_edge(u, v)) {
								if (!bi_idx[uj][vj]->injective_match())
								{
									idx->erase_candidate(uj, vj);
									//delete_and_check(uj, vj);
									q.push(make_pair(uj, vj));
								}
							}
						}
					}
					data_vec.clear();
				}
			}
			else
			{
				for (auto iter = bi_idx[u][v]->cand_map.begin(); iter != bi_idx[u][v]->cand_map.end(); iter++)
				{
					unsigned int uj = iter->first;
					for (int i = 0; i < iter->second.size(); i++)
					{
						unsigned int vj = iter->second[i];
						if (idx->check_candidate(uj, vj))
						{
							if (bi_idx[uj][vj]->delete_edge(u, v))
							{
								if (!bi_idx[uj][vj]->injective_match())
								{
									idx->erase_candidate(uj, vj);
									//delete_and_check(uj, vj);
									q.push(make_pair(uj, vj));
								}
							}
						}
					}
				}
			}
		}
	//	cout << "delete and check finished " << endl;
	}

	void construct_index()
	{
	//	cout << "constructing index" << endl;
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u = i;
			vector<neighbor_info> neighbors;
			qg->get_neighbor(u, neighbors);
			unsigned int degree = neighbors.size();
			vector<unsigned int> vec;
			dg->get_labeled_node(qg->g[u].label, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int v = vec[j];
				idx->set_candidate(u, v);
				bi_idx[u][v] = new bipartite_graph(degree);
				for (int k = 0; k < neighbors.size(); k++)
				{
					unsigned int uj = neighbors[k].id;
					unsigned long long edge = merge_long_long(u, uj);
					unsigned long long reverse_edge = merge_long_long(uj, u);
					vector<unsigned int> data_vec;
					dg->get_labeled_neighbor(v, neighbors[k].v_label, neighbors[k].e_label, data_vec);
					for (int m = 0; m < data_vec.size(); m++)
					{
						unsigned int vj = data_vec[m];
						bi_idx[u][v]->add_edge(uj, vj);
						if (use_edge_view) {
							idx->edges[edge].insert_edge(v, vj);
							in_edges[reverse_edge].insert_edge(vj, v);
						}
					}
				}
			}
		}
	//	cout << "initialization finished" << endl;

		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u = i;
			vector<unsigned int> vec;
			idx->get_candidates(u, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int v = vec[j];
				if (!idx->check_candidate(u, v))
					continue;
				if (!bi_idx[u][v]->injective_match()) {
				//	cout << "erasing candidate " << u << ' ' << v << endl;
					idx->erase_candidate(u, v);
					delete_and_check(u, v);
				//	cout << "erasing finished " << endl;
				}
			}
		}
	//	cout << "construct finished" << endl;
	}

	void delete_edge(unsigned int src, unsigned int src_label, unsigned int dst, unsigned int dst_label, unsigned int edge_label)
	{
		unsigned int vi = src;
		unsigned int vj = dst;
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		vector<unsigned int> src_check;
		vector<unsigned int> dst_check;
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int ui = vec[i].first;
			unsigned int uj = vec[i].second;
			if (bi_idx[ui][vi]->delete_edge(uj, vj))
				src_check.push_back(ui);
			if (bi_idx[uj][vj]->delete_edge(ui, vi))
				dst_check.push_back(uj);
		}
		for (int i = 0; i < src_check.size(); i++)
		{
			unsigned int ui = src_check[i];
			if (idx->check_candidate(ui, vi)&&!bi_idx[ui][vi]->injective_match()) {
				idx->erase_candidate(ui, vi);
				delete_and_check(ui, vi);
			}
		}
		for (int i = 0; i < dst_check.size(); i++)
		{
			unsigned int uj = dst_check[i];
			if (idx->check_candidate(uj, vj) && !bi_idx[uj][vj]->injective_match()) {
				idx->erase_candidate(uj, vj);
				delete_and_check(uj, vj);
			}
		}
	}

	void add_and_check(unsigned int query_node, unsigned int data_node, vector<pair<unsigned int, unsigned int>>& failed)
	{
		queue<pair<unsigned int, unsigned int>> q;
		q.push(make_pair(query_node, data_node));
		while (!q.empty()) {
			pair<unsigned int, unsigned int> p = q.front();
			q.pop();
			unsigned int u = p.first;
			unsigned int v = p.second;;
			vector<unsigned int> vec;
			qg->get_unlabeled_neighbor(u, vec);
			if (use_edge_view)
			{
				for (int i = 0; i < vec.size(); i++)
				{
					unsigned int uj = vec[i];
					unsigned long long edge = merge_long_long(u, uj);
					unsigned long long reverse_edge = merge_long_long(uj, u);
					vector<unsigned int> data_vec;
					in_edges[edge].get_neighbor(v, data_vec);
					for (int j = 0; j < data_vec.size(); j++)
					{
						unsigned int vj = data_vec[j];
						idx->edges[edge].insert_edge(v, vj);
						in_edges[reverse_edge].insert_edge(vj, v);
						bi_idx[uj][vj]->add_edge(u, v);
						if (!idx->check_candidate(uj, vj))
						{
							if (bi_idx[uj][vj]->injective_match())
							{
								idx->set_candidate(uj, vj);
								q.push(make_pair(uj, vj));
							}
							else
								failed.push_back(make_pair(uj, vj));
						}
					}
					data_vec.clear();
				}
			}
			else
			{
				for (auto iter = bi_idx[u][v]->cand_map.begin(); iter != bi_idx[u][v]->cand_map.end(); iter++)
				{
					unsigned int uj = iter->first;

					for (int i = 0; i < iter->second.size(); i++)
					{
						unsigned int vj = iter->second[i];
						bi_idx[uj][vj]->add_edge(u, v);
						if (!idx->check_candidate(uj, vj))
						{
							if (bi_idx[uj][vj]->injective_match())
							{
								idx->set_candidate(uj, vj);
								q.push(make_pair(uj, vj));
							}
							else
								failed.push_back(make_pair(uj, vj));
						}
					}
				}
			}
			vec.clear();
		}
	}
	

	int get_neighbor(unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
		unsigned int dst_qid = (edge_pair & 0xFFFFFFFF);
		unsigned int src_qid = (edge_pair >> 32);
		if (bi_idx[src_qid][candidate_id] && bi_idx[src_qid][candidate_id]->cand_map.find(dst_qid) != bi_idx[src_qid][candidate_id]->cand_map.end()) {
			
			neighbor_list = bi_idx[src_qid][candidate_id]->cand_map[dst_qid];
		/*	for (int i = 0; i < bi_idx[src_qid][candidate_id]->cand_map[dst_qid].size(); i++)
			{
				unsigned int v = bi_idx[src_qid][candidate_id]->cand_map[dst_qid][i];
				if (idx->check_candidate(dst_qid, v))
					neighbor_list.push_back(v);
				else
					cout << "error neighbor" << endl;
			}*/
		}
		return neighbor_list.size();
	}

	int get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
		if (bi_idx[src_qid][candidate_id] && bi_idx[src_qid][candidate_id]->cand_map.find(dst_qid) != bi_idx[src_qid][candidate_id]->cand_map.end())
		{
			neighbor_list = bi_idx[src_qid][candidate_id]->cand_map[dst_qid];
		/*	for (int i = 0; i < bi_idx[src_qid][candidate_id]->cand_map[dst_qid].size(); i++)
			{
				unsigned int v = bi_idx[src_qid][candidate_id]->cand_map[dst_qid][i];
				if (idx->check_candidate(dst_qid, v))
					neighbor_list.push_back(v);
				else
					cout << "error neighbor" << endl;
			}*/
		}
		return neighbor_list.size();
	}
	void insert_edge(unsigned int src, unsigned int src_label, unsigned int dst, unsigned int dst_label, unsigned int edge_label)
	{
		unsigned int vi = src;
		unsigned int vj = dst;
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int ui = vec[i].first;
			unsigned int uj = vec[i].second;
			if (!bi_idx[ui][vi])
				bi_idx[ui][vi] = new bipartite_graph(qg->get_total_degree(ui));
			if (!bi_idx[uj][vj])
				bi_idx[uj][vj] = new bipartite_graph(qg->get_total_degree(uj));
		}
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int ui = vec[i].first;
			unsigned int uj = vec[i].second;
			unsigned long long edge = merge_long_long(ui, uj);
			unsigned long long reverse_edge = merge_long_long(uj, ui);
		//	if (!bi_idx[ui][vi])
	//			bi_idx[ui][vi] = new bipartite_graph(qg->get_total_degree(ui));
		//	if (!bi_idx[uj][vj])
		//		bi_idx[uj][vj] = new bipartite_graph(qg->get_total_degree(uj));
			
			if (use_edge_view)
			{
				idx->edges[edge].insert_edge(vi, vj);
				in_edges[reverse_edge].insert_edge(vj, vi);
			}
			bi_idx[uj][vj]->add_edge(ui, vi);
			if (!idx->check_candidate(uj, vj)) {
				if (bi_idx[uj][vj]->injective_match())
					idx->set_candidate(uj, vj);
			}
			if (idx->check_candidate(uj, vj)) {
				vector<pair<unsigned int, unsigned int>> failed_set;
				add_and_check(uj, vj, failed_set);
				for (int j = 0; j < failed_set.size(); j++) {
					if(!idx->check_candidate(failed_set[j].first, failed_set[j].second))
						delete_and_check(failed_set[j].first, failed_set[j].second);
				}
			}
			else // check anoter direction, otherwise there will be error in the case of circle query 
			{
				if (use_edge_view)
				{
					idx->edges[reverse_edge].insert_edge(vj, vi);
					in_edges[edge].insert_edge(vi, vj);
				}
				bi_idx[ui][vi]->add_edge(uj, vj);
				if (!idx->check_candidate(ui, vi)) {
					if (bi_idx[ui][vi]->injective_match())
						idx->set_candidate(ui, vi);
				}
				if (idx->check_candidate(ui, vi)) {
					vector<pair<unsigned int, unsigned int>> failed_set;
					add_and_check(ui, vi, failed_set);
					for (int j = 0; j < failed_set.size(); j++) {
						if (!idx->check_candidate(failed_set[j].first, failed_set[j].second))
							delete_and_check(failed_set[j].first, failed_set[j].second);
					}
				}
			}
		}
		vec.clear();
	//	cout << "inserting finished" << endl;
	}
	
	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		idx->get_candidates(id, vec);
	}
	bool check_candidate(unsigned int id, unsigned int candidate_id)
	{
		return idx->check_candidate(id, candidate_id);
	}
	pair<unsigned int, unsigned int> get_size()
	{
		unsigned int vertex_sum = 0;
		vector<unsigned int> tmp_vec;
		for (int i = 0; i < qg->g.size(); i++)
		{
			idx->get_candidates(i, tmp_vec);
			vertex_sum += tmp_vec.size();
			tmp_vec.clear();
		}
		unsigned int edge_sum = 0;
		for (auto iter = idx->edges.begin(); iter != idx->edges.end(); iter++) {
			edge_sum += iter->second.get_size();
		//	cout << "edge view " << (iter->first >> 32) << ' ' << (iter->first & 0xFFFFFFFF) << " size " << iter->second.get_size() << endl;
		}
		edge_sum = edge_sum / 2; // CaliG store directed edges. Each undirect edge in the data graph will be tranformed into 2 directed edges, and some of them will be removed later, we 
		// divide the number of edges by 2, so that we can directly compare it with Eq*Eg to evaluate the filtering ability of the edge view 
		return make_pair(vertex_sum, edge_sum);
	}
};
