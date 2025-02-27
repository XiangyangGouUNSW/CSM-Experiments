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
#include<assert.h>
#include<algorithm>
using namespace std;


class NLF_index
{
public:
	data_graph* dg;
	query_graph* qg;
	global_index_template* idx;
	bool use_edge_view;
	unsigned int node_num;


	NLF_index(data_graph* d, query_graph* q, bool use_edge_view_ = false)
	{
		dg = d;
		qg = q;
		node_num = d->vertice_number;
		idx = new global_index_template;
		idx->candidate_map.resize(node_num);
		for (int i = 0; i < node_num; i++)
			idx->candidate_map[i] = 0;
		use_edge_view = use_edge_view_;
	}

	~NLF_index()
	{
		delete idx;
	}

	unsigned int memory_compute()
	{
		unsigned int mem = byte_align(sizeof(dg) + sizeof(qg) + sizeof(idx) + sizeof(use_edge_view) + sizeof(node_num), 8);
		mem += idx->memory_compute(); // memory of query graph and data graph will not be included.
		return mem;
	}
	void resize()
	{
		node_num = dg->vertice_number;
		idx->candidate_map.resize(node_num);
	}
	void add_edge_with_cand(unsigned int query_node, unsigned int new_candidate)
	{
		unsigned int id = query_node;
		vector<neighbor_info> vec;
		qg->get_neighbor(id, vec);
		for (int j = 0; j < vec.size(); j++)
		{
			int neighbor_label = vec[j].v_label;
			unsigned int neighbor_id = vec[j].id;
			int edge_label = vec[j].e_label;
			unsigned long long edge_pair = merge_long_long(id, neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(neighbor_id, id);
			unsigned int candidate_id = new_candidate;
			vector<unsigned int> neigbor_vec;
			dg->get_labeled_neighbor(candidate_id, neighbor_label, edge_label, neigbor_vec);
			for (int m = 0; m < neigbor_vec.size(); m++)
			{
				if(idx->check_candidate(neighbor_id, neigbor_vec[m]))
				{
					idx->edges[edge_pair].insert_edge(candidate_id, neigbor_vec[m]);
					idx->edges[reverse_edge_pair].insert_edge(neigbor_vec[m], candidate_id);
				}
			}
		}
	}

	void delete_edge_with_cand(unsigned int query_node, unsigned int del_candidate)
	{
		unsigned int id = query_node;
		vector<neighbor_info> vec;
		qg->get_neighbor(id, vec);
		for (int j = 0; j < vec.size(); j++)
		{
			int neighbor_label = vec[j].v_label;
			unsigned int neighbor_id = vec[j].id;
			int edge_label = vec[j].e_label;
			unsigned long long edge_pair = merge_long_long(id, neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(neighbor_id, id);
			unsigned int candidate_id = del_candidate;
			vector<unsigned int> neigbor_vec;
			unordered_map<unsigned long long, edge_view>::iterator iter = idx->edges.find(edge_pair);
			if (iter == idx->edges.end())
				continue;
			iter->second.get_neighbor(candidate_id, neigbor_vec);
			for (int m = 0; m < neigbor_vec.size(); m++)
					idx->edges[reverse_edge_pair].delete_edge(neigbor_vec[m], candidate_id);
			idx->edges[edge_pair].delete_node(candidate_id);
		}
	}

	bool NLF_check(unsigned int query_node, unsigned int data_node)
	{
		for (map<unsigned long long, int>::iterator nlf_iter = qg->NLF[query_node].begin(); nlf_iter != qg->NLF[query_node].end(); nlf_iter++)
		{
			if (dg->NLF[data_node].find(nlf_iter->first) == dg->NLF[data_node].end() || dg->NLF[data_node][nlf_iter->first] < nlf_iter->second)
			{
				return false;
			}
		}
		return true;
	}

	bool check_candidate(unsigned int qid, unsigned int cand)
	{
		return idx->check_candidate(qid, cand);
	}
	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);

		for (int i = 0; i < vec.size(); i++)  // if src or dst is not candidates of a query node before and it is now, this query node must be endpoint of an edge with (src_label, edge_label, dst_label), otherwise inserting this edge will not influece the NLF
		{
			unsigned int src_qid = vec[i].first;
			unsigned int dst_qid = vec[i].second;
			bool src_cand = idx->check_candidate(src_qid, src);
			bool dst_cand = idx->check_candidate(dst_qid, dst);
			if (src_cand && dst_cand && use_edge_view)
			{
				unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
				unsigned long long reverse_edge_pair = merge_long_long(dst_qid, src_qid);
				idx->edges[edge_pair].insert_edge(src, dst);
				idx->edges[reverse_edge_pair].insert_edge(dst, src);
			}
			if (!src_cand){
				if (NLF_check(src_qid, src))
				{
					idx->set_candidate(src_qid, src);
					if (use_edge_view)
						add_edge_with_cand(src_qid, src);
				}
			}

			if(!dst_cand)
			{
				if (NLF_check(dst_qid, dst))
				{
					idx->set_candidate(dst_qid, dst);
					if(use_edge_view)
						add_edge_with_cand(dst_qid, dst);
				}
			}
		}
	}

	void delete_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);

		for (int i = 0; i < vec.size(); i++)  // if src or dst is not candidates of a query node before and it is now, this query node must be endpoint of an edge with (src_label, edge_label, dst_label), otherwise inserting this edge will not influece the NLF
		{
			unsigned int src_qid = vec[i].first;
			unsigned int dst_qid = vec[i].second;
			bool src_cand = idx->check_candidate(src_qid, src);
			bool dst_cand = idx->check_candidate(dst_qid, dst);
			if (src_cand && dst_cand && use_edge_view)
			{
				unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
				unsigned long long reverse_edge_pair = merge_long_long(dst_qid, src_qid);
				idx->edges[edge_pair].delete_edge(src, dst);
				idx->edges[reverse_edge_pair].delete_edge(dst, src);
			}
			if (src_cand) {
				if (!NLF_check(src_qid, src))
				{
					idx->erase_candidate(src_qid, src);
					if (use_edge_view)
						delete_edge_with_cand(src_qid, src);
				}
			}

			if (dst_cand)
			{
				if (!NLF_check(dst_qid, dst))
				{
					idx->erase_candidate(dst_qid, dst);
					if (use_edge_view)
						delete_edge_with_cand(dst_qid, dst);
				}
			}
		}
	}
	void construct_index()
	{
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int id = i;
			int label = qg->g[id].label;
			vector<unsigned int> candidates;
			dg->get_labeled_node(label, candidates);
			for (int j = 0; j < candidates.size(); j++)
			{
				unsigned int candidate_id = candidates[j];
				if (NLF_check(id, candidate_id))
					idx->set_candidate(id, candidate_id);
			}
		}
		if (use_edge_view) {
			for (int i = 0; i < qg->g.size(); i++)
			{
				unsigned int id = i;
				vector<neighbor_info> neighbor_vec;
				qg->get_neighbor(id, neighbor_vec);
				vector<unsigned int> cand_vec;
				idx->get_candidates(id, cand_vec);
				for (int j = 0; j < neighbor_vec.size(); j++)
				{
					int neighbor_label = neighbor_vec[j].v_label;
					unsigned int neighbor_id = neighbor_vec[j].id;
					int edge_label = neighbor_vec[j].e_label;
					unsigned long long label_pair = merge_long_long(neighbor_label, edge_label);
					unsigned long long edge_pair = merge_long_long(id, neighbor_id);
					for (int k = 0; k < cand_vec.size(); k++)
					{
						unsigned int candidate_id = cand_vec[k];
						vector<unsigned int> cand_neigbor_vec;
						dg->get_labeled_neighbor(candidate_id, neighbor_label, edge_label, cand_neigbor_vec);
						for (int m = 0; m < cand_neigbor_vec.size(); m++)
						{
							if (idx->check_candidate(neighbor_id, cand_neigbor_vec[m]))
								idx->edges[edge_pair].insert_edge(candidate_id, cand_neigbor_vec[m]);
						}
						cand_neigbor_vec.clear();
					}
				}
			}
		}
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
		edge_sum = edge_sum / 2; //each edge has two copy of views, here we only return the number of edges;
		return make_pair(vertex_sum, edge_sum);
	}

		int get_neighbor(unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
			if (use_edge_view) {
				unordered_map<unsigned long long, edge_view>::iterator iter = idx->edges.find(edge_pair);
				if (iter != idx->edges.end())
					iter->second.get_neighbor(candidate_id, neighbor_list);
			}
			else
			{
				unsigned int dst_qid = (edge_pair & 0xFFFFFFFF);
				unsigned int src_qid = (edge_pair >> 32);
				int e_label = qg->check_edge_label(src_qid, dst_qid);
				assert(e_label != -1);
				int dst_label = qg->g[dst_qid].label;
				vector<unsigned int> vec;
				dg->get_labeled_neighbor(candidate_id, dst_label, e_label, vec);
				for (int i = 0; i < vec.size(); i++)
				{
					unsigned int neighbor_cand = vec[i];
					if (check_candidate(dst_qid, neighbor_cand))
						neighbor_list.push_back(neighbor_cand);
				}
				vec.clear();
			}
			return neighbor_list.size();
	}

	int get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		if (use_edge_view) {
			unordered_map<unsigned long long, edge_view>::iterator iter = idx->edges.find(edge_pair);
			if (iter != idx->edges.end())
				iter->second.get_neighbor(candidate_id, neighbor_list);
		}
		else
		{
			int e_label = qg->check_edge_label(src_qid, dst_qid);
			assert(e_label != -1);
			int dst_label = qg->g[dst_qid].label;
			vector<unsigned int> vec;
			dg->get_labeled_neighbor(candidate_id, dst_label, e_label, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int neighbor_cand = vec[i];
				if (check_candidate(dst_qid, neighbor_cand))
					neighbor_list.push_back(neighbor_cand);
			}
			vec.clear();
		}
		return neighbor_list.size();
	}
	bool check_edge(unsigned int src_qid, unsigned int dst_qid, unsigned int src, unsigned int dst)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		if (use_edge_view) {
			if (idx->edges.find(edge_pair) != idx->edges.end() && idx->edges[edge_pair].check_edge(src, dst))
				return true;
			else
				return false;
		}
		else
		{
			int e_label = qg->check_edge_label(src_qid, dst_qid);
			if (idx->check_candidate(src_qid, src) && idx->check_candidate(dst_qid, dst) && dg->check_edge_label(src, dst, e_label))
				return true;
			else
				return false;
		}
	}

	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		idx->get_candidates(id, vec);
	}


	void set_candidate(unsigned int id, unsigned int candidate)
	{
		idx->set_candidate(id, candidate);
	}
	void erase_candidate(unsigned int id, unsigned int candidate)
	{
		idx->erase_candidate(id, candidate);
	}

};