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
#include<assert.h>
#include<set>
#include<algorithm>
#include"../setting.h"
using namespace std;
bool print = false;

struct nei_map
{
	map<unsigned int, unsigned int> edge_cnt;
	unsigned int matched_parent;
	unsigned int matched_child;
	nei_map()
	{
		matched_parent = 0;
		matched_child = 0;
	}
	~nei_map()
	{
		edge_cnt.clear();
	}
};

class symbi_index
{
public:
	global_index_template* idx;
	global_index_template* down_idx; // in the original paper the down_idx and idx are called D1 and D2 respectively, which indicates if the ancestors of a data vertex matches the ancestors o
	// of query vertex, or if both ancetor and successor match ( in perspective of the DAG).
	data_graph* dg;
	query_graph* qg;
	vector<vector<nei_map*>> nei_maps;
	DAG_query* dag;
	unsigned int node_num;
	bool use_edge_view;


	symbi_index(data_graph* d, query_graph* q, bool use_edge_view_ = false)
	{
		dg = d;
		qg = q;
		node_num = d->vertice_number;
		idx = new global_index_template;
		idx->candidate_map.resize(node_num);
		down_idx = new global_index_template;
		down_idx->candidate_map.resize(node_num);
		for (int i = 0; i < node_num; i++) {
			idx->candidate_map[i] = 0;
			down_idx->candidate_map[i] = 0;
		}
		nei_maps.resize(q->node_number);
		for (int i = 0; i < q->node_number; i++)
		{
			nei_maps[i].resize(node_num);
			for (int j = 0; j < node_num; j++)
				nei_maps[i][j] = NULL;
		}

		dag = new DAG_query;
		dag->symbi_initiate(q);
		use_edge_view = use_edge_view_;
	}

	~symbi_index()
	{
		delete idx;
		delete down_idx;
		delete dag;
		for (int i = 0; i < qg->node_number; i++)
		{
			for (int j = 0; j < node_num; j++)
			{
				if (nei_maps[i][j])
					delete nei_maps[i][j];
			}
		}
	}

	unsigned int memory_compute()
	{
		unsigned int mem = 0;
		mem += sizeof(void*) * 5;
		mem += sizeof(nei_maps);
		mem += sizeof(node_num) + sizeof(use_edge_view);
		mem = byte_align(mem, 8);
		mem += idx->memory_compute();
		mem += down_idx->memory_compute();
		
		mem += nei_maps.capacity() * sizeof(vector<void*>);
		for (int i = 0; i < nei_maps.size(); i++)
		{
			mem += nei_maps[i].capacity() * sizeof(void*);
			for (int j = 0; j < nei_maps[i].size(); j++)
			{
				if (nei_maps[i][j])
				{
					mem += sizeof(unsigned int) * 2 + sizeof(map<unsigned int, unsigned int>);
					mem += nei_maps[i][j]->edge_cnt.size() * (sizeof(unsigned int) * 2 + sizeof(void*) * 2);
				}
			}
		}
		mem += dag->memory_compute();
		return mem;
	}


	void resize()
	{
		node_num = dg->vertice_number;
		idx->candidate_map.resize(node_num);
		down_idx->candidate_map.resize(node_num);
		for (int i = 0; i < node_num; i++) {
			idx->candidate_map[i] = 0;
			down_idx->candidate_map[i] = 0;
		}

		nei_maps.resize(qg->node_number);
		for (int i = 0; i < qg->node_number; i++)
		{
			nei_maps[i].resize(node_num);
			for (int j = 0; j < node_num; j++)
				nei_maps[i][j] = NULL;
		}
	}

	
	bool NLF_check(unsigned int query_node, unsigned int data_node)
	{
		if (index_NLF) {
			for (map<unsigned long long, int>::iterator nlf_iter = qg->NLF[query_node].begin(); nlf_iter != qg->NLF[query_node].end(); nlf_iter++)
			{
				if (dg->NLF[data_node].find(nlf_iter->first) == dg->NLF[data_node].end() || dg->NLF[data_node][nlf_iter->first] < nlf_iter->second)
					return false;
			}
		}
		return true;
	}

	bool NLF_newly_updated(unsigned int query_node, unsigned int data_node, unsigned long long updated_label) // check if the inserted label result into the pass of NLF
	{
		if (index_NLF) {
			unsigned int current_val = dg->NLF[data_node][updated_label];
			if (current_val != qg->NLF[query_node][updated_label]) // either NLF check still fails, or NLF check already passes before update. If the latter happens, we know the candidate is filtered not due to NLF
				return false;
			for (map<unsigned long long, int>::iterator nlf_iter = qg->NLF[query_node].begin(); nlf_iter != qg->NLF[query_node].end(); nlf_iter++)
			{
				if (dg->NLF[data_node].find(nlf_iter->first) == dg->NLF[data_node].end() || dg->NLF[data_node][nlf_iter->first] < nlf_iter->second)
					return false;
			}
		}
		return true;

	}


	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		idx->get_candidates(id, vec);
	}
	bool check_candidate(unsigned int id, unsigned int candidate_id)
	{
		return idx->check_candidate(id, candidate_id);
	}

	int get_neighbor(unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
		if (use_edge_view) {
			unordered_map<unsigned long long, edge_view>::iterator iter = idx->edges.find(edge_pair);
			if (iter != idx->edges.end())
			{
				vector<unsigned int> vec;
				iter->second.get_neighbor(candidate_id, vec);
				unsigned int dst_qid = (edge_pair & 0xFFFFFFFF);
				for (int i = 0; i < vec.size(); i++)
				{
					unsigned int neighbor_cand = vec[i];
					if (idx->check_candidate(dst_qid, neighbor_cand))
						neighbor_list.push_back(neighbor_cand);
				}
				vec.clear();
			}
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
				if (idx->check_candidate(dst_qid, neighbor_cand))
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
			{
				vector<unsigned int> vec;
				iter->second.get_neighbor(candidate_id, vec);
				for (int i = 0; i < vec.size(); i++)
				{
					unsigned int neighbor_cand = vec[i];
					if (idx->check_candidate(dst_qid, neighbor_cand))
						neighbor_list.push_back(neighbor_cand);
				}
				vec.clear();
			}
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
				if (idx->check_candidate(dst_qid, neighbor_cand))
					neighbor_list.push_back(neighbor_cand);
			}
			vec.clear();
		}
		return neighbor_list.size();
	}
	bool check_edge(unsigned int src_qid, unsigned int dst_qid, unsigned int src, unsigned int dst)
	{
		if (idx->check_candidate(src_qid, src) && idx->check_candidate(dst_qid, dst)) {
			int e_label = qg->check_edge_label(src_qid, dst_qid);
			if(dg->check_edge_label(src, dst, e_label))
				return true;
		}
		
		return false;
	}

	
	void top_down_search()
	{
		vector<unordered_set<unsigned int>> candidates;
		vector<unsigned int> visited_degree;
		candidates.resize(qg->g.size());
		visited_degree.resize(qg->g.size());
		for (int i = 0; i < visited_degree.size(); i++)
			visited_degree[i] = 0;
		queue<int> q;
		q.push(dag->root);
		vector<unsigned int> vec;
		dg->get_labeled_node(qg->g[dag->root].label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			if (NLF_check(dag->root, vec[i]))
				candidates[dag->root].insert(vec[i]);
		}
		vec.clear();
		while (!q.empty())
		{
			unsigned int p = q.front();
			q.pop();
			for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				down_idx->set_candidate(p, *iter);
			vector<neighbor_info> vec;
			dag->get_out_neighbor(p, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int neighbor = vec[i].id;
				unsigned long long reverse_edge = merge_long_long(neighbor, p);
			//	cout << "visiting edge " << neighbor << ' ' << p << endl;
				visited_degree[neighbor]++;
				bool to_check = false;
				unsigned int in_degree = dag->get_in_degree(neighbor);
				if (visited_degree[neighbor] == in_degree) {
					to_check = true;
					q.push(neighbor);
				}
				for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				{
					unsigned int cand = *iter;
					vector<unsigned int> neighbor_vec;
					if (use_edge_view)
						idx->edges[merge_long_long(p, neighbor)].get_neighbor(cand, neighbor_vec);
					else
						dg->get_labeled_neighbor(cand, vec[i].v_label, vec[i].e_label, neighbor_vec);
					for (int k = 0; k < neighbor_vec.size(); k++) {
						if (nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.find(p) == nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.end() || nei_maps[neighbor][neighbor_vec[k]]->edge_cnt[p] == 0)
							nei_maps[neighbor][neighbor_vec[k]]->matched_parent++;
						if (nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.find(p) == nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.end())
							nei_maps[neighbor][neighbor_vec[k]]->edge_cnt[p] = 1;
						else
							nei_maps[neighbor][neighbor_vec[k]]->edge_cnt[p]++;
						if (to_check&& nei_maps[neighbor][neighbor_vec[k]]->matched_parent ==in_degree) {
								if(NLF_check(neighbor, neighbor_vec[k]))
									candidates[neighbor].insert(neighbor_vec[k]);
						}
					}
				}
			}
		}
		for (int i = 0; i < candidates.size(); i++)
			candidates[i].clear();
		candidates.clear();
		visited_degree.clear();
	}

	void bottom_up_search()
	{
		vector<unordered_set<unsigned int>> candidates;
		vector<unsigned int> visited_degree;
		candidates.resize(qg->g.size());
		visited_degree.resize(qg->g.size());
		for (int i = 0; i < visited_degree.size(); i++)
			visited_degree[i] = 0;
		vector<unsigned int> leafs;
		dag->get_leaf_nodes(leafs);
		queue<int> q;
		for (int i = 0; i < leafs.size(); i++)
		{
			vector<unsigned int> vec;
			down_idx->get_candidates(leafs[i], vec);
			for (int j = 0; j < vec.size(); j++)
				candidates[leafs[i]].insert(vec[j]);
			q.push(leafs[i]);
		}
		while (!q.empty())
		{
			unsigned int p = q.front();
			q.pop();
			for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				idx->set_candidate(p, *iter);
			vector<neighbor_info> vec;
			dag->get_in_neighbor(p, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int neighbor = vec[i].id;
				unsigned long long edge = merge_long_long(neighbor, p);
				visited_degree[neighbor]++;
				bool to_check = false;
				unsigned int out_degree = dag->get_out_degree(neighbor);
				if (visited_degree[neighbor] == out_degree) {
					to_check = true;
					q.push(neighbor);
				}
				for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				{
					unsigned int cand = *iter;
					vector<unsigned int> neighbor_vec;
					if (use_edge_view)
						idx->edges[merge_long_long(p, neighbor)].get_neighbor(cand, neighbor_vec);
					else
						dg->get_labeled_neighbor(cand, vec[i].v_label, vec[i].e_label, neighbor_vec);
					for (int k = 0; k < neighbor_vec.size(); k++) {
						if (nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.find(p) == nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.end() || nei_maps[neighbor][neighbor_vec[k]]->edge_cnt[p] == 0)
							nei_maps[neighbor][neighbor_vec[k]]->matched_child++;
						if (nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.find(p) == nei_maps[neighbor][neighbor_vec[k]]->edge_cnt.end())
							nei_maps[neighbor][neighbor_vec[k]]->edge_cnt[p] = 1;
						else
							nei_maps[neighbor][neighbor_vec[k]]->edge_cnt[p]++;
						if (to_check)
						{
							if (down_idx->check_candidate(neighbor, neighbor_vec[k])&& nei_maps[neighbor][neighbor_vec[k]]->matched_child == out_degree)
								candidates[neighbor].insert(neighbor_vec[k]);
						}
					}
				}
			}
		}
		for (int i = 0; i < candidates.size(); i++)
			candidates[i].clear();
		candidates.clear();
		visited_degree.clear();
	}


	void build_edge_view()
	{
			for (int i = 0; i < qg->g.size(); i++)
			{
				unsigned int id = i;
				if (dag->get_out_degree(id) == 0)
					continue;
				vector<unsigned int> cand_vec;
				dg->get_labeled_node(qg->g[i].label, cand_vec);
				vector<neighbor_info> neighbors;
				dag->get_out_neighbor(id, neighbors);
				for (int j = 0; j < neighbors.size(); j++)
				{
					unsigned int neighbor_id = neighbors[j].id;
					unsigned long long edge_pair = merge_long_long(id, neighbor_id);
					unsigned long long reverse_edge_pair = merge_long_long(neighbor_id, id);
					for (int k = 0; k < cand_vec.size(); k++)
					{
						vector<unsigned int> vec;
						dg->get_labeled_neighbor(cand_vec[k], neighbors[j].v_label, neighbors[j].e_label, vec);
						for (int m = 0; m < vec.size(); m++)
						{
							idx->edges[edge_pair].insert_edge(cand_vec[k], vec[m]);
							idx->edges[reverse_edge_pair].insert_edge(vec[m], cand_vec[k]);
						}
					}
				}
			}
	}

	void construct_index()
	{
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u = i;
			vector<unsigned int> vec;
			dg->get_labeled_node(qg->g[u].label, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int v = vec[j];
				nei_maps[u][v] = new nei_map;
			}
		}
		if (use_edge_view)
			build_edge_view();
		top_down_search();
		bottom_up_search();
	}

	void down_insert(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int> > &updated_1, queue<pair<unsigned int, unsigned int>>&updated_2)
	{
		if (down_idx->check_candidate(qid, id))
			return;
	//	if (print)
	//		cout << "check " << qid << ' ' << id<< endl;
		down_idx->set_candidate(qid, id);
		vector<neighbor_info> neighbor_vec;
		dag->get_out_neighbor(qid, neighbor_vec);
		for (int i = 0; i < neighbor_vec.size(); i++) {
			unsigned int neighbor_id = neighbor_vec[i].id;
			unsigned long long edge_pair = merge_long_long(qid, neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(neighbor_id, qid);
			vector<unsigned int> cand_neighbors;
	//		if (print)
	//			cout << "check neighbor " << neighbor_id << endl;
			if (use_edge_view)
				idx->edges[edge_pair].get_neighbor(id, cand_neighbors);
			else
				dg->get_labeled_neighbor(id, neighbor_vec[i].v_label, neighbor_vec[i].e_label, cand_neighbors);
	//		if (print)
	//			cout << "get neighbor finished " << cand_neighbors.size()<<endl;
			for (int j = 0; j < cand_neighbors.size(); j++) {
		//		if (!nei_maps[neighbor_id][cand_neighbors[j]])
		//		{
		//			cout << "NULL nei_map " << neighbor_id << ' ' << cand_neighbors[j] << endl;
		//		}
				if (nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.find(qid)== nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.end()|| nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]==0)
				{
					nei_maps[neighbor_id][cand_neighbors[j]]->matched_parent++;
					if (nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.find(qid) == nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.end())
						nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid] = 1;
					else
						nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]++;
					if (nei_maps[neighbor_id][cand_neighbors[j]]->matched_parent == dag->get_in_degree(neighbor_id)&&NLF_check(neighbor_id, cand_neighbors[j])){
						updated_1.push(make_pair(neighbor_id, cand_neighbors[j]));
						if(nei_maps[neighbor_id][cand_neighbors[j]]->matched_child == dag->get_out_degree(neighbor_id))
							updated_2.push(make_pair(neighbor_id, cand_neighbors[j]));
					}
				}
				else
				{
					nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]++;
				}
			}
		}
	}

	void up_insert(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (idx->check_candidate(qid, id))
			return;
		idx->set_candidate(qid, id);
		vector<neighbor_info> neighbor_vec;
		dag->get_in_neighbor(qid, neighbor_vec);
		for (int i = 0; i < neighbor_vec.size(); i++) {
			unsigned int neighbor_id = neighbor_vec[i].id;
			unsigned long long edge_pair = merge_long_long(neighbor_id, qid);
			unsigned long long reverse_edge_pair = merge_long_long(qid, neighbor_id);
			vector<unsigned int> cand_neighbors;
			if (use_edge_view)
				idx->edges[reverse_edge_pair].get_neighbor(id, cand_neighbors);
			else
				dg->get_labeled_neighbor(id, neighbor_vec[i].v_label, neighbor_vec[i].e_label, cand_neighbors);
			for (int j = 0; j < cand_neighbors.size(); j++) {
				if (nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.find(qid) == nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.end() || nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid] == 0)
				{
					nei_maps[neighbor_id][cand_neighbors[j]]->matched_child++;
					if (nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.find(qid) == nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt.end())
						nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid] = 1;
					else
						nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]++;
					if (nei_maps[neighbor_id][cand_neighbors[j]]->matched_child == dag->get_out_degree(neighbor_id)&&down_idx->check_candidate(neighbor_id, cand_neighbors[j]))
							updated_2.push(make_pair(neighbor_id, cand_neighbors[j]));
				}
				else
					nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]++;
			}
			cand_neighbors.clear();
		}
		neighbor_vec.clear();
	}



	void expand_tables(unsigned int id)
	{
		//cout <<"expanding "<< id << ' ' << node_num << endl;
		assert(id == node_num);
		down_idx->candidate_map.push_back(0);
		idx->candidate_map.push_back(0);
		for (int i = 0; i < qg->g.size(); i++) {
			nei_maps[i].push_back(NULL);
		}
		node_num++;
	}

	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		print = false;
		if (src == 50672 && dst == 964075)
			print = true;
	//	cout << "inserting edge " << src << ' ' << src_label << ' ' << dst << ' ' << dst_label <<" in index "<< endl;
		if (src >= down_idx->candidate_map.size())
			expand_tables(src);
		if (dst >= down_idx->candidate_map.size())
			expand_tables(dst);
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			if (!nei_maps[vec[i].first][src])
				nei_maps[vec[i].first][src] = new nei_map;
			if (!nei_maps[vec[i].second][dst])
				nei_maps[vec[i].second][dst] = new nei_map;
		}
		for (int i = 0; i < vec.size(); i++)  // if src or dst is not candidates of a query node before and it is now, this query node must be endpoint of an edge with (src_label, edge_label, dst_label), otherwise inserting this edge will not influece the NLF
		{
		//	if(print)
	//			cout << "mapped to " << vec[i].first << ' ' << vec[i].second << endl;
			int direction = dag->get_direction(vec[i].first, vec[i].second);
			unsigned int pnode, cnode, pcand, ccand, plabel, clabel;
			if (direction == 1)
			{
				pnode = vec[i].first;
				cnode = vec[i].second;
				pcand = src;
				ccand = dst;
				plabel = src_label;
				clabel = dst_label;
			}
			else if (direction == -1)
			{
				pnode = vec[i].second;
				cnode = vec[i].first;
				pcand = dst;
				ccand = src;
				plabel = dst_label;
				clabel = src_label;
			}
			else
			{
				cout << "error, inconsistent dag and query graph" << endl;
				continue;
			}
		//	if (print)
		//		cout << "direction settled" << endl;
		/*	if (!nei_maps[pnode][pcand])
				nei_maps[pnode][pcand] = new nei_map;
			if (!nei_maps[cnode][ccand])
				nei_maps[cnode][ccand] = new nei_map;*/

			unsigned long long edge_pair = merge_long_long(pnode, cnode);
			unsigned long long reverse_edge_pair = merge_long_long(cnode, pnode);
			unsigned long long label_pair = merge_long_long(clabel, edge_label);
			unsigned long long reverse_label_pair = merge_long_long(plabel, edge_label);
			if (use_edge_view) {
				idx->edges[edge_pair].insert_edge(pcand, ccand);
				idx->edges[reverse_edge_pair].insert_edge(ccand, pcand);
			}
			queue<pair<unsigned int, unsigned int>> updated_1; // updated 1 for the record of the new down candidates, and updated_2 for the newly added candidates.
			queue<pair<unsigned int, unsigned int>> updated_2;
			
			bool parent_set = false;

			if (!down_idx->check_candidate(pnode, pcand)) // there may ba the case that all parents of pcand have been matched, but the NLF check fails, and it is not record as a candidate,
				// in this case, we need to recheck it, since the NLF check may be met after inserting the new edge;
			{
				if (nei_maps[pnode][pcand]->matched_parent == dag->get_in_degree(pnode) && NLF_newly_updated(pnode, pcand, label_pair)) {
					updated_1.push(make_pair(pnode, pcand));
					parent_set = true;
					if (nei_maps[pnode][pcand]->matched_child == dag->get_out_degree(pnode))
						updated_2.push(make_pair(pnode, pcand));
				}
			}
			if (!down_idx->check_candidate(cnode, ccand))
			{
				if (nei_maps[cnode][ccand]->matched_parent == dag->get_in_degree(cnode) && NLF_newly_updated(cnode, ccand, reverse_label_pair)) {
					updated_1.push(make_pair(cnode, ccand));
					if (nei_maps[cnode][ccand]->matched_child == dag->get_out_degree(cnode))
						updated_2.push(make_pair(cnode, ccand));
				}
			}

			if(down_idx->check_candidate(pnode, pcand))
			{
				if (nei_maps[cnode][ccand]->edge_cnt.find(pnode) == nei_maps[cnode][ccand]->edge_cnt.end() || nei_maps[cnode][ccand]->edge_cnt[pnode] == 0) {
					nei_maps[cnode][ccand]->matched_parent++;
					if (nei_maps[cnode][ccand]->edge_cnt.find(pnode) == nei_maps[cnode][ccand]->edge_cnt.end())
						nei_maps[cnode][ccand]->edge_cnt[pnode] = 1;
					else
						nei_maps[cnode][ccand]->edge_cnt[pnode]++;
					if (nei_maps[cnode][ccand]->matched_parent == dag->get_in_degree(cnode) && NLF_check(cnode, ccand)) {
						updated_1.push(make_pair(cnode, ccand));
						if (nei_maps[cnode][ccand]->matched_child == dag->get_out_degree(cnode))
							updated_2.push(make_pair(cnode, ccand));
					}
				}
				else
					nei_maps[cnode][ccand]->edge_cnt[pnode]++;
			}

			if(idx->check_candidate(cnode, ccand)){
				if (nei_maps[pnode][pcand]->edge_cnt.find(cnode) == nei_maps[pnode][pcand]->edge_cnt.end() || nei_maps[pnode][pcand]->edge_cnt[cnode] == 0) {
					nei_maps[pnode][pcand]->matched_child++;
					if (nei_maps[pnode][pcand]->edge_cnt.find(cnode) == nei_maps[pnode][pcand]->edge_cnt.end())
						nei_maps[pnode][pcand]->edge_cnt[cnode] = 1;
					else
						nei_maps[pnode][pcand]->edge_cnt[cnode]++;
					if (nei_maps[pnode][pcand]->matched_child ==dag->get_out_degree(pnode) && (down_idx->check_candidate(pnode, pcand)||parent_set))
						updated_2.push(make_pair(pnode, pcand));
				}
				else
					nei_maps[pnode][pcand]->edge_cnt[cnode]++;
			}
		//	if(print)
	//			cout << "queue pop begin" << endl;
			while (!updated_1.empty())
			{
				pair<unsigned int, unsigned int> p = updated_1.front();
				updated_1.pop();
				down_insert(p.first, p.second, updated_1, updated_2);
			}
		//	if(print)
	//			cout << "down insert end" << endl;
			while (!updated_2.empty())
			{
				pair<unsigned int, unsigned int> p = updated_2.front();
				updated_2.pop();
				up_insert(p.first, p.second, updated_2);
			}
		}
	}



	void down_erase(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int> >& updated_1, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (!down_idx->check_candidate(qid, id))
			return;
		down_idx->erase_candidate(qid, id);
		vector<neighbor_info> neighbor_vec;
		dag->get_out_neighbor(qid, neighbor_vec);
		for (int i = 0; i < neighbor_vec.size(); i++) {
			unsigned int neighbor_id = neighbor_vec[i].id;
			unsigned long long edge_pair = merge_long_long(qid, neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(neighbor_id, qid);
			vector<unsigned int> cand_neighbors;
			if (use_edge_view)
				idx->edges[edge_pair].get_neighbor(id, cand_neighbors);
			else
				dg->get_labeled_neighbor(id, neighbor_vec[i].v_label, neighbor_vec[i].e_label, cand_neighbors);
			for (int j = 0; j < cand_neighbors.size(); j++) {
				nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]--;
				if (nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid] == 0)
				{
					nei_maps[neighbor_id][cand_neighbors[j]]->matched_parent--;
					if (down_idx->check_candidate(neighbor_id, cand_neighbors[j]))
					{
						updated_1.push(make_pair(neighbor_id, cand_neighbors[j]));
						if (idx->check_candidate(neighbor_id, cand_neighbors[j])) {
							updated_2.push(make_pair(neighbor_id, cand_neighbors[j]));
						}
					}
				}
			}
		}
	}

	void up_erase(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (!idx->check_candidate(qid, id))
			return;
		idx->erase_candidate(qid, id);
		vector<neighbor_info> neighbor_vec;
		dag->get_in_neighbor(qid, neighbor_vec);
		for (int i = 0; i < neighbor_vec.size(); i++) {
			unsigned int neighbor_id = neighbor_vec[i].id;
			unsigned long long edge_pair = merge_long_long(neighbor_id, qid);
			unsigned long long reverse_edge_pair = merge_long_long(qid, neighbor_id);
			vector<unsigned int> cand_neighbors;
			if (use_edge_view)
				idx->edges[reverse_edge_pair].get_neighbor(id, cand_neighbors);
			else
				dg->get_labeled_neighbor(id, neighbor_vec[i].v_label, neighbor_vec[i].e_label, cand_neighbors);
			for (int j = 0; j < cand_neighbors.size(); j++) {
				nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid]--;
				if (nei_maps[neighbor_id][cand_neighbors[j]]->edge_cnt[qid] == 0)
				{
					nei_maps[neighbor_id][cand_neighbors[j]]->matched_child--;
					if (idx->check_candidate(neighbor_id, cand_neighbors[j])) {
						updated_2.push(make_pair(neighbor_id, cand_neighbors[j]));
					}
				}
			}
			cand_neighbors.clear();
		}
		neighbor_vec.clear();
	}
	void delete_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		for (int i = 0; i < vec.size(); i++)  // if src or dst is not candidates of a query node before and it is now, this query node must be endpoint of an edge with (src_label, edge_label, dst_label), otherwise inserting this edge will not influece the NLF
		{
			int direction = dag->get_direction(vec[i].first, vec[i].second);
			unsigned int pnode, cnode, pcand, ccand;
			if (direction == 1)
			{
				pnode = vec[i].first;
				cnode = vec[i].second;
				pcand = src;
				ccand = dst;
			}
			else if (direction == -1)
			{
				pnode = vec[i].second;
				cnode = vec[i].first;
				pcand = dst;
				ccand = src;
			}
			else
			{
				cout << "error, inconsistent dag and query graph" << endl;
				continue;
			}

			unsigned long long edge_pair = merge_long_long(pnode, cnode);
			unsigned long long reverse_edge_pair = merge_long_long(cnode, pnode);
			if (use_edge_view) {
				idx->edges[edge_pair].delete_edge(pcand, ccand);
				idx->edges[reverse_edge_pair].delete_edge(ccand, pcand);
			}
			queue<pair<unsigned int, unsigned int>> updated_1; // updated 1 for the record of the deleted down candidates, and updated_2 for the deleted candidates.
			queue<pair<unsigned int, unsigned int>> updated_2;
		//	bool parent_deleted = false;

			if (down_idx->check_candidate(pnode, pcand)) 
			{
				if (!NLF_check(pnode, pcand)) {
					updated_1.push(make_pair(pnode, pcand));
					if (idx->check_candidate(pnode, pcand))
						updated_2.push(make_pair(pnode, pcand));
				}
			}
			if (down_idx->check_candidate(cnode, ccand))
			{
				if (!NLF_check(cnode, ccand)) {
					updated_1.push(make_pair(cnode, ccand));
					if (idx->check_candidate(cnode, ccand))
						updated_2.push(make_pair(cnode, ccand));
				}
			}
			if (down_idx->check_candidate(pnode, pcand)){
				nei_maps[cnode][ccand]->edge_cnt[pnode]--;
				if (nei_maps[cnode][ccand]->edge_cnt[pnode] == 0) {
					nei_maps[cnode][ccand]->matched_parent--;
					if (down_idx->check_candidate(cnode, ccand))
						updated_1.push(make_pair(cnode, ccand));
					if (idx->check_candidate(cnode, ccand)) {
						updated_2.push(make_pair(cnode, ccand));
					}
				}
			}
			if (idx->check_candidate(cnode, ccand)) 
			{
				nei_maps[pnode][pcand]->edge_cnt[cnode]--;
				if (nei_maps[pnode][pcand]->edge_cnt[cnode] == 0) {
					nei_maps[pnode][pcand]->matched_child--;
					if (idx->check_candidate(pnode, pcand)) {
						updated_2.push(make_pair(pnode, pcand));
					}
				}
			}
			
			while (!updated_1.empty())
			{
				pair<unsigned int, unsigned int> p = updated_1.front();
				updated_1.pop();
				down_erase(p.first, p.second, updated_1, updated_2);
			}
			while (!updated_2.empty())
			{
				pair<unsigned int, unsigned int> p = updated_2.front();
				updated_2.pop();
				up_erase(p.first, p.second, updated_2);
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
		//	cout << "node " << i << " " << tmp_vec.size() << endl;
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

};
