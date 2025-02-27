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

struct tree_node
{
	unsigned int node_id;
	int parent;
	unsigned int p_elabel;
	vector<unsigned int> child;
	vector<unsigned int> c_elabel;

	tree_node(unsigned int node_id_ = 0)
	{
		node_id = node_id_;
		parent = -1;
	}
	~tree_node()
	{
		child.clear();
		c_elabel.clear();
	}
	void set_parent(unsigned int p, unsigned int edge_lable)
	{
		parent = p;
		p_elabel = edge_lable;
	}
	void add_child(unsigned int c, unsigned int child_lable)
	{
		child.push_back(c);
		c_elabel.push_back(child_lable);
	}
};
class query_tree
{
public:
	unsigned int root;
	vector<tree_node*> tree;
	unsigned int node_num;
	set<unsigned long long> tree_edges;
	query_tree(unsigned int node_num_=0)
	{
		node_num = node_num_;
		tree.resize(node_num);
		for (int i = 0; i < node_num; i++) {
			tree[i] = new tree_node(i);
		}
	}
	~query_tree()
	{
		for (int i = 0; i < tree.size(); i++) {
			if(tree[i])
				delete tree[i];
		}
		tree_edges.clear();
	}
	void add_tree_edge(unsigned int p, unsigned int c, unsigned int elabel)
	{
		tree[c]->set_parent(p, elabel);
		tree[p]->add_child(c, elabel);
		tree_edges.insert(merge_long_long(p, c));

	}
	void get_leaf_nodes(vector<unsigned int>& vec)
	{
		for (int i = 0; i < tree.size(); i++)
		{
			if (tree[i]&&tree[i]->child.empty())
				vec.push_back(i);
		}
	}
	void delete_tree_node(unsigned int u)
	{
		int parent = tree[u]->parent;
		if (parent != -1) {
			for (int i = 0; i < tree[parent]->child.size(); i++)
			{
				if (tree[parent]->child[i] == u)
				{
					tree[parent]->child.erase(tree[parent]->child.begin() + i);
					tree[parent]->c_elabel.erase(tree[parent]->c_elabel.begin() + i);
					break;
				}
			}
		}
		delete tree[u];
		tree[u] = NULL;
	}
	unsigned int memory_compute()
	{
		unsigned int mem = 0;
		mem += sizeof(unsigned int)*2 + sizeof(set<unsigned long long>) + sizeof(vector<tree_node*>);
		mem += tree_edges.size() * (sizeof(unsigned long long) + sizeof(void*) * 2);
		mem += tree.capacity() * sizeof(void*);
		for (int i = 0; i < tree.size(); i++)
		{
			if (!tree[i])
				continue;
			mem += sizeof(unsigned int) * 3;
			mem += sizeof(vector<unsigned int>) * 2;
			mem += (tree[i]->child.size() + tree[i]->c_elabel.size()) * sizeof(unsigned int);
		}
		return mem;
	} 
	void copy(query_tree* qt)
	{
		root = qt->root;
		node_num = qt->node_num;
		for (int i = 0; i < qt->tree.size(); i++)
		{
			tree_node* t = new tree_node;
			t->node_id = qt->tree[i]->node_id;
			t->parent = qt->tree[i]->parent;
			t->p_elabel = qt->tree[i]->p_elabel;
			for (int j = 0; j < qt->tree[i]->child.size(); j++)
			{
				t->child.push_back(qt->tree[i]->child[j]);
				t->c_elabel.push_back(qt->tree[i]->c_elabel[j]);
			}
			tree.push_back(t);
		}
		for (auto iter = qt->tree_edges.begin(); iter != qt->tree_edges.end(); iter++)
			tree_edges.insert(*iter);
		return;
	}

};
bool selectivity_sort(pair<unsigned long long, unsigned int>& p1, pair<unsigned long long, unsigned int>& p2)
{
	return p1.second < p2.second;
}
class TF_index
{
public:
	global_index_template* explicit_idx;
	global_index_template* implicit_idx; 
	data_graph* dg;
	query_graph* qg;
	query_tree* qt;

	unsigned int node_num;
	bool use_edge_view;

	unsigned int get_freq(unsigned int src_label, unsigned int dst_label, unsigned int e_label)
	{
		vector<unsigned int> vec;
		dg->get_labeled_node(src_label, vec);
		unsigned int sum = 0;
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int v = vec[i];
			unsigned long long label_pair = merge_long_long(dst_label, e_label);
			if (dg->NLF[v].find(label_pair) != dg->NLF[v].end())
				sum += dg->NLF[v][label_pair];
		}
		if (src_label == dst_label) // in this case the same edge has been counted twice;
			sum = sum / 2;
		return sum;
	}
	unsigned int choose_root()
	{
		unsigned int root;
		vector<pair<unsigned long long, unsigned int> > edge_selectivity;
		for (unsigned int i = 0; i < qg->g.size(); i++)
		{
			unsigned int u1 = i;
			unsigned int l1 = qg->g[i].label;
			vector<neighbor_info> vec;
			qg->get_neighbor(u1, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				unsigned int u2 = vec[j].id;
				unsigned int l2 = vec[j].v_label;
				unsigned int elabel = vec[j].e_label;
				unsigned int f = get_freq(l1, l2, elabel);
				edge_selectivity.push_back(make_pair(merge_long_long(u1, u2), f));
			}
			vec.clear();
		}
		sort(edge_selectivity.begin(), edge_selectivity.end(), selectivity_sort);
		unsigned long long root_edge = edge_selectivity[0].first;
		unsigned int root_src = (root_edge >> 32);
		unsigned int root_dst = (root_edge & 0xFFFFFFFF);
		unsigned int src_label = qg->g[root_src].label;
		unsigned int dst_label = qg->g[root_dst].label;
		if (src_label != dst_label)
		{
			vector<unsigned int> candidates;
			dg->get_labeled_node(src_label, candidates);
			unsigned int src_freq = candidates.size();
			candidates.clear();
			dg->get_labeled_node(dst_label, candidates);
			unsigned int dst_freq = candidates.size();
			candidates.clear();
			if (src_freq > dst_freq)
				root = root_src;
			else if (dst_freq > src_freq)
				root = root_dst;
			else
			{
				unsigned int src_degree = qg->get_total_degree(root_src);
				unsigned int dst_degree = qg->get_total_degree(root_dst);
				if(src_degree>=dst_degree)
					root = root_src;
				else
					root = root_dst;
			}
		}
		else
		{
			unsigned int src_degree = qg->get_total_degree(root_src);
			unsigned int dst_degree = qg->get_total_degree(root_dst);
			if (src_degree >= dst_degree)
				root = root_src;
			else
				root = root_dst;
		}
		return root;
		
	}
	void choose_edge(unordered_set<unsigned int> &remaining) // here I noticed two versions of implementation. In TurboFlux 2018 paper, it is said that we should choose the edge which minimizes the "intermediate result size", while in 
					  // the new paper "In-depth Analysis of Continuous Subgraph Matching in a Common Delta Query Compilation Framework", it is said that we should choose the edge with the minimum caidinality. I choose the second method.
	{
		if (remaining.empty())
			return;
		unsigned int min_card = 0xFFFFFFFF;
		unsigned int min_src = 0;
		unsigned int min_dst = 0;
		unsigned int min_label = 0;
		for (auto iter = remaining.begin(); iter != remaining.end(); iter++)
		{
			unsigned int  u1 = *iter;
			unsigned int l1 = qg->g[u1].label;;
			vector<neighbor_info> vec;
			qg->get_neighbor(u1, vec);

			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int u2 = vec[i].id;
				if (remaining.find(u2) != remaining.end())
					continue;
				unsigned int card = get_freq(vec[i].v_label, l1, vec[i].e_label);
				if (card < min_card)
				{
					min_card = card;
					min_src = u2;
					min_dst = u1;
					min_label = vec[i].e_label;
				}
			}
		}
		qt->add_tree_edge(min_src, min_dst, min_label);
		remaining.erase(min_dst);

	}
	void generate_query_tree()
	{
		unsigned int r = choose_root();
		qt->root = r;
		unordered_set<unsigned int> remain;
		for (int i = 0; i < qg->g.size(); i++)
		{
			if (i != r)
				remain.insert(i);
		}
		while (!remain.empty())
			choose_edge(remain);
	}
	TF_index(data_graph* d, query_graph* q, bool use_edge_view_ = false)
	{
		dg = d;
		qg = q;
		node_num = d->vertice_number;
		use_edge_view = use_edge_view_;
		explicit_idx = new global_index_template;
		explicit_idx->candidate_map.resize(node_num);
		implicit_idx = new global_index_template;
		implicit_idx->candidate_map.resize(node_num);
		for (int i = 0; i < node_num; i++)
		{
			explicit_idx->candidate_map[i] = 0;
			implicit_idx->candidate_map[i] = 0;
		}
		qt = new query_tree(qg->g.size());
		generate_query_tree();
	}

	~TF_index()
	{
		delete explicit_idx;
		delete implicit_idx;
		delete qt;
	}

	unsigned int memory_compute()
	{
		unsigned int mem = 0;
		mem += sizeof(void*) * 5;
		mem += sizeof(node_num) + sizeof(use_edge_view);
		mem = byte_align(mem, 8);
		mem += explicit_idx->memory_compute();
		mem += implicit_idx->memory_compute();
		mem += qt->memory_compute();
		return mem;
	}


	void resize()
	{
		node_num = dg->vertice_number;
		explicit_idx->candidate_map.resize(node_num);
		implicit_idx->candidate_map.resize(node_num);
		for (int i = 0; i < node_num; i++) {
			explicit_idx->candidate_map[i] = 0;
			implicit_idx->candidate_map[i] = 0;
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



	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		explicit_idx->get_candidates(id, vec);
	}
	bool check_candidate(unsigned int id, unsigned int candidate_id)
	{
		return explicit_idx->check_candidate(id, candidate_id);
	}

	int get_neighbor(unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
		/*unsigned int dst_qid = (edge_pair & 0xFFFFFFFF);
		unsigned int src_qid = (edge_pair >> 32);
		int e_label = qg->check_edge_label(src_qid, dst_qid);
		assert(e_label != -1);
		int dst_label = qg->g[dst_qid].label;
		vector<unsigned int> vec;
		vector<unsigned int> debug;
		dg->get_labeled_neighbor(candidate_id, dst_label, e_label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int neighbor_cand = vec[i];
			if (explicit_idx->check_candidate(dst_qid, neighbor_cand))
				debug.push_back(neighbor_cand);
		}
		vec.clear();*/
		if (use_edge_view) {
			unsigned int src_qid = (edge_pair >> 32);
			unsigned int dst_qid = (edge_pair & 0xFFFFFFFF);
			unsigned long long reverse_edge_pair = merge_long_long(dst_qid, src_qid);
			if (qt->tree_edges.find(edge_pair) != qt->tree_edges.end()) {
				unordered_map<unsigned long long, edge_view>::iterator iter = explicit_idx->edges.find(edge_pair);
				if (iter != explicit_idx->edges.end())
					iter->second.get_neighbor(candidate_id, neighbor_list);
			//	if (neighbor_list.size() != debug.size())
			//		cout << "error, wrong neighbor get by edge view" << endl;
				return neighbor_list.size(); 
			}
			else if (qt->tree_edges.find(reverse_edge_pair) != qt->tree_edges.end()) // a reverse tree edge, in this case we still have to check the explicit_idx
			{
				unordered_map<unsigned long long, edge_view>::iterator iter = explicit_idx->edges.find(edge_pair);
				if (iter != explicit_idx->edges.end())
				{
					vector<unsigned int> vec;
					iter->second.get_neighbor(candidate_id, vec);
					for (int i = 0; i < vec.size(); i++)
					{
						if (explicit_idx->check_candidate(dst_qid, vec[i]))
							neighbor_list.push_back(vec[i]);
					}
					vec.clear();
				}
			//	if (neighbor_list.size() != debug.size())
		//			cout << "error, wrong neighbor get by edge view" << endl;
				return neighbor_list.size();
			}
		}
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
			if (explicit_idx->check_candidate(dst_qid, neighbor_cand))
				neighbor_list.push_back(neighbor_cand);
		}
		vec.clear();
		return neighbor_list.size();
	}

	int get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
	/*	int e_label = qg->check_edge_label(src_qid, dst_qid);
		assert(e_label != -1);
		int dst_label = qg->g[dst_qid].label;
		vector<unsigned int> vec;
		vector<unsigned int> debug;
		dg->get_labeled_neighbor(candidate_id, dst_label, e_label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int neighbor_cand = vec[i];
			if (explicit_idx->check_candidate(dst_qid, neighbor_cand))
				debug.push_back(neighbor_cand);
		}
		vec.clear();*/
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		unsigned long long reverse_edge_pair = merge_long_long(dst_qid, src_qid);
		if (use_edge_view) {
			if (qt->tree_edges.find(edge_pair) != qt->tree_edges.end()) {
				unordered_map<unsigned long long, edge_view>::iterator iter = explicit_idx->edges.find(edge_pair);
				if (iter != explicit_idx->edges.end())
					iter->second.get_neighbor(candidate_id, neighbor_list);
			//	if (neighbor_list.size() != debug.size())
		//			cout << "error, wrong neighbor get by edge view" << endl;
				return neighbor_list.size(); 
			}
			else if (qt->tree_edges.find(reverse_edge_pair) != qt->tree_edges.end()) // a reverse tree edge, in this case we still have to check the explicit_idx
			{
				unordered_map<unsigned long long, edge_view>::iterator iter = explicit_idx->edges.find(edge_pair);
				if (iter != explicit_idx->edges.end())
				{
					vector<unsigned int> vec;
					iter->second.get_neighbor(candidate_id, vec);
					for (int i = 0; i < vec.size(); i++)
					{
						if (explicit_idx->check_candidate(dst_qid, vec[i]))
							neighbor_list.push_back(vec[i]);
					}
					vec.clear();
				}
		//		if (neighbor_list.size() != debug.size())
		//			cout << "error, wrong neighbor get by edge view" << endl;
				return neighbor_list.size();
			}
		}
		int e_label = qg->check_edge_label(src_qid, dst_qid);
		assert(e_label != -1);
		int dst_label = qg->g[dst_qid].label;
		vector<unsigned int> vec;
		dg->get_labeled_neighbor(candidate_id, dst_label, e_label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int neighbor_cand = vec[i];
			if (explicit_idx->check_candidate(dst_qid, neighbor_cand))
				neighbor_list.push_back(neighbor_cand);
		}
		vec.clear();
		return neighbor_list.size();
	}
	bool check_edge(unsigned int src_qid, unsigned int dst_qid, unsigned int src, unsigned int dst)
	{
		if (explicit_idx->check_candidate(src_qid, src) && explicit_idx->check_candidate(dst_qid, dst)) {
			int e_label = qg->check_edge_label(src_qid, dst_qid);
			if (dg->check_edge_label(src, dst, e_label))
				return true;
		}

		return false;
	}

	bool check_child_map(unsigned int qid, unsigned int id)
	{
		bool filtered = false;
	//	bool tmp_filtered = false;

		for (int i = 0; i < qt->tree[qid]->child.size(); i++)
		{
			unsigned int c = qt->tree[qid]->child[i];
		/*	vector<unsigned int> vec;
			dg->get_labeled_neighbor(id, qg->g[c].label, qt->tree[qid]->c_elabel[i], vec);
			bool mapped = false;
			for (int j = 0; j < vec.size(); j++)
			{
				if (explicit_idx->check_candidate(c, vec[j]))
				{
					mapped = true;
					break;
				}
			}
			if (!mapped)
			{
				tmp_filtered = true;
			}*/

			if (use_edge_view) {
				if (!explicit_idx->edges[merge_long_long(qid, c)].check_node(id))
				{
					filtered = true;
			//		if (filtered != tmp_filtered)
				//		cout << "error check" << endl;
					break;
				}
			}
			else
			{
				vector<unsigned int> vec;
				dg->get_labeled_neighbor(id, qg->g[c].label, qt->tree[qid]->c_elabel[i], vec);
				bool mapped = false;
				for (int j = 0; j < vec.size(); j++)
				{
					if (explicit_idx->check_candidate(c, vec[j]))
					{
						mapped = true;
						break;
					}
				}
				if (!mapped)
				{
					filtered = true;
					break;
				}
			}
		}
		return !filtered;
	}

	void top_down_search()
	{
		vector<unordered_set<unsigned int>> candidates;
		candidates.resize(qg->g.size());
		queue<int> q;
		q.push(qt->root);
		vector<unsigned int> vec;
		dg->get_labeled_node(qg->g[qt->root].label, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			if (NLF_check(qt->root, vec[i]))
				candidates[qt->root].insert(vec[i]);
		}
		vec.clear();
		while (!q.empty())
		{
			unsigned int p = q.front();
			q.pop();
			for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				implicit_idx->set_candidate(p, *iter);
			for (int i = 0; i < qt->tree[p]->child.size(); i++)
			{
				unsigned int neighbor = qt->tree[p]->child[i];
				unsigned long long edge_pair = merge_long_long(p, neighbor);
				unsigned long long reverse_edge = merge_long_long(neighbor, p);
				q.push(neighbor);

				for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				{
					unsigned int cand = *iter;
					vector<unsigned int> neighbor_vec;
					dg->get_labeled_neighbor(cand, qg->g[neighbor].label, qt->tree[p]->c_elabel[i], neighbor_vec);
					for (int k = 0; k < neighbor_vec.size(); k++) {
						if (NLF_check(neighbor, neighbor_vec[k])) {
							candidates[neighbor].insert(neighbor_vec[k]);
							if (use_edge_view) {
								implicit_idx->edges[edge_pair].insert_edge(cand, neighbor_vec[k]);
								implicit_idx->edges[reverse_edge].insert_edge(neighbor_vec[k], cand);
							}
						}
					}
				}
			}
		}
		for (int i = 0; i < candidates.size(); i++)
			candidates[i].clear();
		candidates.clear();
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
		qt->get_leaf_nodes(leafs);
		queue<int> q;
		for (int i = 0; i < leafs.size(); i++)
		{
			vector<unsigned int> vec;
			implicit_idx->get_candidates(leafs[i], vec);
			for (int j = 0; j < vec.size(); j++)
				candidates[leafs[i]].insert(vec[j]);
			q.push(leafs[i]);
		}
		while (!q.empty())
		{
			unsigned int p = q.front();
			q.pop();
			for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
				explicit_idx->set_candidate(p, *iter);
			int parent = qt->tree[p]->parent;
			if (parent == -1)
				continue;
			unsigned long long edge_pair = merge_long_long(parent, p);
			unsigned long long reverse_edge = merge_long_long(p, parent);
			visited_degree[parent]++;
			bool to_check = false;
			unsigned int child_num = qt->tree[parent]->child.size();
			if (visited_degree[parent] == child_num) {
					to_check = true;
					q.push(parent);
			}
			for (auto iter = candidates[p].begin(); iter != candidates[p].end(); iter++)
			{
				unsigned int cand = *iter;
				vector<unsigned int> neighbor_vec;
			/*	vector<unsigned int> debug;
				vector<unsigned int> tmp_vec;
				dg->get_labeled_neighbor(cand, qg->g[parent].label, qt->tree[p]->p_elabel, tmp_vec);
				for (int i = 0; i < tmp_vec.size(); i++)
				{
					if (implicit_idx->check_candidate(parent, tmp_vec[i]))
						debug.push_back(tmp_vec[i]);
				}
				tmp_vec.clear();*/

				if (use_edge_view) {
					implicit_idx->edges[reverse_edge].get_neighbor(cand, neighbor_vec);
				//	if (neighbor_vec.size() != debug.size())
				//		cout << "error, wrong neighbor size" << endl;
				}
				else
				{
					vector<unsigned int> tmp_vec;
					dg->get_labeled_neighbor(cand, qg->g[parent].label, qt->tree[p]->p_elabel, tmp_vec);
					for (int i = 0; i < tmp_vec.size(); i++)
					{
						if (implicit_idx->check_candidate(parent, tmp_vec[i]))
							neighbor_vec.push_back(tmp_vec[i]);
					}
					tmp_vec.clear();
				}
				for (int k = 0; k < neighbor_vec.size(); k++) {
					if (use_edge_view) {
						explicit_idx->edges[edge_pair].insert_edge(neighbor_vec[k], cand);
						explicit_idx->edges[reverse_edge].insert_edge(cand, neighbor_vec[k]);
					}
					if (to_check) {
						if (check_child_map(parent, neighbor_vec[k]))
							candidates[parent].insert(neighbor_vec[k]);
					}
				}
		}
		}
		for (int i = 0; i < candidates.size(); i++)
			candidates[i].clear();
		candidates.clear();
		visited_degree.clear();
	}


	void construct_index()
	{
		top_down_search();
		bottom_up_search();
	}


	void down_insert(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int> >& updated_1, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (implicit_idx->check_candidate(qid, id))
			return;
		implicit_idx->set_candidate(qid, id);
		bool filtered = false;
		for (int i = 0; i < qt->tree[qid]->child.size(); i++) {
			unsigned int neighbor_id = qt->tree[qid]->child[i];
			unsigned long long edge_pair = merge_long_long(qid, neighbor_id);
			unsigned long long reverse_edge = merge_long_long(neighbor_id, qid);
			bool child_found = false; 
			vector<unsigned int> cand_neighbors;
			dg->get_labeled_neighbor(id, qg->g[neighbor_id].label, qt->tree[qid]->c_elabel[i], cand_neighbors);
//			if (cand_neighbors.empty())
//				filtered = true;
			for (int j = 0; j < cand_neighbors.size(); j++) {
				if (!NLF_check(neighbor_id, cand_neighbors[j]))
					continue;
				implicit_idx->edges[edge_pair].insert_edge(id, cand_neighbors[j]);
				implicit_idx->edges[reverse_edge].insert_edge(cand_neighbors[j], id);
				if (explicit_idx->check_candidate(neighbor_id, cand_neighbors[j]))
				{
					child_found = true; 
					if (use_edge_view) {
						explicit_idx->edges[edge_pair].insert_edge(id, cand_neighbors[j]);
						explicit_idx->edges[reverse_edge].insert_edge(cand_neighbors[j], id);
					}
				}
				else
					updated_1.push(make_pair(neighbor_id, cand_neighbors[j]));
			}
				if(!child_found)
					filtered = true;
		}
		if (!filtered) // If filtered = false, then either qid is a leaf whit no child, or for each child, this is a explicit candidate connected to id, in either case, id becomes an explicit candidate.
			updated_2.push(make_pair(qid, id));
	}

	void up_insert(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (explicit_idx->check_candidate(qid, id))
			return;
		explicit_idx->set_candidate(qid, id);

		int parent_id = qt->tree[qid]->parent;
		if (parent_id == -1)
			return;
		unsigned long long edge_pair = merge_long_long(parent_id, qid);
		unsigned long long reverse_edge = merge_long_long(qid, parent_id);
		vector<unsigned int> cand_neighbors;
	/*	vector<unsigned int> debug;
		vector<unsigned int> tmp_vec;
		dg->get_labeled_neighbor(id, qg->g[parent_id].label, qt->tree[qid]->p_elabel, tmp_vec);
		for (int i = 0; i < tmp_vec.size(); i++)
		{
			if (implicit_idx->check_candidate(parent_id, tmp_vec[i]))
				debug.push_back(tmp_vec[i]);
		}
		tmp_vec.clear();*/
		if (use_edge_view) {
			implicit_idx->edges[reverse_edge].get_neighbor(id, cand_neighbors);
		//	if (cand_neighbors.size() != debug.size())
	//			cout << "error, wrong neighbor size" << endl;
		}
		else
		{
			vector<unsigned int> tmp_vec;
			dg->get_labeled_neighbor(id, qg->g[parent_id].label, qt->tree[qid]->p_elabel, tmp_vec);
			for (int i = 0; i < tmp_vec.size(); i++)
			{
				if (implicit_idx->check_candidate(parent_id, tmp_vec[i]))
					cand_neighbors.push_back(tmp_vec[i]);
			}
			tmp_vec.clear();
		}
		for (int i = 0; i < cand_neighbors.size(); i++)
		{
			if (use_edge_view) {
				explicit_idx->edges[edge_pair].insert_edge(cand_neighbors[i], id);
				explicit_idx->edges[reverse_edge].insert_edge(id, cand_neighbors[i]);
			}
			if (check_child_map(parent_id, cand_neighbors[i]))
				updated_2.push(make_pair(parent_id, cand_neighbors[i]));
		}
		cand_neighbors.clear();
	}



	void expand_tables(unsigned int id)
	{
		//cout <<"expanding "<< id << ' ' << node_num << endl;
		assert(id == node_num);
		implicit_idx->candidate_map.push_back(0);
		explicit_idx->candidate_map.push_back(0);
		node_num++;
	}

	void check_parent(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int> >& updated_1)
	{
		int parent = qt->tree[qid]->parent;
		if (parent == -1)
			updated_1.push(make_pair(qid, id));
		else
		{
			vector<unsigned int> vec;
			dg->get_labeled_neighbor(id, qg->g[parent].label, qt->tree[qid]->p_elabel, vec);
			bool find_p = false;
			for (int i = 0; i < vec.size(); i++)
			{
				if (implicit_idx->check_candidate(parent, vec[i]))
				{
					implicit_idx->edges[merge_long_long(parent, qid)].insert_edge(vec[i], id);
					implicit_idx->edges[merge_long_long(qid, parent)].insert_edge(id, vec[i]);
					find_p = true;
				}
			}
			if (find_p)
				updated_1.push(make_pair(qid, id));
		}
	}

	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		if (src >= node_num)
			expand_tables(src);
		if (dst >= node_num)
			expand_tables(dst);
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		unordered_set<unsigned int> checked_src;
		unordered_set<unsigned int> checked_dst;
		queue<pair<unsigned int, unsigned int>> updated_1;
		queue<pair<unsigned int, unsigned int>> updated_2;
		for (int i = 0; i < vec.size(); i++)  // if src or dst is not candidates of a query node before and it is now, this query node must be endpoint of an edge with (src_label, edge_label, dst_label), otherwise inserting this edge will not influece the NLF
		{
			unsigned int src_qid = vec[i].first;
			unsigned int dst_qid = vec[i].second;
			// even if the matched query edge is not a tree edge, we still need to check if the NLF check has been meet, since the NLF check is not based on the query 
			// tree, but the entire query graph;
			if (index_NLF) {
				if (checked_src.find(src_qid) == checked_src.end()) {
					checked_src.insert(src_qid);
					if (!implicit_idx->check_candidate(src_qid, src))
					{
						if (NLF_check(src_qid, src))
							check_parent(src_qid, src, updated_1);
					}
				}

				if (checked_dst.find(dst_qid) == checked_dst.end()) {
					checked_dst.insert(dst_qid);
					if (!implicit_idx->check_candidate(dst_qid, dst))
					{
						if (NLF_check(dst_qid, dst))
							check_parent(dst_qid, dst, updated_1);
					}
				}
			}
			unsigned int pnode, cnode, pcand, ccand, plabel, clabel;
			if (qt->tree_edges.find(merge_long_long(vec[i].first, vec[i].second)) != qt->tree_edges.end())
			{
				pnode = vec[i].first;
				cnode = vec[i].second;
				pcand = src;
				ccand = dst;
				plabel = src_label;
				clabel = dst_label;
			}
			else if (qt->tree_edges.find(merge_long_long(vec[i].second, vec[i].first)) != qt->tree_edges.end())
			{
				pnode = vec[i].second;
				cnode = vec[i].first;
				pcand = dst;
				ccand = src;
				plabel = dst_label;
				clabel = src_label;
			}
			else
				continue;

			unsigned long long edge_pair = merge_long_long(pnode, cnode);
			unsigned long long reverse_edge_pair = merge_long_long(cnode, pnode);
			unsigned long long label_pair = merge_long_long(clabel, edge_label);
			unsigned long long reverse_label_pair = merge_long_long(plabel, edge_label);



			if (implicit_idx->check_candidate(pnode, pcand))
			{
				if (NLF_check(cnode, ccand))
				{
					if (use_edge_view) {
						implicit_idx->edges[edge_pair].insert_edge(pcand, ccand);
						implicit_idx->edges[reverse_edge_pair].insert_edge(ccand, pcand);
					}
					updated_1.push(make_pair(cnode, ccand));
				}
				if (explicit_idx->check_candidate(cnode, ccand))
				{
					if (use_edge_view) {
						explicit_idx->edges[edge_pair].insert_edge(pcand, ccand);
						explicit_idx->edges[reverse_edge_pair].insert_edge(ccand, pcand);
					}
					if (check_child_map(pnode, pcand))
						updated_2.push(make_pair(pnode, pcand));
				}
			}
		}

			while (!updated_1.empty())
			{
				pair<unsigned int, unsigned int> p = updated_1.front();
				updated_1.pop();
				down_insert(p.first, p.second, updated_1, updated_2);
			}
			while (!updated_2.empty())
			{
				pair<unsigned int, unsigned int> p = updated_2.front();
				updated_2.pop();
				up_insert(p.first, p.second, updated_2);
			}
	}



	void down_erase(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int> >& updated_1, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (!implicit_idx->check_candidate(qid, id))
			return;
		implicit_idx->erase_candidate(qid, id);
		if (explicit_idx->check_candidate(qid, id))
			updated_2.push(make_pair(qid, id));
		for (int i = 0; i < qt->tree[qid]->child.size(); i++) {
			unsigned int neighbor_id = qt->tree[qid]->child[i];
			unsigned long long edge_pair = merge_long_long(qid, neighbor_id);
			unsigned long long reverse_edge_pair = merge_long_long(neighbor_id, qid);
			vector<unsigned int> cand_neighbors;
			if(use_edge_view)
				implicit_idx->edges[edge_pair].get_neighbor(id, cand_neighbors);
			else
			{
				vector<unsigned int> tmp_vec;
				dg->get_labeled_neighbor(id, qg->g[neighbor_id].label, qt->tree[qid]->c_elabel[i], tmp_vec);
				for (int j = 0; j < tmp_vec.size(); j++)
				{
					if (implicit_idx->check_candidate(neighbor_id, tmp_vec[j]))
						cand_neighbors.push_back(tmp_vec[j]);
				}
				tmp_vec.clear();
			}
			for (int j = 0; j < cand_neighbors.size(); j++) {
				bool parent_deleted = true;
				if (use_edge_view) {
					implicit_idx->edges[edge_pair].delete_edge(id, cand_neighbors[j]);
					parent_deleted = implicit_idx->edges[reverse_edge_pair].delete_edge(cand_neighbors[j], id);
				}
				else
				{
					vector<unsigned int> tmp_vec;
					dg->get_labeled_neighbor(cand_neighbors[j], qg->g[qid].label, qt->tree[qid]->c_elabel[i], tmp_vec);
					for (int k = 0; k < tmp_vec.size(); k++)
					{
						if (implicit_idx->check_candidate(qid, tmp_vec[k]))
						{
							parent_deleted = false;
							break;
						}
					}
				}
				if (parent_deleted)
					updated_1.push(make_pair(neighbor_id, cand_neighbors[j]));
			}
		}
	}

	void up_erase(unsigned int qid, unsigned int id, queue<pair<unsigned int, unsigned int>>& updated_2)
	{
		if (!explicit_idx->check_candidate(qid, id))
			return;
		explicit_idx->erase_candidate(qid, id);
		vector<neighbor_info> neighbor_vec;
		int parent = qt->tree[qid]->parent;
		if (parent == -1)
			return;
		unsigned long long edge_pair = merge_long_long(parent, qid);
		unsigned long long reverse_edge = merge_long_long(qid, parent);
		vector<unsigned int> vec;
		if(use_edge_view)
			implicit_idx->edges[reverse_edge].get_neighbor(id, vec);
		else
		{
			vector<unsigned int> tmp_vec;
			dg->get_labeled_neighbor(id, qg->g[parent].label, qt->tree[qid]->p_elabel, tmp_vec);
			for (int j = 0; j < tmp_vec.size(); j++)
			{
				if (implicit_idx->check_candidate(parent, tmp_vec[j]))
					vec.push_back(tmp_vec[j]);
			}
			tmp_vec.clear();
		}
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int v = vec[i];
			bool child_deleted = true;
			if (use_edge_view) {
				explicit_idx->edges[reverse_edge].delete_edge(id, v);
				child_deleted = explicit_idx->edges[edge_pair].delete_edge(v, id);
			}
			else
			{
				vector<unsigned int> tmp_vec;
				dg->get_labeled_neighbor(vec[i], qg->g[qid].label, qt->tree[qid]->p_elabel, tmp_vec);
				for (int k = 0; k < tmp_vec.size(); k++)
				{
					if (explicit_idx->check_candidate(qid, tmp_vec[k]))
					{
						child_deleted = false;
						break;
					}
				}
			}
			if (child_deleted && explicit_idx->check_candidate(parent, v))
				updated_2.push(make_pair(parent, v));
		}
	}
	void delete_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		vector<pair<unsigned int, unsigned int> > vec;
		qg->get_edge_by_label(src_label, dst_label, edge_label, vec);
		queue<pair<unsigned int, unsigned int>> updated_1;
		queue<pair<unsigned int, unsigned int>> updated_2;
		for (int i = 0; i < vec.size(); i++)  // if src or dst is not candidates of a query node before and it is now, this query node must be endpoint of an edge with (src_label, edge_label, dst_label), otherwise inserting this edge will not influece the NLF
		{
			unsigned int src_qid = vec[i].first;
			unsigned int dst_qid = vec[i].second;
			if (implicit_idx->check_candidate(src_qid, src))
			{
				if (!NLF_check(src_qid, src))
					updated_1.push(make_pair(src_qid, src));
			}

			if (implicit_idx->check_candidate(dst_qid, dst))
			{
				if (!NLF_check(dst_qid, dst))
					updated_1.push(make_pair(dst_qid, dst));
			}

			unsigned int pnode, cnode, pcand, ccand, plabel, clabel;
			if (qt->tree_edges.find(merge_long_long(vec[i].first, vec[i].second)) != qt->tree_edges.end())
			{
				pnode = vec[i].first;
				cnode = vec[i].second;
				pcand = src;
				ccand = dst;
				plabel = src_label;
				clabel = dst_label;
			}
			else if (qt->tree_edges.find(merge_long_long(vec[i].second, vec[i].first)) != qt->tree_edges.end())
			{
				pnode = vec[i].second;
				cnode = vec[i].first;
				pcand = dst;
				ccand = src;
				plabel = dst_label;
				clabel = src_label;
			}
			else
				continue;
			unsigned long long edge_pair = merge_long_long(pnode, cnode);
			unsigned long long reverse_edge_pair = merge_long_long(cnode, pnode);

			if (implicit_idx->check_candidate(pnode, pcand)) {
				if (use_edge_view)
				{
					implicit_idx->edges[edge_pair].delete_edge(pcand, ccand);
					bool parent_deleted = implicit_idx->edges[reverse_edge_pair].delete_edge(ccand, pcand);
					if (parent_deleted)
						updated_1.push(make_pair(cnode, ccand));
				}
				else
				{
					bool parent_deleted = true;
					vector<unsigned int> tmp_vec;
					dg->get_labeled_neighbor(ccand, qg->g[pnode].label, qt->tree[cnode]->p_elabel, tmp_vec);
					for (int i = 0; i < tmp_vec.size(); i++)
					{
						if (implicit_idx->check_candidate(pnode, tmp_vec[i]))
						{
							parent_deleted = false;
							break;
						}
					}
					if (parent_deleted)
						updated_1.push(make_pair(cnode, ccand));
				}
			}
			if (explicit_idx->check_candidate(cnode, ccand))
			{
				if (use_edge_view)
				{
					explicit_idx->edges[reverse_edge_pair].delete_edge(ccand, pcand);
					bool child_deleted = explicit_idx->edges[edge_pair].delete_edge(pcand, ccand);
					if (child_deleted)
						updated_2.push(make_pair(pnode, pcand));
				}
				else
				{
					bool child_deleted = true;
					vector<unsigned int> tmp_vec;
					dg->get_labeled_neighbor(pcand, qg->g[cnode].label, qt->tree[cnode]->p_elabel, tmp_vec);
					for (int i = 0; i < tmp_vec.size(); i++)
					{
						if (explicit_idx->check_candidate(cnode, tmp_vec[i]))
						{
							child_deleted = false;
							break;
						}
					}
					if (child_deleted)
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

	pair<unsigned int, unsigned int> get_size()
	{
		unsigned int vertex_sum = 0;
		vector<unsigned int> tmp_vec;
		for (int i = 0; i < qg->g.size(); i++)
		{
			explicit_idx->get_candidates(i, tmp_vec);
			vertex_sum += tmp_vec.size();
			cout << "node " << i << ' ' << tmp_vec.size() << endl;
			tmp_vec.clear();
		}
		unsigned int edge_sum = 0;
		for (auto iter = explicit_idx->edges.begin(); iter != explicit_idx->edges.end(); iter++) {
			edge_sum += iter->second.get_size();
			//	cout << "edge view " << (iter->first >> 32) << ' ' << (iter->first & 0xFFFFFFFF) << " size " << iter->second.get_size() << endl;
		}
		for (auto iter = implicit_idx->edges.begin(); iter != implicit_idx->edges.end(); iter++) {
			edge_sum += iter->second.get_size();
			//	cout << "edge view " << (iter->first >> 32) << ' ' << (iter->first & 0xFFFFFFFF) << " size " << iter->second.get_size() << endl;
		}
		edge_sum = edge_sum / 2; //each edge has two copy of views, here we only return the number of edges;
		return make_pair(vertex_sum, edge_sum);
	}

};
