#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include<map>
#include "../utility/struct.h"
#include<assert.h>
#include<queue>
#define dense_degree 3
using namespace std;

class query_graph
{
public:
	vector<node> g; // as the node IDs in query_graph are continuous and start from 0, the index in the vector is the corresponding node ID;
	unsigned int edge_number;
	unsigned int node_number;
	int query_type; // 0 means tree, 1 means sparse cyclic, 2 means dense cyclic.
	unordered_map<unsigned long long, vector<int>> edge_label_map; // map from id pair of end points to label of the edge. Here we support parallel edges with different edge labels between the same vertex pair.
	vector<vector<unsigned int> > unlabeled_edges;
	vector<map<unsigned long long, int> > NLF; // NLF of each node.
	unordered_map<unsigned int, vector<unsigned int> > label_node_map; // map from node label to nodes with this label. 
	query_graph(unsigned int vertice_num = 0) {
		edge_number = 0;
		node_number = vertice_num;
		g.resize(node_number);
		NLF.resize(node_number);
		unlabeled_edges.resize(node_number);
	//	cout << g.size() << endl;
		for (int i = 0; i < node_number; i++)
		{
			g[i] = node(i, 0); 
			NLF[i] = map<unsigned long long, int>();
		}
	}
	void clear()
	{
		edge_number = 0;
		node_number = 0;
		query_type = -1;
		for (int i = 0; i < g.size(); i++)
		{
			for (map<unsigned long long, vector<unsigned int>>::iterator iter = g[i].label_edge_map.begin(); iter != g[i].label_edge_map.end(); iter++)
				iter->second.clear();
			g[i].label_edge_map.clear();
		}
		g.clear();
		for (unordered_map<unsigned int, vector<unsigned int>>::iterator iter = label_node_map.begin(); iter != label_node_map.end(); iter++)
			iter->second.clear();
		label_node_map.clear();
		for (int i = 0; i < NLF.size(); i++)
			NLF[i].clear();
		NLF.clear();
		for (int i = 0; i < unlabeled_edges.size(); i++)
			unlabeled_edges[i].clear();
		unlabeled_edges.clear();
		for (unordered_map<unsigned long long, vector<int>>::iterator iter = edge_label_map.begin(); iter != edge_label_map.end(); iter++)
			iter->second.clear();
		edge_label_map.clear();
	}

	~query_graph()
	{
		clear();
	}

	void set_invalid_node(vector<unsigned int>& vec) // mark some nodes as invalid. As in rapidflow we need to delete some nodes from the query graph
	{
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int v = vec[i];
			assert(g[v].label_edge_map.empty());
			g[v].id = 0xFFFFFFFF;
			g[v].label = -1;
		}
	}

	void copy(query_graph* q)
	{
		q->edge_number = edge_number;
		q->node_number = node_number;
		q->resize(node_number);
		for (int i = 0; i < g.size(); i++)
		{
			q->g[i].id = g[i].id;
			q->g[i].label = g[i].label;
			for (map<unsigned long long, vector<unsigned int> >::iterator iter = g[i].label_edge_map.begin(); iter != g[i].label_edge_map.end(); iter++)
				q->g[i].label_edge_map[iter->first] = iter->second;
		}
		for (unordered_map<unsigned long long, vector<int>>::iterator iter = edge_label_map.begin(); iter != edge_label_map.end(); iter++)
			q->edge_label_map[iter->first] = iter->second;
		for (int i = 0; i < unlabeled_edges.size(); i++)
			q->unlabeled_edges[i] = unlabeled_edges[i];
		for (int i = 0; i < NLF.size(); i++)
			q->NLF[i] = NLF[i];
		for (unordered_map<unsigned int, vector<unsigned int> >::iterator iter = label_node_map.begin(); iter != label_node_map.end(); iter++)
			q->label_node_map[iter->first] = iter->second;
	}
	void copy_with_node_delete(unsigned int src, unsigned int dst, query_graph* q)
	{
		q->resize(node_number);
		vector<unsigned int> vec = { src, dst };
		q->set_invalid_node(vec);
		for (int i = 0; i < g.size(); i++)
		{
			if (g[i].id == src || g[i].id == dst)
				continue;
			vector<neighbor_info> vec;
			get_neighbor(g[i].id, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				if (vec[j].id != src && vec[j].id != dst&&vec[j].id<g[i].id)
				{
					q->insert_edge(g[i].id, g[i].label, vec[j].id, vec[j].v_label, vec[j].e_label);
				}
			}
		}
	}
	void resize(unsigned int vertice_num)
	{
		node_number = vertice_num;
		g.resize(node_number);
		NLF.resize(node_number);
		unlabeled_edges.resize(node_number);
		for (int i = 0; i < node_number; i++)
		{
			g[i] = node(i, 0);
			NLF[i] = map<unsigned long long, int>();
		}
	}
	int get_label(unsigned int id)
	{
		return g[id].label;
	}
	int get_total_degree(unsigned int node_id)
	{
		int sum = 0;
		for (map<unsigned long long, vector<unsigned int> >::iterator iter = g[node_id].label_edge_map.begin(); iter != g[node_id].label_edge_map.end(); iter++)
			sum += iter->second.size();
		return sum;
	}
	void insert_edge(unsigned int src, int src_label, unsigned int dst, int dst_label, int edge_label)
	{
		edge_number++;
		g[src].label = src_label;
		g[dst].label = dst_label;
		unsigned long long label_pair = merge_long_long(dst_label, edge_label);
		g[src].label_edge_map[label_pair].push_back(dst);

		unsigned long long reverse_label_pair = merge_long_long(src_label, edge_label);
		g[dst].label_edge_map[reverse_label_pair].push_back(src);

		vector<unsigned int>::iterator vec_iter = lower_bound(label_node_map[src_label].begin(), label_node_map[src_label].end(), src);
		if (vec_iter == label_node_map[src_label].end() || *vec_iter != src)
			label_node_map[src_label].emplace(vec_iter, src);
		if (NLF[src].find(label_pair) != NLF[src].end())
			NLF[src][label_pair] ++;
		else
			NLF[src][label_pair] = 1;

		vec_iter = lower_bound(label_node_map[dst_label].begin(), label_node_map[dst_label].end(), dst);
		if(vec_iter==label_node_map[dst_label].end()||*vec_iter!=dst)
			label_node_map[dst_label].emplace(vec_iter, dst);
		if (NLF[dst].find(reverse_label_pair) != NLF[dst].end())
			NLF[dst][reverse_label_pair]++;
		else
			NLF[dst][reverse_label_pair] = 1;

		unsigned long long id_pair = merge_long_long(src, dst);
		unsigned long long reverse_id_pair = merge_long_long(dst, src);
		edge_label_map[id_pair].push_back(edge_label);
		edge_label_map[reverse_id_pair].push_back(edge_label);

		unlabeled_edges[src].push_back(dst);
		unlabeled_edges[dst].push_back(src);

	}
	void get_neighbor(unsigned int id, vector<neighbor_info>& neighbor_list)
	{
		for (map<unsigned long long, vector<unsigned int> >::iterator iter = g[id].label_edge_map.begin(); iter != g[id].label_edge_map.end(); iter++)
		{
			for (int i = 0; i < iter->second.size(); i++)
			{
				neighbor_list.push_back(neighbor_info(iter->second[i], (iter->first >> 32), (iter->first & 0xFFFFFFFF)));
			}
		}
	}

	void get_labeled_neighbor(unsigned int id, int target_neighbor_label, int target_edge_label, vector<unsigned int>& neighbor_list)
	{
		unsigned long long label_pair = merge_long_long(target_neighbor_label, target_edge_label);
		if (g[id].label_edge_map.find(label_pair) != g[id].label_edge_map.end())
		{
			neighbor_list = g[id].label_edge_map[label_pair];
		}
	}

	void get_unlabeled_neighbor(unsigned int id, vector<unsigned int>& vec)
	{
		vec = unlabeled_edges[id];
	}
	void get_edge_by_label(int src_label, int dst_label, int edge_label, vector<pair<unsigned int, unsigned int> >& query_edges)
	{
		if (label_node_map.find(src_label) != label_node_map.end())
		{
			for (int i = 0; i < label_node_map[src_label].size(); i++)
			{
				unsigned int id = label_node_map[src_label][i];
				unsigned long long label_pair = merge_long_long(dst_label, edge_label);
				if (g[id].label_edge_map.find(label_pair) != g[id].label_edge_map.end())
				{
					for (int j = 0; j < g[id].label_edge_map[label_pair].size(); j++)
						query_edges.push_back(make_pair(id, g[id].label_edge_map[label_pair][j]));
				}
			}
		}
	}

	void check_type()
	{
		unordered_map<unsigned int, int> parent;
		unsigned int v = g[0].id;
		queue<unsigned int> q;
		q.push(v);
		parent[v] = -1;
		bool tree = true;
		while (!q.empty())
		{
			v = q.front();
			q.pop();
			vector<unsigned int> vec;
			get_unlabeled_neighbor(v, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int u = vec[i];
				if (parent.find(u) != parent.end() && u != parent[v])
				{
					tree = false;
					break;
				}
				else if(u!=parent[v])
				{
					parent[u] = v;
					q.push(u);
				}
			}
		}
		if (tree)
			query_type = 0;
		else if (edge_number / node_number <= dense_degree)
			query_type = 1;
		else
			query_type = 2;

	}

	void print_graph()
	{
		for (int i = 0; i < g.size(); i++)
		{
			if (g[i].id != 0xFFFFFFFF) {
				cout << "node " << g[i].id << ' ' << g[i].label << " edge list" << endl;
				for (auto iter = g[i].label_edge_map.begin(); iter != g[i].label_edge_map.end(); iter++)
				{
					unsigned int dst_label = (iter->first >> 32);
					unsigned int edge_label = (iter->first & 0xFFFFFFFF);
					for (int j = 0; j < iter->second.size(); j++)
					{
						cout << iter->second[j] << ' ' << edge_label << ' ' << dst_label << endl;
					}
				}
			}
		}
	}
	void get_vertex_by_label(unsigned int node_label, vector<unsigned int>& node_id)
	{
		if (label_node_map.find(node_label) != label_node_map.end())
		{
			node_id = label_node_map[node_label];
		}
	}

	void check_edge_label(unsigned int src_id, unsigned int dst_id, vector<int>& edge_labels)
	{
		unsigned long long id_pair = merge_long_long(src_id, dst_id);
		if (edge_label_map.find(id_pair) != edge_label_map.end())
			edge_labels = edge_label_map[id_pair];
		else
			return;
	}

	int check_edge_label(unsigned int src_id, unsigned int dst_id)
	{
		unsigned long long id_pair = merge_long_long(src_id, dst_id);
		if (edge_label_map.find(id_pair) != edge_label_map.end())
			return edge_label_map[id_pair][0];
		else
			return -1;

	}


};
