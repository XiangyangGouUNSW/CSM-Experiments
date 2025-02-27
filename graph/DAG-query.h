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
#include<set>
#include<queue>
#include<stack>
#include<algorithm>
using namespace std;

class DAG_query
{
public:
	vector<node> edges;
	vector<node> reverse_edges;
	vector<vector<unsigned int>> unlabeled_edges;
	vector<vector<unsigned int>> unlabeled_reverse_edges;
	unordered_map<int, int> node_depth;
	int root;
	unsigned int node_num;
	DAG_query(unsigned int node_number = 0)
	{
		node_num = node_number;
		if (node_num != 0) {
			edges.resize(node_num);
			reverse_edges.resize(node_num);
			unlabeled_edges.resize(node_num);
			unlabeled_reverse_edges.resize(node_num);
		}
	}
	~DAG_query()
	{
		for (int i = 0; i < edges.size(); i++)
			edges[i].clear();
		for (int i = 0; i < edges.size(); i++)
			reverse_edges[i].clear();
	}
	unsigned int memory_compute()
	{
		unsigned int mem = 0;
		mem += sizeof(vector<node>) * 2 + sizeof(vector<vector<unsigned int>>) * 2;
		mem += sizeof(unordered_map<int, int>);
		mem += sizeof(int) + sizeof(unsigned int);
		mem = byte_align(mem, 8);
		mem += edges.capacity() * sizeof(node);
		for (int i = 0; i < edges.size(); i++)
			mem += edges[i].memory_compute();
		mem += reverse_edges.capacity() * sizeof(node);
		for (int i = 0; i < reverse_edges.size(); i++)
			mem += reverse_edges[i].memory_compute();
		mem += unlabeled_edges.capacity() * sizeof(vector<unsigned int>);
		for (int i = 0; i < unlabeled_edges.size(); i++)
			mem += unlabeled_edges[i].capacity() * sizeof(unsigned int);
		mem += unlabeled_reverse_edges.capacity() * sizeof(vector<unsigned int>);
		for (int i = 0; i < unlabeled_reverse_edges.size(); i++)
			mem += unlabeled_reverse_edges[i].capacity() * sizeof(unsigned int);
		unsigned int kvsize = sizeof(int) * 2;
		mem += node_depth.bucket_count() * sizeof(void*) + node_depth.size() * (kvsize + sizeof(void*));
		return mem;
	}
	int get_direction(unsigned int src, unsigned int dst)
	{
		auto iter = lower_bound(unlabeled_edges[src].begin(), unlabeled_edges[src].end(), dst);
		if (iter != unlabeled_edges[src].end() && *iter == dst)
			return 1;
		else
		{
			iter = lower_bound(unlabeled_edges[dst].begin(), unlabeled_edges[dst].end(), src);
			if (iter != unlabeled_edges[dst].end() && *iter == src)
				return -1;
			else
				return 0;
		}
	}
	int get_in_degree(unsigned int id)
	{
		return unlabeled_reverse_edges[id].size();
	}
	int get_out_degree(unsigned int id)
	{
		return unlabeled_edges[id].size();
	}
	int height_compute(unsigned int root_v, query_graph* qg)
	{
		root = root_v;
		node_depth.clear();
		queue<int> q;
		q.push(root);
		node_depth[root] = 0;
		unsigned int max_level = 0;
		while (!q.empty())
		{
			int v = q.front();
			unsigned int next_level = node_depth[v]+ 1;
			q.pop();
			vector<unsigned int> vec;
			qg->get_unlabeled_neighbor(v, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				if (node_depth.find(vec[i]) == node_depth.end())
				{
					node_depth[vec[i]] = next_level;
					q.push(vec[i]);
					if (next_level > max_level)
						max_level = next_level;
				}
			}
		}
		node_depth.clear();
		return max_level;
	}

	void dfs(unsigned int v, unordered_set<unsigned int>& visited_node, unordered_set<unsigned long long> visited_edges, query_graph* qg, unsigned int root, unsigned int &cycle_num)
	{
		vector<unsigned int> vec;
		qg->get_unlabeled_neighbor(v, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			unsigned int u = vec[i];
			unsigned long long edge;
			if (u > v)
				edge = merge_long_long(u, v);
			else
				edge = merge_long_long(v, u);
			if (visited_edges.find(edge) != visited_edges.end())
				continue;
			if (u == root)
			{
				cycle_num++;
				continue;
			}
			if (visited_node.find(u) != visited_node.end())
				continue;
			visited_node.insert(u);
			visited_edges.insert(edge);
			dfs(u, visited_node, visited_edges, qg, root, cycle_num);
			visited_edges.erase(edge);
			visited_node.erase(u);
		}
	}
	unsigned int cycle_compute(unsigned int root_v, query_graph* qg)
	{

		unsigned int cycle_num = 0;
		unordered_set<unsigned int> visited_node;
		unordered_set<unsigned long long> visited_edges;
		dfs(root_v, visited_node, visited_edges, qg, root_v, cycle_num);
		return cycle_num/2;

	}
	void insert_edge(int src, int dst, int src_label, int dst_label, int edge_label)
	{
		unsigned long long label_pair = merge_long_long(dst_label, edge_label);
		unsigned long long reverse_label_pair = merge_long_long(src_label, edge_label);
		edges[src].id = src;
		edges[src].label = src_label;
		edges[dst].id = dst;
		edges[dst].label = dst_label;
		reverse_edges[src].id = src;
		reverse_edges[src].label = src_label;
		reverse_edges[dst].id = dst;
		reverse_edges[dst].label = dst_label;
		if (edges[src].label_edge_map.find(label_pair) == edges[src].label_edge_map.end())
			edges[src].label_edge_map[label_pair] = vector<unsigned int>();
		if (reverse_edges[dst].label_edge_map.find(reverse_label_pair) == reverse_edges[dst].label_edge_map.end())
			reverse_edges[dst].label_edge_map[reverse_label_pair] = vector<unsigned int>();
		auto iter = lower_bound(edges[src].label_edge_map[label_pair].begin(), edges[src].label_edge_map[label_pair].end(), dst);
		if (iter == edges[src].label_edge_map[label_pair].end() || *iter != dst)
			edges[src].label_edge_map[label_pair].emplace(iter, dst);
		iter = lower_bound(reverse_edges[dst].label_edge_map[reverse_label_pair].begin(), reverse_edges[dst].label_edge_map[reverse_label_pair].end(), src);
		if (iter == reverse_edges[dst].label_edge_map[reverse_label_pair].end() || *iter != src)
			reverse_edges[dst].label_edge_map[reverse_label_pair].emplace(iter, src);
		iter = lower_bound(unlabeled_edges[src].begin(), unlabeled_edges[src].end(), dst);
		if (iter == unlabeled_edges[src].end() || *iter != dst)
			unlabeled_edges[src].emplace(iter, dst);
		iter = lower_bound(unlabeled_reverse_edges[dst].begin(), unlabeled_reverse_edges[dst].end(), src);
		if (iter == unlabeled_reverse_edges[dst].end() || *iter != src)
			unlabeled_reverse_edges[dst].emplace(iter, src);
	}
	void build_DAG(unsigned int root_v, query_graph* qg)
	{
		root = root_v;
		unordered_set<unsigned int> visited_nodes;
		queue<int> q;
		q.push(root);
		node_depth[root] = 0;
		unsigned int max_level = 0;
		while (!q.empty())
		{
			int v = q.front();
			visited_nodes.insert(v);
			unsigned int next_level = node_depth[v]++;
			q.pop();
			vector<neighbor_info> vec;
			qg->get_neighbor(v, vec);
			for (int i = 0; i < vec.size(); i++)
			{
				unsigned int neighbor = vec[i].id;
				if (node_depth.find(neighbor) == node_depth.end() || node_depth[neighbor] == node_depth[neighbor]) // only linked to not visited node or node with the same level;
				{
					if (node_depth.find(neighbor) == node_depth.end()){
					node_depth[neighbor] = next_level;
					q.push(neighbor);
					}
				}
				if (visited_nodes.find(neighbor) == visited_nodes.end())
					insert_edge(v, neighbor, qg->g[v].label, vec[i].v_label, vec[i].e_label);
			}
		}
	}
	void symbi_initiate(query_graph* qg)
	{
		node_num = qg->node_number;
		edges.resize(node_num);
		reverse_edges.resize(node_num);
		unlabeled_edges.resize(node_num);
		unlabeled_reverse_edges.resize(node_num);

		unsigned int max_height = 0;
		unsigned int max_root = 0;
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int height = height_compute(i, qg);
			if (height > max_height)
			{
				max_height = height;
				max_root = i;
			}
		}
		build_DAG(max_root, qg);
	}

	void improved_symbi_initiate(query_graph* qg)
	{
		node_num = qg->node_number;
		edges.resize(node_num);
		reverse_edges.resize(node_num);
		unlabeled_edges.resize(node_num);
		unlabeled_reverse_edges.resize(node_num);

		unsigned int max_height = 0;
		unsigned int max_cycle = 0;
		unsigned int max_root = 0;
		for (int i = 0; i < qg->g.size(); i++)
		{
			unsigned int cycle_num = cycle_compute(i, qg);
			unsigned int height = height_compute(i, qg);
		//	cout << i << ' ' << cycle_num << ' ' << height << endl;
			if (cycle_num > max_cycle)
			{
				max_cycle = cycle_num;
				max_height = height;
				max_root = i;
			}
			else if (cycle_num == max_cycle && height > max_height)
			{
				max_height = height;
				max_root = i;
			}
		}
		build_DAG(max_root, qg);
	}

	void get_out_neighbor(unsigned int id, vector<neighbor_info>& neighbor_list)
	{
		for (map<unsigned long long, vector<unsigned int> >::iterator iter = edges[id].label_edge_map.begin(); iter != edges[id].label_edge_map.end(); iter++)
		{
			for (int i = 0; i < iter->second.size(); i++)
			{
				neighbor_list.push_back(neighbor_info(iter->second[i], (iter->first >> 32), (iter->first & 0xFFFFFFFF)));
			}
		}
	}

	void get_in_neighbor(unsigned int id, vector<neighbor_info>& neighbor_list)
	{
		for (map<unsigned long long, vector<unsigned int> >::iterator iter = reverse_edges[id].label_edge_map.begin(); iter != reverse_edges[id].label_edge_map.end(); iter++)
		{
			for (int i = 0; i < iter->second.size(); i++)
			{
				neighbor_list.push_back(neighbor_info(iter->second[i], (iter->first >> 32), (iter->first & 0xFFFFFFFF)));
			}
		}
	}

	void get_labeled_out_neighbor(unsigned int id, int target_neighbor_label, int target_edge_label, vector<unsigned int>& neighbor_list)
	{
		unsigned long long label_pair = merge_long_long(target_neighbor_label, target_edge_label);
		if (edges[id].label_edge_map.find(label_pair) != edges[id].label_edge_map.end())
		{
			neighbor_list = edges[id].label_edge_map[label_pair];
		}
	}

	void get_labeled_in_neighbor(unsigned int id, int target_neighbor_label, int target_edge_label, vector<unsigned int>& neighbor_list)
	{
		unsigned long long label_pair = merge_long_long(target_neighbor_label, target_edge_label);
		if (reverse_edges[id].label_edge_map.find(label_pair) != reverse_edges[id].label_edge_map.end())
		{
			neighbor_list = reverse_edges[id].label_edge_map[label_pair];
		}
	}

	void get_unlabeled_out_neighbor(unsigned int id, vector<unsigned int>& vec)
	{
		vec = unlabeled_edges[id];
	}
	void get_unlabeled_in_neighbor(unsigned int id, vector<unsigned int>& vec)
	{
		vec = unlabeled_reverse_edges[id];
	}
	void get_leaf_nodes(vector<unsigned int>& vec)
	{
		for (int i = 0; i < unlabeled_edges.size(); i++)
		{
			if (unlabeled_edges[i].empty())
				vec.push_back(i);
		}
	}

};
