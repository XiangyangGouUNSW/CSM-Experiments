#pragma once
#include<iostream>
#include<fstream>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<queue>
#include"../graph/query_graph.h"
#include"../utility/struct.h"
#include"OrderGeneration.h"
using namespace std;


void DFS_match(matching_order* order,unordered_map<unsigned int, unsigned int>& tmp_maping, unordered_set<unsigned int>& used, unsigned int index, query_graph* g, unordered_map<unsigned int, vector<unsigned int>>& match)
{
	unsigned int node_id = order->order[index];
	unsigned int parent_id = order->parent[index];

	unsigned int parent_candidate = tmp_maping[parent_id];
	int candidate_label = g->g[node_id].label;
	int edge_label = g->check_edge_label(parent_id, node_id);

	vector<unsigned int> vec;
	g->get_labeled_neighbor(parent_candidate, candidate_label, edge_label, vec);
	for (int i = 0;i<vec.size();i++)
	{
		unsigned int candidate_id = vec[i];
		if (used.find(candidate_id) != used.end())
			continue;

		bool filtered = false;;
		for (map<unsigned long long, int>::iterator NLF_iter = g->NLF[node_id].begin(); NLF_iter != g->NLF[node_id].end(); NLF_iter++)
		{
			if (g->NLF[candidate_id].find(NLF_iter->first) == g->NLF[candidate_id].end() || g->NLF[candidate_id][NLF_iter->first] < NLF_iter->second)
			{
				filtered = true;
				break;
			}
		}
		if (filtered)
			continue;

		//neighbor check;

		vector<neighbor_info> neighbor_vec;
		g->get_neighbor(node_id, neighbor_vec);
		for (unsigned int i = 0; i < neighbor_vec.size(); i++)
		{
			unsigned int neigbor_id = neighbor_vec[i].id;
			if (tmp_maping.find(neigbor_id) != tmp_maping.end())
			{
				if (g->check_edge_label(candidate_id, tmp_maping[neigbor_id]) != neighbor_vec[i].e_label) {
					filtered = true;
					break;
				}
			}
		}
		
		if (!filtered)
		{
			if (index == order->order.size() - 1)
			{
				for (unordered_map<unsigned int, unsigned int>::iterator map_iter = tmp_maping.begin(); map_iter != tmp_maping.end(); map_iter++)
					match[map_iter->first].push_back(map_iter->second);
				match[node_id].push_back(candidate_id);
			}
			else
			{
				tmp_maping[node_id] = candidate_id;
				used.insert(candidate_id);
				DFS_match(order, tmp_maping, used, index + 1, g, match);
				used.erase(candidate_id);
				tmp_maping.erase(node_id);
			}
		}
	}
	return;
}

void automorphism(query_graph* g, unordered_map<unsigned int, vector<unsigned int>> &match) // algorithm to compute automorphism of the query graph, so that we can group the edges in the query graph into equal sets.
{
	matching_order* order = new matching_order;
	unordered_map<unsigned int, unsigned int> tmp_maping;
	unordered_set<unsigned int> used;
	RI_order_generation(g, order);
	unsigned int root = order->order[0];
	vector<unsigned int> cand;
	for (int i = 0; i < g->g.size(); i++)
	{
		unsigned int v = i;
		if (g->g[i].label!=g->g[root].label)
			continue;
		bool filtered = false;
		for (auto iter = g->NLF[root].begin(); iter != g->NLF[root].end(); iter++)
		{
			if (g->NLF[v].find(iter->first) == g->NLF[v].end() || g->NLF[v][iter->first] < iter->second)
			{
				filtered = true;
				break;
			}
		}
		if (!filtered)
			cand.push_back(v);
	}
	for (int i = 0; i < cand.size(); i++) {
		tmp_maping[root] = cand[i];
		used.insert(cand[i]);
		DFS_match(order, tmp_maping, used, 1, g, match);
		tmp_maping.erase(root);
		used.erase(cand[i]);
	}
}

void generate_autoset(query_graph* g, unordered_map<unsigned int, vector<unsigned int>>& match, map<unsigned long long, vector<pair<unsigned long long, automatch*>>>& auto_set, vector<automatch> &automatches)
{
	unsigned int mapping_number = match.begin()->second.size();
	for (int i = 0; i < mapping_number; i++)
	{
		automatch au;
		for (unordered_map<unsigned int, vector<unsigned int> >::iterator iter = match.begin(); iter != match.end(); iter++)
			au.mapping[iter->first] = iter->second[i];
		automatches.push_back(au);
	}
	
	unordered_set<unsigned long long> selected;
	for (int i =0;i<g->g.size();i++)
	{
		vector<neighbor_info> neigbor_vec;
		g->get_neighbor(i, neigbor_vec);
		for (int j = 0; j < neigbor_vec.size(); j++)
		{
			unsigned int src = i;
			unsigned int dst = neigbor_vec[j].id;
			unsigned long long edge = merge_long_long(src, dst);
			if (selected.find(edge) != selected.end())
				continue;
			else
			{
				vector<pair<unsigned long long, automatch*> > vec;
				for (int k = 0; k < mapping_number; k++)
				{
					unsigned int candidate_s = match[src][k];
					unsigned int candidate_d = match[dst][k];
					unsigned long long matched_edge = merge_long_long(candidate_s, candidate_d);
					if (selected.find(matched_edge) == selected.end()) {
						vec.push_back(make_pair(matched_edge, &automatches[k]));
						selected.insert(matched_edge);
					}
				}
				auto_set[edge] = vec;
			}
		}
	}
}
