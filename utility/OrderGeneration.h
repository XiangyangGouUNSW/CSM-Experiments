#pragma once
#include<iostream>
#include<fstream>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<queue>
#include"../graph/query_graph.h"
using namespace std;

struct matching_order // a general structure of matching order, used for algorithms like Graphflow, Rapidflow and New_SP;
{
	vector<int> order;
	vector<int> parent; // parent is used when we use edge verification in candidate computation for a query vertex
	vector<int> pos;
	map<int, vector<pair<int, int>>> left_neighbors;
};

struct indexing_order // this order is used only for constructing local_indexes in Rapidflow
{
	unsigned int excluded_src;
	unsigned int excluded_dst;
	vector<unsigned int> src_neighbors; 
	vector<unsigned int> dst_neighbors;
	vector<unsigned int> order_of_remain;
	vector<vector<unsigned int>> left_neighbors;
};

struct operation_sequence
{
	vector<int> operation_id;
	vector<int> parent;
	vector<int> operation_type; //the bits from low bit to high bit represents 1. EXP(1) or CPT(0) 2. if the EXP is postponed (1 or 0) 3. If a visited check need to be carried before the EXP 4. If multi-expansion is appliable and we need to pre-compute 
	//the next CPT. 5. If multi-expansion has been applied and this CPT has been computed. 6. If the cover test can be omitted; 
	map<int, vector<pair<int, int>>> left_neighbors;
	map<int, int> CPT_reuse; // this array incidates if the CPT_reuse technique is appliable to a CPT, and if so, it stores the index of the reuse opt operation, amapping from a 
	map<int, vector<pair<int, int>>> compensate_edges;// If CPT reuse can be applied to an operation at a position, but the neighborset of the reused CPT is only a subset of the target CPT, thus we need to carry out a further filter. The vector at the position stored the edges that are absent in the reused CPU, and its label (for quick NLF).
	unsigned int tailEXP;
};


void left_neighbor_compute(matching_order* mo, query_graph* qg)
{
	for (int i = 0; i < mo->order.size(); i++)
	{
		unsigned int id = mo->order[i];
		vector<neighbor_info> vec;
		qg->get_neighbor(id, vec);
		for (int j = 0; j < vec.size(); j++)
		{
			if (mo->pos[vec[j].id] < i)
				mo->left_neighbors[id].push_back(make_pair(vec[j].id, vec[j].e_label));
		}
	}
}


// order computation with RI algorithm, used in Rapidflow and NewSP
void RI_increase_neighbor_score(unsigned int u, query_graph* g, int* neighbor_score, int* two_hop_score, int* remain, matching_order* mo)
{
	vector<neighbor_info> vec;
	g->get_neighbor(u, vec);
	unordered_set<unsigned int> increased;
	for (int i = 0; i < vec.size(); i++)
	{
		if (mo->pos[vec[i].id] == -1) {
			bool newly_connected = false;
			if (neighbor_score[vec[i].id] == 0)
				newly_connected = true;
			neighbor_score[vec[i].id]++;
			if (neighbor_score[u] == 0)
				remain[vec[i].id]--;
			vector<neighbor_info> neighbor_vec;
			g->get_neighbor(vec[i].id, neighbor_vec);
			for (int j = 0; j < neighbor_vec.size(); j++)
			{
				if (mo->pos[neighbor_vec[j].id] > -1)
					continue;
				if (increased.find(neighbor_vec[j].id) == increased.end()) {
					two_hop_score[neighbor_vec[j].id]++;
					increased.insert(neighbor_vec[j].id);//the two hop score of each vertex will only be increased once.
				}
				if (newly_connected)
					remain[neighbor_vec[j].id]--;
			}
		}
	}
}

pair<int, int> RI_choose_node(query_graph* g, int* neighbor_score, int* two_hop_score, int* remain, matching_order* mo)
{
	int max_neighbor = -1;
	int max_two_hop = -1;
	int max_remain = -1;
	int u = 0;
	for (int i = 0; i < g->g.size(); i++)
	{
		if (mo->pos[i] == -1 && g->g[i].id != 0xFFFFFFFF) {
			if (neighbor_score[i] > max_neighbor || (neighbor_score[i] == max_neighbor && two_hop_score[i] > max_two_hop) || (neighbor_score[i] == max_neighbor && two_hop_score[i] == max_two_hop && remain[i] > max_remain))
			{
				max_neighbor = neighbor_score[i];
				max_two_hop = two_hop_score[i];
				max_remain = remain[i];
				u = i;
			}
		}
	}
	vector<neighbor_info> vec;
	g->get_neighbor(u, vec);
	int min_pos = g->g.size();
	int parent = -1;
	for (int i = 0; i < vec.size(); i++)
	{
		if (mo->pos[vec[i].id] > -1 && mo->pos[vec[i].id] < min_pos)
		{
			min_pos = mo->pos[vec[i].id];
			parent = vec[i].id;
		}
	}
	return make_pair(u, parent);
}

void RI_order_generation(query_graph* g, matching_order* mo)
{
	//	cout << "generate order " << endl;
		//g->print_graph();
	int max_degree = -1;
	int max_id = -1;
	//int* pos = new int[g->g.size()];
	int* neighbor_score = new int[g->g.size()];
	int* two_hop_score = new int[g->g.size()];
	int* remain = new int[g->g.size()];
	mo->pos.resize(g->g.size());
	int arranged = 0;
	int invalid = 0;
	for (unsigned int i = 0; i < g->g.size(); i++)
	{
		if (g->g[i].id == 0xFFFFFFFF)
		{
			invalid++;
			remain[i] = -1;
			neighbor_score[i] = -1;
			two_hop_score[i] = -1;
		}
		else {
			remain[i] = g->get_total_degree(i);
			mo->pos[i] = -1;
			neighbor_score[i] = 0;
			two_hop_score[i] = 0;
			if (remain[i] > max_degree)
			{
				max_degree = remain[i];
				max_id = i;
			}
		}
	}

	unsigned int u = max_id;
	unsigned int parent = -1;

	while (arranged < g->g.size() - invalid)
	{
		mo->pos[u] = arranged;
		mo->order.push_back(u);
		mo->parent.push_back(parent);
		RI_increase_neighbor_score(u, g, neighbor_score, two_hop_score, remain, mo);
		pair<int, int> p = RI_choose_node(g, neighbor_score, two_hop_score, remain, mo);
		u = p.first;
		parent = p.second;
		arranged++;
	}
	left_neighbor_compute(mo, g);
} // given a query graph, genere order with RI algorithm.

void edge_oriented_RI_order(query_graph* g, matching_order* mo, unsigned int src, unsigned int dst)
{
	int* neighbor_score = new int[g->g.size()]; // number of 1-hop neighbors for each un-arranged node in the arranged node set.
	int* two_hop_score = new int[g->g.size()]; // number of 2-hop neighbors for each un-arranged node in the arranged node set.
	int* remain = new int[g->g.size()]; // number of 1-hop neighbors which are neither in the arranged node set nor neighbors of arranged nodes for each un-arranged node.s
	mo->pos.resize(g->g.size());
	unsigned int invalid = 0;
	for (unsigned int i = 0; i < g->g.size(); i++)
	{
		if (g->g[i].id == 0xFFFFFFFF) {
			remain[i] = -1;
			neighbor_score[i] = -1;
			two_hop_score[i] = -1;
			invalid++;
		}
		else {
			remain[i] = g->get_total_degree(i);
			mo->pos[i] = -1;
			neighbor_score[i] = 0;
			two_hop_score[i] = 0;
		}
	}
	mo->order.push_back(src);
	mo->order.push_back(dst);
	mo->parent.push_back(-1);
	mo->parent.push_back(src);
	mo->pos[src] = 0;
	mo->pos[dst] = 1;
	unsigned int arranged = 2;

	RI_increase_neighbor_score(src, g, neighbor_score, two_hop_score, remain, mo);
	RI_increase_neighbor_score(dst, g, neighbor_score, two_hop_score, remain, mo);


	while (arranged < g->g.size() - invalid)
	{
		pair<int, int> p = RI_choose_node(g, neighbor_score, two_hop_score, remain, mo);
		mo->pos[p.first] = arranged;
		mo->order.push_back(p.first);
		mo->parent.push_back(p.second);
		RI_increase_neighbor_score(p.first, g, neighbor_score, two_hop_score, remain, mo);
		arranged++;
	}
	left_neighbor_compute(mo, g);
} //RI order with 

// function that further transform RI order into operation sequence in NewSP, including 
void operation_sequence_generation(query_graph* g, matching_order* mo, operation_sequence* os, int qsrc = -1, int qdst = -1, bool use_index = false)
// here we need qsrc and qdst and used_index to distinguish different use cases, simple NewSp, NewSP+index, and NewSP+RapidFlow+Index. 
{
	os->operation_id = mo->order;
	int* right_neighbor = new int[mo->order.size()];
	int* degree = new int[mo->order.size()];
	unordered_map<int, int> parent_map;

	vector<vector<pair<int, int> > >buffered_neighbor_set;
	buffered_neighbor_set.resize(mo->order.size());

	//arrange CPT and EXP operations.
	for (int i = 0; i < mo->order.size(); i++)
	{
		unsigned int id = mo->order[i];
		parent_map[id] = mo->parent[i];
		vector<neighbor_info> vec;
		g->get_neighbor(id, vec);
		degree[id] = vec.size();
		unsigned int min_pos = mo->order.size();
		for (int j = 0; j < vec.size(); j++)
		{
			unsigned int neighbor_id = vec[j].id;
			if (mo->pos[neighbor_id] < i)
				os->left_neighbors[id].push_back(make_pair(neighbor_id, vec[j].e_label));
			if (mo->pos[neighbor_id] > i && mo->pos[neighbor_id] < min_pos) // compute the right most neighbor.
				min_pos = mo->pos[neighbor_id];
			buffered_neighbor_set[id].push_back(make_pair(neighbor_id, vec[j].e_label));
		}
		if (min_pos == mo->order.size()) {
			os->operation_id.push_back(id); // if no right neighbor, postpone the expansion of id to the end.
			right_neighbor[id] = -1;
		}
		else {
			right_neighbor[id] = mo->order[min_pos];
			for (vector<int>::iterator iter = os->operation_id.begin(); iter != os->operation_id.end(); iter++)
			{
				if (*iter == right_neighbor[id])
				{
					os->operation_id.emplace(iter--, id);
					break;
				}
			}
		}
	}

	assert(os->operation_id.size() == 2 * mo->order.size());
	os->parent.resize(os->operation_id.size());
	os->operation_type.resize(os->operation_id.size());

	// set operation type, including CPT/EXP and whether it is a postponed EXP.
	unordered_set<int> visited;
	for (int i = 0; i < os->operation_id.size(); i++)
	{
		int id = os->operation_id[i];
		int operation = 0;
		if (visited.find(id) != visited.end()) { // if the vertex has shown up in the operation sequence, this time the operation is EXP;
			operation = operation | 1;
			if (os->operation_id[i - 1] != id) { // this means the EXP is not just behind CPT, it is a postponed EXP
				operation |= (1 << 1);
				for (int j = i - 1; j >= 0; j--)
				{
					if ((os->operation_type[j] & 1) && g->g[os->operation_id[j]].label == g->g[id].label) { // this means that we need to recheck the visited status of each node in this EXP, as we have expanded other vertices with the same label after CPT.
						operation |= (1 << 2);
						break;
					}
					else if (os->operation_id[j] == id)
						break;
				}
			}
		}
		else {
			visited.insert(id);
		}
		os->operation_type[i] = operation;
		os->parent[i] = parent_map[id];
	}

	//set operation type, including multi-EXP information and cover test omit.
	for (int i = 0; i < os->operation_id.size(); i++)
	{
		if ((os->operation_type[i] & 1) == 0) // a CPT
		{
			if (i + 2 < os->operation_id.size() && (os->operation_type[i + 1] & 1) && (os->operation_type[i + 2] & 1)) // the following 2 operations are EXP, then we can carry out multi-expansion here.
			{
				for (int j = i + 3; j < os->operation_type.size(); j++)
				{
					if ((os->operation_type[j] & 1) == 0) // find next CPT
					{
						unsigned int me_parent = os->parent[j];
						bool parent_expanded = false;
						for (int k = 0; k < i; k++)
						{
							if (os->operation_id[k] == me_parent && (os->operation_type[k] & 1) == 1)
							{
								parent_expanded = true;
								break;
							}
						}
						if(parent_expanded)// multi-expansion can be supported only when there is at least one expanded neighbor of the next CPT vertex
						{
							os->operation_type[i] = (os->operation_type[i] | (1 << 3));
							os->operation_type[j] = (os->operation_type[j] | (1 << 4));
						}
						break;
					}
				}
			}
		}
		if (right_neighbor[os->operation_id[i]] == -1) { // if no right neighbor, we donot need to carry out the cover test, as all it neighbor has been matched.
			os->operation_type[i] = (os->operation_type[i] | (1 << 5));
		}
	}

	// set the index to begin tail EXP, tail EXP can be used to expand the final EXPs together. These EXPs are an independent set of vertices and can be expanded by a simple cartesian product (need to check visited status) 
	os->tailEXP = -1;
	for (int i = os->operation_type.size() - 1; i >= 0; i--)
	{
		if ((os->operation_type[i] & 1) == 0)
		{
			os->tailEXP = i + 1;
			break;
		}
	}

	// compute the CPT reuse status, this is where use_index flag is used. If we use index, it means the candidate set has been filtered by at least NLF rule, or more complex rules. 
	// in this case ,CPT reuse cannot applied only when the the entire neighborset is a subset, otherwise we only need to check the left neighbor.
	for (int i = 0; i < mo->order.size(); i++) {
		int id = mo->order[i];
		unordered_map<int, int> neighbor_set;
		unordered_map<int, int> local_left_neighbors;
		for (int j = 0; j < os->left_neighbors[id].size(); j++)
			local_left_neighbors[os->left_neighbors[id][j].first] = os->left_neighbors[id][j].second;
		for (int j = 0; j < buffered_neighbor_set[id].size(); j++)
			neighbor_set[buffered_neighbor_set[id][j].first] = buffered_neighbor_set[id][j].second;
		int max_common_neighbor = 0;
		int mcn_id = 0;
		int mcn_pos = 0;
		for (int j = 0; j < i; j++)
		{
			unsigned int candidate = mo->order[j];
			if (candidate == qsrc || candidate == qdst) // the CPT os qsrc and qdst cannot be reused.
				continue;
			if (g->g[candidate].label != g->g[id].label)
				continue;
			//	if (right_neighbor[candidate] <= i)
				//	continue; // if the left-most right neighbor of this candiadte is before i, it means this candidate hash been expanded upon CPT of i. in this case, we cannot reuse CPT of this candidate? 
			if (os->left_neighbors[candidate].size() > local_left_neighbors.size())
				continue;
			bool filtered = false;
			for (int k = 0; k < os->left_neighbors[candidate].size(); k++)
			{
				if (local_left_neighbors.find(os->left_neighbors[candidate][k].first) == local_left_neighbors.end() || local_left_neighbors[os->left_neighbors[candidate][k].first] != os->left_neighbors[candidate][k].second) {
					filtered = true;
					break;
				}
			}
			if (filtered)
				continue;
			if (use_index || ((os->operation_type[j] & (1 << 4)) > 0)) // if index is used (1-hop index), candidates are generated from the index, and NLF check must have been conducted,
				// If CPT[j] is pre-computed in multi-expansion technique, NLF check also has been pre-computed. Otherwise, NLF check of candidate will be postponed to its expansion, which has not been conduced yet
			{
				if (degree[candidate] > degree[id])
					continue;
				for (int k = 0; k < buffered_neighbor_set[candidate].size(); k++)
				{
					if (neighbor_set.find(buffered_neighbor_set[candidate][k].first) == neighbor_set.end() || neighbor_set[buffered_neighbor_set[candidate][k].first] != buffered_neighbor_set[candidate][k].second) {
						filtered = true;
						break;
					}
				}
			}
			if (filtered)
				continue;
			if (os->left_neighbors[candidate].size() > max_common_neighbor)
			{
				max_common_neighbor = os->left_neighbors[candidate].size();
				mcn_id = candidate;
				mcn_pos = j;
				if (max_common_neighbor == local_left_neighbors.size())
					break;
			}
		}

		if (max_common_neighbor != 0)
		{
			os->CPT_reuse[id] = mcn_id;
			if (max_common_neighbor < local_left_neighbors.size())
			{
				for (int j = 0; j < os->left_neighbors[mcn_id].size(); j++)
					local_left_neighbors.erase(os->left_neighbors[mcn_id][j].first);
				for (unordered_map<int, int>::iterator iter = local_left_neighbors.begin(); iter != local_left_neighbors.end(); iter++) {
					os->compensate_edges[id].push_back(make_pair(iter->first, iter->second));
				}
			}
			if ((os->operation_type[mcn_pos] & (1 << 4)) > 0 && degree[mcn_id] == degree[id]) // here is a special case if 1. the target of CPT_reuse is a pre-computed CPT, namely its NLF test has been carried out, 2. the neighborhood of these two node are the same,in this case, the NLF check can be omitted.
				os->operation_type[i] = (os->operation_type[i] | (1 << 5));
		}
	}
}

// function to compute the order of building order to 
void RF_indexing_order_generation(query_graph* g, indexing_order* io, unsigned int excluded_src, unsigned int excluded_dst) // the order to compute local index in rapidflow. 
{
	io->excluded_src = excluded_src;
	io->excluded_dst = excluded_dst;
	io->left_neighbors.resize(g->g.size());
	unordered_set<unsigned int> arranged;
	arranged.insert(excluded_src);
	arranged.insert(excluded_dst);
	g->get_unlabeled_neighbor(excluded_src, io->src_neighbors);
	g->get_unlabeled_neighbor(excluded_dst, io->dst_neighbors);
	for (int i = 0; i < io->src_neighbors.size(); i++) {
		if (io->src_neighbors[i] == excluded_dst)
		{
			io->src_neighbors.erase(io->src_neighbors.begin() + i);
			i--;
		}
		else
			arranged.insert(io->src_neighbors[i]);
	}
	for (int i = 0; i < io->dst_neighbors.size(); i++) {

		if (io->dst_neighbors[i] == excluded_src)
		{
			io->dst_neighbors.erase(io->dst_neighbors.begin() + i);
			i--;
		}
		else
			arranged.insert(io->dst_neighbors[i]);
	}
	while (arranged.size() < g->g.size()) {
		unsigned int max_neighbor_num = 0;
		unsigned int max_id = 0;
		for (int i = 0; i < g->g.size(); i++)
		{
			unsigned int id = g->g[i].id;
			if (arranged.find(id) != arranged.end())
				continue;
			unsigned int neighbor_num = 0;
			vector<unsigned int> vec;
			g->get_unlabeled_neighbor(id, vec);
			for (int j = 0; j < vec.size(); j++)
			{
				if (arranged.find(vec[j]) != arranged.end())
					neighbor_num++;
			}
			vec.clear();
			if (neighbor_num > max_neighbor_num)
			{
				max_neighbor_num = neighbor_num;
				max_id = id;
			}
		}

		io->order_of_remain.push_back(max_id);
		arranged.insert(max_id);
		vector<unsigned int> vec;
		g->get_unlabeled_neighbor(max_id, vec);
		for (int i = 0; i < vec.size(); i++)
		{
			if (arranged.find(vec[i]) != arranged.end())
				io->left_neighbors[max_id].push_back(vec[i]);
		}
	}

}

// order computation for GraphFlow, which is a simplified version of RI.

void GF_increase_score(unsigned int u, query_graph* g, int* neighbor_score, matching_order* mo)
{
	vector<unsigned int> vec;
	g->get_unlabeled_neighbor(u, vec);
	for (int i = 0; i < vec.size(); i++)
	{
		if (mo->pos[vec[i]] == -1) // not ordered yet;
			neighbor_score[vec[i]]++;
	}
	return;
}

pair<int, int> GF_choose_node(query_graph* g, int* neighbor_score, int* degree, matching_order* mo)
{
	unsigned int max_score = 0;
	unsigned int max_degree = 0;
	unsigned int max_id = 0;
	for (int i = 0; i < g->g.size(); i++)
	{
		if (mo->pos[i] == -1)
		{
			if (neighbor_score[i] > max_score)
			{
				max_id = i;
				max_score = neighbor_score[i];
				max_degree = degree[i];
			}
			else if (neighbor_score[i] == max_score && degree[i] > max_degree)
			{
				max_id = i;
				max_degree = degree[i];
			}
		}
	}
	unsigned int u = max_id;
	vector<unsigned int> vec;
	g->get_unlabeled_neighbor(u, vec);
	int min_pos = g->g.size();
	int parent = -1;
	for (int i = 0; i < vec.size(); i++)
	{
		if (mo->pos[vec[i]] > -1 && mo->pos[vec[i]] < min_pos)
		{
			min_pos = mo->pos[vec[i]];
			parent = vec[i];
		}
	}
	return make_pair(u, parent);
}
void GraphFlow_order_generation(query_graph* g, matching_order* mo, unsigned int src, unsigned int dst) // this function is slightly different from RI order, when two vertices that have not been ordered have the same number of ordered neighbors,
// RI will first compare the ordered 2-hop neighbors, and then the remain neighbor number (which are neither ordered or neighbor of ordered vertices). But Graphflow will directly compare the total neighbor number. 
{
	int* neighbor_score = new int[g->g.size()]; // number of 1-hop neighbors for each un-arranged node in the arranged node set.
	int* degree = new int[g->g.size()];
	mo->pos.resize(g->g.size());
	unsigned int invalid = 0;
	for (unsigned int i = 0; i < g->g.size(); i++)
	{
		degree[i] = g->get_total_degree(i);
		mo->pos[i] = -1;
		neighbor_score[i] = 0;
	}
	mo->order.push_back(src);
	mo->order.push_back(dst);
	mo->parent.push_back(-1);
	mo->parent.push_back(src);
	mo->pos[src] = 0;
	mo->pos[dst] = 1;
	unsigned int arranged = 2;

	GF_increase_score(src, g, neighbor_score, mo);
	GF_increase_score(dst, g, neighbor_score, mo);


	while (arranged < g->g.size() - invalid)
	{
		pair<int, int> p = GF_choose_node(g, neighbor_score, degree, mo);
		mo->pos[p.first] = arranged;
		mo->order.push_back(p.first);
		mo->parent.push_back(p.second);
		GF_increase_score(p.first, g, neighbor_score, mo);
		arranged++;
	}
	left_neighbor_compute(mo, g);

}



struct ks_order
{
	matching_order mo; // the matching order of kernel vertices, generated with RI algorithm;
	vector<vector<unsigned int>> shell_to_test;// thie vector record the corresponding shell vertex to test for each pos i, where the i_th kernel vertex in mo is the last neighbor of this shell vertex.
	set<unsigned int> shell_vertices;
	set<unsigned int> kernel_vertices;
	map<unsigned int, vector<pair<int, int>>> shell_neighbors;
};

void find_unarranged(ks_order* ko, vector<unsigned int>& vec, unordered_set<unsigned int> &candidates)
{
	for (int i = 0; i < vec.size(); i++)
	{
		unsigned int u = vec[i];
		if (ko->kernel_vertices.find(u) == ko->kernel_vertices.end() && ko->shell_vertices.find(u) == ko->shell_vertices.end())
			candidates.insert(u);
	}
}
bool check_core_neighbor(ks_order* ko, vector<unsigned int>& vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		if (ko->kernel_vertices.find(vec[i]) == ko->kernel_vertices.end())
			return false;
	}
	return true;
}
int count_core_neighbor(ks_order* ko, vector<unsigned int>& vec)
{
	unsigned int cnt = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		if (ko->kernel_vertices.find(vec[i]) != ko->kernel_vertices.end())
			cnt++;
	}
	return cnt;
}
void CailiG_order_generation(query_graph* g, ks_order* ko, unsigned int src, unsigned int dst) // I use the method in the opensource code of CaliG here. It is slightly different
// from the greedy method described in the paper. But the method in the paper is described too briefly and I cannot figure out how it exactly works. Therefore I use the 
//exactly same method in the opensource code. 
{
	unordered_set<unsigned int> candidates;
	vector<unsigned int> vec;
	ko->kernel_vertices.insert(src);
	ko->kernel_vertices.insert(dst);
	ko->mo.order.push_back(src);
	ko->mo.order.push_back(dst);
	ko->shell_to_test.push_back(vector<unsigned int>());
	ko->shell_to_test.push_back(vector<unsigned int>());
	ko->mo.parent.push_back(-1);
	ko->mo.parent.push_back(-1);
	ko->mo.pos.resize(g->g.size());
	ko->mo.pos[src] = 0;
	ko->mo.pos[dst] = 1;
	g->get_unlabeled_neighbor(src, vec);
	find_unarranged(ko, vec, candidates);
	vec.clear();
	g->get_unlabeled_neighbor(dst, vec);
	find_unarranged(ko, vec, candidates);
	vec.clear();
	while (!candidates.empty()) {
		for (auto iter = candidates.begin(); iter != candidates.end();)
		{
			vec.clear();
			unsigned int u = *iter;
			g->get_unlabeled_neighbor(u, vec);
			if (check_core_neighbor(ko, vec))
			{
				sort(vec.begin(), vec.end());
				unsigned int right_neighbor = 0;
				for (int i = 0; i < ko->mo.order.size(); i++)
				{
					unsigned int u2 = ko->mo.order[i];
					auto temp_iter = lower_bound(vec.begin(), vec.end(), u2);
					if (temp_iter != vec.end() && *temp_iter == u2)
						right_neighbor = i;
				}
				ko->shell_vertices.insert(u);
				ko->shell_to_test[right_neighbor].push_back(u);
				for (int i = 0; i < vec.size(); i++)
					ko->shell_neighbors[u].push_back(make_pair(vec[i], g->check_edge_label(vec[i], u)));
				iter = candidates.erase(iter);
			}
			else
				iter++;
		}
		if (candidates.empty())
			break;
		unsigned int maximum_degree = 0;
		unsigned int maximum_kn = 0;
		unsigned int max_cand = 0;
		for (auto iter = candidates.begin(); iter != candidates.end(); iter++)
		{
			vec.clear();
			unsigned int u = *iter;
			g->get_unlabeled_neighbor(u, vec);
			unsigned int d = check_core_neighbor(ko, vec);
			if (d > maximum_kn)
			{
				maximum_kn = d;
				max_cand = u;
			}
			else if (d == maximum_kn)
			{
				if (g->get_total_degree(u) > maximum_degree)
				{
					maximum_degree = g->get_total_degree(u);
					max_cand = u;
				}
			}
		/*	unsigned int degree = g->get_total_degree(u);
			if (g->check_edge_label(u, src) != -1)
				degree--;
			if (g->check_edge_label(u, dst) != -1)
				degree--;
			if (degree > maximum_degree)
			{
				maximum_degree = degree;
				max_cand = u;
			}*/
		}
		vec.clear();
		unsigned int u = max_cand;
		ko->mo.pos[u] = ko->mo.order.size();
		ko->mo.order.push_back(u);
		ko->shell_to_test.push_back(vector<unsigned int>());
		ko->kernel_vertices.insert(u);
		g->get_unlabeled_neighbor(u, vec);
		sort(vec.begin(), vec.end());
		int p = -1;
		for (int i = 0; i < ko->mo.order.size(); i++)
		{
			unsigned int u2 = ko->mo.order[i];
			auto temp_iter = lower_bound(vec.begin(), vec.end(), u2);
			if (temp_iter != vec.end() && *temp_iter == u2)
			{
				if (p == -1)
					p = u2;
				ko->mo.left_neighbors[u].push_back(make_pair(u2, g->check_edge_label(u, u2)));
			}
		}
		ko->mo.parent.push_back(p);
		find_unarranged(ko, vec, candidates);
		candidates.erase(u);
		vec.clear();
	}

}






