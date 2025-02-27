#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include<map>
#define max_int 0xFFFFFFFF
#define merge_long_long(s, d) (((long long)s<<32)|d)
#include<queue>
#include<algorithm> 
#include<set>
using namespace std;

typedef unsigned int v_bitmap; // we use vertex bitmap + edge veiw as the base of all these indexes, where each data vertex is associated with 
// a bitmap, the i_th bit indicating if this data vertex is a candidate of the i_th query vertex. if the query graph is small, v_bitmap can be changed to
// short or char to save space.

int byte_align(int n, int s)
{
	if (n % s == 0)
		return n;
	else
	{
		int a = n / s;
		return (a + 1) * s;
	}
}
struct neighbor_info
{
	unsigned int id;
	int v_label;
	int e_label;
	neighbor_info(unsigned int v, int vertex_label, int edge_label)
	{
		id = v;
		e_label = edge_label;
		v_label = vertex_label;
	}
};

struct node
{
	int id;
	int label;
	map<unsigned long long, vector<unsigned int>> label_edge_map;
	node(int id_=-1, int label_=-1) {
		id = id_;
		label = label_;
	}
	int memory_compute()
	{
		unsigned int mem = sizeof(int) * 2 + sizeof(label_edge_map);
		unsigned int kv_size = sizeof(unsigned long long) + sizeof(vector<unsigned int>);
		kv_size = byte_align(kv_size, 8);
		mem += label_edge_map.size() * (sizeof(void*) * 2 + kv_size);
		return mem;
	}
	int get_original_id()
	{
		return id;
	}
	void set_original_id(unsigned int val)
	{
		id = val;
	}
	void clear()
	{
		for (map<unsigned long long, vector<unsigned int>>::iterator iter = label_edge_map.begin(); iter != label_edge_map.end(); iter++)
			iter->second.clear();
		label_edge_map.clear();
	}

};

struct edge
{
	unsigned int src_id;
	unsigned int dst_id;
	int label;
	edge(unsigned int s, unsigned int d, int l)
	{
		src_id = s;
		dst_id = d;
		label = l;
	}

};

struct edge_hash
{
	size_t operator()(edge e)
	{
		hash<unsigned long long> hasher1;
		hash<int> hasher2;
		unsigned long long id_pair = merge_long_long(e.src_id, e.dst_id);
		size_t result = (hasher1(id_pair) ^ hasher2(e.label));
		return result;
	}
};

class edge_view
{
public:
	unordered_map<unsigned int, vector<unsigned int> > view; // the key is vertex ID, the value is vector of neighbors where the edge between each neighbor and the vertex is a candidate of thee query edge. 
	edge_view() {};
	~edge_view() {
		for (unordered_map<unsigned int, vector<unsigned int> >::iterator iter = view.begin(); iter != view.end(); iter++)
			iter->second.clear();
		view.clear();
	}
	void clear() {
		for (unordered_map<unsigned int, vector<unsigned int> >::iterator iter = view.begin(); iter != view.end(); iter++)
			iter->second.clear();
		view.clear();
	}
	unsigned int get_size()
	{
		unsigned int sum = 0;
		for (auto iter = view.begin(); iter != view.end(); iter++)
			sum += iter->second.size();
		return sum;
	}
	void insert_edge(unsigned int src, unsigned int dst)
	{
		if (view.find(src) != view.end())
		{
			vector<unsigned int>::iterator iter = lower_bound(view[src].begin(), view[src].end(), dst);
			if (iter == view[src].end() || *iter != dst)
				view[src].emplace(iter, dst);
		}
		else
		{
			view[src] = vector<unsigned int>();
			view[src].push_back(dst);
		}
	}
	bool delete_edge(unsigned int src, unsigned int dst)
	{
		if (view.find(src) != view.end())
		{
			vector<unsigned int>::iterator iter = lower_bound(view[src].begin(), view[src].end(), dst);
			if (iter != view[src].end() && *iter == dst) {
				view[src].erase(iter);
				if (view[src].empty())
				{
					view.erase(src);
					return true;
				}
			}
		}
		return false;
	}
	bool check_edge(unsigned int src, unsigned int dst)
	{
		if (view.find(src) != view.end())
		{
			vector<unsigned int>::iterator iter = lower_bound(view[src].begin(), view[src].end(), dst);
			if (iter != view[src].end() && *iter == dst)
				return true;
			else
				return false;
		}
		else
			return false;
	}
	int get_neighbor(unsigned int src, vector<unsigned int>& vec)
	{
		if (view.find(src) != view.end())
		{
			vec = view[src];
		//	for (int i = 0; i < view[src].size(); i++)
			//	vec.push_back(view[src][i]);
			return view[src].size();
		}
		else
			return 0;
	}
	bool check_node(unsigned int id)
	{
		return view.find(id) != view.end();
	}
	void delete_node(unsigned int id)
	{
		if (view.find(id) != view.end())
		{
			view[id].clear();
			view.erase(id);
		}
	}

	unsigned int memory_compute()
	{
		unsigned int mem = sizeof(view);
		mem += view.bucket_count() * sizeof(void*);
		unsigned int mem_per_kv = byte_align((sizeof(vector<unsigned int>) + sizeof(void*) + sizeof(unsigned int)), 8);
		mem += view.size() * mem_per_kv;
		for (auto iter = view.begin(); iter != view.end(); iter++)
			mem += iter->second.capacity() * sizeof(unsigned int);
		return mem;
	}
};


class global_index_template
{
public:
	vector<v_bitmap> candidate_map; // each vertex in the data graph is associated with a v_bitmap, where the i_th bit is 1 indicates that this vertex is a candidate for query vertex i
	unordered_map<unsigned long long, edge_view> edges; // map between each query edge and an edge view. unsigned long long is merged as src id + dst id;

	global_index_template() {};
	~global_index_template()
	{
		for (unordered_map<unsigned long long, edge_view>::iterator iter = edges.begin(); iter != edges.end(); iter++)
			iter->second.clear();
		candidate_map.clear();
		edges.clear();
	}

	void get_neighbor(unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
		unordered_map<unsigned long long, edge_view>::iterator iter = edges.find(edge_pair);
		if (iter != edges.end())
		{
			iter->second.get_neighbor(candidate_id, neighbor_list);
		}
	}

	void get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		unordered_map<unsigned long long, edge_view>::iterator iter = edges.find(edge_pair);
		if (iter != edges.end())
		{
			iter->second.get_neighbor(candidate_id, neighbor_list);
		}
	}
	bool check_edge(unsigned int src_qid, unsigned int dst_qid, unsigned int src, unsigned int dst)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		if (edges.find(edge_pair) != edges.end() && edges[edge_pair].check_edge(src, dst))
			return true;
		else
			return false;
	}

	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		for (int i = 0; i < candidate_map.size(); i++)
		{
			if (((candidate_map[i] >> id) & 1) == 1)
				vec.push_back(i);
		}
	}

	bool check_candidate(unsigned int id, unsigned int candidate)
	{
		if (((candidate_map[candidate] >> id) & 1) == 1)
			return true;
		else
			return false;
	}

	void set_candidate(unsigned int id, unsigned int candidate)
	{
		if(id<candidate_map.size())
			candidate_map[candidate] = (candidate_map[candidate] | (1 << id));
		else
		{
			v_bitmap val = (1 << id);
			candidate_map.push_back(val);
		}
	}
	void erase_candidate(unsigned int id, unsigned int candidate)
	{
		candidate_map[candidate] = (candidate_map[candidate] & (~(1 << id)));
	}
	int memory_compute()
	{
		unsigned int mem = sizeof(candidate_map) + sizeof(edges);
		mem += candidate_map.capacity() * sizeof(v_bitmap);
		mem += edges.bucket_count() * sizeof(void*);
		mem += edges.size() * (sizeof(void*) + sizeof(unsigned long long));
		for (auto iter = edges.begin(); iter != edges.end(); iter++)
			mem += iter->second.memory_compute();
		return mem;
	}
};

class local_index_template
{
public:
	vector<vector<unsigned int>> candidate_map; // each vertex in the data graph is associated with a v_bitmap, where the i_th bit is 1 indicates that this vertex is a candidate for query vertex i
	unordered_map<unsigned long long, edge_view> edges; // map between each query edge and an edge view. unsigned long long is merged as src id + dst id;

	local_index_template() {};
	~local_index_template()
	{
		for (unordered_map<unsigned long long, edge_view>::iterator iter = edges.begin(); iter != edges.end(); iter++)
			iter->second.clear();
		for (int i = 0; i < candidate_map.size(); i++)
			candidate_map[i].clear();
		candidate_map.clear();
		edges.clear();
	}

	void get_neighbor(unsigned int candidate_id, unsigned long long edge_pair, vector<unsigned int>& neighbor_list)
	{
		unordered_map<unsigned long long, edge_view>::iterator iter = edges.find(edge_pair);
		if (iter != edges.end())
		{
			iter->second.get_neighbor(candidate_id, neighbor_list);
		}
	}

	void get_neighbor(unsigned int src_qid, unsigned int dst_qid, unsigned int candidate_id, vector<unsigned int>& neighbor_list)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		unordered_map<unsigned long long, edge_view>::iterator iter = edges.find(edge_pair);
		if (iter != edges.end())
		{
			iter->second.get_neighbor(candidate_id, neighbor_list);
		}
	}
	bool check_edge(unsigned int src_qid, unsigned int dst_qid, unsigned int src, unsigned int dst)
	{
		unsigned long long edge_pair = merge_long_long(src_qid, dst_qid);
		if (edges.find(edge_pair) != edges.end() && edges[edge_pair].check_edge(src, dst))
			return true;
		else
			return false;
	}

	void get_candidates(unsigned int id, vector<unsigned int>& vec)
	{
		vec = candidate_map[id];
	}

	bool check_candidate(unsigned int id, unsigned int candidate)
	{
		auto iter = lower_bound(candidate_map[id].begin(), candidate_map[id].end(), candidate);
		if(iter!=candidate_map[id].end()&&*iter==candidate)
			return true;
		else
			return false;
	}

	void set_candidate(unsigned int id, unsigned int candidate)
	{
		auto iter = lower_bound(candidate_map[id].begin(), candidate_map[id].end(), candidate);
		if (iter == candidate_map[id].end() || *iter != candidate)
			candidate_map[id].emplace(iter, candidate);
	}
	void erase_candidate(unsigned int id, unsigned int candidate)
	{
		auto iter = lower_bound(candidate_map[id].begin(), candidate_map[id].end(), candidate);
		if (iter != candidate_map[id].end() && *iter == candidate)
			candidate_map[id].erase(iter);
	}
};



struct match_result
{
	vector<unsigned int> auto_header;
	vector<unsigned int> matched_data;
};

struct automatch
{
	unordered_map<unsigned int, unsigned int> mapping;
	~automatch()
	{
		mapping.clear();
	}
};

class nei_cnt
{
public: unordered_map<unsigned int, unsigned int> nmap;
	  unsigned int target_u;

};

class bipartite_graph
{
public:
	map<unsigned int, vector<unsigned int>> cand_map; // map from query node to candidate data node, 
	unsigned int target_u;
	bipartite_graph(unsigned int query_node_number_ = 0) {
		target_u = query_node_number_;
	}
	~bipartite_graph()
	{
		for (auto iter = cand_map.begin(); iter != cand_map.end(); iter++)
			iter->second.clear();
		cand_map.clear();
	}
	void print_match()
	{
		for (auto iter = cand_map.begin(); iter != cand_map.end(); iter++) {
			cout << "match of " << iter->first << endl;
			for (int j = 0; j < iter->second.size(); j++)
				cout << iter->second[j] << ' ';
			cout << endl;
		}
	}
	void add_edge(unsigned int query_node, unsigned int data_node)
	{
		unsigned int  v = data_node;
		unsigned int u = query_node;
		if (cand_map.find(u) == cand_map.end())
			cand_map[u] = vector<unsigned int>();
		auto iter = lower_bound(cand_map[u].begin(), cand_map[u].end(), v);
		if (iter == cand_map[u].end() || *iter != v)
			cand_map[u].emplace(iter, v);
	}

	bool delete_edge(unsigned int query_node, unsigned int data_node)
	{
		unsigned int v = data_node;
		unsigned int u = query_node;
		if (cand_map.find(u) == cand_map.end())
			return false;
		auto iter = lower_bound(cand_map[u].begin(), cand_map[u].end(), v);
		if (iter != cand_map[u].end() && *iter == v) {
			cand_map[u].erase(iter);
			if (cand_map[u].empty())
				cand_map.erase(u);
			return true;
		}
		else
			return false;
	}

	bool bfs(map<unsigned int, unsigned int>& data_dist, map<unsigned int, unsigned int>& query_dist,
		map<unsigned int, unsigned int>& data_match, map<unsigned int, unsigned int>& query_match, unsigned int &min_dist)
	{
	//	cout << "bfs start" << endl;
		data_dist.clear();
		query_dist.clear();
		min_dist = 0xFFFFFFFF;
		queue<unsigned int> q;
		for (auto iter = cand_map.begin(); iter != cand_map.end(); iter++)
		{
			unsigned int u = iter->first;
			if (query_match.find(u) == query_match.end())
			{
				q.push(u);
				query_dist[u] = 0;
			}
		}
		while (!q.empty())
		{
			unsigned int u = q.front();
			q.pop();
			if (query_dist[u] > min_dist)
				break;
			for (int i = 0; i < cand_map[u].size(); i++)
			{
				unsigned int v = cand_map[u][i];
				if (data_dist.find(v) == data_dist.end())
				{
					data_dist[v] = query_dist[u] + 1;
					if (data_match.find(v) == data_match.end()) {
						min_dist = data_dist[v];
					}
					else
					{
						query_dist[data_match[v]] = data_dist[v] + 1;
						q.push(data_match[v]);
					}

				}
			}
		}
	//	cout << "bfs finished" << endl;
		return (min_dist != 0xFFFFFFFF);
	}
	bool traverse_path(unsigned int u, map<unsigned int, unsigned int>& data_dist, map<unsigned int, unsigned int>& query_dist,
		map<unsigned int, unsigned int>& data_match, map<unsigned int, unsigned int>& query_match, set<unsigned int>& visited, unsigned int& min_dist)
	{
		for (int i = 0; i < cand_map[u].size(); i++)
		{
			unsigned int v = cand_map[u][i];
			if (visited.find(v) != visited.end()) continue;
			if (data_dist[v] == query_dist[u] + 1) {
				visited.insert(v);
				if (data_match.find(v) != data_match.end() && data_dist[v] == min_dist)
					continue;
				if (data_match.find(v) == data_match.end() || traverse_path(data_match[v], data_dist, query_dist, data_match, query_match, visited, min_dist))
				{
					data_match[v] = u;
					query_match[u] = v;
					return true;
				}
			}
		}
		return false;
	}

	int HK()
	{
	//	cout << "checking match" << endl;
		int ans = 0;
		map<unsigned int, unsigned int> data_dist;
		map<unsigned int, unsigned int> query_dist;
		map<unsigned int, unsigned int> data_match;
		map<unsigned int, unsigned int> query_match;
		set<unsigned int> visited;
		unsigned int min_dist = 0xFFFFFFFF;
		while (bfs(data_dist, query_dist, data_match, query_match, min_dist))
		{
			visited.clear();
			bool find = false;
			for (auto iter = cand_map.begin(); iter != cand_map.end(); iter++)
			{
				unsigned int u = iter->first;
				if (query_match.find(u) == query_match.end())
				{
					if (traverse_path(u, data_dist, query_dist, data_match, query_match, visited, min_dist)) {
						ans++;
						find = true;
					}
				}
			}
			if (!find)
				cout << "cannot extend match " << min_dist<<endl;
		}
		data_dist.clear();
		query_dist.clear();
		data_match.clear();
		query_match.clear();
		visited.clear();
//		cout << "checking match finished" << endl;
		return ans;
	}
	bool brutal_dfs(set<unsigned int> & used, vector<unsigned int>& qvec, unsigned int index )
	{
		unsigned int u = qvec[index];
		for (int i = 0; i < cand_map[u].size(); i++)
		{
			unsigned int v = cand_map[u][i];
			if (used.find(v) != used.end())
				continue;
			if (index + 1 == qvec.size())
				return true;
			used.insert(v);
			bool result = brutal_dfs(used, qvec, index + 1);
			if (result)
				return true;
			else
				used.erase(v);
		}
		return false;
	}
	bool brutal_match_check()
	{
		set<unsigned int> used;
		vector<unsigned int> qvec;
		for (auto iter = cand_map.begin(); iter != cand_map.end(); iter++)
			qvec.push_back(iter->first);
		bool ans= brutal_dfs(used, qvec, 0);
		used.clear();
		qvec.clear();
		return ans;
	}
	bool injective_match()
	{
		if (cand_map.size() < target_u)
			return false;
		//else if (HK() < target_u)
		else if (!brutal_match_check()) //here we tried both Hopcroft¨CKarp algorithm and a brutal search algorithm to determine if there is an injective match. The comparison result shows a draw.  
			//I checked the open-source code a CaliG, it uses brutal search, thus we also adopt it.
			return false;
		else
			return true;

	}

	unsigned int memory_compute()
	{
		unsigned int mem = byte_align(sizeof(cand_map) + sizeof(unsigned int), 8);
		unsigned int mem_per_kv = byte_align((sizeof(vector<unsigned int>) + sizeof(void*)*2 + sizeof(unsigned int)), 8);
		mem += cand_map.size() * mem_per_kv;
		for (auto iter = cand_map.begin(); iter != cand_map.end(); iter++)
			mem += iter->second.capacity() * sizeof(unsigned int);
		return mem;
	}
};

struct cursor
{
	unsigned int vec_pos;
	unsigned int iter_pos;
	unsigned int val;
	void set(unsigned int p1, unsigned int p2, unsigned int v)
	{
		vec_pos = p1;
		iter_pos = p2;
		val = v;
	}
};
bool lf_sort(cursor&p1, cursor& p2)
{
	return p1.val< p2.val;
}

void leap_frog_join(const vector<vector<unsigned int>>& tmp_candidates, vector<unsigned int>& final_candidates)
{
	if (tmp_candidates.size() == 1)
	{
		final_candidates = tmp_candidates[0];
		return;
	}
//	cout << "begin join" << endl;
	vector<cursor> iters;
	unsigned int size = tmp_candidates.size();
	iters.resize(size);
	for (int i = 0; i < size; i++)
		iters[i].set(i, 0, tmp_candidates[i][0]);
	sort(iters.begin(), iters.end(), lf_sort);
	unsigned int cur = 0;
	unsigned int max_val = tmp_candidates[iters[size-1].vec_pos][iters[size-1].iter_pos];
	while (true)
	{
		unsigned int &a = iters[cur].vec_pos;
		unsigned int& b = iters[cur].iter_pos;
		if (tmp_candidates[a][b] == max_val)
		{
			final_candidates.push_back(max_val);
			iters[cur].iter_pos++;
			if (b == tmp_candidates[a].size())
				break;
			else {
				max_val = tmp_candidates[a][b];
				iters[cur].val = max_val;
			}
		}
		else {
			auto iter = lower_bound(tmp_candidates[a].begin()+ b, tmp_candidates[a].end(), max_val);
			if (iter == tmp_candidates[a].end())
				break;
			else {
				max_val = *iter;
				iters[cur].iter_pos = iter - tmp_candidates[a].begin();
				iters[cur].val = max_val;
			}
		}
		cur++;
		if (cur >= iters.size())
			cur = 0;
	}
	iters.clear();
	//cout << "end join" << endl;

}

bool check_order(vector<unsigned int>& vec)
{
	unsigned int min = 0;
	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] < min) {
			cout << "error, unordered vec" << endl;
			return false;
		}
		min = vec[i];
	}
}

