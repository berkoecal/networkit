/*
 * PrunedLabeling.h
 *
 *  Created on: 11.06.2018
 *      Author: berkoec
 */

#ifndef NETWORKIT_CPP_DISTANCE_PRUNEDLABELING_H_
#define NETWORKIT_CPP_DISTANCE_PRUNEDLABELING_H_

#include <malloc.h>
#include <stdlib.h>
#include <stdint.h>
#include <xmmintrin.h>
#include <sys/time.h>
#include <climits>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stack>
#include <queue>
#include <set>
#include <algorithm>
#include <fstream>
#include <utility>
#include "../graph/Graph.h"

namespace NetworKit {

template<int kNumBitParallelRoots = 50>
class PrunedLabeling {
public:
	// Constructs an index from a graph, given as a list of edges.
	// Vertices should be described by numbers starting from zero.
	// Returns |true| when successful.
	bool ConstructIndex();

	// Returns distance vetween vertices |v| and |w| if they are connected.
	// Otherwise, returns |INT_MAX|.
	inline edgeweight QueryDistance(node v, node w);

	// Loads an index. Returns |true| when successful.
//	bool LoadIndex(std::istream &ifs);
//	bool LoadIndex(const char *filename);

	// Stores the index. Returns |true| when successful.
//	bool StoreIndex(std::ostream &ofs);
//	bool StoreIndex(const char *filename);

	int GetNumVertices() { return num_v_; }
	void Free();
	void PrintStatistics();

	PrunedLabeling(const Graph& G)
	: G(G), num_v_(0), index_(NULL), time_load_(0), time_indexing_(0) {}
	virtual ~PrunedLabeling() {
		Free();
	}

private:
	static const uint8_t INF8;  // For unreachable pairs

	struct index_t {
		uint8_t bpspt_d[kNumBitParallelRoots];
		uint64_t bpspt_s[kNumBitParallelRoots][2];  // [0]: S^{-1}, [1]: S^{0}
		uint32_t *spt_v;
		uint8_t *spt_d;
	} __attribute__((aligned(64)));  // Aligned for cache lines

	const Graph& G;
	count num_v_;
	index_t *index_;

	double GetCurrentTimeSec() {
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return tv.tv_sec + tv.tv_usec * 1e-6;
	}

	// Statistics
	double time_load_, time_indexing_;
};


template<int kNumBitParallelRoots>
const uint8_t PrunedLabeling<kNumBitParallelRoots>::INF8 = 100;



template<int kNumBitParallelRoots>
bool PrunedLabeling<kNumBitParallelRoots>
::ConstructIndex() {
	//
	// Prepare the adjacency list and index space
	//
	Free();

	count E = G.numberOfEdges();
	count &V = num_v_;
	V = G.upperNodeIdBound();
	std::vector<std::vector<node> > adj(V);

	index_ = (index_t*)memalign(64, V * sizeof(index_t));
	if (index_ == NULL) {
		num_v_ = 0;
		return false;
	}
	for (node v = 0; v < V; ++v) {
		index_[v].spt_v = NULL;
		index_[v].spt_d = NULL;
	}
	//
	// Order vertices by decreasing order of degree
	//
	time_indexing_ = -GetCurrentTimeSec();
	std::vector<node> inv(V);  // new label -> old label
	{
		// Order
		std::vector<std::pair<float, node> > deg(V);
		for (node v = 0; v < V; ++v) {
			adj[v] = G.neighbors(v);
			// We add a random value here to diffuse nearby vertices
			deg[v] = std::make_pair(adj[v].size() + float(rand()) / RAND_MAX, v);
		}
		std::sort(deg.rbegin(), deg.rend());
		for (node i = 0; i < V; ++i) inv[i] = deg[i].second;

		// Relabel the vertex IDs
		std::vector<node> rank(V);
		for (node i = 0; i < V; ++i) rank[deg[i].second] = i;
		std::vector<std::vector<node> > new_adj(V);
		for (node v = 0; v < V; ++v) {
			for (size_t i = 0; i < adj[v].size(); ++i) {
				new_adj[rank[v]].push_back(rank[adj[v][i]]);
			}
		}
		adj.swap(new_adj);
	}

	//
	// Bit-parallel labeling
	//
	std::vector<bool> usd(V, false);  // Used as root? (in new label)
	{
		std::vector<uint8_t> tmp_d(V);
		std::vector<std::pair<uint64_t, uint64_t> > tmp_s(V);
		std::vector<node> que(V);
		std::vector<std::pair<node, node> > sibling_es(E);
		std::vector<std::pair<node, node> > child_es(E);

		node r = 0;
		for (int i_bpspt = 0; i_bpspt < kNumBitParallelRoots; ++i_bpspt) {
			while (r < V && usd[r]) ++r;
			if (r == V) {
				for (node v = 0; v < V; ++v) index_[v].bpspt_d[i_bpspt] = INF8;
				continue;
			}
			usd[r] = true;

			fill(tmp_d.begin(), tmp_d.end(), INF8);
			fill(tmp_s.begin(), tmp_s.end(), std::make_pair(0, 0));

			node que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			tmp_d[r] = 0;
			que_t1 = que_h;

			int ns = 0;
			std::vector<node> vs;
			sort(adj[r].begin(), adj[r].end());
			for (size_t i = 0; i < adj[r].size(); ++i) {
				node v = adj[r][i];
				if (!usd[v]) {
					usd[v] = true;
					que[que_h++] = v;
					tmp_d[v] = 1;
					tmp_s[v].first = 1ULL << ns;
					vs.push_back(v);
					if (++ns == 64) break;
				}
			}

			for (int d = 0; que_t0 < que_h; ++d) {
				count num_sibling_es = 0, num_child_es = 0;

				for (node que_i = que_t0; que_i < que_t1; ++que_i) {
					node v = que[que_i];

					for (size_t i = 0; i < adj[v].size(); ++i) {
						node tv = adj[v][i];
						int td = d + 1;

						if (d > tmp_d[tv]);
						else if (d == tmp_d[tv]) {
							if (v < tv) {
								sibling_es[num_sibling_es].first  = v;
								sibling_es[num_sibling_es].second = tv;
								++num_sibling_es;
							}
						} else {
							if (tmp_d[tv] == INF8) {
								que[que_h++] = tv;
								tmp_d[tv] = td;
							}
							child_es[num_child_es].first  = v;
							child_es[num_child_es].second = tv;
							++num_child_es;
						}
					}
				}

				for (node i = 0; i < num_sibling_es; ++i) {
					node v = sibling_es[i].first, w = sibling_es[i].second;
					tmp_s[v].second |= tmp_s[w].first;
					tmp_s[w].second |= tmp_s[v].first;
				}
				for (node i = 0; i < num_child_es; ++i) {
					node v = child_es[i].first, c = child_es[i].second;
					tmp_s[c].first  |= tmp_s[v].first;
					tmp_s[c].second |= tmp_s[v].second;
				}

				que_t0 = que_t1;
				que_t1 = que_h;
			}

			for (node v = 0; v < V; ++v) {
				index_[inv[v]].bpspt_d[i_bpspt] = tmp_d[v];
				index_[inv[v]].bpspt_s[i_bpspt][0] = tmp_s[v].first;
				index_[inv[v]].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
			}
		}
	}

	//
	// Pruned labeling
	//
	{
		// Sentinel (V, INF8) is added to all the vertices
		std::vector<std::pair<std::vector<node>, std::vector<uint8_t> > >
		tmp_idx(V, make_pair(std::vector<node>(1, V),
				std::vector<uint8_t>(1, INF8)));

		std::vector<bool> vis(V);
		std::vector<node> que(V);
		std::vector<uint8_t> dst_r(V + 1, INF8);

		for (node r = 0; r < V; ++r) {
			if (usd[r]) continue;
			index_t &idx_r = index_[inv[r]];
			const std::pair<std::vector<node>, std::vector<uint8_t> >
			&tmp_idx_r = tmp_idx[r];
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = tmp_idx_r.second[i];
			}

			node que_t0 = 0, que_t1 = 0, que_h = 0;
			que[que_h++] = r;
			vis[r] = true;
			que_t1 = que_h;

			for (int d = 0; que_t0 < que_h; ++d) {
				for (node que_i = que_t0; que_i < que_t1; ++que_i) {
					node v = que[que_i];
					std::pair<std::vector<node>, std::vector<uint8_t> >
					&tmp_idx_v = tmp_idx[v];
					index_t &idx_v = index_[inv[v]];

					// Prefetch
					_mm_prefetch(&idx_v.bpspt_d[0], _MM_HINT_T0);
					_mm_prefetch(&idx_v.bpspt_s[0][0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.first[0], _MM_HINT_T0);
					_mm_prefetch(&tmp_idx_v.second[0], _MM_HINT_T0);

					// Prune?
					if (usd[v]) continue;
					for (int i = 0; i < kNumBitParallelRoots; ++i) {
						int td = idx_r.bpspt_d[i] + idx_v.bpspt_d[i];
						if (td - 2 <= d) {
							td +=
									(idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
											((idx_r.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) |
													(idx_r.bpspt_s[i][1] & idx_v.bpspt_s[i][0]))
													? -1 : 0;
							if (td <= d) goto pruned;
						}
					}
					for (size_t i = 0; i < tmp_idx_v.first.size(); ++i) {
						node w = tmp_idx_v.first[i];
						int td = tmp_idx_v.second[i] + dst_r[w];
						if (td <= d) goto pruned;
					}

					// Traverse
					tmp_idx_v.first .back() = r;
					tmp_idx_v.second.back() = d;
					tmp_idx_v.first .push_back(V);
					tmp_idx_v.second.push_back(INF8);
					for (size_t i = 0; i < adj[v].size(); ++i) {
						node w = adj[v][i];
						if (!vis[w]) {
							que[que_h++] = w;
							vis[w] = true;
						}
					}
					pruned:
					{}
				}

				que_t0 = que_t1;
				que_t1 = que_h;
			}

			for (node i = 0; i < que_h; ++i) vis[que[i]] = false;
			for (size_t i = 0; i < tmp_idx_r.first.size(); ++i) {
				dst_r[tmp_idx_r.first[i]] = INF8;
			}
			usd[r] = true;
		}

		for (node v = 0; v < V; ++v) {
			count k = tmp_idx[v].first.size();
			index_[inv[v]].spt_v = (uint32_t*)memalign(64, k * sizeof(uint32_t));
			index_[inv[v]].spt_d = (uint8_t *)memalign(64, k * sizeof(uint8_t ));
			if (!index_[inv[v]].spt_v || !index_[inv[v]].spt_d) {
				Free();
				return false;
			}
			for (count i = 0; i < k; ++i) index_[inv[v]].spt_v[i] = tmp_idx[v].first[i];
			for (count i = 0; i < k; ++i) index_[inv[v]].spt_d[i] = tmp_idx[v].second[i];
			tmp_idx[v].first.clear();
			tmp_idx[v].second.clear();
		}
	}

	time_indexing_ += GetCurrentTimeSec();
	return true;
}


template<int kNumBitParallelRoots>
edgeweight PrunedLabeling<kNumBitParallelRoots>
::QueryDistance(node v, node w) {
  if (v >= num_v_ || w >= num_v_) return v == w ? 0 : INT_MAX;

  const index_t &idx_v = index_[v];
  const index_t &idx_w = index_[w];
  int d = INF8;

  _mm_prefetch(&idx_v.spt_v[0], _MM_HINT_T0);
  _mm_prefetch(&idx_w.spt_v[0], _MM_HINT_T0);
  _mm_prefetch(&idx_v.spt_d[0], _MM_HINT_T0);
  _mm_prefetch(&idx_w.spt_d[0], _MM_HINT_T0);

  for (int i = 0; i < kNumBitParallelRoots; ++i) {
    int td = idx_v.bpspt_d[i] + idx_w.bpspt_d[i];
    if (td - 2 <= d) {
      td +=
          (idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][0]) ? -2 :
          ((idx_v.bpspt_s[i][0] & idx_w.bpspt_s[i][1]) | (idx_v.bpspt_s[i][1] & idx_w.bpspt_s[i][0]))
          ? -1 : 0;

      if (td < d) d = td;
    }
  }
  for (node i1 = 0, i2 = 0; ; ) {
    node v1 = idx_v.spt_v[i1], v2 = idx_w.spt_v[i2];
    if (v1 == v2) {
      if (v1 == num_v_) break;  // Sentinel
      int td = idx_v.spt_d[i1] + idx_w.spt_d[i2];
      if (td < d) d = td;
      ++i1;
      ++i2;
    } else {
      i1 += v1 < v2 ? 1 : 0;
      i2 += v1 > v2 ? 1 : 0;
    }
  }

  if (d >= INF8 - 2) d = INT_MAX;
  return d;
}

template<int kNumBitParallelRoots>
void PrunedLabeling<kNumBitParallelRoots>
::Free() {
  for (node v = 0; v < num_v_; ++v) {
    free(index_[v].spt_v);
    free(index_[v].spt_d);
  }
  free(index_);
  index_ = NULL;
  num_v_ = 0;
}

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_DISTANCE_PRUNEDLABELING_H_ */
