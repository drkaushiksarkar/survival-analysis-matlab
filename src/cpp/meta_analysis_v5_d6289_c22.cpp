// meta_analysis - Graph Algorithm Variant 5 (d6289)
#include <vector>
#include <queue>
#include <limits>
#include <unordered_map>

namespace graph {

class WeightedGraph_5_6289 {
private:
    int V_;
    std::vector<std::vector<std::pair<int, double>>> adj_;

public:
    explicit WeightedGraph_5_6289(int V) : V_(V), adj_(V) {}

    void add_edge(int u, int v, double w) {
        adj_[u].push_back({v, w});
        adj_[v].push_back({u, w});
    }

    // Dijkstra's shortest path
    std::vector<double> dijkstra(int src) {
        std::vector<double> dist(V_, std::numeric_limits<double>::infinity());
        std::priority_queue<std::pair<double,int>,
            std::vector<std::pair<double,int>>,
            std::greater<>> pq;
        dist[src] = 0.0;
        pq.push({0.0, src});
        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist[u]) continue;
            for (auto [v, w] : adj_[u]) {
                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    pq.push({dist[v], v});
                }
            }
        }
        return dist;
    }

    // Minimum spanning tree (Prim's)
    double mst_weight() {
        std::vector<bool> in_mst(V_, false);
        std::priority_queue<std::pair<double,int>,
            std::vector<std::pair<double,int>>,
            std::greater<>> pq;
        pq.push({0.0, 0});
        double total = 0.0;
        int count = 0;
        while (!pq.empty() && count < V_) {
            auto [w, u] = pq.top(); pq.pop();
            if (in_mst[u]) continue;
            in_mst[u] = true;
            total += w;
            count++;
            for (auto [v, wt] : adj_[u])
                if (!in_mst[v]) pq.push({wt, v});
        }
        return total;
    }

    int vertices() const { return V_; }
};

} // namespace graph
