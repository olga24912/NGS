#include "DeBruijnGraph.h"

void DeBruijnGraph::add_edge(string v1, string v2, string val, double mass) {
    assert(val.substr(val.size() - k, k) == v2);
    int id1 = add_vertex(v1), id2 = add_vertex(v2);
    int eid = (int) edge_val.size();

    int fl = 0;
    for (int i = 0; i < (int)g[id1].size() && !fl; ++i) {
        if (edge_val[g[id1][i]] == val) {
            fl = 1;
            avg_cover[g[id1][i]] += mass;
        }
    }

    if (fl == 0) {
        edge_val.push_back(val);
        to.push_back(id2);
        from.push_back(id1);
        avg_cover.push_back(1);
        g[id1].push_back(eid);
        gr[id2].push_back(eid);
    }
}

int DeBruijnGraph::add_vertex(string v1) {
    int id;
    if (vertex_id_from_string.count(v1) == 0) {
        id = (int) vertex_val.size();
        vertex_val.push_back(v1);
        g.push_back(vector<int>());
        gr.push_back(vector<int>());
        vertex_id_from_string[v1] = id;
    } else {
        id = vertex_id_from_string[v1];
    }
    return id;
}

void DeBruijnGraph::print_graph_fasta(string file_name) {
    ofstream fout(file_name);

    fout << "> de Bruijn graph with k = " << k << "\n";
    for (int i = 0; i < (int)g.size(); ++i) {
        for (int ed : g[i]) {
            fout << "> vertex # " << i << " -> vertex # " << to[ed] << "\n";
            fout << edge_val[ed] << "\n";
        }
    }

    fout.close();
}

void DeBruijnGraph::print_graph_dot(string file_name) {
    ofstream fout(file_name);

    fout << "digraph deBruijn {\n";
    for (int i = 0; i < (int)g.size(); ++i) {
        for (int ed : g[i]) {
            fout << vertex_val[i] << " -> " << vertex_val[to[ed]] << " [label = \"len = " << edge_val[ed].size()
            << ", cov = " << avg_cover[ed] << "\"];\n";
        }
    }

    fout << "}\n";

    fout.close();
}

void DeBruijnGraph::delete_tails(int min_len, double min_cover) {
    for (int v = 0; v < g.size(); ++v) {
        for (int i = 0; i < g[v].size(); ++i) {
            int e = g[v][i];
            if (g[to[e]].size() == 0 && avg_cover[e] < min_cover && edge_val[e].size() < min_len) {
                deleteEdge(e);
                --i;
            }
        }

        for (int i = 0; i < gr[v].size(); ++i) {
            int e = gr[v][i];
            if (gr[from[e]].size() == 0 && avg_cover[e] < min_cover && edge_val[e].size() < min_len) {
                deleteEdge(e);
                --i;
            }
        }

        simplify(v);
    }
}

void DeBruijnGraph::delete_low_weight_edges(int min_len, double min_cover) {
    for (int v = 0; v < g.size(); ++v) {
        for (int i = 0; i < g[v].size(); ++i) {
            int e = g[v][i];
            if (avg_cover[e] < min_cover && edge_val[e].size() < min_len) {
                deleteEdge(e);
                --i;
            }
        }

        simplify(v);
    }
}

void DeBruijnGraph::deleteEdge(int e) {
    int v = from[e];
    int u = to[e];

    for (int i = 0; i < g[v].size(); ++i) {
        if (g[v][i] == e) {
            swap(g[v][i], g[v][g[v].size() - 1]);
            g[v].resize(g[v].size() - 1);
            --i;
        }
    }


    for (int i = 0; i < gr[u].size(); ++i) {
        if (gr[u][i] == e) {
            swap(gr[u][i], gr[u][gr[u].size() - 1]);
            gr[u].resize(gr[u].size() - 1);
            --i;
        }
    }
}

void DeBruijnGraph::simplify(int v) {
    if (g[v].size() == 1 && gr[v].size() == 1) {
        int fv = from[gr[v][0]];
        int sv = to[g[v][0]];

        int fe = gr[v][0], se = g[v][0];
        to[fe] = sv;
        avg_cover[fe] = avg_cover[fe] * (edge_val[fe].size() - k) + avg_cover[se] * (edge_val[se].size() - k);
        edge_val[fe].append(edge_val[se].substr(k, edge_val[se].size() - k));
        avg_cover[fe] /= (edge_val[fe].size() - k);

        for (int j = 0; j < gr[sv].size(); ++j) {
            if (gr[sv][j] == se) {
                gr[sv][j] = fe;
            }
        }
        g[v].resize(0);
        gr[v].resize(0);
    }
}


