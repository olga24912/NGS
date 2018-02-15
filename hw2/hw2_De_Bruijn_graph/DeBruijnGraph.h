#ifndef INC_04DEBRUIJN_DEBRUIJNGRAPH_H
#define INC_04DEBRUIJN_DEBRUIJNGRAPH_H

#include<bits/stdc++.h>

using namespace std;

/*
 * Class for storage de Bruijn Graph
 */
class DeBruijnGraph {
private:
    int k; //len of vertex string
    vector<string> vertex_val; //from id to string vertex value
    vector<string> edge_val; //full string on edge include vertex value
    vector<double> avg_cover; //edge coverage
    vector<int> to; //from edge id to the end of edge(vertex id)
    vector<int> from; //from edge id to the start of edge(vertex id)

    unordered_map<string, int> vertex_id_from_string; // id from string value of vertex

    vector<vector<int> > g; // graph g[vertex id] = {edge id}
    vector<vector<int> > gr; // reverse graph gr[vertex id] = {edge id}

    int add_vertex(string v1); //add vertex to graph if it is not exists

    void simplify(int v);
    void deleteEdge(int e);

public:
    DeBruijnGraph(int k) {
        this->k = k;
    }

    DeBruijnGraph() {
        this->k = 11;
    }

    //add new edge from v1 to v2 with edge value ed and coverage mass
    void add_edge(string v1, string v2, string ed, double mass);

    //print graph in fasta format to file_name with full edge value
    void print_graph_fasta(string file_name);

    //print graph in dot formt with string value of vertex and len, coverage of edge
    void print_graph_dot(string file_name);

    void delete_tails(int min_len, double min_cover);
    void delete_low_weight_edges(int min_len, double min_cover);
};


#endif //INC_04DEBRUIJN_DEBRUIJNGRAPH_H
