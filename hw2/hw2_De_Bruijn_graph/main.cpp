#include <iostream>
#include "GraphConstructor.h"

using namespace std;

const int def_k = 55;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Need filename for reads");
        return 1;
    }
    int k = def_k;
    if (argc == 3) {
        k = atoi(argv[2]);
        //if (k % 2 == 0) {
        //    printf("k need to be odd");
        //    return 1;
        //}
    }

    GraphConstructor graphConstructor(string(argv[1]), k);
    DeBruijnGraph graph = graphConstructor.read_and_add_all_reads();
    graph.print_graph_fasta("out.fasta");
    graph.print_graph_dot("out.dot");

    graph.delete_tails(k + 10, 10);
    graph.print_graph_dot("delete_tail.dot");

    graph.delete_low_weight_edges(k + 10, 10);
    graph.print_graph_dot("delete_low_weight.dot");
    return 0;
}