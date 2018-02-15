#ifndef INC_04DEBRUIJN_GRAPHCONSTRUCTOR_H
#define INC_04DEBRUIJN_GRAPHCONSTRUCTOR_H

#include "DeBruijnGraph.h"
#include "Hash.h"
/*
 * class for construct de Bruijn graph from reads in fasta/fastq files
 */
class GraphConstructor {
private:
    const int maxEdge = 5000;

    ifstream fin; //stream from fasta/fastq file with reads
    string file_name; //name of file with reads
    DeBruijnGraph graph; //compress de Bruijn graph from reads


    unordered_map<string, char> k_mer_left, k_mer_right; //char that before k_mer and after k_mer
    unordered_set<string> k_mer_not_zip; //value of vertex, that will not be compress
    int k; //len of vertex value

    void process_one_read(string read); //split read on k_mers and add it to graph

    bool correct_read(string read); //check string contain only AGCT

    string create_complementary_read(string read); //chage A<->T, G<->C and reverse

    void find_zip_k_mer(); //find k_mer that will be as vertex

    void analyze_one_read(string read); //find k_mer that will be vertex in one read and save k_mer_left/right

    string find_bg(string k_mer); //if it is first k_mer in some read, find the begin of the edge for this k_mer

    string find_ed(string k_mer); //if it is last k_mer, find the end of the edge.
public:
    GraphConstructor(string file_name, int k) {
        this->file_name = file_name;
        this->k = k;
        fin.open(file_name);
        graph = DeBruijnGraph(k);
    }

    DeBruijnGraph read_and_add_all_reads();//construct compress de Bruijn graph from reads

    ~GraphConstructor() {
        fin.close();
    }
};


#endif //INC_04DEBRUIJN_GRAPHCONSTRUCTOR_H
