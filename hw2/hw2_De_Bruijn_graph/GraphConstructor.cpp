#include "GraphConstructor.h"

DeBruijnGraph GraphConstructor::read_and_add_all_reads() {
    find_zip_k_mer();
    string read;
    while(getline(fin, read)) {
        if (correct_read(read)) {
            process_one_read(read);
            process_one_read(create_complementary_read(read));
        }
    }
    return graph;
}

bool GraphConstructor::correct_read(string read) {
    for (int i = 0; i < (int) read.size(); ++i) {
        if (read[i] != 'A' && read[i] != 'G' && read[i] != 'T' && read[i] != 'C') {
            return 0;
        }
    }
    return 1;
}

void GraphConstructor::process_one_read(string read) {
    string fr = find_bg(read.substr(0, k));
    for (int i = 0; i < (int)read.size() - k; ++i) {
        int j = i + 1;
        while (j < (int)read.size() - k && k_mer_not_zip.count(read.substr(j, k)) == 0) {
            ++j;
        }

        if (j == read.size() - k) {
            string ed = find_ed(read.substr((int)read.size() - k, k));
            if (i != 0) {
                graph.add_edge(read.substr(i, k), ed.substr(ed.size() - k, k),
                               read.substr(i, j - i).append(ed), (read.size() - i - k)/(1.0 * (ed.size() + j - i - k)));
            } else {
                fr = fr.substr(0, fr.size() - k).append(read.substr(0, read.size() - k)).append(ed);
                graph.add_edge(fr.substr(0, k), fr.substr(fr.size() - k, k),
                               fr, (read.size() - k)/(1.0*(fr.size() - k)));
            }
            i = j - 1;
            continue;
        }

        if (i != 0) {
            graph.add_edge(read.substr(i, k), read.substr(j, k), read.substr(i, j - i + k), 1);
        } else {
            graph.add_edge(fr.substr(0, k), read.substr(j, k),
                           fr.append(read.substr(k, j)), j/(1.0 * (fr.size() + (j - k))));
        }
        i = j - 1;

    }
}

string GraphConstructor::create_complementary_read(string read) {
    reverse(read.begin(), read.end());
    for (int i = 0; i < (int)read.size(); ++i) {
        if (read[i] == 'A') {
            read[i] = 'T';
        } else if (read[i] == 'T') {
            read[i] = 'A';
        } else if (read[i] == 'G') {
            read[i] = 'C';
        } else {
            read[i] = 'G';
        }
    }
    return read;
}

void GraphConstructor::find_zip_k_mer() {
    string read;
    while(getline(fin, read)) {
        if (correct_read(read)) {
            analyze_one_read(read);
            analyze_one_read(create_complementary_read(read));
        }
    }

    fin.close();
    fin.open(file_name);
}

void GraphConstructor::analyze_one_read(string read) {
    if (read.size() < k) return;
    for (int i = 0; i <= (int)read.size() - k; ++i) {
        string cur = read.substr(i, k);
        if (i - 1 >= 0 && (k_mer_left.count(cur) == 0 || k_mer_left[cur] == read[i - 1])) {
            k_mer_left[cur] = read[i - 1];
        } else if (i - 1 >= 0) {
            k_mer_not_zip.insert(cur);
        }


        if (i + k < read.size() && (k_mer_right.count(cur) == 0 || k_mer_right[cur] == read[i + k])) {
            k_mer_right[cur] = read[i + k];
        } else if (i + k < read.size()) {
            k_mer_not_zip.insert(cur);
        }
    }
}

string GraphConstructor::find_bg(string k_mer) {
    while(k_mer_not_zip.count(k_mer.substr(0, k)) == 0) {
        if (k_mer_left.count(k_mer.substr(0, k)) == 0) {
            return k_mer;
        } else {
            k_mer = k_mer_left[k_mer.substr(0, k)] + k_mer;
        }
    }
    return k_mer;
}

string GraphConstructor::find_ed(string k_mer) {
    while(k_mer_not_zip.count(k_mer.substr(k_mer.size() - k, k)) == 0) {
        if (k_mer_right.count(k_mer.substr(k_mer.size() - k, k)) == 0) {
            return k_mer;
        } else {
            k_mer = k_mer + k_mer_right[k_mer.substr(k_mer.size() - k, k)];
        }
    }
    return k_mer;
}
