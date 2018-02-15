//
// Created by olga on 08.10.16.
//

#ifndef INC_04DEBRUIJN_HASH_H
#define INC_04DEBRUIJN_HASH_H

#include "bits/stdc++.h"

using namespace std;

class Hash {
private:
    const int p = 179;
    string s;
    int pos = 0;
    int h;
    vector<int> pp;
public:
    Hash(int n);
    int push_char(char c);
    int set_string(string s);
    int pop_char();
    int getHash();
};


#endif //INC_04DEBRUIJN_HASH_H
