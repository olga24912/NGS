//
// Created by olga on 08.10.16.
//

#include "Hash.h"

Hash::Hash(int n) {
    pp.resize(n);
    pp[0] = 1;
    for (int i = 1; i < n; ++i) {
        pp[i] = pp[i - 1] * p;
    }
}

int Hash::push_char(char c) {
    s += c;
    h = h * p + c;
    if (s.size() >= pp.size()) {
        for (int i = 0, j = pos; j < s.size(); ++j, ++i) {
            s[i] = s[j];
        }
        s.resize(s.size() - pos);
    }
    return h;
}

int Hash::set_string(string s) {
    this->s = s;
    h = 0;
    for (int i = 0; i < (int)s.size(); ++i) {
        h = h * p + s[i];
    }
    return h;
}

int Hash::pop_char() {
    h -= s[pos] * pp[(s.size() - pos) - 1];
    ++pos;
    return h;
}

int Hash::getHash() {
    return h;
}
