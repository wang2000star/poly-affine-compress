#include "circuit.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
}

Circuit read_circuit(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "  ERROR: cannot open " << path << "\n"; std::exit(1); }
    std::stringstream buf;
    buf << f.rdbuf();
    std::string text = buf.str();
    for (auto& c : text) if (c == '\n' || c == '\r') c = ' ';
    while (text.find("  ") != std::string::npos)
        text.replace(text.find("  "), 2, " ");

    Circuit circ;
    auto p = text.find("INORDER");
    if (p == std::string::npos) { std::cerr << "  ERROR: no INORDER\n"; std::exit(1); }
    auto eq = text.find('=', p);
    auto semi = text.find(';', p);
    if (eq == std::string::npos || semi == std::string::npos) { std::cerr << "  ERROR: bad INORDER\n"; std::exit(1); }
    {
        std::string list = trim(text.substr(eq + 1, semi - eq - 1));
        std::stringstream ss(list); std::string name;
        while (ss >> name) circ.inputs.push_back(name);
    }
    circ.n_inputs = (int)circ.inputs.size();

    p = text.find("OUTORDER");
    if (p == std::string::npos) { std::cerr << "  ERROR: no OUTORDER\n"; std::exit(1); }
    eq = text.find('=', p);
    semi = text.find(';', p);
    if (eq == std::string::npos || semi == std::string::npos) { std::cerr << "  ERROR: bad OUTORDER\n"; std::exit(1); }
    {
        std::string list = trim(text.substr(eq + 1, semi - eq - 1));
        std::stringstream ss(list); std::string name;
        while (ss >> name) circ.outputs.push_back(name);
    }

    for (auto& inp : circ.inputs)
        circ.name_to_idx[inp] = -(int)circ.name_to_idx.size() - 1;

    auto body_start = text.find(';', std::max(
        text.find("INORDER"), text.find("OUTORDER"))) + 1;
    std::string body = text.substr(body_start);
    std::stringstream ss(body);
    std::string token;
    while (std::getline(ss, token, ';')) {
        token = trim(token);
        if (token.empty()) continue;
        auto m = token.find('=');
        if (m == std::string::npos) continue;
        std::string lhs = trim(token.substr(0, m));
        std::string rhs = trim(token.substr(m + 1));
        if (rhs.empty()) continue;

        Stmt s; s.name = lhs;
        if (rhs == "0" || rhs == "false") { s.op = Op::CONST0; }
        else if (rhs == "1" || rhs == "true") { s.op = Op::CONST1; }
        else if (rhs[0] == '!') { s.op = Op::NOT; s.arg1 = trim(rhs.substr(1)); }
        else if (rhs.find('*') != std::string::npos) {
            s.op = Op::AND;
            auto star = rhs.find('*');
            s.arg1 = trim(rhs.substr(0, star));
            s.arg2 = trim(rhs.substr(star + 1));
        } else if (rhs.find('+') != std::string::npos) {
            s.op = Op::XOR;
            auto plus = rhs.find('+');
            s.arg1 = trim(rhs.substr(0, plus));
            s.arg2 = trim(rhs.substr(plus + 1));
        } else {
            s.op = Op::INPUT;
            s.arg1 = rhs;
        }
        circ.stmts.push_back(s);
        circ.name_to_idx[lhs] = (int)circ.stmts.size();
    }
    return circ;
}

int count_gates(const Circuit& circ, Op op) {
    return (int)std::count_if(circ.stmts.begin(), circ.stmts.end(),
        [op](const Stmt& s) { return s.op == op; });
}

int count_and_gates(const Circuit& circ) { return count_gates(circ, Op::AND); }
int count_xor_gates(const Circuit& circ) { return count_gates(circ, Op::XOR); }
int count_not_gates(const Circuit& circ) { return count_gates(circ, Op::NOT); }
