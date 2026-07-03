#ifndef CIRCUIT_H
#define CIRCUIT_H

#include "types.h"
#include <string>

Circuit read_circuit(const std::string& path);
int count_gates(const Circuit& circ, Op op);
int count_and_gates(const Circuit& circ);
int count_xor_gates(const Circuit& circ);
int count_not_gates(const Circuit& circ);

#endif
