#ifndef NJS_HPP
#define NJS_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

int give_index(int i, int j, int n);
void error(const std::string &message, int OTU1, int OTU2);
int H(double t);
double cnxy(int x, int y, int n, std::vector<double>& D);
int mxy(int x, int y, int n, std::vector<double>& D);
double nxy(int x, int y, int n, std::vector<double>& D);
int cxy(int x, int y, int n, std::vector<double>& D);
void cpp_choosePair( std::vector<double>& D, int n, std::vector<double>& R, std::vector<int>& s, int& sw, int& x, int& y, int fS);
void cpp_njs(std::vector<double>& D, int N, std::vector<int>& edge1, std::vector<int>& edge2, std::vector<double>& edge_length, int fsS);

#endif
