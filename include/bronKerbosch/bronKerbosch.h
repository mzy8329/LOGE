#ifndef _BRONKERBOSCH_H_
#define _BRONKERBOSCH_H_

#include <set>
#include <vector>
#include <eigen3/Eigen/Dense>

std::vector<std::set<int>> getClique(const Eigen::MatrixXf& _adj_matrix);
void bronKerbosch(std::set<int> _R, std::set<int> _P, std::set<int> _X, const Eigen::MatrixXi& _adj_matrix, std::vector<std::set<int>>& _clique_list);

#endif // _BRONKERBOSCH_H_