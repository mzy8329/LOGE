#include "bronKerbosch.h"

std::vector<std::set<int>> getClique(const Eigen::MatrixXf& _adj_matrix)
{
    std::vector<std::set<int>> clique_list;
    Eigen::MatrixXi adj_matrix_no_weight = (_adj_matrix.array() != 0).cast<int>();

    std::set<int> R, P, X;
    for (int i = 0; i < adj_matrix_no_weight.cols(); i++)
    {
        P.insert(i);
    }
    bronKerbosch(R, P, X, adj_matrix_no_weight, clique_list);
    return clique_list;
}

void bronKerbosch(std::set<int> _R, std::set<int> _P, std::set<int> _X, const Eigen::MatrixXi& _adj_matrix, std::vector<std::set<int>>& _clique_list)
{
    if (_P.empty() && _X.empty())
    {
        if (_R.size() > 1)
        {
            _clique_list.push_back(_R);

        }
        return;
    }
    static int depth = 0;

    Eigen::VectorXi adj_nums = _adj_matrix.colwise().sum();
    auto P_X = _P;
    P_X.insert(_X.begin(), _X.end());

    std::pair<int, int> u(-1, -1);
    for (auto i : P_X)
    {
        if (adj_nums[i] > u.second)
        {
            u.first = i;
            u.second = adj_nums[i];
        }
    }

    std::set<int> u_r_neighbor;
    for (auto i : _P)
    {
        if (_adj_matrix(u.first, i) == 0)
        {
            u_r_neighbor.insert(i);
        }
    }

    for (int v : u_r_neighbor)
    {
        auto R_new = _R;
        R_new.insert(v);

        std::set<int> P_new;
        std::set<int> X_new;
        for (int w : _P)
        {
            if (_adj_matrix(v, w) != 0) { P_new.insert(w); }
        }
        for (int w : _X)
        {
            if (_adj_matrix(v, w) != 0) { X_new.insert(w); }
        }
        bronKerbosch(R_new, P_new, X_new, _adj_matrix, _clique_list);

        _P.erase(v);
        _X.insert(v);
    }
}

