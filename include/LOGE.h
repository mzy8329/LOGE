#ifndef _LOGE_H_
#define _LOGE_H_

#include <Eigen/Dense>
#include <vector>
#include <set>
#include <thread>

#include <bronKerbosch/bronKerbosch.h>

template <typename T>
using match = std::pair<const T*, const T*>;
template <typename T>
using matches = std::vector<match<T>>;
using Bone = struct
{
    std::set<int> union_data;
    std::set<int> bone_data;
    std::set<int> c_index;
    int max_length;
};



std::set<int> intersectSets(const std::set<int>& _a, const std::set<int>& _b) {
    std::set<int> result;
    set_intersection(_a.begin(), _a.end(), _b.begin(), _b.end(),
        inserter(result, result.begin()));
    return result;
}

std::set<int> unionSets(const std::set<int>& _a, const std::set<int>& _b)
{
    std::set<int> result;
    std::set_union(_a.begin(), _a.end(), _b.begin(), _b.end(), std::inserter(result, result.begin()));
    return result;
}


template <typename T>
class LOGE
{
public:
    LOGE(float _match_dist_similarity_th, float _noise_th) :
        match_dist_similarity_th_(_match_dist_similarity_th),
        noise_th_(_noise_th) {
    };
    ~LOGE() {};

    std::vector<matches<T>> getLOGEMatches(const std::vector<T>& _data_list_A, const std::vector<T>& _data_list_B);

protected:
    virtual matches<T> getMatchingPairs(const std::vector<T>& _data_list_A, const std::vector<T>& _data_list_B);
    virtual float getDist(const T& _A, const T& _B);

private:
    Eigen::MatrixXf getAdjMatrix(const matches<T>& _matches);
    void getAdjWeightedMatrix(const std::set<int>& _clique, const matches<T>& _origin_matches, Eigen::MatrixXf& _adj_matrix_with_weight);
    std::vector<std::set<int>> denoiseCliques(const std::vector<std::set<int>>& _C3_list_with_noise, const Eigen::MatrixXf& _adj_matrix_with_weight);
    Bone findBone(const std::vector<std::set<int>>& _C3_list_no_noise);

    float match_dist_similarity_th_;
    float noise_th_;
    float bone_th_;
};

template <typename T>
std::vector<matches<T>> LOGE<T>::getLOGEMatches(const std::vector<T>& _data_list_A, const std::vector<T>& _data_list_B)
{
    matches<T> origin_matches = getMatchingPairs(_data_list_A, _data_list_B);
    int u_num = origin_matches.size();
    if (u_num < 2) { return {}; }

    Eigen::MatrixXf adj_matrix = getAdjMatrix(origin_matches);
    std::vector<std::set<int>> C1_list = getClique(adj_matrix);

    Eigen::MatrixXf adj_matrix_with_weight = Eigen::MatrixXf::Zero(u_num, u_num);
    std::vector<std::thread> thread_list;
    for (const auto& clique_set : C1_list)
    {
        thread_list.emplace_back(std::thread(this->getAdjWeightedMatrix, clique_set, origin_matches, std::ref(adj_matrix_with_weight)));
    }
    for (const auto& thread : thread_list)
    {
        thread.join();
    }

    std::vector<std::set<int>> C3_list_with_noise = getClique(adj_matrix_with_weight);
    std::vector<std::set<int>> C3_list_no_noise = denoiseCliques(C3_list_with_noise);

    std::vector<std::set<int>> final_C_list;
    std::vector<matches<T>> all_matches;
    while (!C3_list_no_noise.empty())
    {
        Bone bone = findBone(C3_list_no_noise);
        int index = 0;
        for (int id : bone.c_index)
        {
            C3_list_no_noise.erase(C3_list_with_noise.begin() + id - index);
            index++;
        }
        if (bone.union_data.size() >= 2)
        {
            matches<T> m;
            for (const int& i : bone.union_data)
            {
                m.emplace_back(origin_matches[i]);
            }
            all_matches.emplace_back(m);
        }
    }
    return all_matches;
}

template <typename T>
Eigen::MatrixXf LOGE<T>::getAdjMatrix(const matches<T>& _matches)
{
    Eigen::MatrixXf adjMatrix = Eigen::MatrixXf::Zero(_matches.size(), _matches.size());
    for (int i = 0; i < _matches.size() - 1; i++)
    {
        for (int j = i + 1; j < _matches.size(); j++)
        {
            if (_matches[i].first != _matches[j].first
                && _matches[i].second != _matches[j].second)
            {
                adjMatrix(i, j) = 1;
                adjMatrix(j, i) = 1;
            }
        }
    }
    return adjMatrix;
}

template <typename T>
void LOGE<T>::getAdjWeightedMatrix(const std::set<int>& _clique, const matches<T>& _origin_matches, Eigen::MatrixXf& _adj_matrix_with_weight)
{
    int c_size = clique.size();
    int match_num = 0;

    Eigen::MatrixXf sub_adj_matrix = Eigen::MatrixXf::Zero(c_size, c_size);
    std::vector<int> clique(_clique.begin(), _clique.end());
    for (int i = 0; i < c_size - 1; i++)
    {
        for (int j = i + 1; j < c_size; j++)
        {
            float left_d = getDist(_origin_matches[clique[i]].first, _origin_matches[clique[j]].first);
            float right_d = getDist(_origin_matches[clique[i]].first, _origin_matches[clique[j]].first);
            bool similarity = std::abs(left_d - right_d) / (left_d + right_d) < match_dist_similarity_th_;
            if (similarity)
            {
                sub_adj_matrix(i, j) = 1;
                sub_adj_matrix(j, i) = 1;
                match_num++;
            }
        }
    }
    if (match_num < 1) { return; }

    std::vector<std::set<int>> C2_list = getClique(sub_adj_matrix);
    for (const auto& sub_clique_set : C2_list)
    {
        std::vector<int> sub_clique(sub_clique_set.begin(), sub_clique_set.end());
        for (int i = 0; i < sub_clique.size() - 1; i++)
        {
            for (int j = i + 1; j < sub_clique.size(); j++)
            {
                _adj_matrix_with_weight(clique[sub_clique[i]], clique[sub_clique[j]]) += 1;
                _adj_matrix_with_weight(clique[sub_clique[j]], clique[sub_clique[i]]) += 1;
            }
        }
    }
}

template <typename T>
std::vector<std::set<int>> LOGE<T>::denoiseCliques(const std::vector<std::set<int>>& _C3_list_with_noise, const Eigen::MatrixXf& _adj_matrix_with_weight)
{
    std::vector<std::set<int>> clique_denoised;
    for (const auto& clique : _C3_list_with_noise)
    {
        if (sub_clique.size() < 2) { continue; }

        std::vector<int> num_list;
        int den = 0;
        for (int i : clique)
        {
            int num = 0;
            for (int j : clique)
            {
                num += _adj_matrix_with_weight(i, j);
            }
            den += num;
            num_list.emplace_back(num);
        }

        int n = sub_clique.size();
        float th = noise_th_ / float((n / 2 - 1) + noise_th_);
        std::set<int> clique_af_p;
        if (n > 2)
        {
            for (int i : clique)
            {
                float p = num_list[i] / float(den);
                if (p >= th)
                {
                    clique_af_p.insert(i);
                }
            }
            clique_denoised.emplace_back(clique_af_p);
        }
        else if (n == 2)
        {
            clique_denoised.emplace_back(clique);
        }
        else
        {
            continue;
        }
    }
    return clique_denoised;
}

template <typename T>
Bone LOGE<T>::findBone(const std::vector<std::set<int>>& _C3_list_no_noise)
{
    Bone bone;
    for (int i = 0; i < _C3_list_no_noise.size(); i++)
    {
        if (bone.data.empty())
        {
            bone.union_data = _C3_list_no_noise[i];
            bone.bone_data = _C3_list_no_noise[i];
            bone.c_index.insert(i);
            bone.max_length = bone.union_data.size();
        }
        else
        {
            std::set<int> inter_sec = intersectSets(bone.bone_data, _C3_list_no_noise[i]);
            if (inter_sec.size() > _C3_list_no_noise[i].size() * bone_th_
                && inter_sec.size() > bone.max_length * bone_th_)
            {
                bone.union_data = unionSets(bone.union_data, _C3_list_no_noise[i]);
                bone.bone_data = inter_sec;
                bone.c_index.insert(i);
                bone.max_length = std::max(bone.max_length, _C3_list_no_noise[i].size());
            }
        }
    }
    return bone;
}


#endif //_LOGE_H_