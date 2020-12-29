//
// Created by fenix on 27.12.20.
//

#ifndef INC_2D_SLAM_CUDA_SLAM_INTEGRATE_HEADER_H
#define INC_2D_SLAM_CUDA_SLAM_INTEGRATE_HEADER_H

#include <utility>
#include <vector>

typedef std::pair<double, double> double2;
typedef std::pair<double2, double2> double_double2;

typedef std::pair<double, double> size_t2;
typedef std::pair<size_t2, size_t2> double_size_t2;

typedef std::vector<std::vector<double>> matrix2d_double;
typedef std::vector<std::vector<size_t>> matrix2d_size_t;
typedef std::vector<std::vector<bool>> matrix2d_bool;

namespace np {
    template<typename E /*element*/>
    static std::vector<std::vector<E>> create_2d_with(size_t y, size_t x, E val) {
        std::vector<std::vector<E>> res(y);
        std::vector<E> tmp_vector(x);
        std::fill(tmp_vector.begin(), tmp_vector.end(), val);
        std::fill(res.begin(), res.end(), tmp_vector);
        return res;
    }

    template<typename E /*element*/>
    static std::vector<std::vector<E>> create_1d_with(size_t x, E val) {
        std::vector<E> res(x);
        std::fill(res.begin(), res.end(), val);
        return res;
    }

    template<typename E /*element*/>
    static std::vector<std::vector<E>> ones(size_t y, size_t x) {
        return create_2d_with(y, x, 1u);
    }

    template<typename E /*element*/>
    static std::vector<std::vector<E>> ones(size_t x) {
        return create_1d_with(x, 1u);
    }

    template<typename E /*element*/>
    static std::vector<std::vector<E>> zeros(size_t y, size_t x) {
        return create_2d_with(y, x, 0u);
    }

    template<typename E /*element*/>
    static std::vector<std::vector<E>> zeros(size_t x) {
        return create_1d_with(x, 0u);
    }

    // return empty vector on Error
    template<typename E /*element*/>
    static std::vector<E> linspace(E start, E end, size_t num) {
        std::vector<E> res{};
        if ((end - start) % num != 0)
            return res; // error
        res.reserve(num);
        E dx = (end - start) / num;
        for (E val = start; val < end; val += dx)
            res.emplace_back(val);
        return res;
    }

    template<typename A /*element*/, typename E>
    inline void add_1d(A v, E val) {
        for (auto &el : v)
            el += val;
    }

    template<typename A /*element*/, typename E>
    inline void mul_2d(A v, E val) {
        for (auto &v_1d : v)
            for (auto &el : v_1d)
                el *= val;
    }

    template<typename T>
    std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>>
    meshgrid(const std::vector<T> &v1, const std::vector<T> &v2) {
        std::vector<std::vector<T>> res1, res2;
        for (size_t i = 0; i < v2.size(); ++i)
            res1.emplace_back(v1);

        std::vector<T> temp(v1.size());
        for (size_t j = 0; j < v2.size(); ++j) {
            for (size_t i = 0; i < v1.size(); ++i) {
                temp[i] = v2[j];
            }
            res2.emplace_back(temp);
        }

        return std::pair<std::vector<std::vector<T>>, std::vector<std::vector<T>>>{res1, res2};
    };
}

#endif //INC_2D_SLAM_CUDA_SLAM_INTEGRATE_HEADER_H
