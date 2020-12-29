//
// Created by fenix on 27.12.20.
//

#ifndef INC_2D_SLAM_CUDA_SLAM_INTEGRATE_HEADER_H
#define INC_2D_SLAM_CUDA_SLAM_INTEGRATE_HEADER_H

#include <utility>
#include <vector>
#include <algorithm>
#include <cmath>

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
        return create_2d_with(y, x, (size_t) 1u);
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
        return create_1d_with(x, (size_t) 0u);
    }

    // return empty vector on Error
    static std::vector<size_t> linspace(size_t start, size_t end, size_t num) {
        std::vector<size_t> res{};
        if ((end - start) % num != 0)
            return res; // error
        res.reserve(num);
        size_t dx = (end - start) / num;
        for (size_t val = start; val < end; val += dx)
            res.emplace_back(val);
        return res;
    }

    // return empty vector on Error
    static std::vector<double> linspace(double start, double end, size_t num) {
        std::vector<double> res{};
        res.reserve(num);
        const double dx = (end - start) / num;
        double val = start;
        while (val < end) {
            res.emplace_back(val);
            val += dx;
        }
        return res;
    }

    template<typename A /*element*/, typename E>
    inline void add_1d(A v, E val) {
        for (auto &el : v)
            el += val;
    }

    template<typename E>
    std::vector<std::vector<E>> add_2d(std::vector<std::vector<E>> a, const std::vector<std::vector<E>> &b) {
        for (int j = 0; j < a.size(); ++j) {
            for (int i = 0; i < a[0].size(); ++i) {
                a[j][i] += b[j][i];
            }
        }
        return a;
    }

    template<typename E>
    std::vector<std::vector<E>> add_2d(std::vector<std::vector<E>> &&a, std::vector<std::vector<E>> &&b) {
        for (int j = 0; j < a.size(); ++j) {
            for (int i = 0; i < a[0].size(); ++i) {
                a[j][i] += b[j][i];
            }
        }
        return std::move(a);
    }

    template<typename E>
    std::vector<std::vector<E>> pow_2d(std::vector<std::vector<E>> a) {
        for (int i = 0; i < a.size(); ++i) {
            a[i] *= a[i];
        }
        return a;
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

    template<typename T>
    std::vector<std::vector<T>> &fliplr(std::vector<std::vector<T>> &v) {
        for (auto &el : v)
            std::reverse(el.begin(), el.end());
        return v;
    }

    template<typename T>
    std::vector<std::vector<T>> &flipud(std::vector<std::vector<T>> &v) {
        std::reverse(v.begin(), v.end());
        return v;
    }

    template<typename T>
    std::vector<std::vector<double>> sqrt_2d(std::vector<std::vector<T>> &&v) {
        std::vector<std::vector<double>> v_res{std::move(v)};
        for (auto &el : v)
            std::for_each(el.begin(), el.end(), [](T &el) { el = sqrt(el); });
        return std::move(v);
    }
}

#endif //INC_2D_SLAM_CUDA_SLAM_INTEGRATE_HEADER_H
