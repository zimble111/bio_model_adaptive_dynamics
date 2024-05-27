#pragma once

#include <Eigen/Dense>

using dl = double;
using vec = Eigen::ArrayXd;
using mat = Eigen::MatrixXd;

//! Подсчёт свёртки
vec conv_dim(const vec &u, const vec &v, ssize_t points, dl A, ssize_t dim);
