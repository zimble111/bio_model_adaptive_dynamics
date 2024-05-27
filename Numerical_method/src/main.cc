#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <mutex>
#include <iostream>
#include <memory>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <filesystem>
#include <iostream>
#include <filesystem>
#include <initializer_list>
#include <thread>
#include <string>
#include "clipp.h"
#include "conv.h"

constexpr dl PI = 3.1415926535897932384626433832795;

using dl = double;
using vec = Eigen::Array<dl, Eigen::Dynamic, 1>;
using mat = Eigen::Matrix<dl, Eigen::Dynamic, Eigen::Dynamic>;

using Eigen::pow;
using Eigen::exp;
using std::abs;

dl y_calc(vec w, vec C, dl d, dl h, ssize_t N, dl A, ssize_t dim);

//! Функция возвращающая значение функции Гаусса
/*!
  \f$norm(x, \mu, \sigma)=\frac{1}{\displaystyle{\sigma\sqrt{2\pi}}}e^{\displaystyle{-\frac{(x-\mu)^2}{2\sigma^2}}}\f$
  \param x это \f$x\f$ может быть скаляром или вектором Eigen
  \param mu это \f$\mu\f$, медиана, может быть скаляром или вектором Eigen
  \param sigma это \f$\sigma\f$, среднеквадратичное отклонение, может быть скаляром или вектором Eigen
  \return вектор или скаляр значения функции Гауса
*/
template<class T, class M, class S>
T normpdf(const T& x, const M& mu, const S& sigma) {
    return pow(2. * PI * pow(sigma, 2), -0.5) * exp(-0.5 * pow((x - mu) / sigma, 2));
}

//! Подсчёт функции Гауса для матрицы векторов с домножением на число
/*!
  \param x это вектор, на котором рассчитывается функция Гауса
  \param mu это матрица медиан: \f$M=\pmatrix{\mu_{1,1} & \cdots & \mu_{1,n} \\ & \ddots & \\ \mu_{m,1} & \cdots & \mu_{m,n}}\f$
  \param sigma это матрица среднеквадратичных отклонений: \f$S=\pmatrix{\sigma_{1,1} & \cdots & \sigma_{1,n} \\ & \ddots & \\ \sigma_{m,1} & \cdots & \sigma_{m,n}}\f$
  \param mul это матрица чисел, на которые домножактся результат: \f$V=\pmatrix{\sigma_{1,1} & \cdots & \sigma_{1,n} \\ & \ddots & \\ \sigma_{m,1} & \cdots & \sigma_{m,n}}\f$
  \return матрицу векторов: \f$\pmatrix{V_{1,1}\cdot norm(\overrightarrow{x}, M_{1,1}, S_{1,1}) & \cdots & V_{1,n}\cdot norm(\overrightarrow{x}, M_{1,n}, S_{1,n}) \\  & \ddots & \\ V_{m,1}\cdot norm(\overrightarrow{x}, M_{m,1}, S_{m,1}) & \cdots & V_{m,n}\cdot norm(\overrightarrow{x}, M_{m,n}, S_{m,n})}\f$
*/
Eigen::Array<vec, Eigen::Dynamic, Eigen::Dynamic> normpdf(
        const vec &x,
        const mat &mu,
        const mat &sigma,
        const mat &mul
) {
    Eigen::Array<vec, Eigen::Dynamic, Eigen::Dynamic> res;
    res.resizeLike(mu);
    for (ssize_t i = 0; i < res.rows(); ++i) {
        for (ssize_t j = 0; j < res.cols(); ++j) {
            res(i, j) = mul(i, j) * normpdf(x, mu(i, j), sigma(i, j));
        }
    }
    return res;
}

//! Подсчёт ошибки
/*! Ошибка равна: 
 \f$
    \sum_{i,j}\left(
        \begin{array}{l|r}
            \left|\displaystyle{\frac{y_{i,j} - \text{old}_{i,j}}{y_{i,j}}}\right| & y_{i,j} \neq 0 \\
                1 & y_{i,j} = 0
        \end{array}
    \right)
\f$
 \param y новое значение
 \param old старое значение
 */
dl mis_calc(const mat &y, const mat &old) {
    dl res = 0;
    for (ssize_t i = 0; i < y.rows(); ++i) {
        for (ssize_t j = 0; j < y.cols(); ++j) {
            if (old(i, j) == 0.) {
                if (y(i, j) != 0) {
                    res += 1;
                }
                continue;
            }
            res += abs((y(i, j) - old(i, j)) / y(i, j));
        }
    }
    return res;
}
struct Base {
    ssize_t Species;
    vec N;
    Eigen::Array<vec, Eigen::Dynamic, Eigen::Dynamic> D;
    Eigen::Array<vec, Eigen::Dynamic, Eigen::Dynamic> w;
    mat dmat;
    Eigen::Array<vec, Eigen::Dynamic, 1> m;
    mat y;
    vec b;
    vec dvec;
    vec sm;
    mat sw;
    dl h;
    dl A;
    dl al;
    ssize_t points;
    ssize_t dim;
    size_t iters;

    /*! Сделана, что бы не вызавать каждый раз с параметрами, которые не меняются во время итераций
      \param u вектор \f$u\f$
      \param v вектор \f$v\f$
      \return вектор \f$(u * v)\f$
    */
    vec Conv(const vec& u, const vec& v) const {
        return conv_dim(u, v, points, A, dim);
    }

    //! Распечатывает вектор \ref N в stderr 
    void PrintN() {
        for (ssize_t i = 0; i < Species; ++i) {
            std::cerr << (i == 0 ? "" : "\t") << N(i);
        }
        std::cerr << std::endl;
    }

    //! Подсчёт вектора N
    /*! Для подсчёта решается матричное уравнение \f$X\cdot Y=B-D\f$
        \param Y считается с помощью \ref y_calc()
        \param B \ref Base::b
        \param D \ref Base::dvec
        \return \f$X\f$ вектор значений решения данного уравнения
     */
    vec n_calc() {
        Eigen::VectorXd bd = b - dvec;
        mat y;
        y.resize(Species, Species);
        for (ssize_t i = 0; i < Species; ++i) {
            for (ssize_t j = 0; j < Species; ++j) {
                y(i, j) = y_calc(   w(i, j),
                                    D(i, j),
                                    dmat(i, j),
                                    h,
                                    points,
                                    A,
                                    dim
                );
            }
        }
        vec x = y.colPivHouseholderQr().solve(bd);
        return x;
    }
    //! Подсчёт \f$D(i, j)\f$
    /*! Считает
        \f$D(i,j) = 
            \left\{
                \begin{array}{l|r}
                    \displaystyle{\frac{
                        \left[
                            {({m_i+m_j})*D_{ij}}
                            \right]-w_{ij}-w_{ji}- \sum\left(
                                \displaystyle\frac{\alpha}{2}N_k(
                                    (D_{ij}+2)([w_{ik}*D_{jk}]+[w_{jk}*D_{ik}])
                                ) + [w_{ik}D_{ik}*D_{jk}] + [w_{jk}D_{jk}*D_{ik}]
                            \right)
                    }{
                        \displaystyle\frac{\alpha}{2}\left(
                            d_i+d_j-b_i-b_j+\sum_k\left(
                                N_k\left(d'_{ik}+d'_{jk}\right)
                            \right)
                        \right)+
                        (b_i+b_j)+w_{ij}+w_{ji}
                    }} & i \neq j \\
                    \displaystyle\frac{
                        \displaystyle\frac{m_i}{m_i}+\left[m_i*D_{ii}\right] - \displaystyle\sum_k
                        \left(
                            \displaystyle\frac{\alpha}{2}N_k
                            \left(
                                \left(D_{ii}+2\right)[w_{ik}*D_{ik}]
                            \right)
                            +
                            [W_{ik}D_{ik}*D_{ik}]
                        \right)
                    }{
                        \displaystyle\frac{\alpha}{2}\left(d_i-b_i+\displaystyle\sum_kN_kd'_{ik}\right) + b_i + w_{ii}
                    }
                \end{array}
            \right.
        \f$
        \sa \f$D\f$ D
        \sa \f$d'\f$ dmat
        \sa \f$d\f$ dvec
        \sa \f$w\f$ w
        \sa \f$m\f$ m
        \param i номер вида
        \param j номер вида
        \return вектор значение
    */
    vec d_calc(ssize_t i, ssize_t j) const {
        auto ones = vec::Ones(points);
        vec left, right;

        if (i != j) {
            left = (dvec(i) + dvec(j) - b(i) - b(j)) * ones;
            for (ssize_t k = 0; k != Species; ++k) {
                left += N(k) * (dmat(i, k) + dmat(j, k)) * ones;
            }
            left *= al / 2;
            left += (b(i) + b(j)) * ones + w(i, j) + w(j, i);

            right = Conv(m(i) + m(j), D(i, j)) - w(i, j) - w(j, i);
            for (ssize_t k = 0; k != Species; ++k) {
                right -= (al / 2) * N(k) * (
                        (D(i, j) + 2) * (Conv(w(i, k), D(j, k)) + Conv(w(j, k), D(i, k))) +
                        Conv(w(i, k) * D(i, k), D(j, k)) +
                        Conv(w(j, k) * D(j, k), D(i, k))
                );
            }

            return (right / left);
        } else {
            left = (dvec(i) - b(i)) * ones;
            for (ssize_t k = 0; k != Species; ++k) {
                left += (N(k) * dmat(i, k)) * ones;
            }
            left *= al / 2;
            left += (b(i) * ones) + w(i, i);

            right = (m(i) / N(i)) + Conv(m(i), D(i, i));
            for (ssize_t k = 0; k != Species; ++k) {
                right -= (al / 2) * N(k) * ((D(i, i) + 2) * Conv(w(i, k), D(i, k)) + Conv(w(i, k) * D(i, k), D(i, k)));
            }

            return (right / left);
        }
    }
    size_t solver() {
        N.setConstant(Species, 1, 100);
        y.resize(Species, Species);
        dl eps = 1e-8;
        dl eps2 = 1e-4;
        size_t max_iter = 500;
        mat y_old = mat::Ones(Species, Species);
        auto N_old = N;
        N_old = 100000 * vec::Ones(Species);
        dl mis2 = 1000, mis1 = 1;
        iters = 0;
        PrintN();
        while ((mis1 > eps || mis2 > eps2) && iters < max_iter) {
            for (ssize_t i = 0; i < Species; ++i) {
                for (ssize_t j = (i + 1); j < Species; ++j) {
                    if (i < j) {
                        D(i, j) = D(j, i) = d_calc(i, j);
                    }
                }
            }
            N = n_calc();
            PrintN();

            for (ssize_t i = 0; i < Species; ++i) {
                D(i, i) = d_calc(i, i);
            }

            for (ssize_t i = 0; i < Species; ++i) {
                for (ssize_t j = (i + 1); j < Species; ++j) {
                    y(i, j) = y_calc(w(i, j), D(i, j), dmat(i, j), h, points, A, dim);
                }
            }

            N = n_calc();

            mis1 = mis_calc(y, y_old);
            y_old = y;
            ++iters;
            mis2 = 0;
            for (ssize_t i = 0; i < Species; ++i) {
                mis2 += abs(N(i) - N_old(i));
            }
            N_old = N;
            PrintN();
        }
        return iters;
    }
};


//! Подсчёт \f$y_{i,j}\f$
/*! Считает \f$y_{i,j}=\displaystyle\int_{R^n}w_{ij}(\xi)C_{ij}(\xi)d\xi\f$ с помощью метода прямоугольников
    \param w \ref Base::w ядро конкурентной смертности
    \param D \ref Base::D второй момент
    \param d \ref Base::dmat коэффициент конкурентной смертности
    \param h \ref Base::h расстояние между точками
    \param N \ref Base::points количество точек
    \param A \ref Base::A размер области
    \param dim \ref Base::dim размерность пространства
 */
dl y_calc(
        vec w,
        vec D,
        dl d,
        dl h,
        ssize_t N,
        dl A,
        ssize_t dim
) {
    vec r;
    switch (dim) {
        case 1:
            return  h * (w * D ).head(N - 1).sum() + d;
        case 2:
            r.setLinSpaced(N, 0, A);
            return 2 * PI * h * (w * D * r).head(N - 1).sum() + d;
        case 3:
            r.setLinSpaced(N, 0, A);
            return 4 * PI * h * (w * D * r * r).head(N - 1).sum() + d;
        default:
            throw std::runtime_error("Unspecifinded dimention\n");
    }
}

void plotD(const Base &A, const std::string &fname) {
    using std::ofstream;
    using std::string;
    namespace fs = std::__fs::filesystem;
    std::cout << "Begin write files" << std::endl;
    std::cout << "Begin write N" << std::endl;
    {
        ofstream out(fname + ".N", std::ofstream::trunc);
        for (ssize_t i = 0; i < A.Species; ++i) {
            out << A.N(i) << " ";
        }
        out << "\n";
    }
    std::cout << "End write N" << std::endl;
    std::cout << "Begin write other" << std::endl;
    {
        ofstream out(fname + ".other", std::ofstream::trunc);
        for (ssize_t i = 0; i < A.Species; ++i) {
            for (ssize_t j = 0; j < A.Species; ++j) {
                out << A.y(i, j) << " ";
            }
            out << "\n";
        }
        out << A.iters << " ";
        for (ssize_t i = 1; i < A.Species; ++i) {
            out << "0 ";
        }
        out << "\n";
    }
    std::cout << "End write other" << std::endl;
    std::cout << "Begin write data" << std::endl;
    {
        ofstream out(fname + ".data", std::ofstream::trunc);
        for (ssize_t i = 0; i < A.Species; ++i) {
            for (ssize_t j = 0; j < A.Species; ++j) {
                if (i <= j) {
                    for (ssize_t k = 0; k < A.D(i, j).size(); ++k) {
                        out << (A.D(i, j)(k) + 1.)  << " ";
                    }
                    out << "\n";
                }
            }
        }
    }
    std::cout << "End write data" << std::endl;
    std::cout << "Begin write py" << std::endl;
    ofstream out(fname + ".py", std::ofstream::trunc);
    out << "#!/usr/bin/env python3\n\n" <<
        "import matplotlib.pyplot as plt\n" <<
        "import numpy as np\n\n" <<
        "x = np.linspace(0, " << A.A << ", " << A.points << ")\n" <<
        "D = np.genfromtxt(\'" << (fname + ".data") << "\')\n";
    ssize_t k = 0;
    for (ssize_t i = 0; i < A.Species; ++i) {
        for (ssize_t j = 0; j < A.Species; ++j) {
            if (i <= j) {
                out << "plt.plot(x, D[" << (k++) << "], label=r'$D_{" << i + 1 << ", " << j + 1 << "}$')\n";
            }
        }
    }
    out <<
        "plt.legend()\n" <<
        "plt.show()\n";
    out.close();
    fs::permissions(fname + ".py", fs::perms::owner_all | fs::perms::group_all | fs::perms::others_exec);
    std::cout << "End write py" << std::endl;
    std::cout << "End write files" << std::endl;
}

enum class Examples {
    Solve,
    Help,
};

vec MakeArrayX(const std::vector<dl>& in) {
    vec res;
    ssize_t width = static_cast<ssize_t>(in.size());
    res.resize(width);
    for (auto i = 0; i < width; ++i) {
        res(i) = in[static_cast<size_t>(i)];
    }
    return res;
}

mat MakeArrayXX(const std::vector<dl>& in, ssize_t width) {
    mat res;
    ssize_t height = static_cast<ssize_t>(in.size()) / width;
    res.resize(width, height);
    for (auto i = 0; i < width; ++i) {
        for (ssize_t j = 0; j < height; ++j) {
            res(i, j) = in[static_cast<size_t>(j * width + i)];
        }
    }
    return res;
}

int main(int argc, char* argv[]) {
    using clipp::option;
    using clipp::value;
    using clipp::values;
    using clipp::command;
    using clipp::parse;
    using clipp::make_man_page;
    using std::to_string;

    Base A;
    std::string Title;

    A.al = 0.4;
    A.dim = 1;
    A.points = 512;
    A.A = 2;
    A.Species = 3;
    A.b = vec::Constant(A.Species, 1, 0.4);
    A.dvec = vec::Constant(A.Species, 1, 0.2);
    A.dmat = mat::Constant(A.Species, A.Species, 0.001);
    A.sw = mat::Constant(A.Species, A.Species, 0.04);
    A.sm = vec::Constant(A.Species, 1, 0.04);

    std::vector<dl> b, dvec, dmat, sw, sm;

    auto mode = Examples::Solve;
    auto solve_mode = (
            clipp::required("-t", "--title") & value("Title", Title),
                    option("-species") & value(to_string(A.Species), A.Species),
                    option("-dim") & value(to_string(A.dim), A.dim),
                    option("-al") & value(to_string(A.al), A.al),
                    option("-points") & value(to_string(A.points), A.points),
                    option("-a") & value(to_string(A.A), A.A),
                    option("-b") & values("b [Species]", b),
                    option("-dvec") & values("d [Species]", dvec),
                    option("-dmat") & values("d\' [Species * Species]", dmat),
                    option("-sw") & values("sigma w [Species * Species]", sw),
                    option("-sm") & values("sigma m [Species]", sm),
                    option("-h", "--help").set(mode, Examples::Help)
    );

    auto help_mode = (command("help").set(mode, Examples::Help));
    auto cli = (solve_mode | help_mode);
    vec r;
    if (parse(argc, argv, cli)) {
        switch (mode) {
            case Examples::Solve:
                r.setLinSpaced(A.points, 0, A.A);
                A.h = r(1) - r(0);
                if (!b.empty()) {
                    if (static_cast<ssize_t>(b.size()) == A.Species) {
                        A.b = MakeArrayX(b);
                    } else {
                        std::cout << "b size: " << b.size() << " != Species: " << A.Species << std::endl;
                        return 1;
                    }
                }
                if (!dvec.empty()) {
                    if (static_cast<ssize_t>(dvec.size()) == A.Species) {
                        A.dvec = MakeArrayX(dvec);
                    } else {
                        std::cout << "d size: " << dvec.size() << " != Species: " << A.Species << std::endl;
                        return 1;
                    }
                }
                if (!dmat.empty()) {
                    if (static_cast<ssize_t>(dmat.size()) == A.Species * A.Species) {
                        A.dmat = MakeArrayXX(dmat, A.Species);
                    } else {
                        std::cout << "d' size: " << dmat.size() << " != Species * Species: " << A.Species * A.Species << std::endl;
                        return 1;
                    }
                }
                if (!sw.empty()) {
                    if (static_cast<ssize_t>(sw.size()) == A.Species * A.Species) {
                        A.sw = MakeArrayXX(sw, A.Species);
                    } else {
                        std::cout << "sigma w size: " << sw.size() << " != Species * Species: " << A.Species * A.Species << std::endl;
                        return 1;
                    }
                }
                if (!sm.empty()) {
                    if (static_cast<ssize_t>(sm.size()) == A.Species) {
                        A.sm = MakeArrayXX(sm, A.Species);
                    } else {
                        std::cout << "sigma m size: " << sm.size() << " != Species * Species: " << A.Species * A.Species << std::endl;
                        return 1;
                    }
                }
                A.m = normpdf(r, vec::Zero(A.Species), A.sm, A.b);
                A.w = normpdf(r, mat::Zero(A.Species, A.Species), A.sw, A.dmat);
                A.D.resize(A.Species, A.Species);
                for (ssize_t i = 0; i < A.Species; ++i) {
                    for (ssize_t j = 0; j < A.Species; ++j) {
                        A.D(i, j) = vec::Zero(A.points);
                    }
                }
                A.solver();
                for (auto i = 0; i < A.Species; ++i) {
                    std::cout << "N[" << i << "] = " << A.N(i) << std::endl;
                }
                plotD(A, Title);
                return 0;
            default:
                std::cout << make_man_page(cli, argv[0]) << std::endl;
                return 0;
        }
    } else {
        std::cout << clipp::usage_lines(cli, argv[0]) << std::endl;
        return 0;
    }

}
