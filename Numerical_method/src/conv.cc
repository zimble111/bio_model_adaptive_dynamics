#include <cstddef>
#include <gsl/gsl_dht.h>
#include <mutex>

#include <unsupported/Eigen/FFT>

#include "conv.h"

constexpr dl PI = 3.1415926535897932384626433832795;

//! Выбор типа свёртки
enum class ConvType {
    full, /*!< Полная свёртка */
    same, /*!< Центральная часть свёртки таково же размера как и первый аргумент */
};

//! Свёртка
/*! Свёртка делается с помощью быстрого преобразования фурье
    \param ina первый вектор свёртки
    \param inb второй вектор свёртки
    \param type тип свёртки
    \return вектор
 */
vec conv(const vec &ina, const vec &inb, ConvType type = ConvType::same) {
    assert(ina.size() == inb.size());
    auto size = ina.size() * 2;
    Eigen::VectorXd a(size), b(size);
    a.setZero();
    b.setZero();
    a.head(ina.size()) = ina;
    b.head(inb.size()) = inb;
    Eigen::VectorXcd ca, cb;
    Eigen::VectorXd res;
    static Eigen::FFT<dl> fft;
    fft.fwd(ca, a);
    fft.fwd(cb, b);
    ca = ca.cwiseProduct(cb);
    fft.inv(res, ca);
    if (type == ConvType::same) {
        return res.segment((size - 1 - inb.size()) / 2, ina.size());
    }
    return res.head(size - 1);
}

//! Мьютекс для преобразования Ханкеля
/*! Нужен был для паралельного вычисления, что бы не пересекались. Сейчас возможно уже не нужен */
std::mutex mut_for_2D_hankel;
//! 
vec hankel_2D_gsl(vec f, double A)
{

    volatile static size_t the_size = 0; // /2
    static gsl_dht *the_dht;

    vec res(f.size());
    if (the_size == 0) {
        std::lock_guard<std::mutex> acc(mut_for_2D_hankel);
        if (the_size == 0) {
            size_t the_size_1 = static_cast<size_t>(f.size());
            the_dht = gsl_dht_new(the_size_1 / 2 + the_size_1 % 2, 0, A);
            the_size = the_size_1;
        }
    }
    double the_in[the_size];
    double the_out[the_size];

    for (size_t i = 0; i < the_size / 2 + the_size % 2; i++) {
        the_in[the_size / 2 + the_size % 2 - i - 1] = f(static_cast<Eigen::Index>(i));
    }

    gsl_dht_apply(the_dht, the_in, the_out);

    for (size_t i = 0; i < the_size / 2 + the_size % 2; i++) {
        res(static_cast<ssize_t>(the_size / 2 + the_size % 2 - i - 1)) =
        res(static_cast<ssize_t>(the_size / 2 + i)) = the_out[i];
    }
    return res;
}

// двумерное преобразование Фурье радиально-симметричных функций
// середина входного вектора считается точкой ноль
vec conv_radial_2D(const vec& a, const vec& b, double A) {
    vec res;
    res = hankel_2D_gsl((hankel_2D_gsl(a, A)).cwiseProduct(hankel_2D_gsl(b, A)) * 4 * PI * PI, A) * (1 / (2 * PI));
    return res;
}

vec conv_dim(const vec &u, const vec &v, ssize_t points, dl A, ssize_t dim) {
    vec r;
    switch (dim) {
        case 1:
            return conv(u, v);
        case 2:
            return conv_radial_2D(u, v, A);
        case 3:
            r.setLinSpaced(points, 0, A);
            r = r * u;
            return 4 * PI * conv(r, v);
        default:
            throw std::runtime_error("Unspecifinded dimention " + std::to_string(dim) + "\n");
    }
}
