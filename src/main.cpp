
#include <cmath>
#include <complex>
#include <format>
#include <iostream>
#include <numbers>

#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "gui.h"

namespace {

using std::numbers::pi;
constexpr double kSpeedOfLight = 2.99792458e8;  // Velocidade da luz.
constexpr double kEpsilon0 = 8.85418782e-12;    // Permissividade elétrica vácuo.

template <typename T>
T Square(T n) {
  return n * n;
}

}  // namespace

class Dipole {
 public:
  Dipole(double length, double frequency, double radius, int segments)
      : length_(length),
        frequency_(frequency),
        radius_(radius),
        delta_(length / (segments + 1)),
        wavelength_(kSpeedOfLight / frequency),
        k_((2.0 * pi) / wavelength_),
        half_lenght_(length / 2.0),
        elements_(segments),
        Z_(Eigen::MatrixXcd(segments, segments)),
        V_(Eigen::VectorXcd::Zero(segments)) {
    const double omega = 2 * pi * frequency_;
    v_in = std::complex(0.0, -omega * kEpsilon0);
    z_in = 1.0 / v_in;
    V_.array() = v_in;
    FillImpedances();
    Solve();
  }

  Eigen::VectorXcd GetCurrents() { return I_; }

 private:
  void Solve() {
    Eigen::JacobiSVD<Eigen::MatrixXcd> solver(
        Z_, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV);
    I_ = solver.solve(V_);
  }

  std::complex<double> Psi(double m, double n) {
    const double zm = m * delta_ - half_lenght_;
    const double zn = n * delta_ - half_lenght_;
    const double r = std::sqrt(Square(zm - zn) - Square(radius_));

    const auto m_eq_n = [=] {
      return 1 / (2 * pi * delta_) * std::log(delta_ / radius_) - std::complex(0.0, k_ / (4 * pi));
    }();

    const auto m_neq_n = [=] {
      const auto pow = std::complex(0.0, k_ * r);
      const auto num = std::exp(pow);
      return num / (4 * pi * r);
    }();

    return (m == n) ? m_eq_n : m_neq_n;
  }

  std::complex<double> Impedance(int m, int n) {
    const auto a_mn = Square(delta_) * Psi(m, n);
    const auto phi_mn =
        Psi(m - .5, n - .5) - Psi(m - .5, n + .5) - Psi(m + .5, n - .5) + Psi(m + .5, n + .5);
    return Square(k_) * a_mn - phi_mn;
  }

  void FillImpedances() {
    for (int i = 0; i < elements_; ++i) {
      for (int j = 0; j < elements_; ++j) {
        Z_(i, j) = Impedance(i, j);
      }
    }
  }

  const double length_;
  const double frequency_;
  const double radius_;
  const double delta_;
  const double wavelength_;
  const double k_;
  const double half_lenght_;
  const int elements_;

  std::complex<double> z_in;
  std::complex<double> v_in;
  Eigen::MatrixXcd Z_;
  Eigen::VectorXcd V_;
  Eigen::VectorXcd I_;
};

int main() {
  // Inutil para esse problema, só importa a razão do comprimento de onda.
  double frequency = 2.4 * 1e6;
  double wavelenght = kSpeedOfLight / frequency;
  double len = wavelenght / 2;
  double radius = wavelenght * 1e-4;

  int elements = 3;

  Dipole dipole(len, frequency, radius, elements);
  std::cout << dipole.GetCurrents() << '\n';

  ShowImgui();

  return 0;
}
