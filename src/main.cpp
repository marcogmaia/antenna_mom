
#include <cmath>
#include <complex>
#include <format>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include <imgui.h>
#include <implot.h>
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
        V_(Eigen::VectorXcd::Zero(segments)),
        xs_(segments + 2, 0),
        ys_(segments + 2, 0) {
    const int central_element = elements_ / 2;
    const double omega = 2 * pi * frequency_;
    v_in = std::complex(0.0, -omega * kEpsilon0);
    V_(central_element) = v_in;
    FillImpedances();
    I_ = Solve();
    z_in = 1.0 / I_(central_element);
    std::iota(xs_.begin(), xs_.end(), 0);
    std::ranges::transform(
        xs_, xs_.begin(), [this](double n) { return n * delta_ - length_ / 2.0; });
    std::ranges::transform(
        I_, std::next(ys_.begin()), [](std::complex<double> value) { return std::abs(value); });
  }

  Eigen::VectorXcd GetCurrents() { return I_; }

  const std::vector<double>& get_xs() const { return xs_; }

  const std::vector<double>& get_ys() const { return ys_; }

 private:
  Eigen::VectorXcd Solve() {
    Eigen::JacobiSVD<Eigen::MatrixXcd> solver(
        Z_, Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV);
    return solver.solve(V_);
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

  std::vector<double> xs_;
  std::vector<double> ys_;
};

void PlotDipole(std::string_view legend, Dipole& dipole) {
  const auto& xs = dipole.get_xs();
  const auto& ys = dipole.get_ys();
  ImPlot::PlotLine(legend.data(), xs.data(), ys.data(), xs.size());
}

int main() {
  // Inutil para esse problema, só importa a razão do comprimento de onda.
  // double frequency = 2.4 * 1e6;
  const double frequency = 2.4 * 1e9;
  double wavelenght = kSpeedOfLight / frequency;
  double len = wavelenght / 2;
  double radius = wavelenght * 1e-4;

  int elements = 3;

  Dipole dipole1(len, frequency, radius, 3);
  Dipole dipole2(len, frequency, radius, 7);
  Dipole dipole3(len, frequency, radius, 19);
  Dipole dipole4(len, frequency, radius, 55);

  auto* window = GuiInit();
  
  // Main loop
  while (!glfwWindowShouldClose(window)) {
    GuiNewFrame();
    glfwPollEvents();

    if (ImPlot::BeginPlot("Dipole")) {
      PlotDipole("N=3", dipole1);
      PlotDipole("N=7", dipole2);
      PlotDipole("N=19", dipole3);
      PlotDipole("N=55", dipole4);
    }
    ImPlot::EndPlot();

    ClearBackGround(window);
    Render(window);
  }

  GuiTerminate(window);

  return 0;
}
