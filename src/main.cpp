
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
    Z_ = V_.array() / I_.array();
    z_in = 1.0 / I_(central_element);
    const auto zin_test = V_(central_element) / I_(central_element);
    std::iota(xs_.begin(), xs_.end(), 0);
    std::ranges::transform(
        xs_, xs_.begin(), [this](double n) { return n * delta_ - length_ / 2.0; });
    std::ranges::transform(
        I_, std::next(ys_.begin()), [](std::complex<double> value) { return std::abs(value); });
  }

  // Accessors
  const Eigen::VectorXcd& GetCurrents() const { return I_; }
  const std::vector<double>& get_xs() const { return xs_; }
  const std::vector<double>& get_ys() const { return ys_; }
  const std::vector<double>& get_resistance() const { return resistance_; }
  const std::vector<double>& get_reactance() const { return reactance_; }

  int get_elements() const { return elements_; }

 private:
  Eigen::VectorXcd Solve() const {
    Eigen::FullPivHouseholderQR<Eigen::MatrixXcd> solver(Z_);
    return solver.solve(V_);
  }

  std::complex<double> Psi(double m, double n) const {
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

  std::complex<double> Impedance(int m, int n) const {
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
  std::vector<double> resistance_;
  std::vector<double> reactance_;
};

void PlotDipoleCurrent(const Dipole& dipole) {
  const auto& xs = dipole.get_xs();
  const auto& ys = dipole.get_ys();
  ImPlot::PlotLine(
      std::format("N={}", dipole.get_elements()).c_str(), xs.data(), ys.data(), xs.size());
}

void PlotDipoleReactanceResistance(const Dipole& dipole) {}

void PlotResistances(const Dipole& dipole) {}

void AdjustableDipole() {
  if (ImGui::Begin("Adjustable")) {
    static int segs = 3;
    constexpr double frequency = 2.4 * 1e9;
    constexpr double wavelenght = kSpeedOfLight / frequency;
    constexpr double len = wavelenght / 2;
    constexpr double radius = wavelenght * 1e-4;
    static auto dipole = std::make_unique<Dipole>(len, frequency, radius, segs);

    if (ImGui::SliderInt("segsments", &segs, 3, 101)) {
      segs = segs & 1 ? segs : segs + 1;
      dipole = std::make_unique<Dipole>(len, frequency, radius, segs);
    }

    if (ImPlot::BeginPlot("##AdjustableDipole")) {
      PlotDipoleCurrent(*dipole);
      ImPlot::EndPlot();
    }
  }
  ImGui::End();
}

namespace {

using Dipoles = std::vector<Dipole>;

}

void ShowInterface(const Dipoles& dipoles) {
  auto* window = GuiInit();

  // Main loop
  while (!glfwWindowShouldClose(window)) {
    GuiNewFrame();
    glfwPollEvents();

    if (ImGui::Begin("Dipole")) {
      if (ImPlot::BeginPlot("##Dipole")) {
        for (auto& dipole : dipoles) {
          PlotDipoleCurrent(dipole);
        }
        ImPlot::EndPlot();
      }
    }
    ImGui::End();

    AdjustableDipole();

    ClearBackGround(window);
    Render(window);
  }

  GuiTerminate(window);
}

int main() {
  // Inútil para esse problema, só importa a razão do comprimento de onda.
  // double frequency = 2.4 * 1e6;
  const double frequency = 2.4 * 1e9;
  double wavelenght = kSpeedOfLight / frequency;
  double len = wavelenght / 2;
  double radius = wavelenght * 1e-4;

  int elements = 3;

  std::vector<Dipole> dipoles;
  for (auto& elems : {3, 7, 19, 29}) {
    dipoles.emplace_back(len, frequency, radius, elems);
  }

  ShowInterface(dipoles);

  return 0;
}
