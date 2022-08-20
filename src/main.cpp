
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
        currents_(segments + 2, 0) {
    const int central_element = elements_ / 2;
    const double omega = 2 * pi * frequency_;
    v_in = std::complex(0.0, -omega * kEpsilon0);
    V_(central_element) = v_in;
    FillImpedances();
    I_ = Solve();
    i_in_ = I_(central_element);
    Z_ = V_.array() / I_.array();
    z_in_ = 1.0 / i_in_;
    std::iota(xs_.begin(), xs_.end(), 0);
    std::ranges::transform(
        xs_, xs_.begin(), [this](double n) { return n * delta_ - length_ / 2.0; });
    std::ranges::transform(I_, std::next(currents_.begin()), [](std::complex<double> value) {
      return std::abs(value);
    });
  }

  // Accessors
  const Eigen::VectorXcd& GetCurrents() const { return I_; }
  const std::vector<double>& get_xs() const { return xs_; }
  const std::vector<double>& get_currents() const { return currents_; }

  std::complex<double> get_input_impedance() const { return z_in_; }
  std::complex<double> get_input_current() const { return i_in_; }

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

    if (m == n) {
      return std::log(delta_ / radius_) / (2 * pi * delta_) - std::complex(0.0, k_ / (4 * pi));
    }
    const auto pow = std::complex(0.0, -k_ * r);
    const auto num = std::exp(pow);
    return num / (4 * pi * r);
  }

  std::complex<double> Impedance(int m, int n) const {
    const auto a_mn = Square(delta_) * Psi(m, n);
    const auto phi_mn = Psi(m - 0.5, n - 0.5) - Psi(m - 0.5, n + 0.5) - Psi(m + 0.5, n - 0.5) +
                        Psi(m + 0.5, n + 0.5);
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

  std::complex<double> z_in_;
  std::complex<double> v_in;
  std::complex<double> i_in_;
  Eigen::MatrixXcd Z_;
  Eigen::VectorXcd V_;
  Eigen::VectorXcd I_;

  std::vector<double> xs_;
  std::vector<double> currents_;
};

void PlotDipoleCurrent(const Dipole& dipole) {
  const auto& xs = dipole.get_xs();
  const auto& ys = dipole.get_currents();

  ImPlot::PlotLine(
      std::format("N={}", dipole.get_elements()).c_str(), xs.data(), ys.data(), xs.size());
}

void PlotDipoleReactanceResistance(const Dipole& dipole) {
  const auto& xs = dipole.get_xs();
  const auto& ys = dipole.get_currents();
  ImPlot::PlotLine(
      std::format("N={}", dipole.get_elements()).c_str(), xs.data(), ys.data(), xs.size());
}

void AdjustableDipole() {
  if (ImGui::Begin("Adjustable")) {
    static int segs = 3;
    constexpr double frequency = 2.4 * 1e6;
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

struct ImpedancePlotter {
  ImpedancePlotter(double len, double frequency, double radius) {
    reactances.reserve(128);
    resistances.reserve(128);
    elems.reserve(128);
    for (int segs = 3; segs < 205; segs += 2) {
      Dipole dipole{len, frequency, radius, segs};
      const auto input_impedance = dipole.get_input_impedance();
      elems.emplace_back(segs);
      resistances.emplace_back(input_impedance.real());
      reactances.emplace_back(input_impedance.imag());
    }
  };

  void Plot() const {
    ImPlot::PlotLine("Resistance", elems.data(), resistances.data(), resistances.size());
    ImPlot::PlotLine("Reactance", elems.data(), reactances.data(), reactances.size());
  }

  std::vector<double> reactances;
  std::vector<double> resistances;
  std::vector<double> elems;
};

void ShowInterface(const Dipoles& dipoles, const ImpedancePlotter& impedance_plotter) {
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

    if (ImGui::Begin("Convergência com N -- Impedâncias")) {
      if (ImPlot::BeginPlot("##Impedances")) {
        impedance_plotter.Plot();
      }
      ImPlot::EndPlot();
    }
    ImGui::End();

    ImPlot::SetNextLineStyle(ImVec4{1.0, 1.0, 1.0, 1.0});
    AdjustableDipole();

    ClearBackGround(window);
    Render(window);
  }

  GuiTerminate(window);
}

int main() {
  // Inútil para esse problema, só importa a razão do comprimento de onda.
  // double frequency = 2.4 * 1e9;
  const double frequency = 2.4 * 1e9;
  double wavelenght = kSpeedOfLight / frequency;
  double len = wavelenght / 2;
  double radius = wavelenght * 1e-4;
  int elements = 3;

  std::vector<Dipole> dipoles;
  for (auto& elems : {3, 7, 19, 29}) {
    dipoles.emplace_back(len, frequency, radius, elems);
  }

  ImpedancePlotter impedance_plotter{len, frequency, radius};

  ShowInterface(dipoles, impedance_plotter);

  return 0;
}
