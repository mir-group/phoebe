#include "app.h"
#include "active_bandstructure.h"
#include "bands_app.h"
#include "bandstructure.h"
#include "constants.h"
#include "context.h"
#include "dos_app.h"
#include "exceptions.h"
#include "io.h"
#include "particle.h"
#include "phonon_transport_app.h"
#include "points.h"
#include "polarization_app.h"
#include "qe_input_parser.h"
#include "utilities.h"
#include "window.h"
#include <cmath>
#include <string>

// app factory
std::unique_ptr<App> App::loadApp(const std::string &choice) {
  const std::vector<std::string> choices = {
      "phononTransport",      "phononDos",           "electronWannierDos",
      "electronFourierDos",   "phononBands",         "electronWannierBands",
      "electronFourierBands", "electronPolarization"};

  // check if the app choice is valid, otherwise we stop.
  if (std::find(choices.begin(), choices.end(), choice) == choices.end()) {
    Error e("The app name is not valid, didn't find an app to launch.");
  }

  if (choice == "phononTransport") {
    return std::unique_ptr<App>(new PhononTransportApp);
  } else if (choice == "phononDos") {
    return std::unique_ptr<App>(new PhononDosApp);
  } else if (choice == "electronWannierDos") {
    return std::unique_ptr<App>(new ElectronWannierDosApp);
  } else if (choice == "electronFourierDos") {
    return std::unique_ptr<App>(new ElectronFourierDosApp);
  } else if (choice == "phononBands") {
    return std::unique_ptr<App>(new PhononBandsApp);
  } else if (choice == "electronWannierBands") {
    return std::unique_ptr<App>(new ElectronWannierBandsApp);
  } else if (choice == "electronFourierBands") {
    return std::unique_ptr<App>(new ElectronFourierBandsApp);
  } else if (choice == "electronPolarization") {
    return std::unique_ptr<App>(new ElectronPolarizationApp);
  } else {
    return std::unique_ptr<App>(nullptr);
  }
}

void App::run(Context &context) {
  (void)context; // suppress unused variable compiler warnings
  Error e("Base class app doesn't have a run()");
}

// no requirements for the base class.
void App::checkRequirements(Context &context) {
  (void)context; // suppress unused variable compiler warnings
}

void App::throwErrorIfUnset(const std::string &x, const std::string &name) {
  if (x.empty()) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const std::vector<std::string> &x,
                            const std::string &name) {
  if (x.size() == 0) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const double &x, const std::string &name) {
  if (std::isnan(x)) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::VectorXi &x, const std::string &name) {
  if (x.size() == 0) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::Vector3i &x, const std::string &name) {
  if (x.size() == 0) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::VectorXd &x, const std::string &name) {
  if (x.size() == 0) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::MatrixXd &x, const std::string &name) {
  if (x.rows() == 0) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::Tensor<double, 3> &x,
                            const std::string &name) {
  if (x.dimension(0) == 0) {
    Error e("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwWarningIfUnset(const std::string &x, const std::string &name) {
  if (x.empty()) {
    Warning e("Input variable " + name + " hasn't been found in input");
  }
}
