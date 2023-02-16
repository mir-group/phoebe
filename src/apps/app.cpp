#include "app.h"
#include "bands_app.h"
#include "context.h"
#include "dos_app.h"
#include "electron_wannier_transport_app.h"
#include "elph_plot_app.h"
#include "elel_plot_app.h"
#include "elph_qe_to_phoebe_app.h"
#include "exceptions.h"
#include "lifetimes_app.h"
#include "phonon_transport_app.h"
#include "transport_epa_app.h"
#include "elel_to_phoebe.h"
#include <cmath>
#include <string>

// app factory
std::unique_ptr<App> App::loadApp(const std::string &choice) {
  const std::vector<std::string> choices = {"phononTransport",
                                            "phononDos",
                                            "electronWannierDos",
                                            "electronFourierDos",
                                            "phononBands",
                                            "electronWannierBands",
                                            "electronFourierBands",
                                            "electronWannierTransport",
                                            "elPhQeToPhoebe",
                                            "elElToPhoebe",
                                            "elPhCouplingPlot",
                                            "elElCouplingPlot",
                                            "electronLifetimes",
                                            "phononLifetimes",
                                            "transportEpa"};

  // check if the app choice is valid, otherwise we stop.
  if (std::find(choices.begin(), choices.end(), choice) == choices.end()) {
    Error("The app name is not valid, didn't find an app to launch.");
  }

  if (choice == "phononTransport") {
    return std::unique_ptr<App>(new PhononTransportApp);
  } else if (choice == "electronWannierTransport") {
    return std::unique_ptr<App>(new ElectronWannierTransportApp);
  } else if (choice == "elPhQeToPhoebe") {
    return std::unique_ptr<App>(new ElPhQeToPhoebeApp);
  } else if (choice == "elElToPhoebe") {
    return std::unique_ptr<App>(new ElElToPhoebeApp);
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
  } else if ( choice == "transportEpa" ) {
    return std::unique_ptr<App>(new TransportEpaApp);
  } else if (choice == "elPhCouplingPlot") {
    return std::unique_ptr<App>(new ElPhCouplingPlotApp);
  } else if (choice == "elElCouplingPlot") {
    return std::unique_ptr<App>(new ElElCouplingPlotApp);
  } else if (choice == "electronLifetimes") {
    return std::unique_ptr<App>(new ElectronLifetimesApp);
  } else if (choice == "phononLifetimes") {
    return std::unique_ptr<App>(new PhononLifetimesApp);
  } else {
    return std::unique_ptr<App>(nullptr);
  }
}

void App::run(Context &context) {
  (void)context; // suppress unused variable compiler warnings
  Error("Base class app doesn't have a run()");
}

// no requirements for the base class.
void App::checkRequirements(Context &context) {
  (void)context; // suppress unused variable compiler warnings
}

void App::throwErrorIfUnset(const std::string &x, const std::string &name) {
  if (x.empty()) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const std::vector<std::string> &x,
                            const std::string &name) {
  if (x.empty()) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const double &x, const std::string &name) {
  if (std::isnan(x)) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::VectorXi &x, const std::string &name) {
  if (x.size() == 0) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::Vector3i &x, const std::string &name) {
  if (x.size() == 0) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::VectorXd &x, const std::string &name) {
  if (x.size() == 0) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::MatrixXd &x, const std::string &name) {
  if (x.rows() == 0) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwErrorIfUnset(const Eigen::Tensor<double, 3> &x,
                            const std::string &name) {
  if (x.dimension(0) == 0) {
    Error("Input variable " + name + " hasn't been found in input");
  }
}

void App::throwWarningIfUnset(const std::string &x, const std::string &name) {
  if (x.empty()) {
    Warning("Input variable " + name + " hasn't been found in input");
  }
}
