#ifndef APP_H
#define APP_H

#include "active_bandstructure.h"
#include "bandstructure.h"
#include "context.h"

/** Base class for launching apps.
 * An app is a subprogram, used e.g. to compute phonon transport, or a DoS.
 * The base class is almost just an interface, with the implementation given
 * in the derived classes.
 *
 * To create a new app,
 * 1) create a new subclass of App.
 * 2) implement the run() method, with the main body of the subprogram.
 * 3) Implement checkRequirements(), which checks that the user input has
 * loaded the necessary input flags.
 * 4) Update the base class App header, adding a new app name to choices.
 * 5) Update the base class loadApp method in the app.cpp file, adding the
 * possibility of loading the new application.
 */
class App {
    /** Factory method, used to load the desired subclass.
     * Note: auto_ptr transfers the ownership to the returned pointer itself
     * the returned app pointer is then destructed at the end of main().
     * @param choice: a string witht the app name.
     * @return app: a unique_pointer to the subclass app instance.
     */
    static std::unique_ptr<App> loadApp(const std::string &choice);

    /** Launches the subprogram.
     * @param context: object Context with the user input.
     */
    virtual void run(Context &context);

    /** Checks that the user wrote a complete input file.
     * Calls an error if the requirements are not satisfied.
     * To be used before a call to run().
     * The list of requirements is specified by each App subclass.
     * @param context: object Context with the user input.
     */
    virtual void checkRequirements(Context &context);
protected:
    static const inline std::vector<std::string> choices { "phononTransport",
            "phononDos", "electronWannierDos", "electronFourierDos",
            "phononBands", "electronWannierBands", "electronFourierBands",
            "electronPolarization" };

    void throwErrorIfUnset(const std::string &x, const std::string &name);
    void throwErrorIfUnset(const std::vector<std::string> &x,
            const std::string &name);
    void throwErrorIfUnset(const double &x, const std::string &name);
    void throwErrorIfUnset(const Eigen::VectorXi &x, const std::string &name);
    void throwErrorIfUnset(const Eigen::Vector3i &x, const std::string &name);
    void throwErrorIfUnset(const Eigen::VectorXd &x, const std::string &name);
    void throwErrorIfUnset(const Eigen::MatrixXd &x, const std::string &name);
    void throwErrorIfUnset(const Eigen::Tensor<double, 3> &x,
            const std::string &name);

    void throwWarningIfUnset(const std::string &x, const std::string &name);
};

#endif
