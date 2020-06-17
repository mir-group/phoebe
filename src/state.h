#ifndef STATE_H
#define STATE_H

#include "points.h"
#include "exceptions.h"
#include "utilities.h"

class DetachedState {
public:
    // this class works just like state, but doesn't have a reference to a
    // bandstructure.
    // It's basically a simple container of energies and other stuff

    DetachedState(Eigen::Vector3d &point_, Eigen::VectorXd &energies_,
            long numBands1_, long numBands2_, Eigen::MatrixXcd &eigenvectors_,
            Eigen::Tensor<std::complex<double>, 3> *velocities_ = nullptr);
    DetachedState(const DetachedState &that); // copy constructor
    DetachedState& operator=(const DetachedState &that); // assignment op

    Eigen::Vector3d getCoords(const int &basis);
    double getEnergy(const long &bandIndex, double chemicalPotential = 0.);
    Eigen::VectorXd getEnergies(double chemicalPotential = 0.);
    Eigen::Vector3d getVelocity(const long &bandIndex);
    Eigen::Vector3cd getVelocity(const long &bandIndex1,
            const long &bandIndex2);
    Eigen::MatrixXd getGroupVelocities();
    Eigen::Tensor<std::complex<double>, 3> getVelocities();
    void getEigenvectors(Eigen::Tensor<std::complex<double>, 3> &eigs);
    void getEigenvectors(Eigen::MatrixXcd &eigs);
protected:
    // pointers to the bandstructure, I don't want to duplicate storage here
    Eigen::Vector3d point;
    Eigen::VectorXd energies;
    long numBands1;
    long numBands2;
    Eigen::Tensor<std::complex<double>, 3> velocities;
    Eigen::MatrixXcd eigenvectors;
};

/** Class containing harmonic information for all bands at a given k(q) point.
 * State is a base class. PhState and ElState should be used in the code.
 */
class State {
public:
    /** Class constructor
     * @param point: a Point instance with the coordinates in the brillouin
     * zone.
     * @param energies: a vector with all the energies at a given point
     * @param velocities: a complex tensor vector with the velocity operator at
     * this point. Dimensions(numBands,numBands,3). The diagonal over bands is
     * the group velocity, the off-diagonal are linked to the dipole operator.
     */
    State(Point &point_, double *energies_, long numBands1_, long numBands2_,
            std::complex<double> *velocities_ = nullptr,
            std::complex<double> *eigenvectors_ = nullptr);
    State(const State &that); // copy constructor
    State& operator=(const State &that); // assignment operator

    /** get the wavevector (Point object)
     * @return point: a Point object.
     */
    Point getPoint();

    /** get the cartesian coordinates of a wavevector (Point object)
     * @return point: a Eigen::Vector3d object.
     */
    Eigen::Vector3d getCoords(const int &basis);

    /** get the energy of a single band
     * @param bandIndex: integer from 0 to numBands-1
     * @return energy: Bloch state energy in rydbergs.
     */
    double getEnergy(const long &bandIndex, double chemicalPotential = 0.);

    /** get all energies at a given point
     * @return energies: a vector of energies in rydbergs for all bands present
     */
    Eigen::VectorXd getEnergies(double chemicalPotential = 0.);

    /** get the group velocity of a single Bloch state.
     * @param bandIndex: integer from 0 to numBands-1.
     * @return velocity: the 3d-vector of the group velocity.
     */
    Eigen::Vector3d getVelocity(const long &bandIndex);

    /** get the off-diagonal velocity operator of two Bloch states.
     * @param bandIndex1: bloch index of the bra state,
     * integer from 0 to numBands-1.
     * @param bandIndex2: bloch index of the ket state,
     * integer from 0 to numBands-1.
     * @return velocity: the 3d-vector of the velocity.
     */
    Eigen::Vector3cd getVelocity(const long &bandIndex1,
            const long &bandIndex2);

    /** get the group velocities of all bands for given k/q point
     * @return groupVelocities: double matrix of size (numBands,3) with the
     * group velocities.
     */
    Eigen::MatrixXd getGroupVelocities();

    /** get the velocities of all bands for given k/q point
     * @return Velocities: a complex tensor of dimensions (numBands,numBands,3)
     * with the velocity operator.
     */
    Eigen::Tensor<std::complex<double>, 3> getVelocities();

    /** get the weight of the k/q point. Used for integrations over the
     * brillouin zone with an irreducible mesh of points.
     * @return weight: a vector of size (numBands).
     */
    double getWeight();

    /** get the eigenvectors for the current Point
     * @input/output eigenvectors: a complex tensor of size
     * (3,numAtoms,numBands) with the phonon eigenvectors. Error if
     * the eigenvectors are not set.
     */
    void getEigenvectors(Eigen::Tensor<std::complex<double>, 3> &eigs);

    /** get the eigenvectors for the current Point
     * @input/output eigenvectors: a complex matrix of size (numBands,numBands)
     * with the electron eigenvectors. Error if eigenvectors are not set.
     */
    void getEigenvectors(Eigen::MatrixXcd &eigs);
protected:
    // pointers to the bandstructure, I don't want to duplicate storage here
    Point point;
    double *energies;
    long numBands1; // this is full
    long numBands2; // this is reduced by the window
    std::complex<double> *velocities = nullptr;
    std::complex<double> *eigenvectors = nullptr;
    bool hasVelocities = false;
    bool hasEigenvectors = false;
};

#endif
