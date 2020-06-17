#ifndef BANDSTRUCTURE_H
#define BANDSTRUCTURE_H

#include "points.h"
#include "state.h"
#include "particle.h"
#include "exceptions.h"
#include "utilities.h"

class BaseBandStructure {
public:
    virtual Particle getParticle();

    virtual Points getPoints();
    virtual Point getPoint(const long &pointIndex);
    virtual long getNumPoints();
    virtual long getNumBands(); // this only works in FullBandStructure
    virtual long hasWindow();

    virtual State getState(Point &point);

    /** Returns a State object
     * Same as getState(Point & point), but the wavevector is identified
     * with its integer index
     * @param pointIndex: index of the wavevector, range [0,numPoints[
     * @return State: a State object evaluated at Point.
     */
    virtual State getState(const long &pointIndex);

    // needed in the BTE
    virtual long getIndex(const WavevectorIndex &ik, const BandIndex &ib);
    virtual long getNumStates();
    virtual const double& getEnergy(const long &stateIndex);
    virtual Eigen::Vector3d getGroupVelocity(const long &stateIndex);
    virtual Eigen::Vector3d getWavevector(const long &stateIndex);

    // we need the same logic to flatten indices
    virtual void setEnergies(Point &point, Eigen::VectorXd &energies_);
    virtual void setEigenvectors(Point &point, Eigen::MatrixXcd &eigenvectors_);
    virtual void setVelocities(Point &point,
            Eigen::Tensor<std::complex<double>, 3> &velocities_);
};

class ActiveBandStructure;
// forward declaration of friend class

/** FullBandStructure is the class that stores the energies, velocities and
 * eigenvectors of a quasiparticle computed on a set of wavevectors (as defined
 * by Points() ) in the Brillouin zone.
 */
class FullBandStructure: public BaseBandStructure {
public:
    /** Constructor of the FullBandStructure
     * @param numBands: an integer with the number of bands in the system
     * @param particle: a Particle object that contains the type of
     * quasiparticle
     * @param withVelocities: a boolean to decide whether to store velocities
     * @param withEigenvectors: a boolean to decide whether to store eigenvecs
     * @param points: the underlying mesh of wavevectors.
     */
    FullBandStructure(long numBands_, Particle &particle_, bool withVelocities,
            bool withEigenvectors, Points &points_);

    FullBandStructure();

    /** Copy constructor
     */
    FullBandStructure(const FullBandStructure &that);

    /** Assignment operator
     */
    FullBandStructure& operator =(const FullBandStructure &that);

    /** Get the Particle object associated with this class
     * @return particle: a Particle object, describing e.g. whether this
     * is a phonon or electron bandStructure
     */
    Particle getParticle();

    /** Returns the wavevectors on which the bandstructure is computed.
     * @return Points: the object representing the Brillouin zone wavevectors.
     */
    Points getPoints();

    /** Returns a wavevector, given a wavevector index.
     * The wavevector index runs from 0 to numPoints-1
     */
    Point getPoint(const long &pointIndex);

    /** Returns the total number of k/q-points.
     * @return numPoints: the total number of wavevectors of the bandStructure.
     */
    long getNumPoints();

    /** Returns the number of bands.
     * @return numPoints: the total number of wavevectors of the bandStructure.
     */
    long getNumBands();

    long hasWindow();

    /** Returns a State object.
     * The state object, defined elsewhere, is a container that holds all bands
     * eigenvalues, eigenvectors and velocities (if available), at a fixed
     * wavevector. This is used in particular for the construction of the
     * scattering operator.
     * @param point: a Point object containing the desired wavevector
     * @return State: a State object evaluated at Point.
     */
    State getState(Point &point);

    /** Returns a State object.
     * The state object, defined elsewhere, is a container that holds all bands
     * eigenvalues, eigenvectors and velocities (if available), at a fixed
     * wavevector. This is used in particular for the construction of the
     * scattering operator.
     * @param pointIndex: an integer index of the desired wavevector.
     * @return State: a State object evaluated at Point.
     */
    State getState(const long &pointIndex);

    /** Builds a Bloch state index, which runs on both wavevector index and
     * band index. ik runs from 0 to numPoints-1, ib from 0 to numBands-1.
     * It's used to view the various matrices such as energy as a 1D vector,
     * and can be used in combination with get() methods.
     * @param wavevectorIndex: strong-typed index on wavevector
     * @return stateIndex: integer from 0 to numStates-1=numBands*numPoints-1
     */
    long getIndex(const WavevectorIndex &ik, const BandIndex &ib);

    /** Returns the total number of Bloch states, equal to numPoints*numBands.
     * @return numStates: the total number of Bloch states in the class.
     */
    long getNumStates();

    /** Returns the energy of a quasiparticle from its Bloch index
     * Used for accessing the bandstructure in the BTE.
     * @param stateIndex: an integer index in range [0,numStates[
     * @return energy: the value of the QP energy for that given Bloch index.
     * Phonon energies are referred to zero, with negative energies being
     * actually complex phonon frequencies. Electronic energies are not saved
     * with any particular reference, and should be used together with the
     * chemical potential computed by StatisticsSweep. By policy, it's in
     * rydbergs units.
     */
    const double& getEnergy(const long &stateIndex);

    /** Returns the energy of a quasiparticle from its Bloch index
     * Used for accessing the bandstructure in the BTE.
     * @param stateIndex: an integer index in range [0,numStates[
     * @return velocity: a 3d vector with velocity. By policy, we save it in
     * the cartesian basis and in atomic rydberg units.
     */
    Eigen::Vector3d getGroupVelocity(const long &stateIndex);

    /** Returns the energy of a quasiparticle from its Bloch index
     * Used for accessing the bandstructure in the BTE.
     * @param stateIndex: an integer index in range [0,numStates[
     * @return wavevector: a 3d vector with the wavevector in cartesian
     * coordinates in units of Bohr^-1.
     */
    Eigen::Vector3d getWavevector(const long &stateIndex);

    /** This method overrides setEnergies, but uses a Point object to find the
     * k-point.
     */
    void setEnergies(Point &point, Eigen::VectorXd &energies_);

    /** Method to save quasiparticle eigenvectors inside FullBandStructure().
     * Note that in this case, eigenvectors are passed as a matrix, which is
     * the case e.g. for the Wannier interpolation, where the eigenvectors
     * represent the unitary transformation matrix U for Wannier localization.
     * @param point: a Point object with the coordinates of the wavevector,
     * which should come from the same Point class stored in FullBandStructure
     * @param eigenvectors: a complex matrix of size (numBands,numBands)
     */
    void setEigenvectors(Point &point, Eigen::MatrixXcd &eigenvectors_);

    /** Saves in the class the velocities computed at a particular point.
     * @param point: a Point object representing the wavevector where these
     * velocities have been computed.
     * @param velocities: a rank-3 tensor of size (numBands,numBands,3)
     * containing the matrix elements of the velocity operator. Diagonal
     * elements are the quasiparticle group velocities.
     */
    void setVelocities(Point &point,
            Eigen::Tensor<std::complex<double>, 3> &velocities_);

    // class specific method

    /** Method to save quasiparticle eigenvectors inside FullBandStructure().
     * @param point: a vector of 3 crystal coordinates. The method will look
     * for the wavevector index.
     * @param energies: a vector of size (numBands) with the quasiparticle
     * energies
     */
    void setEnergies(Eigen::Vector3d &point, Eigen::VectorXd &energies_);

    /** Returns all electronic energies for all wavevectors at fixed band index
     * Used by the Fourier interpolation of the band structure.
     * @param bandIndex: index in [0,numBands[ for the quasiparticle band
     * @return energies: a vector of size numPoints with the QP energies.
     * Energies are ordered as the underlying wavevectors in Points.
     */
    Eigen::VectorXd getBandEnergies(long &bandIndex);

protected:
    // stores the quasiparticle kind
    Particle particle;

    // link to the underlying mesh of points
    Points &points; // these may be FullPoints or PathPoints

    // matrices storing the raw data
    Eigen::MatrixXd energies; // size(bands,points)
    Eigen::MatrixXcd velocities; // size(3*bands^2,points)
    Eigen::MatrixXcd eigenvectors; // size(bands^2,points)

    // pointers to the raw data, used to move
    double *rawEnergies = nullptr;
    std::complex<double> *rawVelocities = nullptr;
    std::complex<double> *rawEigenvectors = nullptr;
    // these are integers used to move the pointers in the raw data
    long energiesCols;
    long velocitiesCols;
    long eigenvectorsCols;

    // auxiliary variables
    long numBands = 0;
    long numAtoms = 0;

    bool hasEigenvectors = false;
    bool hasVelocities = false;

    // method to find the index of the kpoint, from its crystal coordinates
    long getIndex(Eigen::Vector3d &pointCoords);

    // ActiveBandStructure, which restricts the FullBandStructure to a subset
    // of wavevectors, needs low-level access to the raw data (for now)
    friend class ActiveBandStructure;
};

#endif
