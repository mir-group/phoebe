#include "points.h"

/** Class containing harmonic information for all bands at a given k-q point.
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
	 * @param dnde: derivative of the Fermi/Bose distribution wrt energy.
	 * @param dndt: derivative of the Fermi/Bose distribution wrt temperature
	 */
	State(Point& point_,
			Eigen::VectorXd& energies_,
			Eigen::Tensor<std::complex<double>,3>& velocities_,
			Eigen::VectorXd& dnde_, Eigen::VectorXd& dndt_);

	/** get the wavevector (Point object)
	 * @return point: a Point object.
	 */
	Point getPoint();

	/** get the energy of a single band
	 * @param bandIndex: integer from 0 to numBands-1
	 * @return energy: Bloch state energy in rydbergs.
	 */
	double getEnergy(const int bandIndex);

	/** get all energies at a given point
	 * @return energies: a vector of energies in rydbergs for all bands present
	 */
	Eigen::VectorXd getEnergies();

	/** get the group velocity of a single Bloch state.
	 * @param bandIndex: integer from 0 to numBands-1.
	 * @return velocity: the 3d-vector of the group velocity.
	 */
	Eigen::Vector3d getGroupVelocity(const int bandIndex);

	/** get the group velocities of all bands for given k/q point
	 * @return groupVelocities: double matrix of size (numBands,3) with the
	 * group velocities.
	 */
	Eigen::MatrixXd getGroupVelocities();

	/** get the velocities of all bands for given k/q point
	 * @return Velocities: a complex tensor of dimensions (numBands,numBands,3)
	 * with the velocity operator.
	 */
	Eigen::Tensor<std::complex<double>,3> getVelocities();

	/** get all values of dn/dT (the derivative of the equilibrium
	 * distribution wrt temperature) for the given k/q point.
	 * @return dndt: a vector of size (numBands)
	 */
	Eigen::VectorXd getDndt();

	/** get all values of dn/de (the derivative of the equilibrium
	 * distribution wrt energy) for the given k/q point.
	 * @return dnde: a vector of size (numBands)
	 */
	Eigen::VectorXd getDnde();

	/** get the weight of the k/q point. Used for integrations over the
	 * brillouin zone with an irreducible mesh of points.
	 * @return weight: a double.
	 */
	double getWeight();
protected:
	// pointers to the bandstructure, I don't want to duplicate storage here
	Point& point; // in cryst coords
	Eigen::VectorXd energies;
	Eigen::VectorXd dndt;
	Eigen::VectorXd dnde;
	Eigen::Tensor<std::complex<double>,3> velocities;
	int numBands;
};

/** Describes electronic Bloch states at fixed kpoint.
 * See the documentation of State for further details
 */
class ElState: public State {
public:
	ElState(Point& point_,
			Eigen::VectorXd& energies_,
			Eigen::Tensor<std::complex<double>,3>& velocities_,
			Eigen::VectorXd& dnde_, Eigen::VectorXd& dndt_);
};

/** Describes phonon Bloch states at fixed qpoint.
 * See the documentation of State for further details.
 * The only difference for public members is the addition of the phonon
 * eigenvectors in the constructor.
 */
class PhState: public State {
public:
	PhState(Point& point_,
			Eigen::VectorXd& energies_,
			Eigen::Tensor<std::complex<double>,3>& eigenvectors_,
			Eigen::Tensor<std::complex<double>,3>& velocities_,
			Eigen::VectorXd& dnde_,
			Eigen::VectorXd& dndt_);
	/** get the eigenvectors for the current Point
	 * @return eigenvectors: a complex tensor of size (3,numAtoms,numBands)
	 * with the phonon eigenvectors
	 */
	Eigen::Tensor<std::complex<double>,3> getEigenvectors();
protected:
	Eigen::Tensor<std::complex<double>,3> eigenvectors;
};
