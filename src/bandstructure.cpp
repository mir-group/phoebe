#include "bandstructure.h"
#include "exceptions.h"

//auxiliary functions to work with the velocity
Eigen::VectorXcd packXcd(Eigen::Tensor<std::complex<double>,3> a,
		int size1, int size2, int size3) {
	Eigen::VectorXcd b(size1*size2*size3);
	for ( int i=0; i<size1; i++ ) {
		for ( int j=0; j<size2; j++ ) {
			for ( int k=0; k<size3; k++ ) {
				b(i*size2*size3 + j*size3 + k) = a(i,j,k);
			}
		}
	}
	return b;
}

Eigen::Tensor<std::complex<double>,3> unpackXcd(Eigen::VectorXcd b,
		int size1, int size2, int size3) {
	Eigen::Tensor<std::complex<double>,3> a(size1,size2,size3);
	for ( int i=0; i<size1; i++ ) {
		for ( int j=0; j<size2; j++ ) {
			for ( int k=0; k<size3; k++ ) {
				a(i,j,k) = b(i*size2*size3 + j*size3 + k);
			}
		}
	}
	return a;
}

BaseBandStructure::BaseBandStructure(int numBands_, FullPoints* fullPoints_,
		IrreduciblePoints* irreduciblePoints_) {
	numBands = numBands_;

	if ( fullPoints != nullptr ) {
		useIrreducible = false;
		fullPoints = fullPoints_;
	} else if ( irreduciblePoints != nullptr ) {
		useIrreducible = true;
		irreduciblePoints = irreduciblePoints_;
	} else {
		Error e("BaseBandStructure must provide one Points mesh.", 1);
	}
	if ( ( fullPoints != nullptr ) && ( irreduciblePoints != nullptr ) ) {
		Error e("BaseBandStructure must provide only one Points mesh.", 1);
	}

	Eigen::MatrixXd energies_(getNumPoints(), numBands);
	Eigen::MatrixXcd velocities_(getNumPoints(), numBands*numBands*3);
	energies.setZero();
	velocities_.setZero();
	energies = energies_;
	velocities = velocities_;

	Eigen::MatrixXd dndt_(getNumPoints(), numBands);
	Eigen::MatrixXd dnde_(getNumPoints(), numBands);
	dndt_.setZero();
	dnde_.setZero();
	dndt = dndt_;
	dnde = dnde_;

}

int BaseBandStructure::getNumPoints() {
	if ( useIrreducible ) {
		return irreduciblePoints->getNumPoints();
	} else {
		return fullPoints->getNumPoints();
	}
}

int BaseBandStructure::getIndex(Eigen::Vector3d& pointCoords) {
	if ( useIrreducible ) {
		return irreduciblePoints->getIndex(pointCoords);
	} else {
		return fullPoints->getIndex(pointCoords);
	}
}

double BaseBandStructure::getWeight(int pointIndex) {
	if ( useIrreducible ) {
		return irreduciblePoints->getWeight(pointIndex);
	} else {
		return fullPoints->getWeight(pointIndex);
	}
}

Eigen::Vector3d BaseBandStructure::getPoint(const int pointIndex) {
	if ( useIrreducible ) {
		return irreduciblePoints->getPoint(pointIndex);
	} else {
		return fullPoints->getPoint(pointIndex);
	}
}

void BaseBandStructure::populate() {
	Error e("populate() not implemented in BaseBandStructure", 1);
}

void BaseBandStructure::setEnergies(Eigen::Vector3d& pointCoords,
		Eigen::VectorXd& energies_) {
	int ik = getIndex(pointCoords);
	energies.row(ik) = energies_;
}

void BaseBandStructure::setVelocities(Eigen::Vector3d& pointCoords,
		Eigen::Tensor<std::complex<double>,3>& velocities_) {
	Eigen::VectorXcd tmpVelocities_(numBands*numBands*3);
	tmpVelocities_ = packXcd(velocities_, numBands, numBands, 3);
	int ik = getIndex(pointCoords);
	velocities.row(ik) = tmpVelocities_;
}

//TODO: maybe combine set chemPot with set Temp?
void BaseBandStructure::setChemicalPotential(double chemPot) {
	chemicalPotential = chemPot;
	BaseBandStructure::setOccupations();
}

void BaseBandStructure::setTemperature(double temp) {
	temperature = temp;
	BaseBandStructure::setOccupations();
}

State BaseBandStructure::getStateFromPointIndex(int pointIndex) {
	Error e("Base getStateFromPointIndex is not implemented", 1);
}

void BaseBandStructure::setOccupations() {
	if ( temperature == 0. ) {
		Error e("setOccupations error: must set temperature first" ,1);
	}
	double arg, x, y;
	if ( statistics=="bose" ) {
		for ( int i=0; i<energies.rows(); i++ ) {
			y = energies(i) - chemicalPotential;
			arg = y / 2. / temperature;
			x = sinh(arg);
			dndt(i) = y / 4. / temperature / temperature / x / x;
			dnde(i) = - y / 4. / temperature / x / x;
		}
	}
	else if ( statistics=="fermi" ) {
		for ( int i=0; i<energies.rows(); i++ ) {
			y = energies(i) - chemicalPotential;
			arg = y / 2. / temperature;
			x = cosh(arg);
			dndt(i) = y / 4. / temperature / temperature / x / x;
			dnde(i) = - y / 4. / temperature / x / x;
		}
	}
}

PhBandStructure::PhBandStructure(int numBands_,
		FullPoints* fullPoints_,
		IrreduciblePoints* irreduciblePoints_) : BaseBandStructure(
				numBands_, fullPoints_, irreduciblePoints_) {
	numAtoms = numBands_ / 3;
	Eigen::MatrixXcd eigenvectors_(getNumPoints(),3*numAtoms*numBands);
	eigenvectors_.setZero();
	eigenvectors = eigenvectors_;
}

void PhBandStructure::setEigenvectors(Eigen::Vector3d& pointCoords,
		Eigen::Tensor<std::complex<double>,3>& eigenvectors_) {
	int ik = getIndex(pointCoords);
	Eigen::VectorXcd tmp = packXcd(eigenvectors_, 3, numAtoms, numBands);
	eigenvectors.row(ik) = tmp;
}

PhState PhBandStructure::getStateFromPointIndex(int pointIndex) {
	Eigen::Tensor<std::complex<double>,3> thisVel = unpackXcd(
			velocities.row(pointIndex), numBands, numBands, 3);
	Eigen::Tensor<std::complex<double>,3> thisEig = unpackXcd(
			eigenvectors.row(pointIndex), 3, numAtoms, numBands);
	Eigen::VectorXd thisDndt = dndt.row(pointIndex);
	Eigen::VectorXd thisDnde = dnde.row(pointIndex);
	Eigen::Vector3d p;
	p = getPoint(pointIndex);
	double weight = getWeight(pointIndex);
	Eigen::VectorXd thisEn = energies.row(pointIndex);
	PhState s(p, weight, thisEn, thisEig,
			thisVel, thisDndt, thisDnde);
	return s;
}

ElState ElBandStructure::getStateFromPointIndex(int pointIndex) {
	Eigen::VectorXd thisEn = energies.row(pointIndex);
	thisEn = thisEn.array() - chemicalPotential;
	Eigen::Tensor<std::complex<double>,3> thisVel = unpackXcd(
			velocities.row(pointIndex), numBands, numBands, 3);
	Eigen::VectorXd thisDndt = dndt.row(pointIndex);
	Eigen::VectorXd thisDnde = dnde.row(pointIndex);
	Eigen::Vector3d p = getPoint(pointIndex);
	double weight = getWeight(pointIndex);
	ElState s(p, weight, thisEn, thisVel, thisDndt, thisDnde);
	return s;
}

void PhBandStructure::populate(PhononH0 phononH0) {
	std::vector<Eigen::Vector3d> points;
	if ( useIrreducible ) {
		points = irreduciblePoints->getPoints();
	} else {
		points = fullPoints->getPoints();
	}

	for ( auto q : points ) {
		auto [energies, eigenvectors] = phononH0.diagonalize(q);
		setEnergies(q, energies);
		setEigenvectors(q, eigenvectors);
	}
}

