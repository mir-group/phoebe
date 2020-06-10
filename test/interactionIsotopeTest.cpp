#include "gtest/gtest.h"
#include "points.h"
#include "state.h"
#include "qe_input_parser.h"
#include "ph_scattering.h"
#include "bandstructure.h"
#include <fstream>

TEST (InteractionIsotope,Wphisoiq4) {
  // Here we compare the phonon-isotope scattering rates 
  // against those from ShengBTE. The test material is cubic
  // silicon and the phonon q-mesh is 8x8x8.

  //Parse espresso ifc2 file
  Context context;
  context.setPhD2FileName("../test/interaction3ph/QEspresso.fc");
  QEParser qeParser;
  auto [crystal,phononH0] = qeParser.parsePhHarmonic(context);
  phononH0.setAcousticSumRule("simple");

  //Total number of wave vectors
  int nq = 8*8*8;
  //For full BZ mesh in crystal coordinates
  Eigen::Vector3i qMesh;
  qMesh << 8, 8, 8;
  FullPoints points(crystal, qMesh);
  
  //Number of atoms
  int numAtoms = crystal.getNumAtoms();
  //Number of phonon branches
  int numBranches = 3*numAtoms;
  //Mass of basis atoms
  auto atomicMasses = crystal.getAtomicMasses();
  //Mass variance for isotopic scattering
  //TODO should get this through a getMassVariance method
  Eigen::VectorXd massVariance(numAtoms);
  massVariance << 0.00020164270215294942, 0.00020164270215294942; //from src/constants/periodic_table.h
  
  //Phonon energies initialized to zero
  Eigen::VectorXd energies(numBranches);
  energies.setConstant(0.0); 

  // Test q-point
  int iq = 4; // = (0.5, 0, 0) in crystal coordinates
  Eigen::VectorXd WIsotopeiq(numBranches);
  WIsotopeiq.setConstant(0.0);

  //Reference value from ShengBTE with 0.01 Ry fixed Gaussian broadening
  Eigen::VectorXd WRef(numBranches);
  WRef << 0.3327465318e-03, 0.3327465318e-03, 0.4151744806e-02, 0.4978611835e-02, 0.6913784969e-02,  0.6913784969e-02;
  
  context.setSmearingWidth(0.01); // in rydbergs
  GaussianDeltaFunction smearing(context);
  
  //Eigenvector and angular frequencies at iq
  auto ip = points.getPoint(iq);
  //cout << ip.getCoords(Points::crystalCoords);
  auto [omegasiq,eviq] = phononH0.diagonalize(ip);
  //auto vsiq = phononH0.diagonalizeVelocity(ip);
  
  for(int jq = 0; jq < nq; jq++){ //integrate over
    auto jp = points.getPoint(jq);
    //Eigenvector and angular frequencies at jq
    auto [omegasjq,evjq] = phononH0.diagonalize(jp);
    //auto vsjq = phononH0.diagonalizeVelocity(jp);
    
    for(int ib = 0; ib < numBranches; ib++){ //phonon branches of the test q-point
      double fac = pow(omegasiq(ib),2)/nq;
      for(int jb = 0; jb < numBranches; jb++){ //integrate over
	//Using fixed gaussian for now
	double deltaWeight = smearing.getSmearing(omegasiq(ib) - omegasjq(jb));

	if(deltaWeight == 0) continue;

	for(int p = 0; p < numAtoms; p++){ //sum over atomic basis
	  double aux = 0.0;
	  for(int kdim = 0; kdim < 3; kdim++){ //inner product over Cartesian space
	    //Recall that the eigenvectors were mass-normalized
	    aux += pow(abs(std::conj(eviq(kdim,p,ib))*evjq(kdim,p,jb)),2)*pow(atomicMasses(p),2);
	  }
	  WIsotopeiq(ib) += aux*deltaWeight*massVariance(p)*fac; //Ry
	  //
	}//p
      }//jb
    }//ib
  }//jq

  cout << "branch " << "ShengBTE    " << "Phoebe    " << "\n";
  for(int ib = 0; ib < numBranches; ib++){
    //Convert from Ry to THz
    WIsotopeiq(ib) *= 32889.83;
    cout << ib << "    " << WRef(ib) << " " << WIsotopeiq(ib) << " THz" << "\n";
  }
  
}

