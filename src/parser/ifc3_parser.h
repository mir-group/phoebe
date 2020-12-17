#include "context.h"
#include "crystal.h"
#include "interaction_3ph.h"

/** Class for the parsing of the 3rd derivative of the total energy with
 * respect to ionic displacements, and the creation of a dedicated object.
 */
class IFC3Parser {
public:
  /** Interface to parse the 3rd derivative files.
   * Will open the file and decide whether it's a ShengBTE or QuantumESPRESSO
   * format, and call the dedicated method.
   * This is the method that should be called in the rest of the code.
   * @param crystal: Crystal object. To be obtained parsing the dynamical
   * matrix first.
   * @param context: the Context object with the user input (to get filename)
   * @return Interaction3Ph: the object that contains the 3rd derivative.
   */
  static Interaction3Ph parse(Context &context, Crystal &crystal);

  /** Parse the 3rd derivative of the total energy wrt ionic displacements.
   * using the ShengBTE format.
   * @param crystal: Crystal object. To be obtained parsing the dynamical
   * matrix first.
   * @param context: the Context object with the user input (to get filename)
   * @return Interaction3Ph: the object that contains the 3rd derivative.
   */
  static Interaction3Ph parseFromShengBTE(Context &context, Crystal &crystal);

  /** Parse the 3rd derivative of the total energy wrt ionic displacements.
   * using the Phono3py format.
   * @param context: the Context object with the user input (to get filename)
   * @return Interaction3Ph: the object that contains the 3rd derivative.
   */
  static Interaction3Ph parseFromPhono3py(Context &context, Crystal &crystal);

  /** In development!
   * TODO: while the parsing is done, and it's easy, QE uses a different
   * convention. Let D3(r1,r2,r3) be the 3rd derivative matrix (cartesian and
   * atomic indices excluded for simplicity). ShengBTE defines D3 such that
   * r1 is always in the origin. QE instead assumes that r3 is in the origin.
   * Thus, one should either (1) modify interaction3Ph to put the suitable
   * phases, or (2) Fourier transform D3(r1,r2,r3) -> D3(q1,q2,q3) and then
   * transform back setting r1=0.
   * @param crystal: Crystal object. To be obtained parsing the dynamical
   * matrix first.
   * @param context: the Context object with the user input (to get filename)
   * @return Interaction3Ph: the object that contains the 3rd derivative.
   */
  static Interaction3Ph parseFromQE(Context &context, Crystal &crystal);

private:
  /** Find the index of a 3d vector in a matrix of 3d vectors.
   * Calls an error if the vector is not found
   *
   * @param cellPositions2: the matrix containing the list of vectors to be
   * searched, of size (numVectors,3)
   * @param position2: the target of the search algorithm
   * @return idx: the row index of position2 in cellPositions2
   */
  long findIndexRow(Eigen::MatrixXd &cellPositions2,
                    Eigen::Vector3d &position2);

  /** This function, used for phonon3py, folds a list of Bravais vectors into
   * the Wigner Seitz zone of a supercell.
   *
   * @param supBravaisVectors: list of Bravais lattice vectors, in cartesian
   * coordinates, which is assumed to be in the form i*a1+j*a2+k*a3, with
   * a1,a2,a3 the primitive lattice vectors, and i,j,k NON-NEGATIVE integers
   * @param crystal: Crystal object describing the primitive unit cell
   * @param grid: identifies the size of the supercell, compared to the
   * primitive unit cell
   * @return a tuple with: a matrix of size (numWSVectors,3) containing the list
   * of Bravais lattice vectors folded in the Wigner Seitz zone, containing also
   * reducible copies; a map as a Eigen::VectorXi such that map(wsIndex) maps
   * to the equivalent lattice vector in the input list.
   */
  std::tuple<Eigen::MatrixXd, Eigen::VectorXi>
  buildWignerSeitz(const Eigen::MatrixXd &supBravaisVectors, Crystal &crystal,
                   Eigen::Vector3i grid);
};
