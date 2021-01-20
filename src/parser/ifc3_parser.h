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

private:
  static Eigen::MatrixXd wsInit(Crystal &crystal, Eigen::Vector3i qCoarseGrid);
  Eigen::MatrixXd rws;
};
