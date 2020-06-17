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
};
