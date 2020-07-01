
#ifndef EPA_PARSER_H
#define EPA_PARSER_H

#include "context.h"
#include "interaction_epa.h"
//#include "interaction_epa.h"

/** Class for parsing the squared electron-phonon matrix elements
 * \f$g_{\nu}^2(\epsilon_i,\epsilon_j)\$f averaged on the electron energy
 * grid using moving-least-squares strategy as well as averaged phonon frequencies
 *\f$ \omega_{\nu} \$f for all phonon modes.
 * After the parsing, the dedicated object is created. At the moment, calculation of \f$g^2\$f
 * and \f$ \omega_{\nu} \$f is  performed using epa2mls.f90 code from D. Wee and G. Samsonidze;
 * Later this class will be extended or new class will be created to read the electron-
 * phonon matrix elements from QE and perform this averaging inside Phoebe.
 * **IMPORTANT**: currently valence and conduction bands are treated separately in EPA
 * which likely will not work for metals. We will probably change this implementation when the
 * averaging procedure will be implemented in phoebe
 */

class EpaParser {
    public:
    
    /** Interface to parse averaged squared electron-phonon matrix elements and phonon
     * frequencies.
     * It parses content of epa.e file.
     * @param[in] context: the Context objext with the user input (to get filename)
     * @return[in] interactionEpa: InteractionEPA object containing the matrix of squared
     * electron-phonon matrix elements and averaged phonon frequencies.
     */
    
    static InteractionEpa parseAvCouplings(Context & context);
};

#endif
