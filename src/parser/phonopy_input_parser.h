#ifndef PHONOPYPARSER_H
#define PHONOPYPARSER_H

#include <string>
#include "phonon_h0.h"
#include "context.h"

/** Class used to parse the raw data from phonopy 
 */
class PhonopyParser {
public:
    /** parsing of force constants.
     * @param context: the object containing the user input.
     * @return Crystal: crystal is the object describing the crystal structure.
     * @return PhononH0: the object containing the force Constants and the
     * functionality to compute the phonon energies.
     */
    static std::tuple<Crystal, PhononH0> parsePhHarmonic(Context &context);

};

#endif
