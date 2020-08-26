#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <assert.h>
#include <string>
#include <iostream>

/** Object used to print an error message, and stop the code.
 */
class Error {
public:
    /** object constructor.
     * @param errorMessage: message to be print to the user.
     * @param errCode: return integer error code, should be different from 0!
     */
    Error(const std::string &errMessage, const int &errCode = 1);
};

/** Object to print a warning to the user, without stopping the code
 */
class Warning {
public:
    Warning(const std::string &errMessage);
};

#endif
