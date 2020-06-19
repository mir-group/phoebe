#ifndef PHONONTRANSPORTAPP_H
#define PHONONTRANSPORTAPP_H

#include <string>
#include "app.h"

/** Main driver for the transport calculation
 *
 */
class PhononTransportApp: public App {
public:
    void run(Context &context);
    void checkRequirements(Context &context);
};

#endif
