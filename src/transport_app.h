#include <string>

/** Main driver for the transport calculation
 *
 */
class TransportApp {
public:
//	TransportApp();
	/** Setup of the calculations
	 * reads input files, sets up parallel environment
	 * @param argc, argv: the parameters provided on the command line
	 */
	TransportApp(int argc, char** argv);
private:
	std::tuple<std::string, std::string> parseInputArgs(int argc, char* argv[]);
};
