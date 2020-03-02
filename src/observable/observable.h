#include <string>
#include <vector>

using namespace std;

class Context {
public:
	std::string units;

	std::vector<double> getObservable(std::vector<double>);
};
