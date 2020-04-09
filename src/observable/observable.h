#include <string>
#include <vector>

class Observable {
public:
	std::string units;
	std::vector<double> getObservable(std::vector<double>);
};
