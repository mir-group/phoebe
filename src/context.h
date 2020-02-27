#include <string>
#include <vector>

using namespace std;

class Context {
public:
	std::vector<int> qCoarseMesh = {0,0,0};
	std::vector<int> qFineMesh = {0,0,0};
	std::vector<int> kCoarseMesh = {0,0,0};
	std::vector<int> kFineMesh = {0,0,0};

	void setQCoarseMesh(int* x);
	std::vector<int> getQCoarseMesh();

	void setupFromInput(string fileName);

};
