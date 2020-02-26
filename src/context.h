#include <string>

using namespace std;

class Context {
public:
	int qCoarseMesh[3];
	int qFineMesh[3];
	int kCoarseMesh[3];
	int kFineMesh[3];

	void setQCoarseMesh(int* x);
	const int* getQCoarseMesh();

	void setupFromInput(string fileName);

};
