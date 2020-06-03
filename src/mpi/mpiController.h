#include <vector>
#include <complex>
#include <chrono>
#ifdef MPI_AVAIL
#include <mpi.h>
#endif

#ifndef MPICONTROLLER_H
#define MPICONTROLLER_H

class MPIcontroller{
	
	//MPI_Comm com; // MPI communicator -- might not need this, can just call default comm "MPI_COMM_WORLD"
	int size; // number of MPI processses
	int rank;

	#ifdef MPI_AVAIL 
	double startTime; // the time for the entire mpi operation
	#else 
	std::chrono::steady_clock::time_point startTime; 
	#endif    

	public:
		// MPIcontroller class constructors -----------------------------------
		/** a constructor which sets up the MPI environment, initializes the communicator, and starts a timer **/
		MPIcontroller();
		~MPIcontroller(); //{ if (!MPI::Is_finalized()) finalize(); }
		
		// Calls finalize and potentially reports statistics, time or handles errors
		void finalize() const; 
	
		// Collective communications functions -----------------------------------
		template<typename T> void bcast(T* dataIn) const;
		template<typename T> void reduceSum(T* dataIn, T* dataOut) const;
		template<typename T> void reduceMax(T* dataIn, T* dataOut) const;
		template<typename T> void reduceMin(T* dataIn, T* dataOut) const;

		// point to point functions -----------------------------------
		//template<typename T> void send(T&& data) const;
		//template<typename T> void recv(T&& data) const;
	
		// Asynchronous functions
		void barrier();
	
		// Utility functions -----------------------------------
		bool mpiHead() const{ return rank==0; } // a function to tell us if this process is the head
		int getRank() const { return rank; }
		int getSize() const { return size; }    

		// Error reporting and statistics
		void errorReport(int errCode) const; // collect errors from processes and reports them, then kills the code
		void time() const; // returns time elapsed since mpi started
	
		// IO functions
		// TODO: implement these functions, if we need them. 
		void mpiWrite();
		void mpiRead();
		void mpiAppend(); 
			
};

// we need to use the concept of a "type traits" object to serialize the standard cpp types
// // and then we can define specific implementations of this for types which are not standard
// // https://stackoverflow.com/questions/42490331/generic-mpi-code
namespace mpi { 
		#ifdef MPI_AVAIL
		// Forward declaration for a basic container type 
		template <typename ...> struct containerType;

		// Define a macro to shorthand define this container for scalar types
		// Size for basic scalar types will always be 1
		#define MPIDataType(cppType,mpiType) \
		template<> struct containerType<cppType> { \
				static inline cppType* getAddress(cppType* data) { return data; } \
				static inline size_t getSize(cppType* data) { return 1; } \
				static inline MPI_Datatype getMPItype() { return mpiType; } \
		};

		// Use definition to generate containers for scalar types
		MPIDataType(int, MPI_INT)
		MPIDataType(unsigned int, MPI_UNSIGNED)
		MPIDataType(float, MPI_FLOAT)
		MPIDataType(double, MPI_DOUBLE)

		MPIDataType(std::complex<double>, MPI_DOUBLE)
		MPIDataType(std::complex<float>, MPI_FLOAT)
		#undef MPIDataType 

		// A container for a std::vector
		template<typename T> struct containerType<std::vector<T>>{
				static inline T* getAddress(std::vector<T>* data) { return data->data(); }
				static inline size_t getSize(std::vector<T>* data) { return data->size(); }
				static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
		};
		// TODO: Define any other useful containers (a matrix? a standard array?)
  
		#endif
}

// default constructor
MPIcontroller::MPIcontroller(){
	
	#ifdef MPI_AVAIL
	// start the MPI environment
	MPI_Init(NULL, NULL);
	
	// get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// get rank of current process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// start a timer
	startTime = MPI_Wtime();
	
	// set this so that MPI returns errors and lets us handle them, rather
	// than using the default, MPI_ERRORS_ARE_FATAL
	MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	#else
 
	startTime = std::chrono::steady_clock::now();
	#endif
}

// default destructor
MPIcontroller::~MPIcontroller(){
	#ifdef MPI_AVAIL
	if (!MPI::Is_finalized()) finalize();
	#endif
}

// TODO: any other stats would like to output here? 
void MPIcontroller::finalize() const {
	#ifdef MPI_AVAIL
	fprintf(stdout, "Final time for rank %3d: %3f\n ", rank, MPI_Wtime() - startTime );
	MPI_Finalize();
	#else
	fprintf(stdout, "Final time for rank %3d: %3f\n ", 0,  - startTime );   
	#endif
}

// Collective communications functions -----------------------------------
template<typename T> void MPIcontroller::bcast(T* dataIn) const{
	using namespace mpi;     
	#ifdef MPI_AVAIL
	if(size==1) return; 
	int errCode;

	errCode = MPI_Bcast( containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn), containerType<T>::getMPItype(), 0, MPI_COMM_WORLD);
	if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
	#endif
}
template<typename T> void MPIcontroller::reduceSum(T* dataIn, T* dataOut) const{
	using namespace mpi;     
	#ifdef MPI_AVAIL
	if(size==1) return; 
	int errCode;

	errCode = MPI_Reduce(containerType<T>::getAddress(dataIn),containerType<T>::getAddress(dataOut),containerType<T>::getSize(dataIn),containerType<T>::getMPItype(), MPI_SUM, 0, MPI_COMM_WORLD);
	if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
	#endif
}
template<typename T> void MPIcontroller::reduceMax(T* dataIn, T* dataOut) const{
	using namespace mpi; 
	#ifdef MPI_AVAIL
	if(size==1) return; 
	int errCode;

	errCode = MPI_Reduce(containerType<T>::getAddress(dataIn),containerType<T>::getAddress(dataOut),containerType<T>::getSize(dataIn),containerType<T>::getMPItype(), MPI_MAX, 0, MPI_COMM_WORLD);
	if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
	#endif
}
template<typename T> void MPIcontroller::reduceMin(T* dataIn, T* dataOut) const{
	using namespace mpi;     
	#ifdef MPI_AVAIL
	if(size==1) return; 
	int errCode;

	errCode = MPI_Reduce(containerType<T>::getAddress(dataIn),containerType<T>::getAddress(dataOut),containerType<T>::getSize(dataIn),containerType<T>::getMPItype(), MPI_MIN, 0, MPI_COMM_WORLD);
	if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
	#endif
}

// Asynchronous support functions -----------------------------------
void MPIcontroller::barrier(){
	#ifdef MPI_AVAIL
	if(size==1) return; 
	int errCode;
	
	errCode = MPI_Barrier(MPI_COMM_WORLD);
	if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
	#endif
}

// Utility functions  -----------------------------------

// get the error string and print it to stderr before returning
void MPIcontroller::errorReport(int errCode) const{
	#ifdef MPI_AVAIL
	char errString[BUFSIZ];
	int lengthOfString;
	
	MPI_Error_string(errCode, errString, &lengthOfString);
	fprintf(stderr, "Error from rank %3d: %s\n", rank, errString);
	MPI_Abort(MPI_COMM_WORLD, errCode);
	#else 
	// TODO: how are we throwing non-mpi errors? 
	#endif
}

void MPIcontroller::time() const{
	#ifdef MPI_AVAIL
	fprintf(stdout, "Time for rank %3d : %3f\n", rank, MPI_Wtime() - startTime );
	#else
	std::cout << "Time for rank 0 :" << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - seconds).count() << " secs" << std::endl;
	#endif
}


#endif
