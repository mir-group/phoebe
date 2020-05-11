
// include statements
// TODO: if we don't compile with MPI, we can't import it, so maybe this should be in an ifdef
#include <vector>
#include <complex>
#include <mpi.h>

#ifndef MPICONTROLLER_H
#define MPICONTROLLER_H

class MPIcontroller{
    
    //MPI_Comm com; // MPI communicator -- might not need this, can just call default comm "MPI_COMM_WORLD"
    int size; // number of MPI processses
    int rank;
    double startTime; // the time for the entire mpi operation
    double lastTime;
    
    public:
        // MPIcontroller class constructors -----------------------------------
        /** a constructor which sets up the MPI environment, initializes the communicator, and starts a timer **/
        MPIcontroller();
        ~MPIcontroller(){ if (!MPI::Is_finalized()) finalize(); }
        
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
        // TODO: upgrade this function so that it also takes a string. Then, append the string to some buffer stored in this object
        //     fprintf(stderr, "Error from rank %3d: %s\n", rank, errString);
        // or at least, something that prints a time statistic for a given block
        void time() const; // returns time elapsed since mpi started
    
        // IO functions
        // do we want some parallel hdf5 options? how does that work?
        // wgropp.cs.illinois.edu/courses/cs598-s16/lectures/lecture32.pdf  -- look at slide 10 here?
        // need to think about how to do this properly. Also, are we doing parallel hdf5?
        //void mpiWrite();
        //void mpiRead();
            
};

// we need to use the concept of a "type traits" object to serialize the standard cpp types
// // and then we can define specific implementations of this for types which are not standard
// // https://stackoverflow.com/questions/42490331/generic-mpi-code
namespace mpi { 
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
        // TODO: A container for an array

        // TODO: Define here a container for a matrix
}

// default constructor
MPIcontroller::MPIcontroller(){
    
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
}

// TODO: any other stats would like to output here? 
void MPIcontroller::finalize() const {
    fprintf(stdout, "Final time for rank %3d: %3f\n ", rank, MPI_Wtime() - startTime );
    MPI_Finalize();
}

// Collective communications functions -----------------------------------
template<typename T> void MPIcontroller::bcast(T* dataIn) const{

    using namespace mpi;     
    if(size==1) return; 
    int errCode;

    errCode = MPI_Bcast( containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn), containerType<T>::getMPItype(), 0, MPI_COMM_WORLD);
    if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
    
}
template<typename T> void MPIcontroller::reduceSum(T* dataIn, T* dataOut) const{

    using namespace mpi;     
    if(size==1) return; 
    int errCode;

    errCode = MPI_Reduce(containerType<T>::getAddress(dataIn),containerType<T>::getAddress(dataOut),containerType<T>::getSize(dataIn),containerType<T>::getMPItype(), MPI_SUM, 0, MPI_COMM_WORLD);
    if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
}
template<typename T> void MPIcontroller::reduceMax(T* dataIn, T* dataOut) const{
    
    using namespace mpi; 
    if(size==1) return; 
    int errCode;

    errCode = MPI_Reduce(containerType<T>::getAddress(dataIn),containerType<T>::getAddress(dataOut),containerType<T>::getSize(dataIn),containerType<T>::getMPItype(), MPI_MAX, 0, MPI_COMM_WORLD);
    if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
}
template<typename T> void MPIcontroller::reduceMin(T* dataIn, T* dataOut) const{

    using namespace mpi;     
    if(size==1) return; 
    int errCode;

    errCode = MPI_Reduce(containerType<T>::getAddress(dataIn),containerType<T>::getAddress(dataOut),containerType<T>::getSize(dataIn),containerType<T>::getMPItype(), MPI_MIN, 0, MPI_COMM_WORLD);
    if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
}

// Asynchronous support functions -----------------------------------
void MPIcontroller::barrier(){
    
    if(size==1) return; 
    int errCode;
    
    errCode = MPI_Barrier(MPI_COMM_WORLD);
    if(errCode != MPI_SUCCESS) {  errorReport(errCode); } 
}

// Utility functions  -----------------------------------

// get the error string and print it to stderr before returning
void MPIcontroller::errorReport(int errCode) const{
    char errString[BUFSIZ];
    int lengthOfString;
    
    MPI_Error_string(errCode, errString, &lengthOfString);
    fprintf(stderr, "Error from rank %3d: %s\n", rank, errString);
    MPI_Abort(MPI_COMM_WORLD, errCode);
}

void MPIcontroller::time() const{
    fprintf(stdout, "Time for rank %3d : %3f\n", rank, MPI_Wtime() - lastTime );
}


#endif
