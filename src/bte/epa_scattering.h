#ifndef EPA_SCATTERING_H
#define EPA_SCATTERING_H
#include "context.h"
#include "statistics_sweep.h"
#include "bandstructure.h"
#include "vector_bte.h"
#include "eigen.h"

class EpaScattering {
public:
    /** EPA scattering vector constructor
     * @param context: object containing user input configuration.
     * @param statisticsSweep: object controlling the loops over temperature
     * and chemical potential.
     * @param bandStructure: BandStructure object.
     */
  //  EpaScattering(Context &context_, StatisticsSweep &statisticsSweep_,
        //             FullBandStructure &bandStructure_);
    
    
    EpaScattering(Context & context_, StatisticsSweep & statisticsSweep_, FullBandStructure & fullBandStructure_, Eigen::VectorXd & energies_);
    
    EpaScattering(const EpaScattering &that);
    
    EpaScattering& operator=(const EpaScattering &that);
    
    ~EpaScattering();
    
    static BaseVectorBTE setup(Context & context, StatisticsSweep & statisticsSweep, FullBandStructure & fullBandStructure, Eigen::VectorXd & energies);
    
    /** Copy constructor
     */
    //EpaScattering(const EpaScattering &that);
    
    /** Copy assignment operator
     */
   // EpaScattering& operator=(const EpaScattering &that);
    
    /** This method needs to be called right after the constructor.
     * It will setup variables of the scattering matrix, and compute at least
     * the linewidths.
     * If the user input ask to store the matrix in memory, the scattering
     * matrix gets allocated and built here, through a call to a virtual member
     * builder(), which is defined in subclasses.
     */
    //void setup();
    
    /** Returns the diagonal matrix elements.
     * @return diagonal: a VectorBTE containing the linewidths.
     */
   // BaseVectorBTE rates();
    
    /** Computes the product A*f - diag(A)*f
     * where A is the scattering matrix, f is the vector of quasiparticle
     * populations, and diag(A) is the diagonal of the matrix.
     */
    
private:
    Context &context;
    StatisticsSweep &statisticsSweep;
    FullBandStructure &fullBandStructure;
    Eigen::VectorXd energies;
};

#endif
