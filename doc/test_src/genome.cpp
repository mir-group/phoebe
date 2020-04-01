#include "genome.hpp"
#include <cstdlib> // drand48()

namespace Penna
{

//declaration of the static variables
age_type Genome::mutation_rate_;

void Genome::set_mutation_rate( age_type m )
{
    mutation_rate_ = m; 
}

age_type Genome::count_bad( age_type n ) const
{
    return (genes_<<(number_of_genes-n-1)).count();
}

void Genome::mutate()
{
    // Mutate a random selection of M genes
    for( size_t i = 0; i < mutation_rate_; ++i )
        genes_.flip( int( drand48()*number_of_genes ) );
}

} // end namespace Penna
