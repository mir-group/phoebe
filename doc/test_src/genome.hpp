/* Programming Techniques for Scientific Simulations, HS 2016
 * Exercise 4.3
 */

#ifndef GENOME_HPP
#define GENOME_HPP

#include <bitset>
#include <limits>

namespace Penna
{

typedef unsigned age_type;

/**
 * Genome class. Each gene is represented by one bit.
 */
class Genome
{
public:
    /// Up to this size bitset is a lot faster
    static const age_type number_of_genes = 
        std::numeric_limits<unsigned long>::digits;
    
    static void set_mutation_rate( age_type );
    static void set_bad_threshold( age_type );
    
    /// Default constructor: Initialize genes to all good.
    /// Genome() {};  // provided by the compiler
    
    /// Count number of bad genes in first n years.
    age_type count_bad( age_type ) const;
    /// mutate the genome by flipping of M genes
    void mutate();

private:
    /// Parameter M in Penna's paper
    static age_type mutation_rate_;

    std::bitset<number_of_genes> genes_;
};

} // end namespace Penna

#endif // !defined GENOME_HPP
