/*
 * This is the test program for genome class
 */

#include "genome.hpp"
#include <cassert>
#include <iostream>


void test(int M)
{
    Penna::Genome::set_mutation_rate(M);
    
    Penna::Genome parent_genome;
    Penna::Genome child_genome;
    
    for (std::size_t index=0; index < 1000; ++index)
    {
        child_genome = Penna::Genome(parent_genome);
        child_genome.mutate();
        
        unsigned diff = (child_genome.count_bad(0)==parent_genome.count_bad(0) ? 0 : 1);
        for (std::size_t age=1; age < Penna::Genome::number_of_genes; ++age) 
        {
    	    assert( (child_genome.count_bad(age-1) == child_genome.count_bad(age)) || (child_genome.count_bad(age-1)+1 == child_genome.count_bad(age)));
            unsigned child_gene = child_genome.count_bad(age)-child_genome.count_bad(age-1);
            unsigned parent_gene = parent_genome.count_bad(age)-parent_genome.count_bad(age-1);
            if (parent_gene!=child_gene) ++diff;
        }
        assert ( diff <= M && (diff%2 == M%2));
        parent_genome = child_genome;  
    }
}


int main(int agrc, char** argv)
{
    //seed the random number generator
    srand(42);
    
    //run the test
    for (int M=1; M<5; ++M) 
        test(M);
    std::cout<<"Genome tests passed." << std::endl;
    return 0;
}
