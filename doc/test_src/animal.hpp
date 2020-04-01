/* Programming Techniques for Scientific Simulations, HS 2016
 * Exercise 4.3
 */

#ifndef ANIMAL_HPP
#define ANIMAL_HPP

#include "genome.hpp"


/** @brief Namespace for variables and function describing animals.
*  
* 	This is detailed description of namespace Penna.
*/
namespace Penna
{   

/** @brief Class describing animal having a genome and age. 
* 
* 	This is detailed description of class Animal.
*/
class Animal
{
public:
    static const age_type maximum_age = Genome::number_of_genes;

	/** @brief Function which opens a file.
	* 
	*   @param pathname The name of the file.
	*   @param flags Opening flags.
	*/
    int open(const char pathname,int flags);
    static void set_maturity_age( age_type );
    static void set_probability_to_get_pregnant ( double );
    /// interface to genome bad_threshold
    static void set_bad_threshold( age_type );
    static age_type get_bad_threshold();

    /// Default constructor: Uses all good genome.
    Animal();
    /// Constructor using a given genome.
    Animal( const Genome& );

    /// Returns age larger than maturity age and bad state
    bool is_dead() const;
    /// Returns that the animal is pregnant
    bool is_pregnant() const;
    /// Age of animal
    age_type age() const;

    /// Make the animal grow older by one year.
    void grow();
    /// Create a baby animal inheriting its genome from this except for some random mutations.
    Animal give_birth();

private:
    /// Bad genes threshold
    static age_type bad_threshold_;
    /// Maturity age. Parameter R in Penna's paper.
    static age_type maturity_age_;
    /// Probability to be pregnant
    static double probability_to_get_pregnant_;

    const Genome gen_;
    age_type age_;
    bool pregnant_;
};

} // end namespace Penna

#endif // !defined ANIMAL_HPP
