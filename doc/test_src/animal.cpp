#include "animal.hpp"
#include <cassert>
#include <cstdlib> // drand48()



namespace Penna
{


// Definition of static data members
age_type Animal::bad_threshold_;
age_type Animal::maturity_age_;
double Animal::probability_to_get_pregnant_;

void Animal::set_maturity_age( age_type r )
{
    maturity_age_ = r;
}

void Animal::set_probability_to_get_pregnant( double p )
{
    probability_to_get_pregnant_ = p;
}

void Animal::set_bad_threshold( age_type t )
{
  bad_threshold_ = t;
}

age_type Animal::get_bad_threshold()
{
  return bad_threshold_;
}

/// Default constructor: Uses all good genome.
Animal::Animal()
:   gen_()
,   age_(0)
{
}

/// Constructor using a given genome.
Animal::Animal( const Genome& gen )
:   gen_(gen)
,   age_(0)
{
}


bool Animal::is_dead() const
{
    return age_ > maximum_age || gen_.count_bad(age_) > bad_threshold_;
}

bool Animal::is_pregnant() const
{
    return pregnant_;
}

age_type Animal::age() const
{
    return age_;
}

void Animal::grow()
{
    assert( !is_dead() );
    ++age_;
    if (age_ > maturity_age_ && !pregnant_)
    {
      if (probability_to_get_pregnant_>drand48())
        pregnant_=true;
    }
}

Animal Animal::give_birth()
{
    assert( pregnant_ == true );
    pregnant_=false;
    Genome child_genome=gen_;
    child_genome.mutate();
    return Animal( child_genome );
}

} // end namespace Penna
