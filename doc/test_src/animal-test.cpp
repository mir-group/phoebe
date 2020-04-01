/**
 * This is the test program for animal class
 */

#include "animal.hpp"
#include <cassert>
#include <iostream>


void test (unsigned bad_threshold, unsigned maturity_age)
{
    Penna::Animal::set_bad_threshold(bad_threshold);
    Penna::Animal::set_maturity_age(maturity_age);
    Penna::Genome genome;

    for (std::size_t index=0; index < 1000; ++index)
    {
        genome.mutate();
        Penna::Animal animal(genome);

        while (!animal.is_dead()) {
            assert(genome.count_bad(animal.age()) <= bad_threshold);

            if (animal.is_pregnant()) {
                Penna::Animal child = animal.give_birth();
                assert(child.age()==0);
                assert(!child.is_pregnant());
            }
            animal.grow();
        }
        assert(animal.age()==Penna::Animal::maximum_age+1 ||
               genome.count_bad(animal.age())==bad_threshold+1);

    }
}

int main(int agrc, char** argv)
{
    //seed the random number generator
    srand(42);
    // normal cases
    test(4,1);
    test(3,4);
    test(2,5);
    // boundary cases
    test(0,5);
    test(4,Penna::Animal::maximum_age);
    test(Penna::Animal::maximum_age,6);
    test(Penna::Animal::maximum_age,Penna::Animal::maximum_age);

    std::cout<<"Animal tests passed." << std::endl;
    return 0;
}
