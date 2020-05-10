#ifndef UTILS_H
#define UTILS_H

// returns the remainder of the division between two integers.
// Note: c++ defines the % operator such that: (a/b)*b + a%b == a    (for b!=0)
// this works as mod() in fortran or % in python when a and b are positive
// But the behavior is not the same for negative integers
// the function below can be used instead
long mod(long a, long b); //{
//	return ( a%b + b ) % b;
//}

// checks if string ends with a suffix
bool hasSuffix(const std::string & str, const std::string & suffix);

#endif
