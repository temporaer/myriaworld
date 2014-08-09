#include <iostream>
#include "cached_function.hpp"

long fib(long i){
    if(i==0)
        return 0;
    if(i==1)
        return 1;
    return fib(i-1) + fib(i-2);
}

int
main(int argc, char **argv)
{
    fscache::function_cache c;
    std::cout << c(fib, atoi(argv[1])) << std::endl;
    std::cout << c(fib, atoi(argv[1])) << std::endl;
    return 0;
}
