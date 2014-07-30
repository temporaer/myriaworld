/**
 *
 * Published under three-clause BSD license.
 *
 * Copyright 2014 Hannes Schulz <schulz@ais.uni-bonn.de>
 * 
 * Implements an almost transparent disk cache for C++ functions.
 * Depends heavily on the C++11 features (auto, decltype, rvalue references,
 * variadic templates)
 * 
 * Another (optional) dependency is boost.log, which is contained 
 * in boost versions * >=1.55.
 * 
 * Preconditions are:
 * - All function arguments must be hashable, side effects are not considered.
 * - Returned object must be serializable
 * - Returned object must be default constructible
 *
 * Usage example below the header code.
 * 
 * Compile using:
 * $ g++ -DBOOST_ALL_DYN_LINK -DCFTEST -std=c++11 -x c++ cached_function.hpp -lboost_system -lboost_filesystem -lboost_serialization -pthread -lboost_log
 * and run with
 * $ ./a.out 42
 */
#ifndef __CACHED_FUNCTION_HPP__
#     define __CACHED_FUNCTION_HPP__
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>

#define CACHED(cache, func, ...) cache(#func, func, __VA_ARGS__)

namespace fscache{
    namespace fs = boost::filesystem;
    namespace detail{
        template <typename T>
            size_t hash_combine(std::size_t seed, const T& t) {
                boost::hash_combine(seed, t);
                return seed;
            }
        template <typename T, typename... Params>
            size_t hash_combine(std::size_t seed, const T& t, const Params&... params) {
                boost::hash_combine(seed, t);
                return hash_combine(seed, params...);
            }
    }
    struct function_cache{
        fs::path m_path;
        function_cache(std::string path = fs::current_path().string())
        :m_path(fs::path(path) / "cache"){
            fs::create_directories(m_path);
        }

        template<typename Func, typename... Params>
            auto operator()(Func f, Params&&... params) -> decltype(f(params...)) {
                return (*this)("anonymous", f, std::forward<Params>(params)...);
            }
        template<typename Func, typename... Params>
            auto operator()(std::string descr, Func f, Params&&... params) -> decltype(f(params...)) {
                std::size_t seed = detail::hash_combine(0, descr, params...);
                return (*this)(descr, seed, f, std::forward<Params>(params)...);
            }
        template<typename Func, typename... Params>
            auto operator()(std::string descr, std::size_t seed, Func f, Params&&... params) -> decltype(f(params...)) {
                typedef decltype(f(params...)) retval_t;
                std::string fn = descr + "-" + boost::lexical_cast<std::string>(seed);
                fn = (m_path / fn).string();
                if(fs::exists(fn)){
                    std::ifstream ifs(fn);
                    boost::archive::binary_iarchive ia(ifs);
                    retval_t ret;
                    ia >> ret;
                    BOOST_LOG_TRIVIAL(info) << "Cached access from file "<<fn;
                    return ret;
                }
                retval_t ret = f(std::forward<Params>(params)...);
                BOOST_LOG_TRIVIAL(info) << "Non-cached access, file "<<fn;
                std::ofstream ofs(fn);
                boost::archive::binary_oarchive oa(ofs);
                oa << ret;
                return ret;
            }
    };

}
#endif /* __CACHED_FUNCTION_HPP__ */


#ifdef CFTEST
#include <iostream>
#include <boost/serialization/vector.hpp>

long fib(long i){
    if(i==0)
        return 0;
    if(i==1)
        return 1;
    return fib(i-1) + fib(i-2);
}

std::vector<int> times(const std::vector<int>& v, int factor){
    std::vector<int> v2 = v;
    for(int& i : v2)
        i *= factor;
    return v2;
}

int
main(int argc, char **argv)
{
    fscache::function_cache c;
    if(argc != 2){
        std::cout << "Usage: " << argv[0] << " N" << std::endl;
        std::cout << " where N is the index of the fibonacci number to compute" << std::endl;
        exit(1);
    }
    // name cache same as function name
    std::cout << CACHED(c, fib, atoi(argv[1])) << std::endl;
    std::cout << CACHED(c, fib, atoi(argv[1])) << std::endl;

    // most convenient
    std::cout << c(fib, atoi(argv[1])+1) << std::endl;
    std::cout << c(fib, atoi(argv[1])+1) << std::endl;
 
    // unhashable arg, use own seed
    std::cout << c("fib", 28725, fib, atoi(argv[1])+2) << std::endl;
    std::cout << c("fib", 28725, fib, atoi(argv[1])+2) << std::endl;

    std::vector<int> v(10000, atoi(argv[1])), v2, v3;
    v2 = CACHED(c, times, v, 5);
    v3 = CACHED(c, times, v, 5);
    assert(v2 == v3);

    return 0;
}
#endif
