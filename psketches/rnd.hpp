#ifndef RND_HPP_
#define RND_HPP_
#include <algorithm>
#include <vector>
#include <random>
#include <memory>

template<typename T>
class rnd_gen {
public:
    rnd_gen(const typename T::param_type &params, unsigned int seed = 0);

    std::shared_ptr<std::vector<typename T::result_type>> get_stream(long n);
    typename T::result_type get_value();

private:
    T pdis;
    std::mt19937 gen;
};

#include "rnd_impl.hpp"

#endif
