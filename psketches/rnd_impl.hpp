#ifndef RND_IMPL_HPP_
#define RND_IMPL_HPP_

#include <random>

template <typename T>
rnd_gen<T>::rnd_gen(const typename T::param_type &params, unsigned int seed) {
    if (!seed) {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        seed = rd(); // seeded with rd()
    }
    gen = std::mt19937(seed); // Standard mersenne_twister_engine 
    pdis = T(params);
}

template <typename T>
std::shared_ptr<std::vector<typename T::result_type>> rnd_gen<T>::get_stream(long n) {
    auto buffer = std::make_shared<std::vector<typename T::result_type>>();
    for ( int i = 0; i < n; i++) {
        buffer->push_back(pdis(gen));
    }
    return buffer;
}

template <typename T>
typename T::result_type rnd_gen<T>::get_value() {
    return pdis(gen);
}

#endif