#ifndef _UDDSKETCH_H_
#define _UDDSKETCH_H_

#include <cstddef>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <vector>
#include <limits>

class UDDSketch {
    public:
        UDDSketch(double alpha, int m = 0): alpha(alpha), m(m) {
                this->gamma = (1+alpha)/(1-alpha);
                this->log_gamma = std::log(this->gamma);
                this->inverse_log_gamma = 1.0 / std::log(this->gamma);
                this->initial_alpha = this->alpha;
                this->zero_bucket = 0;
                this->min_addressable_value = std::fmax(exp((std::numeric_limits<long>::min()) * this->log_gamma), std::numeric_limits<double>::min());
        }
        UDDSketch(std::vector<uint8_t> &s);
        void add(double value);
        int remove(double value);
        double get_quantile(double q);
        std::vector<uint8_t> serialize();
        static UDDSketch deserialize(void* bytes, std::size_t size);
        long get_bucket_key(double value);
        long get_bucket_count(long key);
        long get_zero_bucket_count();
        void sketch_dump();
        double get_bucket_value(long bkey);
        long get_sketch_count();
        int get_number_buckets();
        double get_initial_alpha();
        double get_alpha();
        double get_gamma();
        int get_max_number_buckets();
        std::size_t get_max_size_bytes(double min = 0, double max = 0);
        std::size_t get_size_bytes();
        void merge(UDDSketch& other);
    private:
        void  collapse(int n = 1);
        double get_min_addressable_value();

        std::unordered_map<long, long> store;
        double alpha;
        double gamma;
        double log_gamma;
        double inverse_log_gamma;
        double initial_alpha;
        double min_addressable_value;
        int m;
        long zero_bucket;
        
};

#endif // _UDDSKETCH_H_
