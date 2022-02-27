#include "uddsketch.h"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <stdexcept>
#include <map>
#include <iostream>
#include <sstream>




inline float natural_log(float x) {

  // ASSUMING:
  // - non-denormalized numbers i.e. x > 2^−126
  // - integer is 32 bit. float is IEEE 32 bit.

  // INSPIRED BY:
  // - https://stackoverflow.com/a/44232045
  // - http://mathonweb.com/help_ebook/html/algorithms.htm#ln
  // - https://en.wikipedia.org/wiki/Fast_inverse_square_root

  // FORMULA:
  // x = m * 2^p =>
  //   ln(x) = ln(m) + ln(2)p,

  // first normalize the value to between 1.0 and 2.0
  // assuming normalized IEEE float
  //    sign  exp       frac
  // 0b 0    [00000000] 00000000000000000000000
  // value = (-1)^s * M * 2 ^ (exp-127)
  //
  // exp = 127 for x = 1,
  // so 2^(exp-127) is the multiplier

  // evil floating point bit level hacking
  unsigned int bx = * (unsigned int *) (&x);

  // extract exp, since x>0, sign bit must be 0
  unsigned int ex = bx >> 23;
  signed int t = (signed int)ex-(signed int)127;
  unsigned int s = (t < 0) ? (-t) : t;

  // reinterpret back to float
  //   127 << 23 = 1065353216
  //   0b11111111111111111111111 = 8388607
  bx = 1065353216 | (bx & 8388607);
  x = * (float *) (&bx);


  // use remez algorithm to find approximation between [1,2]
  // - see this answer https://stackoverflow.com/a/44232045
  // - or this usage of C++/boost's remez implementation
  //   https://computingandrecording.wordpress.com/2017/04/24/
  // e.g.
  // boost::math::tools::remez_minimax<double> approx(
  //    [](const double& x) { return log(x); },
  // 4, 0, 1, 2, false, 0, 0, 64);
  //
  // 4th order is:
  // { -1.74178, 2.82117, -1.46994, 0.447178, -0.0565717 }
  //
  // 3rd order is:
  // { -1.49278, 2.11263, -0.729104, 0.10969 }

  return

  /* less accurate */
    -1.49278+(2.11263+(-0.729104+0.10969*x)*x)*x

  /* OR more accurate */
  // -1.7417939+(2.8212026+(-1.4699568+(0.44717955-0.056570851*x)*x)*x)*x

  /* compensate for the ln(2)s. ln(2)=0.6931471806 */
    + 0.6931471806*t;
}


inline double UDDSketch::get_min_addressable_value() {
    return std::fmax(exp((std::numeric_limits<long>::min()) * this->log_gamma), std::numeric_limits<double>::min());
}

int UDDSketch::get_number_buckets() {
   return this->store.size();
}

std::size_t UDDSketch::get_size_bytes() {
    // init_alpha|alpha|m|zero_bucket|store_size|[bucket_key|bucket_count|bucket_key|bucket_count|...]
    size_t bytes_size = 2*sizeof(double) + sizeof(int) + sizeof(long) + sizeof(size_t) + this->store.size() * (sizeof(long) + sizeof(long));
    return bytes_size;
}

std::size_t UDDSketch::get_max_size_bytes(double min, double max) {
    // alpha|m|zero_bucket|store_size|[bucket_key|bucket_count|bucket_key|bucket_count|...]
    size_t bytes_size;

    if (this->m > 0) {
        bytes_size = 2*sizeof(double) + sizeof(int) + sizeof(long) + sizeof(size_t) + this->m * (sizeof(long) + sizeof(long));
    } else {
        if (min == 0.0 && max == 0.0) {
            throw std::logic_error("Bad domain range");
        }
        if (!(min > 0)) {
            min = this->get_min_addressable_value();
        }
        int max_m = std::ceil(std::log(max/min)/this->log_gamma);
        bytes_size = 2*sizeof(double) + sizeof(int) + sizeof(long) + sizeof(size_t) + max_m * (sizeof(long) + sizeof(long));
    }

    return bytes_size;
}

int UDDSketch::get_max_number_buckets() {
    return this->m;
}

double UDDSketch::get_alpha() {
    return this->alpha;
}

double UDDSketch::get_initial_alpha() {
    return this->initial_alpha;
}

double UDDSketch::get_gamma() {
    return this->gamma;
}

long UDDSketch::get_bucket_count(long key) {
    auto bucket = this->store.find(key);
    return (bucket != this->store.end()? bucket->second : -1);
}

long UDDSketch::get_zero_bucket_count() {
    return this->zero_bucket;
}

long UDDSketch::get_sketch_count(){
    long tcount = 0;
    for (auto e : this->store) {
        tcount += e.second;
    }
    return tcount+this->zero_bucket;
}

inline long UDDSketch::get_bucket_key(double value) {
    long bkey = (int)std::ceil( natural_log(value) * this->inverse_log_gamma );
    return bkey;
}

inline double UDDSketch::get_bucket_value(long bkey) {
    //double value = exp(bkey * this->log_gamma) * (1 - this->alpha);
    double value = pow(this->gamma, bkey) * (1 - this->alpha); //  2gamma^bkey / (gamma + 1)
    return value;
}

void UDDSketch::sketch_dump() {
    std::cout << "initial alpha: " << this->initial_alpha
              << "alpha: " << this->alpha << "\n"
              << "gamma: " << this->gamma << "\n"
              << "log_gamma: " << this->log_gamma << "\n"
              << "zero bucket: " << this->zero_bucket << "\n"
              << "m: " << this->m << "\n"
              << "final m: " << this->store.size() << "\n";
    std::map<long, long> tmpstore;
    for (auto b : this->store) {
        tmpstore[b.first] = b.second;
    }
    for (auto b : tmpstore) {
        std::cout << b.first << ":" << b.second << "\n";
    }
}

void UDDSketch::add(double value) {
    if (value < this->min_addressable_value) {
        this->zero_bucket++;
    } else {
        long bkey = (long)std::ceil( natural_log(value) * this->inverse_log_gamma );
        this->store[bkey]++;
        //std::cout << value << ":" << std::log(value) << ":" << bkey << "\n";
        while((this->m > 0) && (this->store.size() > this->m)) {
            this->collapse();
        }
    }
}

int UDDSketch::remove(double value) {
    if (value < get_min_addressable_value()) {
        if (this->zero_bucket > 0){
            this->zero_bucket--;
            return 0;
        }
        return 1;
    } else {
        long bkey = (long)std::ceil( natural_log(value) * this->inverse_log_gamma );
        auto search = this->store.find(bkey);
        if (search != this->store.end() && search->second > 0) {
            search->second--;
            if (search->second == 0) {
                this->store.erase(search);
            }
            return 0;
        }
        return 1;
    }
}

void UDDSketch::collapse(int n) {
    if (n <= 0) return;
    std::unordered_map<long, long> new_store;
    int twopow = std::pow(2, n);
    for (const auto& el : this->store) {
        long bkey = (long)std::ceil((double)el.first / twopow);
        new_store[bkey] += el.second;
    }
    this->store.swap(new_store);
    this->gamma = std::pow(this->gamma, twopow);
    this->log_gamma = twopow * this->log_gamma;
    this->inverse_log_gamma = 1.0 / this->log_gamma;
    this->min_addressable_value = get_min_addressable_value();
    this->alpha = (gamma - 1)/(gamma + 1);
}

double UDDSketch::get_quantile(double q) {
    long rank = std::floor(q * (get_sketch_count() - 1));
    if (rank < this->zero_bucket) {
        return 0.0;
    }
    std::map<long, long> ordered_store;

    for (auto e : this->store) {
        ordered_store[e.first] = e.second;
    }

    auto it = ordered_store.begin();
    long count = this->zero_bucket + it->second;
    long bkey = it->first;
    while (count <= rank && it != ordered_store.end()) {
        it++;
        count += it->second;
        bkey = it->first;
    }
    return get_bucket_value(bkey);
}

template<typename T>
static inline size_t copy_to_mem(const T& item, void* dst) {
  memcpy(dst, &item, sizeof(T));
  return sizeof(T);
}

template<typename T>
static inline size_t copy_from_mem(void* src, T& item) {
  memcpy(&item, src, sizeof(T));
  return sizeof(T);
}

std::vector<uint8_t> UDDSketch::serialize() {
    size_t byte_size = get_size_bytes();
    std::vector<uint8_t> bytes(byte_size, 0);
    uint8_t* bptr = bytes.data();
    bptr += copy_to_mem<double>(this->initial_alpha, bptr);
    bptr += copy_to_mem<double>(this->alpha, bptr);
    bptr += copy_to_mem<int>(this->m, bptr);
    bptr += copy_to_mem<long>(this->zero_bucket, bptr);
    size_t store_size = this->store.size();
    bptr += copy_to_mem<size_t>(store_size, bptr);
    for (const auto &b : this->store) {
        bptr += copy_to_mem<long>(b.first, bptr);
        bptr += copy_to_mem<long>(b.second, bptr);
    }
    return bytes;
}

static inline void check_serialized_size(size_t size, size_t ssize) {
    size_t size_bytes = 2*sizeof(double) + sizeof(int) + sizeof(long) + sizeof(size_t) + ssize * (sizeof(long) + sizeof(long));
    if(size_bytes != size) throw std::logic_error("Serialized size error");
}

UDDSketch UDDSketch::deserialize(void *bytes, std::size_t size){
    double initial_alpha = 0.0;
    double alpha = 0.0;
    uint8_t *bptr = (uint8_t*)bytes;
    bptr += copy_from_mem<double>(bptr, initial_alpha);
    bptr += copy_from_mem<double>(bptr, alpha);
    UDDSketch sketch(alpha);
    sketch.initial_alpha = initial_alpha;
    bptr += copy_from_mem<int>(bptr, sketch.m);
    bptr += copy_from_mem<long>(bptr, sketch.zero_bucket);
    size_t store_size = 0;
    bptr += copy_from_mem<size_t>(bptr, store_size);
    check_serialized_size(size, store_size);
    for (int i = 0; i < store_size; i++) {
        long bkey;
        long bcount;
        bptr += copy_from_mem<long>(bptr, bkey);
        bptr += copy_from_mem<long>(bptr, bcount);
        sketch.store[bkey] = bcount;
    }

    return sketch;
}

void UDDSketch::merge(UDDSketch& other ) {

    if (this->initial_alpha != other.initial_alpha) {
        throw std::logic_error("Merging sketches have incompatible alpha values");
    }
    double gammafrac = this->log_gamma > other.log_gamma? this->log_gamma/other.log_gamma : other.log_gamma / this->log_gamma;
    int twopow =  round(gammafrac);
    double tmp = std::log2(twopow);
    double c;
    if (modf(tmp, &c) != 0) {
	std::cout << twopow << " " << tmp << ": " << this->initial_alpha << " " << this->log_gamma << " " << other.log_gamma << "\n";
	throw std::logic_error("Merging sketches have incompatible alpha values: expected power of two value here");
    }
    this->log_gamma > other.log_gamma ?  other.collapse(c) : this->collapse(c);

    this->zero_bucket += other.zero_bucket;
    for (auto b : other.store ) {
        this->store[b.first] += b.second;
    }

    while((this->m > 0) && (this->store.size() > this->m)) {
        this->collapse();
    }
}
