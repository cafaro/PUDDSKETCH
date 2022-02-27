#include <cstdlib>
#include <mpi.h>
#include <functional>
#include <limits>
#include <cstdint>
#include <cassert>
#include <chrono>
#include "DataSketches/kll_sketch.hpp"
#include "rnd.hpp"
#include "unistd.h"

struct opts {
    double k;         // value of parameter k for the algorithm's sketch
    long loc_n;       // number of items per process
    std::string dist; // random distribution of generated items
    double p;         // first parameter of input stream random distribution
    double q;         // second parameter of input stream random distribution
    unsigned int s;   // seed
    bool csv;         // csv output
};

struct results {
    int num_procs;
    int sketch_size;
    int theo_max_size;
    int64_t seq_time;
    int64_t reduce_time;
    int64_t tot_time;
    std::vector<double> qs;
};

using sketch_type = datasketches::kll_sketch<double>;

using myclock=std::chrono::high_resolution_clock;
using unifint=std::uniform_int_distribution<unsigned>;
using unifreal=std::uniform_real_distribution<double>;

void print_stats(const opts& o, const results& r);

void merge_op( void *in, void *inout, int *len, MPI_Datatype *dptr )
{
    uint8_t *iin = (uint8_t*)in;
    uint8_t *iinout = (uint8_t *)inout;
    sketch_type sketch_in = sketch_type::deserialize(&iin[sizeof(size_t)], *(size_t*)iin);
    sketch_type sketch_inout = sketch_type::deserialize(&iinout[sizeof(size_t)], *(size_t*)iinout);
    sketch_inout.merge(sketch_in);

    auto tmpsk = sketch_inout.serialize();
    size_t s = tmpsk.size();
    std::memcpy(iinout, &s, sizeof(size_t));
    std::memcpy(&iinout[sizeof(size_t)], tmpsk.data(), tmpsk.size());
}

void print_usage(char *prog_name) {
  printf("\nUsage: %s -f <input file name>\n"
         "    Options\n"
         "      -k <double>: parameter K for sketch. Default 200.0\n"
         "      -n <int>: number of items per process. Default 1000\n"
         "      -d <string>: input distribution: unif|norm|exp. Default unif\n"
         "      -p <double>: random distribution first parameter. Default 1.0\n"
         "      -q <double>: random distribution second parameter. Default 1000.0\n"
         "      -s <uint>: random distribution seed. Default 0 (random seed)\n"
         "      -o : csv output.\n\n",
         prog_name);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int pid, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    opts opts{200.0,1000,"unif",1.0,1000.0,0,false};
    int opt;
    while ((opt = getopt(argc, argv, "hk:n:d:p:q:s:o")) != -1) {
      switch (opt) {
      case 'k':
        opts.k = std::stod(optarg);
        break;
      case 'n':
        opts.loc_n = std::stol(optarg);
        break;
      case 'd':
        opts.dist = std::string(optarg);
        break;
      case 'p':
        opts.p = std::stod(optarg);
        break;
      case 'q':
        opts.q = std::stod(optarg);
        break;
      case 's':
        opts.s = std::stoi(optarg);
        break;
      case 'o':
        opts.csv = true;
        break;
      case 'h':
        print_usage(argv[0]);
        exit(EXIT_SUCCESS);
      default: /* '?' */
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
    }

    auto u = rnd_gen<unifint>(unifint::param_type{0, std::numeric_limits<unsigned>::max()}, opts.s);
    unsigned local_seed = 0;
    for (int i = 0; i <= pid; i++){
        local_seed = u.get_value();
    }

    std::shared_ptr<std::vector<double>> in_stream;
    if (opts.dist == "unif") {
        auto rg = rnd_gen<unifreal>(unifreal::param_type{opts.p, opts.q}, local_seed);
        in_stream = rg.get_stream(opts.loc_n);
    } else if (opts.dist == "norm") {
        auto rg = rnd_gen<std::normal_distribution<>>(std::normal_distribution<>::param_type{opts.p, opts.q}, local_seed);
        in_stream = rg.get_stream(opts.loc_n);
    } else if (opts.dist == "exp") {
        auto rg = rnd_gen<std::exponential_distribution<>>(std::exponential_distribution<>::param_type{opts.p}, local_seed);
        in_stream = rg.get_stream(opts.loc_n);
    }

    int theo_max_size = 0;
    double qs_rank[99];
    for (int i = 1; i < 100; i++){
        qs_rank[i-1] = i/100.0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto start_time = myclock::now();

    sketch_type sketch(opts.k);
    for (int i = 0; i < opts.loc_n; i++) {
        sketch.update(in_stream->at(i));
    }

    auto stop_time = myclock::now();
    auto seq_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

    int64_t upd_time_ms = seq_duration.count();
    int64_t max_upd_time_ms;
    MPI_Allreduce(&upd_time_ms, &max_upd_time_ms, 1, MPI_INT64_T, MPI_MAX, MPI_COMM_WORLD);

    start_time = myclock::now();

    theo_max_size = sketch_type::get_max_serialized_size_bytes(opts.k, opts.loc_n);
    int len = theo_max_size+sizeof(size_t);

    uint8_t sketch_bytes[len];
    uint8_t merged_sketch_bytes[len];

    std::vector<uint8_t> tmp_sketch = sketch.serialize();
    size_t s1 = tmp_sketch.size();

    std::memcpy(sketch_bytes, &s1, sizeof(size_t));
    std::memcpy(&sketch_bytes[sizeof(size_t)], tmp_sketch.data(), tmp_sketch.size());

    MPI_Datatype SKETCH_BYTES_T;
    MPI_Type_contiguous (len, MPI_UINT8_T, &SKETCH_BYTES_T);
    MPI_Type_commit(&SKETCH_BYTES_T);

    MPI_Op op;
    MPI_Op_create((MPI_User_function*)merge_op, 1, &op);

    MPI_Reduce(sketch_bytes, merged_sketch_bytes, 1, SKETCH_BYTES_T, op, 0, MPI_COMM_WORLD);

    if (!pid) {
        auto merged_sketch = sketch_type::deserialize(&merged_sketch_bytes[sizeof(size_t)], *(size_t*)merged_sketch_bytes);

        stop_time = myclock::now();
        auto red_duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

        results res;
        res.num_procs = np;
        res.qs = merged_sketch.get_quantiles(qs_rank, 99);
        res.sketch_size = merged_sketch.get_serialized_size_bytes();
        res.theo_max_size = theo_max_size;
        res.seq_time = max_upd_time_ms;
        res.reduce_time = red_duration.count();
        res.tot_time = max_upd_time_ms + red_duration.count();
        print_stats(opts, res);
    }

    MPI_Type_free(&SKETCH_BYTES_T);
    MPI_Finalize();
}

void print_stats(const opts& o, const results& r) {
    if (!o.csv) {
        std::cout << "Parallel KLL algorithm:\n"
                  << "K: " << o.k << "\n"
                  << "Final sketch size: " << r.sketch_size << "\n"
                  << "Theoretical max sketch size: " << r.theo_max_size << "\n"
                  << "Number of elements: " << o.loc_n * r.num_procs << ", " << o.loc_n << " per process\n"
                  << "Dataset distribution: " << o.dist << "\n"
                  << "Random distribution params: " << o.p << ", " << o.q << "\n"
                  << "Number of processes: " << r.num_procs << "\n"
                  << "Sequential processing time: " << r.seq_time << "ms \n"
                  << "Parallel reduction time: " << r.reduce_time << "ms \n"
                  << "Total running time: " << r.tot_time << "ms \n"
                  << "Estimated median: " << r.qs[49] << "\n";
    } else {
        std::cout << o.s << ","
                  << o.loc_n * r.num_procs << ","
                  << o.dist << ","
                  << o.p << ","
                  << o.q << ","
                  << o.k << ","
                  << r.num_procs << ","
                  << r.seq_time << ","
                  << r.reduce_time << ","
                  << r.tot_time << ","
                  << r.sketch_size << ","
                  << r.theo_max_size;
        for (auto q : r.qs) {
            std::printf(",%.10lf", q);
        }
        std::cout << "\n";
    }
}
