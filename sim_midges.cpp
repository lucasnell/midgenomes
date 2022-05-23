#include <Rcpp.h>

#include <random>
#include <vector>
#include <cstdint>

// [[Rcpp::plugins(openmp, cpp11)]]
// [[Rcpp::depends(RcppProgress)]]


#include <progress.hpp>         // for the progress bar
#ifdef _OPENMP
#include <omp.h>                // OpenMP
#endif
#include "pcg/pcg_random.hpp"   // pcg prng
#include "pcg/pcg_extras.hpp"  // pcg 128-bit integer type


using namespace Rcpp;



typedef uint_fast8_t uint8;
typedef uint_fast32_t uint32;
typedef int_fast32_t sint32;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;
typedef pcg_extras::pcg128_t uint128;

template<typename T> using vector2d = std::vector<std::vector<T>>;





/*
 For consistent results from multi-thread operations, you should call `mt_seeds`
 when outside multi-thread mode, then create a PRNG (`pcg64` object) once
 inside multi-thread mode, one PRNG per thread.
 For each simulation rep, use a `std::vector<uint64>` vector (from inside
 the object output from `mt_seeds`) to seed the PRNG using `seed_pcg`.
 Using seeds for each rep (not thread) means that if you run `set.seed(999)`
 in R before running a process that can be run across multiple threads,
 you'll get the same output no matter how many threads you use.
 */

// To sample for seeds before multi-thread operations
inline std::vector<std::vector<uint64>> mt_seeds(const uint32& n_reps) {

    std::vector<std::vector<uint64>> sub_seeds(n_reps, std::vector<uint64>(8));

    for (uint32 i = 0; i < n_reps; i++) {
        // These are 32-bit integers cast as 64-bit for downstream compatibility
        sub_seeds[i] = as<std::vector<uint64>>(Rcpp::runif(8,0,4294967296));
    }

    return sub_seeds;
}


// Change seed of existing rng from 8 32-bit seeds (casted to 64-bit)
inline void seed_pcg(pcg64& eng, const std::vector<uint64>& sub_seeds) {

    uint128 seed1;
    uint128 seed2;
    // Converting to two 128-bit seeds for pcg64
    uint128 seed64_1 = (sub_seeds[0]<<32) + sub_seeds[1];
    uint128 seed64_2 = (sub_seeds[2]<<32) + sub_seeds[3];
    uint128 seed64_3 = (sub_seeds[4]<<32) + sub_seeds[5];
    uint128 seed64_4 = (sub_seeds[6]<<32) + sub_seeds[7];
    seed1 = (seed64_1<<64) + seed64_2;
    seed2 = (seed64_3<<64) + seed64_4;

    eng.seed(seed1, seed2);

    return;
}


/*
 Allows you to verify that you're able to use multiple threads.
 */
//[[Rcpp::export]]
bool using_openmp() {
    bool out = false;
#ifdef _OPENMP
    out = true;
#endif
    return out;
}


//' Check that the number of threads doesn't exceed the number available, and change
//' to 1 if OpenMP isn't enabled.
//'
//' @noRd
//'
inline void thread_check(uint32& n_threads) {

#ifdef _OPENMP
    if (n_threads == 0) n_threads = 1;
    if (n_threads > omp_get_max_threads()) {
        std::string max_threads = std::to_string(omp_get_max_threads());
        std::string err_msg = "\nThe number of requested threads (" +
            std::to_string(n_threads) +
            ") exceeds the max available on the system (" +
            max_threads + ").";
        stop(err_msg.c_str());
    }
#else
    n_threads = 1;
#endif

    return;
}

// For checking for user interrupts every N iterations:
inline bool interrupt_check(uint32& iters,
                            Progress& prog_bar,
                            const uint32& N = 100) {
    ++iters;
    if (iters > N) {
        if (prog_bar.is_aborted() || prog_bar.check_abort()) return true;
        prog_bar.increment(iters);
        iters = 0;
    }
    return false;
}





// Population info for a single set of simulations.
struct Population {

    uint32 n_alleles;           // Number of alleles
    std::vector<double> P;      // Allele frequencies
    double mu;                  // Mutation rate

    Population() : n_alleles(), P(), mu() {}
    Population(const std::vector<double>& P_,
               const double& mu_)
        : n_alleles(P_.size()), P(P_), mu(mu_) {}
    Population(const Population& other)
        : n_alleles(other.n_alleles), P(other.P), mu(other.mu) {}
    Population& operator=(const Population& other) {
        n_alleles = other.n_alleles;
        P = other.P;
        mu = other.mu;
        return *this;
    }

    // Update allele frequencies and counts, return average heterozygosity
    double update(const double& N, pcg64& eng) {
        /*
         Use binomial if hapN < 1000 or if, for allele i:
            hapN * P[i] < 100 OR hapN * (1 - P[i]) < 100
         */
        double hapN = N * 2;  // <-- haploid N
        bool low_N = hapN < 1000;
        // Don't want to be anywhere near integer limits:
        bool high_N = hapN >= 2147483647.0;
        double S, M1, M2, X;
        uint32 rS;
        // If you need to use this for binomial:
        uint32 rhapN = static_cast<uint32>(std::round(hapN));
        double p_sum = 0;
        for (uint32 i = 0; i < n_alleles; i++) {
            S = one_sample(high_N, low_N, hapN, rhapN, P[i], eng);
            rS = static_cast<uint32>(std::round(S));
            M1 = one_sample(high_N, low_N, hapN - S, rhapN - rS, mu, eng);
            M2 = one_sample(high_N, low_N, S, rS, mu, eng);
            X = S + M1 - M2;
            P[i] = X / hapN;
            p_sum += P[i] * P[i];
            p_sum += (1 - P[i]) * (1 - P[i]);
        }
        double het = 1 - p_sum / static_cast<double>(n_alleles);
        return het;
    }

private:


    //' This does all the checks the below functions don't do,
    //' and always returns a double.
    //' Note that `n` here does NOT indicate census pop. size (`N`).
    //' It's the number of trials.
    //'
    double one_sample(const bool& high_N,
                      const bool& low_N,
                      const double& n,
                      const uint32& rn,
                      const double& p,
                      pcg64& eng) {

        if (p == 0 || n <= 0) return 0.0;
        if (p == 1) return n;

        bool use_binom = !high_N && (low_N || (n * p < 100) ||
                                     (n * (1 - p) < 100));
        double x;
        if (use_binom) {
            uint32 rx = one_sample_internal(rn, p, eng);
            x = static_cast<double>(rx);
        } else {
            x = one_sample_internal(n, p, eng);
        }
        return x;
    }


    //' Below two functions are for a single sampling event.
    //' Checks for p == 1 or p == 0 should happen beforehand.

    //' using normal approximation (happens when `n` is a double)
    double one_sample_internal(const double& n,
                               const double& p,
                               pcg64& eng) {
        double sig = std::sqrt(n * p * (1 - p));
        double x = n*p + norm(eng) * sig;
        // This should never happen, but just in case:
        if (x < 0) x = 0;
        if (x > n) x = n;
        return x;
    }
    //' using binomial distribution (when `n` is integer):
    uint32 one_sample_internal(const uint32& n,
                               const double& p,
                               pcg64& eng) {
        binom.param(std::binomial_distribution<uint32>::param_type(n, p));
        uint32 x = binom(eng);
        return x;
    }

    std::binomial_distribution<uint32> binom;
    std::normal_distribution<double> norm =
        std::normal_distribution<double>(0, 1);

};




//'
//' To calculate average heterozygosity.
//'
void avg_het__(double& het, const std::vector<double>& P) {
    double p_sum = 0;
    for (uint32 i = 0; i < P.size(); i++) {
        p_sum += P[i] * P[i];
        p_sum += (1 - P[i]) * (1 - P[i]);
    }
    het = 1 - p_sum / static_cast<double>(P.size());
    return;
}
//' Same thing, but exported to R.
// [[Rcpp::export]]
double avg_het(const std::vector<double>& P) {
    double het;
    avg_het__(het, P);
    return het;
}

//' For a given vector of census population sizes and starting allele
//' frequencies (some can and should be zeros/ones), simulate changes in these
//' counts though time, then calculate average heterozygosity.
//'
// [[Rcpp::export]]
vector2d<double> sim_het(const std::vector<double>& N,
                         const std::vector<double>& P0,
                         const double& mu,
                         const uint32& n_reps,
                         uint32 n_threads,
                         const bool& show_progress = false) {

    thread_check(n_threads);
    Progress prog_bar(N.size() * n_reps, show_progress);
    std::vector<int> status_codes(n_threads, 0);

    vector2d<double> out(n_reps, std::vector<double>(N.size()));
    // Fill in first row of average heterozygosities
    double het;
    avg_het__(het, P0);
    for (uint32 i = 0; i < n_reps; i++) out[i][0] = het;

    const std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);


#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(n_threads) if (n_threads > 1)
{
#endif
    pcg64 eng;

    // Write the active seed per thread or just write one of the seeds.
#ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
#else
    uint32 active_thread = 0;
#endif
    int& status_code(status_codes[active_thread]);

    // Parallelize the Loop
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (uint32 i = 0; i < n_reps; i++) {
        if (status_code != 0) continue;
        seed_pcg(eng, seeds[i]);
        Population pop(P0, mu);
        for (uint32 j = 1; j < N.size(); j++) {
            out[i][j] = pop.update(N[j], eng);
        }
    }

#ifdef _OPENMP
}
#endif

    return out;
}

