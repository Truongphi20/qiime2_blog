#include <iostream>
#include <vector>
#include <string>
#include "src/dada.h"


// Forward declaration of the function in Rmain.cpp
Rcpp::List dada_uniques(std::vector< std::string > seqs, std::vector<int> abundances, std::vector<bool> priors,
                        Rcpp::NumericMatrix err,
                        Rcpp::NumericMatrix quals,
                        int match, int mismatch, int gap,
                        bool use_kmers, double kdist_cutoff,
                        int band_size,
                        double omegaA, double omegaP, double omegaC, bool detect_singletons,
                        int max_clust,
                        double min_fold, int min_hamming, int min_abund,
                        bool use_quals,
                        bool final_consensus,
                        bool vectorized_alignment,
                        int homo_gap,
                        bool multithread,
                        bool verbose,
                        int SSE,
                        bool gapless,
                        bool greedy);

// Mocking Rcpp types if running in a standalone C++ test harness 
// Note: If you are using Rcpp directly, you must initialize an R environment
int main() {
    // 1. DATA PREPARATION (Based on your ls.str() output)
    
    // seqs : chr [1:2490]
    std::vector<std::string> seqs = {
        "TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGATGGATGTTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGATACTGGATGTCTT",
        // ... add more representative sequences here for testing
    };

    // abundances : int [1:2490]
    std::vector<int> abundances = {1497, 693, 519, 357, 269}; 
    // Truncate vectors for manual C++ testing if needed, or load from CSV
    abundances.resize(seqs.size());

    // priors : logi [1:2490]
    std::vector<bool> priors(seqs.size(), false);

    // err : num [1:16, 1:41]
    // Initializing a 16x41 matrix with a dummy value (0.25) or your specific rates
    Rcpp::NumericMatrix err(16, 41);
    std::fill(err.begin(), err.end(), 0.05); 

    // quals : num [1:120, 1:2490] (120 positions, 2490 sequences)
    int max_len = 120;
    Rcpp::NumericMatrix quals(max_len, seqs.size());
    std::fill(quals.begin(), quals.end(), 39.0); // High quality default

    // 2. SCALAR PARAMETERS (From your opts list)
    int match = 5;
    int mismatch = -4;
    int gap = -8;
    int homo_gap = -8;
    bool use_kmers = true;
    double kdist_cutoff = 0.42;
    int band_size = 16;
    double omegaA = 1e-40;
    double omegaP = 1e-04;
    double omegaC = 1e-40;
    bool detect_singletons = false;
    int max_clust = 0; // 0 means it will be set to nraw in the function
    double min_fold = 1.0;
    int min_hamming = 1;
    int min_abund = 1;
    bool use_quals = true;
    bool final_consensus = false;
    bool vectorized_alignment = true;
    bool multithread = false;
    bool verbose = true; // Turn on for debugging
    int SSE = 2;
    bool gapless = true;
    bool greedy = true;

    // 3. EXECUTION
    std::cout << "Starting dada_uniques test with " << seqs.size() << " sequences..." << std::endl;

    try {
        Rcpp::List result = dada_uniques(
            seqs, abundances, priors, err, quals, 
            match, mismatch, gap, use_kmers, kdist_cutoff, 
            band_size, omegaA, omegaP, omegaC, detect_singletons, 
            max_clust, min_fold, min_hamming, min_abund, 
            use_quals, final_consensus, vectorized_alignment, 
            homo_gap, multithread, verbose, SSE, gapless, greedy
        );

        std::cout << "Success! Results generated." << std::endl;
        // Accessing map result as an example
        Rcpp::IntegerVector cluster_map = result["map"];
        std::cout << "First sequence assigned to cluster: " << cluster_map[0] << std::endl;

    } catch (std::exception &e) {
        std::cerr << "Error during DADA2 execution: " << e.what() << std::endl;
    }

    return 0;
}