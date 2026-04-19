#include <iostream>
#include <vector>
#include <string>
#include <cstdint>

// Rcpp and Embedded R Headers
#include <Rcpp.h>
#include <Rembedded.h>

/**
 * FORWARD DECLARATION
 * Must match the signature in Rmain.cpp exactly.
 */
extern Rcpp::List dada_uniques(
    std::vector<std::string> seqs, 
    std::vector<int> abundances, 
    std::vector<bool> priors,
    Rcpp::NumericMatrix err, 
    Rcpp::NumericMatrix quals,
    int match, int mismatch, int gap,
    bool use_kmers, double kdist_cutoff,
    int band_size,
    double omegaA, double omegaP, double omegaC, 
    bool detect_singletons,
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
    bool greedy
);

int main(int argc, char** argv) {
    // 1. Initialize Embedded R
    char *r_args[] = {(char*)"R", (char*)"--silent", (char*)"--vanilla"};
    Rf_initEmbeddedR(sizeof(r_args)/sizeof(r_args[0]), r_args);

    std::cout << "--- DADA2 Standalone Debug Harness (Pure Rcpp) ---" << std::endl;

    try {
        // 2. Load Rcpp namespace to prevent JIT/Precious list errors
        Rcpp::Function library = Rcpp::Environment::base_env()["library"];
        library("Rcpp");

        // 3. Use R's native readRDS to load the data
        std::cout << "Reading RDS file directly using R..." << std::endl;
        Rcpp::Function read_rds = Rcpp::Environment::base_env()["readRDS"];
        
        // This natively returns an Rcpp::List containing your environment variables
        Rcpp::List data = read_rds("/workspaces/qiime2_blog/commands/tmp_data/dada_uniques_env.rds");

        // 4. Effortless Rcpp Data Extraction
        std::vector<std::string> seqs        = Rcpp::as<std::vector<std::string>>(data["seqs"]);
        std::vector<int> abundances          = Rcpp::as<std::vector<int>>(data["abundances"]);
        std::vector<bool> priors             = Rcpp::as<std::vector<bool>>(data["priors"]);
        
        Rcpp::NumericMatrix err              = data["err"];
        Rcpp::NumericMatrix quals            = data["quals"];

        std::cout << "Data loaded natively. Seqs count: " << seqs.size() << std::endl;

        // 5. Invoke DADA2 Core
        std::cout << "Invoking dada_uniques..." << std::endl;
        Rcpp::List result = dada_uniques(
            seqs, abundances, priors, err, quals,
            data["match"], data["mismatch"], data["gap"],
            data["use_kmers"], data["kdist_cutoff"],
            data["band_size"], data["omegaA"], data["omegaP"],
            data["omegaC"], data["detect_singletons"],
            data["max_clust"], data["min_fold"], data["min_hamming"],
            data["min_abund"], data["use_quals"], data["final_consensus"],
            data["vectorized_alignment"], data["homo_gap"],
            data["multithread"], data["verbose"], data["SSE"],
            data["gapless"], data["greedy"]
        );

        std::cout << "--- DADA2 finished successfully ---" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "FATAL ERROR: " << e.what() << std::endl;
    }

    // 6. Cleanup
    Rf_endEmbeddedR(0);
    return 0;
}