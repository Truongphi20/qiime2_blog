#include <iostream>
#include <vector>
#include <string>
#include <dlfcn.h>      // Required for dlopen fix
#include <stdlib.h>     // Required for setenv
#include <Rcpp.h>
#include <Rembedded.h>  // Required for Rf_initEmbeddedR
#include "src/dada.h"

/**
 * FORWARD DECLARATION
 * This must match the entry point in your DADA2 C++ source (typically Rmain.cpp).
 */
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

int main(int argc, char *argv[]) {
    std::cout << "--- DADA2 Standalone Debug Harness ---" << std::endl;

    // 1. FORCE GLOBAL SYMBOL VISIBILITY
    // This resolves 'Rcpp_precious_remove' by making Rcpp symbols global before R starts.
    const char* rcpp_so = "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R/library/Rcpp/libs/Rcpp.so";
    void* rcpp_handle = dlopen(rcpp_so, RTLD_NOW | RTLD_GLOBAL);
    if (!rcpp_handle) {
        std::cerr << "CRITICAL: Could not manually load Rcpp.so: " << dlerror() << std::endl;
        return 1;
    }
    std::cout << "Rcpp symbols exported globally." << std::endl;

    // 2. ENVIRONMENT OVERRIDES
    const char* r_home_path = "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R";
    setenv("R_HOME", r_home_path, 1);
    setenv("R_ENABLE_JIT", "0", 1); 
    setenv("R_LIBS_SITE", "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R/library", 1);

    // 3. INITIALIZE EMBEDDED R
    std::cout << "Initializing Embedded R engine..." << std::endl;
    char *r_args[] = {(char*)"R", (char*)"--silent", (char*)"--vanilla"};
    Rf_initEmbeddedR(sizeof(r_args)/sizeof(r_args[0]), r_args);

    try {
        Rcpp::Environment global = Rcpp::Environment::global_env();

        // 4. LOAD DATA
        std::cout << "Loading debug_state.RData..." << std::endl;
        Rcpp::Function load("load");
        load("debug_state.RData");

        // 5. EXTRACT DATA
        auto seqs        = Rcpp::as<std::vector<std::string>>(global["seqs"]);
        auto abundances  = Rcpp::as<std::vector<int>>(global["abundances"]);
        auto priors      = Rcpp::as<std::vector<bool>>(global["priors"]);
        Rcpp::NumericMatrix err   = global["err"];
        Rcpp::NumericMatrix quals = global["quals"];

        std::cout << "Data loaded: " << seqs.size() << " sequences ready." << std::endl;

        // 6. EXECUTION
        std::cout << "Calling dada_uniques..." << std::endl;
        Rcpp::List result = dada_uniques(
            seqs, abundances, priors, err, quals, 
            global["match"], global["mismatch"], global["gap"], 
            global["use_kmers"], global["kdist_cutoff"], 
            global["band_size"], global["omegaA"], global["omegaP"], 
            global["omegaC"], global["detect_singletons"], 
            global["max_clust"], global["min_fold"], global["min_hamming"], 
            global["min_abund"], global["use_quals"], global["final_consensus"], 
            global["vectorized_alignment"], global["homo_gap"], 
            global["multithread"], global["verbose"], global["SSE"], 
            global["gapless"], global["greedy"]
        );

        std::cout << "--- DADA2 Execution Finished ---" << std::endl;
        
        if (result.containsElementNamed("clustering")) {
            Rcpp::DataFrame clustering = result["clustering"];
            std::cout << "Result: " << clustering.nrow() << " clusters found." << std::endl;
        }

    } catch (std::exception &e) {
        std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
    }

    // 7. CLEANUP
    if (rcpp_handle) dlclose(rcpp_handle);
    Rf_endEmbeddedR(0);
    return 0;
}