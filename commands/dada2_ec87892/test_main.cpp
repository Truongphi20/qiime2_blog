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
 * This must exactly match the signature in src/Rmain.cpp
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

    // 1. GLOBAL SYMBOL VISIBILITY FIX
    // We force-load the Rcpp shared object into the global namespace to 
    // resolve the 'Rcpp_precious_remove' error before the R engine starts.
    const char* rcpp_so_path = "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R/library/Rcpp/libs/Rcpp.so";
    void* rcpp_handle = dlopen(rcpp_so_path, RTLD_NOW | RTLD_GLOBAL);
    if (!rcpp_handle) {
        std::cerr << "CRITICAL: Could not dlopen Rcpp.so: " << dlerror() << std::endl;
        return 1;
    }

    // 2. ENVIRONMENT OVERRIDES
    // Force R to use the Conda environment and disable JIT to prevent init crashes.
    const char* conda_r_home = "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R";
    const char* conda_r_libs = "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R/library";
    
    setenv("R_HOME", conda_r_home, 1);
    setenv("R_LIBS_SITE", conda_r_libs, 1);
    setenv("R_LIBS_USER", "", 1);
    setenv("R_ENABLE_JIT", "0", 1); 

    // 3. INITIALIZE EMBEDDED R
    std::cout << "Initializing Embedded R engine..." << std::endl;
    char *r_args[] = {(char*)"R", (char*)"--silent", (char*)"--vanilla"};
    Rf_initEmbeddedR(sizeof(r_args)/sizeof(r_args[0]), r_args);

    try {
        Rcpp::Environment global = Rcpp::Environment::global_env();

        // 4. LOAD SAVED STATE
        // This assumes you ran save(list=target_vars, file="debug_state.RData") in R
        std::cout << "Loading debug_state.RData..." << std::endl;
        Rcpp::Function load("load");
        load("debug_state.RData");

        // 5. EXTRACT VARIABLES
        // We cast them from R types to the C++ types expected by dada_uniques
        std::cout << "Extracting variables from R environment..." << std::endl;
        auto seqs        = Rcpp::as<std::vector<std::string>>(global["seqs"]);
        auto abundances  = Rcpp::as<std::vector<int>>(global["abundances"]);
        auto priors      = Rcpp::as<std::vector<bool>>(global["priors"]);
        Rcpp::NumericMatrix err   = global["err"];
        Rcpp::NumericMatrix quals = global["quals"];

        std::cout << "Data loaded successfully. (" << seqs.size() << " sequences)" << std::endl;

        // 6. EXECUTION
        // Note: Using global["var"] directly for scalars to keep things clean.
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
        
        // Quick results summary
        if (result.containsElementNamed("clustering")) {
            Rcpp::DataFrame clustering = result["clustering"];
            std::cout << "Final cluster count: " << clustering.nrow() << std::endl;
        }

    } catch (std::exception &e) {
        std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
    }

    // 7. CLEANUP
    std::cout << "Shutting down R engine." << std::endl;
    if (rcpp_handle) dlclose(rcpp_handle);
    Rf_endEmbeddedR(0);
    
    return 0;
}