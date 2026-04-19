#include <iostream>
#include <vector>
#include <string>
#include <cstdint>

// Rcpp and Embedded R Headers
#include <Rcpp.h>
#include <Rembedded.h>

// rds2cpp Headers
#include "rds2cpp/rds2cpp.hpp"

/**
 * FORWARD DECLARATION
 * This signature must match Rmain.cpp exactly for the linker to resolve it.
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

/**
 * Helper to find list elements by name in rds2cpp GenericVector.
 */
const rds2cpp::RObject* find_element(const rds2cpp::GenericVector* vptr, 
                                    const std::vector<std::string>& names, 
                                    const std::string& key) {
    for (size_t i = 0; i < names.size(); ++i) {
        if (names[i] == key) return vptr->data[i].get();
    }
    return nullptr;
}

int main(int argc, char** argv) {
    // 1. Initialize Embedded R
    // Mandatory for Rcpp types (NumericMatrix, List).
    char *r_args[] = {(char*)"R", (char*)"--silent", (char*)"--vanilla"};
    Rf_initEmbeddedR(sizeof(r_args)/sizeof(r_args[0]), r_args);

    std::cout << "--- DADA2 Standalone Debug Harness ---" << std::endl;

    try {
        // 2. Load RDS via rds2cpp
        std::string rds_path = "/workspaces/qiime2_blog/commands/tmp_data/dada_uniques_env.rds";
        auto file_info = rds2cpp::parse_rds(rds_path, rds2cpp::ParseRdsOptions());
        auto vptr = static_cast<const rds2cpp::GenericVector*>(file_info.object.get());

        // 3. Map names
        std::vector<std::string> element_names;
        for (const auto& attr : vptr->attributes) {
            if (file_info.symbols[attr.name.index].name == "names") {
                auto nptr = static_cast<const rds2cpp::StringVector*>(attr.value.get());
                for (const auto& s : nptr->data) {
                    element_names.push_back(s.value.has_value() ? *(s.value) : "");
                }
            }
        }

        // Extraction lambdas
        auto get_int = [&](const std::string& k) {
            return static_cast<const rds2cpp::IntegerVector*>(find_element(vptr, element_names, k))->data[0];
        };
        auto get_dbl = [&](const std::string& k) {
            return static_cast<const rds2cpp::DoubleVector*>(find_element(vptr, element_names, k))->data[0];
        };
        auto get_log = [&](const std::string& k) {
            return (bool)static_cast<const rds2cpp::LogicalVector*>(find_element(vptr, element_names, k))->data[0];
        };

        // 4. Data Conversion
        // Sequences
        auto s_obj = static_cast<const rds2cpp::StringVector*>(find_element(vptr, element_names, "seqs"));
        std::vector<std::string> seqs;
        for(const auto& s : s_obj->data) {
            if(s.value.has_value()) seqs.push_back(*(s.value));
        }

        // Abundances & Priors
        std::vector<int> abundances = static_cast<const rds2cpp::IntegerVector*>(find_element(vptr, element_names, "abundances"))->data;
        auto p_raw = static_cast<const rds2cpp::LogicalVector*>(find_element(vptr, element_names, "priors"))->data;
        std::vector<bool> priors(p_raw.begin(), p_raw.end());

        // Matrix Bridge (rds2cpp flat vector -> Rcpp::NumericMatrix)
        auto err_data = static_cast<const rds2cpp::DoubleVector*>(find_element(vptr, element_names, "err"))->data;
        Rcpp::NumericMatrix err(16, 41, err_data.begin());

        auto quals_data = static_cast<const rds2cpp::DoubleVector*>(find_element(vptr, element_names, "quals"))->data;
        int nraw = seqs.size();
        int maxlen = quals_data.size() / nraw;
        Rcpp::NumericMatrix quals(maxlen, nraw, quals_data.begin());

        std::cout << "Data loaded. Calling dada_uniques..." << std::endl;

        // 5. Execution
        Rcpp::List result = dada_uniques(
            seqs, abundances, priors, err, quals,
            get_int("match"), get_int("mismatch"), get_int("gap"),
            get_log("use_kmers"), get_dbl("kdist_cutoff"),
            get_int("band_size"), get_dbl("omegaA"), get_dbl("omegaP"),
            get_dbl("omegaC"), get_log("detect_singletons"),
            get_int("max_clust"), get_dbl("min_fold"), get_int("min_hamming"),
            get_int("min_abund"), get_log("use_quals"), get_log("final_consensus"),
            get_log("vectorized_alignment"), get_int("homo_gap"),
            get_log("multithread"), get_log("verbose"), get_int("SSE"),
            get_log("gapless"), get_log("greedy")
        );

        std::cout << "--- DADA2 Finished Successfully ---" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "FATAL: " << e.what() << std::endl;
    }

    Rf_endEmbeddedR(0);
    return 0;
}