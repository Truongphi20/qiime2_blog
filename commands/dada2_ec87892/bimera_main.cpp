#include <iostream>
#include <vector>
#include <string>

// Rcpp and Embedded R Headers
#include <Rcpp.h>
#include <Rembedded.h>

/**
 * FORWARD DECLARATION
 * This must match the signature of the function you are testing.
 */
extern Rcpp::DataFrame C_table_bimera2(
    Rcpp::IntegerMatrix mat, 
    std::vector<std::string> seqs, 
    double min_fold, 
    int min_abund, 
    bool allow_one_off, 
    int min_one_off_par_dist, 
    int match, 
    int mismatch, 
    int gap_p, 
    int max_shift
);

int main(int argc, char** argv) {
    // 1. Initialize Embedded R
    char *r_args[] = {(char*)"R", (char*)"--silent", (char*)"--vanilla"};
    Rf_initEmbeddedR(sizeof(r_args)/sizeof(r_args[0]), r_args);

    std::cout << "--- DADA2 Bimera Debug Harness ---" << std::endl;

    try {
        // 2. Load Rcpp namespace
        Rcpp::Environment base = Rcpp::Environment::base_env();
        Rcpp::Function library = base["library"];
        library("Rcpp");

        // 3. Load the data from RDS
        std::cout << "Loading workspace data..." << std::endl;
        Rcpp::Function read_rds = base["readRDS"];
        Rcpp::List data = read_rds("/workspaces/qiime2_blog/commands/tmp_data/debug_bimera_data.rds"); 

        // 4. Extract Data 
        // Note: Rcpp::as<> handles the conversion from R objects to C++ types
        Rcpp::IntegerMatrix mat      = data["mat"];
        std::vector<std::string> sqs = Rcpp::as<std::vector<std::string>>(data["seqs"]);
        
        double min_fold              = data["min_fold"];
        int min_abund                = data["min_abund"];
        bool allow_one_off           = data["allow_one_off"];
        int min_one_off_par_dist     = data["min_one_off_par_dist"];
        int match                    = data["match"];
        int mismatch                 = data["mismatch"];
        int gap_p                    = data["gap_p"];
        int max_shift                = data["max_shift"];

        std::cout << "Data Loaded. Matrix: " << mat.nrow() << "x" << mat.ncol() << std::endl;
        std::cout << "Sequences: " << sqs.size() << std::endl;

        // 5. Invoke the Bimera algorithm
        std::cout << "Running C_table_bimera2..." << std::endl;
        Rcpp::DataFrame result = C_table_bimera2(
            mat, sqs, min_fold, min_abund, allow_one_off, 
            min_one_off_par_dist, match, mismatch, gap_p, max_shift
        );

        // 6. Inspect Result
        Rcpp::IntegerVector nflag = result["nflag"];
        std::cout << "Success! Processed " << nflag.size() << " columns." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "FATAL ERROR: " << e.what() << std::endl;
    }

    // 7. Cleanup
    Rf_endEmbeddedR(0);
    return 0;
}