#include <Rcpp.h>
#include <Rembedded.h> // Required for standalone execution
#include <stdlib.h>

using namespace Rcpp;

// [[Rcpp::export]]
void read_r_data(std::string file_path) {
    Environment env = Environment::global_env();
    
    Function load("load");
    load(file_path);

    // Cast the binding to a specific type (e.g., int) so it can be printed
    int match_val = Rcpp::as<int>(env["match"]);
    Rcout << "Read object 'match'. Value: " << match_val << std::endl;
}

int main() {
    // Standard setup for embedded R
    setenv("R_HOME", "/opt/conda/envs/qiime2-amplicon-2026.1/lib/R", 1);
    char *r_args[] = {(char*)"R", (char*)"--silent", (char*)"--vanilla"};
    Rf_initEmbeddedR(sizeof(r_args)/sizeof(r_args[0]), r_args);

    std::string file_path = "/workspaces/qiime2_blog/debug_state.RData";
    
    try {
        read_r_data(file_path);
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    Rf_endEmbeddedR(0);
    return 0;
}