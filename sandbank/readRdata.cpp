#include <iostream>
#include <string>
#include <vector>
#include "rds2cpp/rds2cpp.hpp"

int main() {
    std::string file_path = "/workspaces/qiime2_blog/commands/tmp_data/dada_uniques_env.rds";

    auto file_info = rds2cpp::parse_rds(file_path, rds2cpp::ParseRdsOptions());
    const auto& ptr = file_info.object;

    if (ptr->type() == rds2cpp::SEXPType::VEC) {
        auto vptr = static_cast<const rds2cpp::GenericVector*>(ptr.get());

        // 1. Extract the names safely
        std::vector<std::string> element_names;
        for (const auto& attr : vptr->attributes) {
            if (file_info.symbols[attr.name.index].name == "names") {
                auto nptr = static_cast<const rds2cpp::StringVector*>(attr.value.get());
                
                // Use .data if .value is not recognized
                for (const auto& r_str : nptr->data) {
                    if (r_str.value.has_value()) {
                        element_names.push_back(*(r_str.value));
                    } else {
                        element_names.push_back(""); 
                    }
                }
            }
        }

        // 2. Loop through elements
        for (size_t i = 0; i < vptr->data.size(); ++i) {
            std::string current_name = (i < element_names.size()) ? element_names[i] : "";

            if (current_name == "match") {
                auto match_ptr = static_cast<const rds2cpp::IntegerVector*>(vptr->data[i].get());
                std::cout << "Found match: " << match_ptr->data[0] << std::endl;
            }
            
            if (current_name == "seqs") {
                auto seq_ptr = static_cast<const rds2cpp::StringVector*>(vptr->data[i].get());
                // Use .data here as well
                std::cout << "Found " << seq_ptr->data.size() << " sequences." << std::endl;
            }
        }
    }

    return 0;
}