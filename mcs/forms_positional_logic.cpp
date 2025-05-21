#include "main.h"    // For Forms class declaration
#include <vector>
#include <string>
#include <set>
#include <algorithm> // For std::remove for cleaning up working_mcs in older C++ versions if needed
#include <iostream>
#include <fstream>
#include <random>    // For std::mt19937, std::uniform_int_distribution

// --- Phase 1: Pattern Reduction for Positional MCS ---
// Reduces this->Vars_Glob (full C(N,K) patterns) based on filters from mcs_i_const_ref.
// Stores the result in this->reducedVars_Glob and this->NReducedVars_Glob.
// Saves the reduced patterns to output_reduced_patterns_filename.

void Forms::perform_phase1_pattern_reduction(
    const Forms& mcs_i_const_ref, // MCS_i (filters to check against)
    const Forms& global_pattern_source, // Source of full base patterns (e.g., f_global_pattern_holder)
    const std::string& output_reduced_patterns_filename)
{
    // Log using global_pattern_source.NVars_Glob for the count of basic patterns being processed
    std::cout << "  Phase 1 (Positional): Reducing " << global_pattern_source.NVars_Glob << " basic patterns against MCS_i ("
        << mcs_i_const_ref.Nform1_Glob << " filters)..." << std::endl;
    std::cout << "    Outputting reduced patterns to: " << output_reduced_patterns_filename << std::endl;

    // Clear any previous reducedVars_Glob in 'this' object (this part is correct for 'this' object)
    if (this->reducedVars_Glob != nullptr) {
        for (int i = 0; i < this->NReducedVars_Glob; ++i) delete[] this->reducedVars_Glob[i];
        delete[] this->reducedVars_Glob;
        this->reducedVars_Glob = nullptr;
        this->NReducedVars_Glob = 0;
    }

    // Check the global_pattern_source for available basic patterns
    if (global_pattern_source.Vars_Glob == nullptr || global_pattern_source.NVars_Glob == 0) {
        std::cerr << "    ERROR (perform_phase1): Full basic patterns (global_pattern_source.Vars_Glob) not available to reduce." << std::endl;
        // 'this' object's reducedVars_Glob remains null, NReducedVars_Glob remains 0.
        return;
    }

    std::vector<bool*> kept_pattern_pointers;

    // Iterate ALL original basic patterns from global_pattern_source
    for (int pat_idx = 0; pat_idx < global_pattern_source.NVars_Glob; ++pat_idx) {
        const bool* current_basic_pattern_bits = global_pattern_source.Vars_Glob[pat_idx];
        bool pattern_starts_with_an_mcs_i_filter = false;

        if (mcs_i_const_ref.Nform1_Glob > 0 && mcs_i_const_ref.forms1Glob != nullptr) {
            for (int f_idx = 0; f_idx < mcs_i_const_ref.Nform1_Glob; ++f_idx) {
                const int* filter_from_mcs_i_bits = mcs_i_const_ref.forms1Glob[f_idx];
                // Actual span of the filter from MCS_i. mcs_i_const_ref.N_glob is correct here.
                int filter_actual_span = mcs_i_const_ref.forms1Glob[f_idx][mcs_i_const_ref.N_glob];

                // The pattern length from global_pattern_source must be at least as long as the filter span.
                // global_pattern_source.N_glob is the length of current_basic_pattern_bits.
                if (filter_actual_span <= 0 || filter_actual_span > global_pattern_source.N_glob) continue;

                bool current_filter_matches_prefix = true;
                if (current_basic_pattern_bits[0] == false && filter_from_mcs_i_bits[0] == 1) {
                    current_filter_matches_prefix = false;
                }
                else {
                    for (int k = 0; k < filter_actual_span; ++k) {
                        if (filter_from_mcs_i_bits[k] == 1 && current_basic_pattern_bits[k] == false) {
                            current_filter_matches_prefix = false;
                            break;
                        }
                    }
                }

                if (current_filter_matches_prefix) {
                    pattern_starts_with_an_mcs_i_filter = true;
                    break;
                }
            }
        }

        if (!pattern_starts_with_an_mcs_i_filter) {
            // Keep this pattern: allocate new memory (size from global_pattern_source.N_glob)
            // and copy. The reduced patterns will have the same length as the original base patterns.
            // 'this->N_glob' should be set to global_pattern_source.N_glob for consistency
            // when these reduced patterns are used by 'this' object later.
            bool* pattern_copy = new bool[global_pattern_source.N_glob];
            for (int k = 0; k < global_pattern_source.N_glob; ++k) {
                pattern_copy[k] = current_basic_pattern_bits[k];
            }
            kept_pattern_pointers.push_back(pattern_copy);
        }
    }

    // Populate this->reducedVars_Glob (this is correct for 'this' object)
    this->NReducedVars_Glob = kept_pattern_pointers.size();
    if (this->NReducedVars_Glob > 0) {
        this->reducedVars_Glob = new bool* [this->NReducedVars_Glob];
        for (size_t i = 0; i < kept_pattern_pointers.size(); ++i) {
            this->reducedVars_Glob[i] = kept_pattern_pointers[i];
        }
        // The N_glob for 'this' object should reflect the length of the patterns it now holds.
        // This line assumes that the N_glob of the object being modified ('this') should
        // match the N_glob of the source patterns. This is generally true in your setup.
        this->N_glob = global_pattern_source.N_glob;
    }
    else {
        this->reducedVars_Glob = nullptr;
        // If no patterns are kept, this->N_glob might not be strictly necessary to set,
        // but it doesn't hurt to keep it consistent or 0.
        // this->N_glob = 0; or this->N_glob = global_pattern_source.N_glob;
    }

    std::cout << "    Phase 1 (Positional): Reduction complete. " << this->NReducedVars_Glob << " patterns remain." << std::endl;

    std::ofstream fout_reduced(output_reduced_patterns_filename);
    if (!fout_reduced.is_open()) {
        std::cerr << "    ERROR (perform_phase1): Could not open file to write reduced patterns: "
            << output_reduced_patterns_filename << std::endl;
    }
    else {
        // Header information:
        // this->N_glob now reflects the length of the reduced patterns (should be same as original).
        // this->N_2_glob (K) is a property of the original C(N,K) patterns. It's good to log the
        // original K value from global_pattern_source, or ensure this->N_2_glob was set appropriately.
        fout_reduced << "# Reduced patterns for Positional MCS step. Based on MCS_i (M="
            << mcs_i_const_ref.Nsovp1_Glob << "). Original N=" << global_pattern_source.N_glob
            << ", K=" << global_pattern_source.N_2_glob << std::endl;
        fout_reduced << "NReducedVars_Glob = " << this->NReducedVars_Glob << std::endl;

        for (int i = 0; i < this->NReducedVars_Glob; ++i) {
            // Print bits up to this->N_glob (which should now be pattern length)
            for (int j = 0; j < this->N_glob; ++j) {
                fout_reduced << (this->reducedVars_Glob[i][j] ? '1' : '0');
            }
            fout_reduced << std::endl;
        }
        fout_reduced.close();
        std::cout << "    Reduced patterns saved to " << output_reduced_patterns_filename << std::endl;
    }
}


// --- Phase 2: Single Filter Pruning Step for Positional MCS ---
// Operates on this->forms1Glob (which is a candidate MCS_{i+1}, initially copy of MCS_i filters)
// Uses this->reducedVars_Glob (populated by Phase 1) to guide protection.
// Removes ONE unprotected filter randomly.
// Returns true if a filter was removed, false if stable or no filters left/no unprotected found.
// Implements one iteration of PDF Page 9-10, Steps 2.2, 2.3.
bool Forms::perform_phase2_single_filter_prune_step(std::mt19937& rng_engine) {
    // std::cout << "    Phase 2 (Positional): Pruning step. Nform1_in=" << this->Nform1_Glob
    //           << ", NReducedVars=" << this->NReducedVars_Glob << "..." << std::endl;

    if (this->Nform1_Glob == 0) {
        // std::cerr << "      No filters in current MCS candidate to prune. Stable." << std::endl;
        return false; // No filters to remove
    }

    // If no reduced patterns, all current filters are considered unprotected (unless Nform1_Glob is 0)
    bool all_filters_unprotected_due_to_no_reduced_patterns =
        (this->reducedVars_Glob == nullptr || this->NReducedVars_Glob == 0);

    std::vector<bool> is_filter_protected(this->Nform1_Glob, false);

    if (!all_filters_unprotected_due_to_no_reduced_patterns) {
        for (int i = 0; i < this->NReducedVars_Glob; ++i) { // Iterate through each reduced_pattern
            // Convert reduced_pattern (bool*) to std::string for substring search
            std::string reduced_pattern_str_version;
            reduced_pattern_str_version.reserve(this->N_glob);
            for (int bit_idx = 0; bit_idx < this->N_glob; ++bit_idx) {
                reduced_pattern_str_version += (this->reducedVars_Glob[i][bit_idx] ? '1' : '0');
            }

            std::vector<int> covering_filter_indices_for_this_pattern;
            for (int f_idx = 0; f_idx < this->Nform1_Glob; ++f_idx) { // Iterate current filters
                if (this->forms1Glob[f_idx] == nullptr) continue; // Should not happen if Nform1_Glob is correct

                int filter_actual_span = this->forms1Glob[f_idx][this->N_glob];
                if (filter_actual_span <= 0 || filter_actual_span > this->N_glob) continue;

                // Convert filter (int*) to std::string
                std::string filter_str_version;
                filter_str_version.reserve(filter_actual_span);
                for (int bit_idx = 0; bit_idx < filter_actual_span; ++bit_idx) {
                    filter_str_version += (this->forms1Glob[f_idx][bit_idx] == 1 ? '1' : '0');
                }

                // Check if filter_str_version is a SUBSTRING of reduced_pattern_str_version
                if (reduced_pattern_str_version.find(filter_str_version) != std::string::npos) {
                    covering_filter_indices_for_this_pattern.push_back(f_idx);
                }
            }

            if (covering_filter_indices_for_this_pattern.size() == 1) {
                is_filter_protected[covering_filter_indices_for_this_pattern[0]] = true;
            }
        }
    }
    // else: all filters remain unprotected if no reduced patterns to guide protection.

    std::vector<int> unprotected_filter_actual_indices;
    for (int i = 0; i < this->Nform1_Glob; ++i) {
        if (!is_filter_protected[i]) {
            unprotected_filter_actual_indices.push_back(i);
        }
    }

    if (!unprotected_filter_actual_indices.empty()) {
        // Randomly choose one unprotected filter to remove
        std::uniform_int_distribution<> dist(0, unprotected_filter_actual_indices.size() - 1);
        int list_idx_to_remove = dist(rng_engine);
        int actual_filter_idx_in_forms1Glob = unprotected_filter_actual_indices[list_idx_to_remove];

        // std::cout << "      Removing unprotected filter at current list index " << actual_filter_idx_in_forms1Glob << ". (Filter: ";
        // for(int bit_idx=0; bit_idx < this->forms1Glob[actual_filter_idx_in_forms1Glob][this->N_glob]; ++bit_idx) {
        //     std::cout << this->forms1Glob[actual_filter_idx_in_forms1Glob][bit_idx];
        // }
        // std::cout << ")" << std::endl;

        // Perform removal from this->forms1Glob
        delete[] this->forms1Glob[actual_filter_idx_in_forms1Glob]; // Delete the filter's data
        // Shift subsequent pointers down
        for (int k = actual_filter_idx_in_forms1Glob; k < this->Nform1_Glob - 1; ++k) {
            this->forms1Glob[k] = this->forms1Glob[k + 1];
        }
        this->Nform1_Glob--; // Decrement count

        // If Nform1_Glob becomes 0, forms1Glob array itself might be deleted by the calling Forms object's destructor eventually.
        // No need to reallocate forms1Glob smaller here, just manage the count Nform1_Glob.
        // Setting the "last" (now out of bounds) pointer to nullptr can be good hygiene if array size fixed.
        if (this->Nform1_Glob > 0 && this->forms1Glob) {
            //this->forms1Glob[this->Nform1_Glob] = nullptr; // If forms1Glob was a fixed-size buffer
        }
        else if (this->Nform1_Glob == 0) {
            // If all filters were removed, forms1Glob might become an array of 0 pointers.
            // Its actual deletion (delete[] this->forms1Glob) handled by destructor or if reallocated.
            // For now, if Nform1_Glob is 0, the array of pointers might still exist but be empty.
        }

        return true; // A filter was removed
    }
    else {
        // std::cout << "      All " << this->Nform1_Glob << " filters are protected or Nform1_Glob is 0. Stable." << std::endl;
        return false; // No unprotected filter found (or no filters at all). Stable.
    }
}