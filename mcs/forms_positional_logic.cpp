// make_positions_phase1_2.cpp
#include "main.h"    // For Forms class declaration
#include <vector>
#include <string>
#include <set>
#include <algorithm> 
#include <iostream>
#include <fstream>
#include <random>    

// --- Helper to read patterns from file for perform_phase1_pattern_reduction ---
// This is now a private member of Forms, declared in main.h
// Refined Read_Patterns_FromFile logic
bool Forms::Read_Patterns_FromFile(const std::string& filename,
    bool**& out_vars, int& out_nvars,
    int& out_n_val, int& out_k_val) {
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "  ERROR (Read_Patterns_FromFile): Could not open pattern file: " << filename << std::endl;
        out_vars = nullptr; out_nvars = 0; out_n_val = 0; out_k_val = -1;
        return false;
    }

    std::string line;
    std::vector<std::string> pattern_lines_data;
    out_nvars = 0; out_n_val = -1; out_k_val = -1;
    int count_from_header = -1;
    bool in_pattern_block = false; // Flag to indicate we are past headers

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        if (!in_pattern_block && line[0] == '#') { // Comment or header line
            size_t n_pos = line.find("N_glob =");
            if (n_pos == std::string::npos) n_pos = line.find("N=");
            if (n_pos != std::string::npos) {
                size_t val_start = line.find("=", n_pos) + 1;
                try { out_n_val = std::stoi(line.substr(val_start)); }
                catch (...) {}
            }
            size_t k_pos = line.find("N_2_glob (K) =");
            if (k_pos != std::string::npos) {
                size_t val_start = k_pos + strlen("N_2_glob (K) =");
                try { out_k_val = std::stoi(line.substr(val_start)); }
                catch (...) {}
            }
            continue;
        }

        if (!in_pattern_block) { // Look for the count line only if not already in pattern block
            size_t count_pos = line.find("NVars_Glob =");
            if (count_pos == std::string::npos) count_pos = line.find("NReducedVars_Glob =");
            if (count_pos != std::string::npos) {
                size_t val_start = line.find("=", count_pos) + 1;
                try { count_from_header = std::stoi(line.substr(val_start)); }
                catch (...) {}
                in_pattern_block = true; // After count line, expect patterns or EOF
                continue;
            }
        }

        // If we are past the header/count lines, or if no count line was found yet,
        // and it's not a comment, treat as a pattern.
        // However, to be safer, once 'in_pattern_block' is true, ONLY accept patterns.
        if (in_pattern_block && line[0] == '#') {
            std::cerr << "  Warning (Read_Patterns_FromFile): Encountered comment line after NReducedVars_Glob line. Stopping pattern read." << std::endl;
            break; // Stop if comments appear after data count
        }

        if (in_pattern_block || (line[0] != '#' && line.find("Vars_Glob =") == std::string::npos)) { // Add check for Vars_Glob
            pattern_lines_data.push_back(line);
        }
        else if (!in_pattern_block && line[0] != '#') {
            // This case means we haven't seen a count line yet, and this line is not a comment.
            // It could be an early pattern if the file format is loose.
            // Or it could be an unexpected line. For now, let's assume it's a pattern.
            // This is risky if the file format isn't strict.
            std::cerr << "  Warning (Read_Patterns_FromFile): Reading line as pattern before count header: " << line << std::endl;
            pattern_lines_data.push_back(line);
        }


    } // end while getline
    fin.close();

    // ... (rest of the function: process pattern_lines_data, check counts, allocate) ...
    // The logic for determining out_n_val and checking consistency remains important.

    if (pattern_lines_data.empty()) {
        out_nvars = (count_from_header > 0) ? count_from_header : 0;
        if (out_nvars > 0) { // Header says patterns, but none found.
            std::cerr << "  ERROR (Read_Patterns_FromFile): Header indicated " << out_nvars
                << " patterns, but no pattern lines were read from " << filename << std::endl;
            out_vars = nullptr; out_n_val = (out_n_val != -1 ? out_n_val : 0); return false; // Indicate error
        }
        std::cout << "  Read_Patterns_FromFile: No pattern data lines found or expected in " << filename << std::endl;
        out_vars = nullptr; out_n_val = (out_n_val != -1 ? out_n_val : 0);
        return true;
    }

    if (out_n_val == -1) { // N not found in header
        out_n_val = pattern_lines_data[0].length();
    }
    else { // N found in header, verify
        if (static_cast<int>(pattern_lines_data[0].length()) != out_n_val) {
            std::cerr << "  Warning (Read_Patterns_FromFile): Header N " << out_n_val
                << " mismatches pattern length " << pattern_lines_data[0].length()
                << " in " << filename << ". Using actual pattern line length: " << pattern_lines_data[0].length() << std::endl;
            out_n_val = pattern_lines_data[0].length();
        }
    }

    out_nvars = pattern_lines_data.size();
    if (count_from_header != -1 && count_from_header != out_nvars) {
        std::cerr << "  Warning (Read_Patterns_FromFile): Header count " << count_from_header
            << " mismatches actual pattern lines " << out_nvars << " in " << filename << ". Using line count." << std::endl;
    }

    out_vars = new bool* [out_nvars];
    for (int i = 0; i < out_nvars; ++i) {
        if (pattern_lines_data[i].length() != static_cast<size_t>(out_n_val)) {
            std::cerr << "  ERROR (Read_Patterns_FromFile): Inconsistent pattern length in " << filename
                << ". Expected " << out_n_val << ", got " << pattern_lines_data[i].length() << " for pattern " << i
                << " Content: [" << pattern_lines_data[i] << "]" << std::endl;
            for (int j = 0; j < i; ++j) delete[] out_vars[j];
            delete[] out_vars; out_vars = nullptr; out_nvars = 0; return false;
        }
        out_vars[i] = new bool[out_n_val];
        for (int j = 0; j < out_n_val; ++j) {
            out_vars[i][j] = (pattern_lines_data[i][j] == '1');
        }
    }
    std::cout << "  Read_Patterns_FromFile: Read " << out_nvars << " patterns of length " << out_n_val << " from " << filename << ". K from header (if any): " << out_k_val << std::endl;
    return true;
}


// --- Phase 1: Pattern Reduction for Positional MCS (File-based cumulative) ---
void Forms::perform_phase1_pattern_reduction(
    const Forms& mcs_i_ref,
    const std::string& input_patterns_filename,
    const std::string& output_reduced_patterns_filename)
{
    bool** source_patterns_data = nullptr;
    int source_nvars = 0;
    int source_n_val = 0; // N of the patterns being read
    int source_k_val = 0; // K of the patterns being read (if available)

    if (!Read_Patterns_FromFile(input_patterns_filename, source_patterns_data, source_nvars, source_n_val, source_k_val)) {
        std::cerr << "  ERROR (perform_phase1): Failed to read source patterns from " << input_patterns_filename << std::endl;
        clear_reducedVars_Glob(); return;
    }

    this->N_glob = source_n_val; // 'this' object's N_glob should reflect length of patterns it will hold in reducedVars_Glob
    // this->N_2_glob could be source_k_val or context from mcs_i_ref or initial setup.

    std::cout << "  Phase 1 (Positional): Reducing " << source_nvars << " patterns from " << input_patterns_filename
        << " against MCS_i (" << mcs_i_ref.Nform1_Glob << " filters)..." << std::endl;
    clear_reducedVars_Glob();

    if (source_patterns_data == nullptr || source_nvars == 0) {
        std::cout << "    Phase 1: No source patterns from file to reduce. NReducedVars_Glob will be 0." << std::endl;
        if (source_patterns_data) delete[] source_patterns_data; // Should be null if nvars is 0
        // Save empty output file
        std::ofstream fout_empty(output_reduced_patterns_filename);
        if (fout_empty.is_open()) {
            fout_empty << "# No source patterns provided or read from " << input_patterns_filename << std::endl;
            fout_empty << "NReducedVars_Glob = 0" << std::endl;
            fout_empty.close();
        }
        return;
    }

    std::vector<bool*> kept_pattern_pointers_vec;
    for (int pat_idx = 0; pat_idx < source_nvars; ++pat_idx) {
        const bool* current_basic_pattern_bits = source_patterns_data[pat_idx];
        bool pattern_starts_with_an_mcs_i_filter = false;
        if (mcs_i_ref.Nform1_Glob > 0 && mcs_i_ref.forms1Glob != nullptr) {
            for (int f_idx = 0; f_idx < mcs_i_ref.Nform1_Glob; ++f_idx) {
                const int* filter_from_mcs_i_bits = mcs_i_ref.forms1Glob[f_idx];
                int filter_actual_span = mcs_i_ref.forms1Glob[f_idx][mcs_i_ref.N_glob];
                if (filter_actual_span <= 0 || filter_actual_span > source_n_val) continue;
                bool current_filter_matches_prefix = true;
                if (source_n_val > 0 && filter_from_mcs_i_bits[0] == 1 && !current_basic_pattern_bits[0]) {
                    current_filter_matches_prefix = false;
                }
                else if (source_n_val > 0) {
                    for (int k = 0; k < filter_actual_span; ++k) {
                        if (filter_from_mcs_i_bits[k] == 1 && !current_basic_pattern_bits[k]) {
                            current_filter_matches_prefix = false; break;
                        }
                    }
                }
                else { current_filter_matches_prefix = false; }
                if (current_filter_matches_prefix) { pattern_starts_with_an_mcs_i_filter = true; break; }
            }
        }
        if (!pattern_starts_with_an_mcs_i_filter) {
            bool* pattern_copy = new bool[source_n_val];
            for (int k = 0; k < source_n_val; ++k) pattern_copy[k] = current_basic_pattern_bits[k];
            kept_pattern_pointers_vec.push_back(pattern_copy);
        }
    }

    this->NReducedVars_Glob = kept_pattern_pointers_vec.size();
    if (this->NReducedVars_Glob > 0) {
        this->reducedVars_Glob = new bool* [this->NReducedVars_Glob];
        for (size_t i = 0; i < this->NReducedVars_Glob; ++i) {
            this->reducedVars_Glob[i] = kept_pattern_pointers_vec[i];
        }
    } // else reducedVars_Glob is already nullptr from clear_reducedVars_Glob

    std::cout << "    Phase 1 (Positional): Reduction complete. " << this->NReducedVars_Glob << " patterns remain." << std::endl;
    std::ofstream fout_reduced(output_reduced_patterns_filename);
    if (!fout_reduced.is_open()) {
        std::cerr << "    ERROR (perform_phase1): Could not open " << output_reduced_patterns_filename << std::endl;
    }
    else {
        fout_reduced << "# Reduced patterns. Based on MCS_i (M=" << mcs_i_ref.Nsovp1_Glob
            << "). Input: " << input_patterns_filename << std::endl;
        fout_reduced << "# N_glob (length of these patterns) = " << this->N_glob << std::endl;
        fout_reduced << "NReducedVars_Glob = " << this->NReducedVars_Glob << std::endl;
        for (int i = 0; i < this->NReducedVars_Glob; ++i) {
            for (int j = 0; j < this->N_glob; ++j) {
                fout_reduced << (this->reducedVars_Glob[i][j] ? '1' : '0');
            }
            fout_reduced << std::endl;
        }
        fout_reduced.close();
        std::cout << "    Reduced patterns saved to " << output_reduced_patterns_filename << std::endl;
    }
    if (source_patterns_data != nullptr) {
        for (int i = 0; i < source_nvars; ++i) delete[] source_patterns_data[i];
        delete[] source_patterns_data;
    }
}

// --- Phase 2: Single Filter Pruning Step for Positional MCS ---
bool Forms::perform_phase2_single_filter_prune_step(std::mt19937& rng_engine) {
    if (this->Nform1_Glob == 0) return false;

    bool all_filters_unprotected_due_to_no_reduced_patterns =
        (this->reducedVars_Glob == nullptr || this->NReducedVars_Glob == 0);

    std::vector<bool> is_filter_protected(this->Nform1_Glob, false);

    if (!all_filters_unprotected_due_to_no_reduced_patterns) {
        for (int i = 0; i < this->NReducedVars_Glob; ++i) {
            std::string reduced_pattern_str_version;
            reduced_pattern_str_version.reserve(this->N_glob); // N_glob of reduced patterns
            for (int bit_idx = 0; bit_idx < this->N_glob; ++bit_idx) {
                reduced_pattern_str_version += (this->reducedVars_Glob[i][bit_idx] ? '1' : '0');
            }
            std::vector<int> covering_filter_indices_for_this_pattern;
            for (int f_idx = 0; f_idx < this->Nform1_Glob; ++f_idx) {
                if (this->forms1Glob[f_idx] == nullptr) continue;
                int filter_actual_span = this->forms1Glob[f_idx][this->N_glob]; // N_glob of filter storage
                if (filter_actual_span <= 0 || filter_actual_span > this->N_glob) continue; // N_glob of reduced patterns
                std::string filter_str_version;
                filter_str_version.reserve(filter_actual_span);
                for (int bit_idx = 0; bit_idx < filter_actual_span; ++bit_idx) {
                    filter_str_version += (this->forms1Glob[f_idx][bit_idx] == 1 ? '1' : '0');
                }
                if (reduced_pattern_str_version.find(filter_str_version) != std::string::npos) {
                    covering_filter_indices_for_this_pattern.push_back(f_idx);
                }
            }
            if (covering_filter_indices_for_this_pattern.size() == 1) {
                is_filter_protected[covering_filter_indices_for_this_pattern[0]] = true;
            }
        }
    }

    std::vector<int> unprotected_filter_actual_indices;
    for (int i = 0; i < this->Nform1_Glob; ++i) {
        if (!is_filter_protected[i]) {
            unprotected_filter_actual_indices.push_back(i);
        }
    }

    if (!unprotected_filter_actual_indices.empty()) {
        std::uniform_int_distribution<> dist(0, unprotected_filter_actual_indices.size() - 1);
        int list_idx_to_remove = dist(rng_engine);
        int actual_filter_idx_in_forms1Glob = unprotected_filter_actual_indices[list_idx_to_remove];

        delete[] this->forms1Glob[actual_filter_idx_in_forms1Glob];
        for (int k = actual_filter_idx_in_forms1Glob; k < this->Nform1_Glob - 1; ++k) {
            this->forms1Glob[k] = this->forms1Glob[k + 1];
        }
        this->Nform1_Glob--;
        if (this->Nform1_Glob == 0) { // All filters removed
            delete[] this->forms1Glob; // Delete the outer array too
            this->forms1Glob = nullptr;
        }
        else {
            // Optional: Set the last (now unused) pointer to nullptr if not reallocating smaller
            // this->forms1Glob[this->Nform1_Glob] = nullptr; 
        }
        return true;
    }
    else {
        return false;
    }
}