// forms1.cpp
#include "main.h"    // For Forms class declaration
#include <fstream>   // For std::ofstream
#include <iostream>  // For std::cerr, std::cout
#include <vector>    // For std::vector (used internally for permutations if that method is chosen)
#include <algorithm> // For std::next_permutation / std::prev_permutation
#include <cmath>     // For pow
#include <iomanip>   // For std::setw, std::flush

// Helper function to calculate base_val^n_power (integer power)
double Forms::degInt(int n_power, double base_val) {
    double result = 1.0;
    for (int i = 0; i < n_power; i++) {
        result *= base_val;
    }
    return result;
}

// FindWord
bool Forms::FindWord(std::ifstream& fin, std::string target_word) {
    char char_from_file;
    size_t match_idx;

    while (fin.get(char_from_file)) {
        if (char_from_file == target_word[0]) {
            match_idx = 1;
            while (match_idx < target_word.length() && fin.get(char_from_file)) {
                if (char_from_file == target_word[match_idx]) {
                    match_idx++;
                }
                else {
                    if (char_from_file == target_word[0]) { // Check if the mismatch char is start of new target
                        std::streambuf* pbuf = fin.rdbuf();
                        pbuf->sungetc(); //fin.unget(); is C++11, for older use rdbuf
                    }
                    break;
                }
            }
            if (match_idx == target_word.length()) {
                return true;
            }
        }
    }
    return false;
}

// Builds the basic binary patterns (C(N_glob, N_2_glob)) starting with '1'.
// Populates this->Vars_Glob and this->NVars_Glob.
// Saves patterns to "tavnits.txt".
void Forms::create_mas1() {
    std::cout << "create_mas1: Generating basic patterns (N=" << N_glob << ", K=" << N_2_glob << ")..." << std::endl;

    if (N_glob <= 0 || N_2_glob <= 0 || N_2_glob > N_glob) {
        std::cerr << "  ERROR (create_mas1): Invalid N_glob (" << N_glob
            << ") or N_2_glob (" << N_2_glob << "). Cannot generate patterns." << std::endl;
        this->NVars_Glob = 0;
        this->Vars_Glob = nullptr;
        return;
    }
    if (N_2_glob == 0) { // No ones to place
        std::cerr << "  Warning (create_mas1): N_2_glob is 0. No patterns with ones will be generated." << std::endl;
        // If N_2_glob is 0, but patterns must start with '1', this is impossible.
        // The prompt specified "start with a '1'", so K must be at least 1.
        // The N_2_glob > 0 check above handles this.
    }


    clear_Vars_Glob(); // Clear any existing Vars_Glob
    if (this->VarsDescr_Glob) {
        delete[] this->VarsDescr_Glob;
        this->VarsDescr_Glob = nullptr;
    }

    std::vector<bool*> collected_patterns_vec;

    // --- Method: Using std::prev_permutation (more direct for C(N,K) starting with '1') ---
    // This is generally more efficient than iterating 2^N.
    std::cout << "  Generating patterns using permutations (N=" << N_glob << ", K=" << N_2_glob << ")..." << std::endl;

    if (N_2_glob < 1) { // Must have at least one '1' to fix the first bit
        std::cerr << "  ERROR (create_mas1): N_2_glob (K) must be at least 1 for patterns starting with '1'." << std::endl;
        this->NVars_Glob = 0; this->Vars_Glob = nullptr; return;
    }
    if (N_glob - 1 < N_2_glob - 1) { // Not enough remaining spots for remaining K-1 ones
        std::cerr << "  ERROR (create_mas1): Not enough spots (N-1) for remaining K-1 ones." << std::endl;
        this->NVars_Glob = 0; this->Vars_Glob = nullptr; return;
    }


    std::string p_str_template(N_glob - 1, '0'); // Template for N-1 bits
    if (N_2_glob - 1 > 0) { // If more than one '1' is needed overall
        for (int i = 0; i < N_2_glob - 1; ++i) {
            p_str_template[i] = '1'; // Place K-1 ones at the start of the N-1 template
        }
    }
    // Example: N=5, K=3. N-1=4, K-1=2. p_str_template = "1100"
    // Example: N=5, K=1. N-1=4, K-1=0. p_str_template = "0000"

    std::sort(p_str_template.begin(), p_str_template.end(), std::greater<char>()); // Start with "11...00.."

    do {
        bool* new_pattern_alloc = new bool[N_glob];
        new_pattern_alloc[0] = true; // First bit is always '1'
        for (int j = 0; j < N_glob - 1; ++j) {
            new_pattern_alloc[j + 1] = (p_str_template[j] == '1');
        }
        collected_patterns_vec.push_back(new_pattern_alloc);

        if (collected_patterns_vec.size() % 50000 == 0) { // Log progress
            std::cout << "\r    Collected " << collected_patterns_vec.size() << " patterns..." << std::flush;
        }
    } while (std::prev_permutation(p_str_template.begin(), p_str_template.end()));

    std::cout << "\n  Permutation generation complete. Total patterns: " << collected_patterns_vec.size() << std::endl;
    // END OF PERMUTATION METHOD

    this->NVars_Glob = collected_patterns_vec.size();
    if (this->NVars_Glob > 0) {
        this->Vars_Glob = new bool* [this->NVars_Glob];
        for (size_t i = 0; i < this->NVars_Glob; ++i) { // Use size_t for vector index
            this->Vars_Glob[i] = collected_patterns_vec[i];
        }
    }
    else {
        this->Vars_Glob = nullptr;
    }
    // collected_patterns_vec goes out of scope, its memory is managed (vector itself).
    // The bool* it contained are now owned by Vars_Glob.

    std::ofstream fout("tavnits.txt");
    if (!fout.is_open()) {
        std::cerr << "  ERROR (create_mas1): Could not open tavnits.txt for writing." << std::endl;
    }
    else {
        fout << "# Basic patterns: N_glob = " << N_glob
             << ", N_2_glob (K) = " << N_2_glob 
             << ", NVars_Glob = " << this->NVars_Glob << std::endl;
        for (int i = 0; i < this->NVars_Glob; ++i) {
            // fout << "\n"; // Original had newline before each pattern
            for (int j = 0; j < N_glob; ++j) {
                fout << (this->Vars_Glob[i][j] ? '1' : '0');
            }
            fout << std::endl; // Newline after each pattern
        }
        fout.close();
        std::cout << "  Basic patterns saved to tavnits.txt" << std::endl;
    }
    std::cout << "create_mas1: Finished. Generated " << this->NVars_Glob << " basic patterns." << std::endl;
}