#ifndef MAIN_H_SDF
#define MAIN_H_SDF

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <random> // For std::mt19937
#include <map>    // For FilterMapCollection

// Structures for the filter map (used by build_filter_map and main.cpp for searches)
struct MapEntryData {
    std::string masked_word;
    std::vector<int> occurrences;
};

struct FilterMapCollection {
    std::map<std::string, int> word_to_index_map; // Maps masked_word to index in 'entries_vector'
    std::vector<MapEntryData> entries_vector;      // Stores the actual MapEntryData
};

class Forms {
public:
    // --- Parameters & Counts ---
    int N_glob;         // N, typically N_PARAM_VAL. For an object holding filters, it's their max length. For reduced patterns, it's their length.
    int N_2_glob;       // K, target 1s in C(N,K) base patterns.
    int Nsovp1_Glob;    // M, target 1s in filters.

    int NVars_Glob;     // Count of patterns in Vars_Glob (full base patterns C(N,K))
    bool** Vars_Glob;   // Full base patterns (e.g., from create_mas1)

    int Nform1_Glob;    // Count of filters in forms1Glob (current MCS)
    int form1Size;      // N_glob + 1, full allocation size for each filter in forms1Glob
    int** forms1Glob;   // Filters (MCS)

    int NReducedVars_Glob; // Count of patterns in reducedVars_Glob
    bool** reducedVars_Glob; // Reduced base patterns (e.g., for Positional MCS Phase 2)

    // Original members from prompt, ensure initialized if used beyond Forms1_anal1
    int* VarsDescr_Glob; // Purpose not fully clear from all contexts, manage if used
    int* NForms_pos_Glob;
    int*** forms1_pos_Glob;
    int* SeqPos;

    // --- Constructors / Destructor / Assignment ---
    Forms();
    ~Forms();
    Forms(const Forms& other);            // Copy Constructor
    Forms& operator=(const Forms& other); // Assignment Operator

    // --- Methods ---
    // from forms1.cpp
    void create_mas1(); // Generates Vars_Glob
    double degInt(int n_power, double base_val);
    bool FindWord(std::ifstream& fin, std::string target_word);

    // from forms2.cpp
    void chetv_struct_Generation(const std::string& output_initial_forms_filename, const Forms& pattern_source);
    void Forms1_anal1_Global_Refinement(const std::string& input_initial_forms_filename, const Forms& global_pattern_source, const std::string& output_refined_forms_filename);
    void Read_Forms1(std::string fname); // Reads filters into this->forms1Glob

    // from forms_positional_logic.cpp
    void perform_phase1_pattern_reduction(const Forms& mcs_i_ref, const std::string& input_patterns_filename, const std::string& output_reduced_patterns_filename);
    bool perform_phase2_single_filter_prune_step(std::mt19937& rng_engine);

    // Utility
    void save_forms_to_file(const std::string& filename, const std::string& comment_prefix) const;
    std::string pattern_to_string(const bool* pattern, int length) const;
    void clear_reducedVars_Glob();

private:
    void clear_all_memory();
    void deep_copy_members(const Forms& other);

    void clear_Vars_Glob();
    void clear_forms1Glob();
    void clear_VarsDescr_Glob();
    // Add clear methods for NForms_pos_Glob, forms1_pos_Glob, SeqPos if they are dynamic

    // Helper for perform_phase1_pattern_reduction if reading from file internally
    bool Read_Patterns_FromFile(const std::string& filename, bool**& out_vars, int& out_nvars, int& out_n_val, int& out_k_val);
};

#endif // MAIN_H_SDF