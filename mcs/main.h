// main.h
#ifndef MAIN_H_SDF
#define MAIN_H_SDF

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <random> // For std::mt19937

// Forward declaration not needed if Forms is fully defined here
// class Forms;

// ... (FilterMapCollection structs if they stay here, or move to make_positions.h)

class Forms {
public:
    // --- Parameters & Counts ---
    int N_glob;
    int N_2_glob;
    int Nsovp1_Glob;

    int NVars_Glob;
    bool** Vars_Glob;

    int Nform1_Glob;
    int form1Size;
    int** forms1Glob;

    int NReducedVars_Glob;
    bool** reducedVars_Glob;

    int* VarsDescr_Glob; // Added based on your original code
    // These were in your original main.h, ensure they are initialized if used
    int* NForms_pos_Glob;
    int*** forms1_pos_Glob;
    int* SeqPos;


    // --- Constructors / Destructor / Assignment ---
    Forms();
    ~Forms();
    Forms(const Forms& other);
    Forms& operator=(const Forms& other);

    // --- Methods ---
    // from forms1.cpp
    void create_mas1();
    double degInt(int n, double q); // Helper for create_mas1
    bool FindWord(std::ifstream& fin, std::string aaa); // Helper for Read_Forms1

    // from forms2.cpp (or new files)
    void chetv_struct_Generation(const std::string& output_initial_forms_filename,
                                 const Forms& pattern_source);
    void Forms1_anal1_Global_Refinement(const std::string& input_forms_filename_to_read,
                                        const Forms& global_pattern_source,
                                        const std::string& output_refined_forms_filename);
    void Read_Forms1(std::string fname); // Helper for Forms1_anal1_Global_Refinement

    // For Positional MCS
    void perform_phase1_pattern_reduction(const Forms& mcs_i_const_ref,
                                          const Forms& global_pattern_source,
                                          const std::string& output_reduced_patterns_filename);
    bool perform_phase2_single_filter_prune_step(std::mt19937& rng_engine);

    // Utility
    void save_forms_to_file(const std::string& filename, const std::string& comment_prefix) const;
    // load_forms_from_file might be redundant if Read_Forms1 is the main loader

private:
    void clear_all_memory();
    void deep_copy_members(const Forms& other);

    // Optional allocation helpers (could be inlined or used internally)
    // void allocate_Vars_Glob(int count, int length);
    // void allocate_forms1Glob(int count, int full_length);
    // void allocate_reducedVars_Glob(int count, int length);

    // Specific clear functions (called by clear_all_memory)
    void clear_Vars_Glob();
    void clear_forms1Glob();
    void clear_reducedVars_Glob();
    void clear_VarsDescr_Glob(); // If dynamically allocated
    // Add clear methods for NForms_pos_Glob, forms1_pos_Glob, SeqPos if they are dynamic
};

#endif // MAIN_H_SDF