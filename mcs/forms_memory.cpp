// forms_memory.cpp
#include "main.h"    // Crucial: Includes the Forms class declaration
#include <iostream>  // For std::cerr debugging (optional)
#include <cstring>   // For std::memcpy (though not ideal for bool**)

// --- Private Helper Methods Implementation ---
void Forms::clear_Vars_Glob() {
    if (Vars_Glob != nullptr) { // Check before dereferencing/deleting
        for (int i = 0; i < NVars_Glob; ++i) {
            if (Vars_Glob[i] != nullptr) { // Check inner pointers too
                delete[] Vars_Glob[i];
                // Vars_Glob[i] = nullptr; // Good practice, though not strictly needed before outer delete
            }
        }
        delete[] Vars_Glob;
        Vars_Glob = nullptr;
    }
    NVars_Glob = 0;
}

void Forms::clear_forms1Glob() {
    if (forms1Glob != nullptr) {
        for (int i = 0; i < Nform1_Glob; ++i) {
            if (forms1Glob[i] != nullptr) {
                delete[] forms1Glob[i];
            }
        }
        delete[] forms1Glob;
        forms1Glob = nullptr;
    }
    Nform1_Glob = 0;
    // form1Size is a parameter
}

void Forms::clear_reducedVars_Glob() {
    if (reducedVars_Glob != nullptr) {
        for (int i = 0; i < NReducedVars_Glob; ++i) {
            if (reducedVars_Glob[i] != nullptr) {
                delete[] reducedVars_Glob[i];
            }
        }
        delete[] reducedVars_Glob;
        reducedVars_Glob = nullptr;
    }
    NReducedVars_Glob = 0;
}

void Forms::clear_VarsDescr_Glob() {
    if (VarsDescr_Glob != nullptr) {
        delete[] VarsDescr_Glob;
        VarsDescr_Glob = nullptr;
    }
    // NVars_Glob is the count, reset when Vars_Glob is cleared
}


void Forms::clear_all_memory() {
    clear_Vars_Glob();
    clear_forms1Glob();
    clear_reducedVars_Glob();
    clear_VarsDescr_Glob();

    // Placeholder for other dynamic members if they are added and managed by Forms
    if (NForms_pos_Glob) { delete[] NForms_pos_Glob; NForms_pos_Glob = nullptr; }
    if (SeqPos) { delete[] SeqPos; SeqPos = nullptr; }
    if (forms1_pos_Glob) {
        // This is int***, needs careful nested deletion based on its structure
        // Example: if forms1_pos_Glob[i] is an array of int**
        // for (int i = 0; i < some_count_for_forms1_pos_Glob; ++i) {
        //     if (forms1_pos_Glob[i]) {
        //         for (int j = 0; j < NForms_pos_Glob[i] (hypothetical); ++j) {
        //             delete[] forms1_pos_Glob[i][j];
        //         }
        //         delete[] forms1_pos_Glob[i];
        //     }
        // }
        // delete[] forms1_pos_Glob;
        // forms1_pos_Glob = nullptr;
        // For now, since its allocation is not defined, just nulling:
        forms1_pos_Glob = nullptr; // Needs proper deallocation if used
    }
}

void Forms::deep_copy_members(const Forms& other) {
    // Copy parameters first, as they might be needed for allocations
    N_glob = other.N_glob;
    N_2_glob = other.N_2_glob;
    Nsovp1_Glob = other.Nsovp1_Glob;
    form1Size = other.form1Size; // Crucial for forms1Glob allocation size

    // Deep copy Vars_Glob
    NVars_Glob = other.NVars_Glob; // Copy count first
    if (other.Vars_Glob != nullptr && NVars_Glob > 0 && N_glob > 0) {
        Vars_Glob = new bool* [NVars_Glob];
        for (int i = 0; i < NVars_Glob; ++i) {
            Vars_Glob[i] = new bool[N_glob];
            for (int j = 0; j < N_glob; ++j) {
                Vars_Glob[i][j] = other.Vars_Glob[i][j];
            }
        }
    }
    else {
        Vars_Glob = nullptr;
        NVars_Glob = 0; // Ensure count is 0 if pointer is null
    }

    // Deep copy forms1Glob
    Nform1_Glob = other.Nform1_Glob; // Copy count first
    if (other.forms1Glob != nullptr && Nform1_Glob > 0 && form1Size > 0) {
        forms1Glob = new int* [Nform1_Glob];
        for (int i = 0; i < Nform1_Glob; ++i) {
            forms1Glob[i] = new int[form1Size]; // Use copied form1Size
            for (int j = 0; j < form1Size; ++j) {
                forms1Glob[i][j] = other.forms1Glob[i][j];
            }
        }
    }
    else {
        forms1Glob = nullptr;
        Nform1_Glob = 0;
    }

    // Deep copy reducedVars_Glob
    NReducedVars_Glob = other.NReducedVars_Glob; // Copy count first
    if (other.reducedVars_Glob != nullptr && NReducedVars_Glob > 0 && N_glob > 0) {
        reducedVars_Glob = new bool* [NReducedVars_Glob];
        for (int i = 0; i < NReducedVars_Glob; ++i) {
            reducedVars_Glob[i] = new bool[N_glob];
            for (int j = 0; j < N_glob; ++j) {
                reducedVars_Glob[i][j] = other.reducedVars_Glob[i][j];
            }
        }
    }
    else {
        reducedVars_Glob = nullptr;
        NReducedVars_Glob = 0;
    }

    // Deep copy VarsDescr_Glob
    if (other.VarsDescr_Glob != nullptr && NVars_Glob > 0) { // Assuming size matches NVars_Glob
        VarsDescr_Glob = new int[NVars_Glob];
        for (int i = 0; i < NVars_Glob; ++i) {
            VarsDescr_Glob[i] = other.VarsDescr_Glob[i];
        }
    }
    else {
        VarsDescr_Glob = nullptr;
    }

    // Initialize other pointers to nullptr, their deep copy would be more complex
    // and depends on how NForms_pos_Glob (if it's a count array) is structured.
    NForms_pos_Glob = nullptr;
    forms1_pos_Glob = nullptr;
    SeqPos = nullptr;
}


// --- Public Constructors / Destructor / Assignment Operator Implementation ---
Forms::Forms() :
    N_glob(0), N_2_glob(0), Nsovp1_Glob(0),
    NVars_Glob(0), Vars_Glob(nullptr),
    Nform1_Glob(0), form1Size(0), forms1Glob(nullptr),
    NReducedVars_Glob(0), reducedVars_Glob(nullptr),
    VarsDescr_Glob(nullptr),
    NForms_pos_Glob(nullptr), forms1_pos_Glob(nullptr), SeqPos(nullptr) // Initialize new members
{
    // std::cout << "Forms Default Constructor: All members initialized." << std::endl;
}

Forms::~Forms() {
    // std::cout << "Forms Destructor: Clearing all memory." << std::endl;
    clear_all_memory();
}

Forms::Forms(const Forms& other) : // Copy Constructor
    // Initialize all members to default before copying, crucial for RAII if deep_copy_members throws
    N_glob(0), N_2_glob(0), Nsovp1_Glob(0),
    NVars_Glob(0), Vars_Glob(nullptr),
    Nform1_Glob(0), form1Size(0), forms1Glob(nullptr),
    NReducedVars_Glob(0), reducedVars_Glob(nullptr),
    VarsDescr_Glob(nullptr),
    NForms_pos_Glob(nullptr), forms1_pos_Glob(nullptr), SeqPos(nullptr)
{
    // std::cout << "Forms Copy Constructor: Deep copying from other." << std::endl;
    deep_copy_members(other);
}

Forms& Forms::operator=(const Forms& other) { // Assignment Operator
    // std::cout << "Forms Assignment Operator: Handling assignment." << std::endl;
    if (this != &other) {       // Self-assignment check
        clear_all_memory();      // Clear existing resources in 'this' object
        deep_copy_members(other); // Copy resources from 'other' to 'this'
    }
    return *this;
}