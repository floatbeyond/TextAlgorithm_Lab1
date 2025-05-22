// forms_memory.cpp
#include "main.h"    // Crucial: Includes the Forms class declaration
#include <iostream>
#include <cstring>   // For std::memcpy in some deep copy scenarios (though not ideal for bool**)

// --- Private Helper Methods Implementation ---
void Forms::clear_Vars_Glob() {
    if (Vars_Glob != nullptr) {
        for (int i = 0; i < NVars_Glob; ++i) {
            if (Vars_Glob[i] != nullptr) {
                delete[] Vars_Glob[i];
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
    // form1Size is a parameter, not cleared here
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
}

void Forms::clear_all_memory() {
    clear_Vars_Glob();
    clear_forms1Glob();
    clear_reducedVars_Glob();
    clear_VarsDescr_Glob();

    if (NForms_pos_Glob) { delete[] NForms_pos_Glob; NForms_pos_Glob = nullptr; }
    if (SeqPos) { delete[] SeqPos; SeqPos = nullptr; }
    if (forms1_pos_Glob) {
        // Proper deallocation for int*** forms1_pos_Glob depends on its allocation structure.
        // Assuming it's not deeply used or managed elsewhere for now.
        // If it were Forms*** forms1_pos_Glob[outer][inner_count_from_NForms_pos_Glob]
        // then a nested loop deleting Forms objects then the arrays of pointers.
        // For now, just to prevent a trivial leak if it were simply new int**:
        // delete[] forms1_pos_Glob; // This is wrong if it's int*** and not just int**
        forms1_pos_Glob = nullptr; // Placeholder - requires actual deallocation logic if used
    }
}

void Forms::deep_copy_members(const Forms& other) {
    N_glob = other.N_glob;
    N_2_glob = other.N_2_glob;
    Nsovp1_Glob = other.Nsovp1_Glob;
    form1Size = other.form1Size;

    NVars_Glob = other.NVars_Glob;
    if (other.Vars_Glob != nullptr && NVars_Glob > 0 && other.N_glob > 0) { // Use other.N_glob for source dimension
        Vars_Glob = new bool* [NVars_Glob];
        for (int i = 0; i < NVars_Glob; ++i) {
            Vars_Glob[i] = new bool[other.N_glob]; // Allocate with other.N_glob
            for (int j = 0; j < other.N_glob; ++j) {
                Vars_Glob[i][j] = other.Vars_Glob[i][j];
            }
        }
    }
    else {
        Vars_Glob = nullptr;
        NVars_Glob = 0;
    }

    Nform1_Glob = other.Nform1_Glob;
    if (other.forms1Glob != nullptr && Nform1_Glob > 0 && other.form1Size > 0) { // Use other.form1Size
        forms1Glob = new int* [Nform1_Glob];
        for (int i = 0; i < Nform1_Glob; ++i) {
            forms1Glob[i] = new int[other.form1Size]; // Allocate with other.form1Size
            for (int j = 0; j < other.form1Size; ++j) {
                forms1Glob[i][j] = other.forms1Glob[i][j];
            }
        }
    }
    else {
        forms1Glob = nullptr;
        Nform1_Glob = 0;
    }

    NReducedVars_Glob = other.NReducedVars_Glob;
    // N_glob for reducedVars_Glob should be consistent with other.N_glob if source was Vars_Glob,
    // or specific N_glob of other if reducedVars_Glob had different length patterns.
    // Assuming reduced patterns have same length as original N_glob from 'other' if 'other' holds them.
    int reduced_pattern_length = other.N_glob; // Or a specific N for reducedVars if different
    // If 'other' itself has reducedVars, use other.N_glob related to its reducedVars

    if (other.reducedVars_Glob != nullptr && NReducedVars_Glob > 0 && reduced_pattern_length > 0) {
        reducedVars_Glob = new bool* [NReducedVars_Glob];
        for (int i = 0; i < NReducedVars_Glob; ++i) {
            reducedVars_Glob[i] = new bool[reduced_pattern_length];
            for (int j = 0; j < reduced_pattern_length; ++j) {
                reducedVars_Glob[i][j] = other.reducedVars_Glob[i][j];
            }
        }
    }
    else {
        reducedVars_Glob = nullptr;
        NReducedVars_Glob = 0;
    }

    if (other.VarsDescr_Glob != nullptr && NVars_Glob > 0) { // Assuming size matches NVars_Glob
        VarsDescr_Glob = new int[NVars_Glob];
        for (int i = 0; i < NVars_Glob; ++i) {
            VarsDescr_Glob[i] = other.VarsDescr_Glob[i];
        }
    }
    else {
        VarsDescr_Glob = nullptr;
    }

    NForms_pos_Glob = nullptr; // Deep copy not implemented for these complex members
    forms1_pos_Glob = nullptr;
    SeqPos = nullptr;
}

Forms::Forms() :
    N_glob(0), N_2_glob(0), Nsovp1_Glob(0),
    NVars_Glob(0), Vars_Glob(nullptr),
    Nform1_Glob(0), form1Size(0), forms1Glob(nullptr),
    NReducedVars_Glob(0), reducedVars_Glob(nullptr),
    VarsDescr_Glob(nullptr),
    NForms_pos_Glob(nullptr), forms1_pos_Glob(nullptr), SeqPos(nullptr) {
}

Forms::~Forms() {
    clear_all_memory();
}

Forms::Forms(const Forms& other) :
    N_glob(0), N_2_glob(0), Nsovp1_Glob(0), // Initialize primitives
    NVars_Glob(0), Vars_Glob(nullptr),
    Nform1_Glob(0), form1Size(0), forms1Glob(nullptr),
    NReducedVars_Glob(0), reducedVars_Glob(nullptr),
    VarsDescr_Glob(nullptr),
    NForms_pos_Glob(nullptr), forms1_pos_Glob(nullptr), SeqPos(nullptr)
{
    deep_copy_members(other);
}

Forms& Forms::operator=(const Forms& other) {
    if (this != &other) {
        clear_all_memory();
        deep_copy_members(other);
    }
    return *this;
}