# Advanced Text Searching Algorithms: MCS Comparison Framework

## 1. Overview

This C++ project implements and compares several text searching algorithms, with a primary focus on variations of the Minimal Cover Set (MCS) approach. It provides a framework for:
*   Generating synthetic text data and query sets.
*   Implementing a naive (brute-force) search algorithm as a baseline.
*   Implementing a "Usual" MCS-based search algorithm.
*   Implementing a more sophisticated "Positional-Associated" MCS algorithm that iteratively refines filter sets based on their performance at specific alignments.
*   Evaluating these algorithms based on preprocessing time, search time, and the number of matches found.

## 3. Algorithms Implemented

### 3.1. Naive Search
A straightforward brute-force search that slides the query across the entire text, calculating similarity at each position. Serves as a correctness baseline and a performance benchmark for worst-case scenarios.

### 3.2. "Usual" MCS (MCS_0)
*   **Filter Generation (`chetv_struct_Generation`):**
    *   Generates an initial set of filters (MCS_0) from a source of "basic patterns" (typically all C(N,K) binary strings starting with '1').
    *   A basic pattern is considered "covered" if any existing filter in the MCS is a **substring** of it.
    *   If not covered, a candidate filter (prefix of the basic pattern with `M` ones, starting with '1') is derived and added if unique.
*   **Optional Global Refinement (`Forms1_anal1_Global_Refinement`):**
    *   An optional, intensive two-pass refinement of MCS_0 using all basic patterns.
    *   Pass 1: Marks filters as essential if they are the sole **substring cover** for any basic pattern.
    *   Pass 2: Iteratively removes non-essential filters, promoting others if they become sole covers, until a stable, minimal set is achieved.
*   **Filter Map Building (`build_filter_map`):**
    *   Slides each filter from the (refined) MCS_0 over the entire text.
    *   For each alignment, a "masked text segment" is created (text characters are kept for filter '1's, '_' for filter '0's).
    *   These masked segments and their start positions in the text are stored in a hash map for fast lookup.
*   **Search (`search_usual_mcs`):**
    *   For a given query, each filter from MCS_0 is slid across the query.
    *   At each alignment, a "masked query segment" is generated.
    *   This masked query segment is looked up in the precomputed filter map.
    *   Hits from the map provide candidate text positions, which are then fully verified against the original query using similarity calculation.

### 3.3. "Positional-Associated" MCS
This is an iterative algorithm that generates a sequence of MCS sets (`MCS_0, MCS_1, ..., MCS_p_final`), where `MCS_p` is intended to be more effective for features/patterns starting around position `p`.

For each position `p` (from 0 to `N-1`):
*   **`mcs_iter_current`** starts as the refined `MCS_p-1` (or the globally refined `MCS_0` for `p=0`).
*   **Step 1: Identify "Problem Patterns" for `MCS_p` Refinement (Image's Stage 2.1 inspired)**
    *   Function: `perform_positional_pattern_pruning_phase_a` (called with `position_index = 0` against `f_global_pattern_holder`).
    *   Identifies patterns from the *full global set* that are **NOT prefixed (at index 0)** by any filter in the current `mcs_iter_current`. These are `problem_patterns_for_mcs_p_refinement`.
*   **Step 2: Refine `mcs_iter_current` to become `MCS_p` (Image's Stages 2.2-2.4 inspired)**
    *   Function: `mcs_iter_current.refine_mcs_iteratively_phase_b`.
    *   This function iteratively refines `mcs_iter_current` using the `problem_patterns_for_mcs_p_refinement` as guiding patterns (or all global patterns if the problem set is empty but `mcs_iter_current` has filters).
    *   Protection Criterion: A filter is protected if it's the sole **substring cover** (filter '0's are wildcards) for any pattern in the guiding set. Unprotected filters are removed one by one (randomly if multiple) until all remaining are protected.
    *   The refined `mcs_iter_current` becomes the final `MCS_p` and is stored.
*   **Step 3: Prepare "Flowing Patterns" (`P_p+1`) for the *next* position's cycle**
    *   Function: `perform_positional_pattern_pruning_phase_a` (called with `position_index = p`).
    *   Uses the *just refined* `MCS_p` to prune the set of patterns `P_p` (that flowed into the current iteration `p`).
    *   Pruning is based on **positional prefix match at index `p`**.
    *   The output `P_p+1` flows to the next iteration.
*   **Search (`search_positional_mcs` - Advanced Version):**
    *   For a query, and for each internal query position `q_pos` from `0` to `W-1`:
        *   Selects `MCS_q_pos` from the generated sequence.
        *   Applies filters from `MCS_q_pos` to the query segment starting at `q_pos` (i.e., `query.substr(q_pos, filter_span)`).
        *   Looks up the masked query sub-segment in the `map_from_mcs0`.
        *   Candidate text positions are derived and then verified.

## 4. Key Parameters (Configured in `main.cpp`)

*   `X_PARAM`: A parameter influencing `N` and `K`.
*   `N_PARAM_VAL` (W): The length of queries and text windows.
*   `K_PARAM_VAL` (K): The number of '1's in the C(N,K) basic binary patterns.
*   `MATCH_SIMILARITY_THRESHOLD`: The minimum similarity (0.0 to 1.0) for a query and text segment to be considered a match.
*   `M_VALUES_PARAM`: A list of `M` values to test (number of '1's in filters derived by `chetv_struct_Generation`).
*   `TEXT_SIZE_PARAM`: Length of the randomly generated main text.
*   `NUM_QUERIES_PARAM`: Number of queries to extract from the main text.
*   `ALPHABET_SIZE_Y_PARAM`: Size of the alphabet for text generation (e.g., 4 for DNA-like, 20 for protein-like). Automatically determined or can be overridden.

## 5. Project Structure (Key Files)

*   `main.cpp`: Orchestrates the entire experiment, including parameter setup, text/query generation, algorithm execution, and results summarization.
*   `main.h`: Header file for `Forms` class, `FilterMapCollection` structs, and declarations for shared utility functions or global parameters.
*   `forms1.cpp`: Contains `Forms::create_mas1` (basic C(N,K) pattern generation) and other utility methods of the `Forms` class.
*   `forms2.cpp`: Contains `Forms::chetv_struct_Generation` (initial MCS filter generation), `Forms::Forms1_anal1_Global_Refinement` (optional global MCS refinement), and I/O methods like `Forms::Read_Forms1`, `Forms::save_forms_to_file`.
*   `forms_memory.cpp`: Contains memory management for the `Forms` class (constructors, destructor, `clear_*` methods, `deep_copy_members`).
*   `make_positions.cpp` / `make_positions.h`: Contains `build_filter_map` for creating the text filter map. Also likely contains `perform_positional_pattern_pruning_phase_a` (or it might be in a `forms_positional_logic.cpp`).
*   `forms_positional_logic.cpp` (Assumed or integrate into `forms2.cpp`): Would contain `Forms::refine_mcs_iteratively_phase_b`.
*   `tavnits.txt`: Output file for basic C(N,K) patterns.
*   `main_text_Y*.txt`: The generated main text file.
*   `*_matches_*.txt`: Output files detailing match positions for each query by each algorithm.
*   `Usual_MCS_*.txt`, `Pos_MCS_*.txt`: Intermediate files storing generated filters for different MCS stages.
*   `Usual_Map_*.txt`: The filter map built from MCS_0.
*   `experiment_conclusions.txt`: Summary of performance results and parameters.

## 6. Prerequisites

*   A C++11 (or later) compatible compiler (e.g., g++, Clang, MSVC).
*   Standard C++ Library.
*   The project uses `<random>` for `std::mt19937`.

## 7. Building the Project

*   **Using an IDE (e.g., Visual Studio, CLion):**
    *   Create a new C++ project.
    *   Add all `.cpp` files (`main.cpp`, `forms1.cpp`, `forms2.cpp`, `forms_memory.cpp`, `make_positions.cpp`, and any new ones for positional logic) to the project.
    *   Add all `.h` files (`main.h`, `make_positions.h`) to the project or ensure they are in the include path.
    *   Build the project.
*   **Using g++ (or similar command-line compiler):**
    ```bash
    g++ -std=c++11 -o mcs_experiment main.cpp forms1.cpp forms2.cpp forms_memory.cpp make_positions.cpp [other_cpp_files.cpp] -Wall -O2
    ```
    (Add `-I.` if headers are in the same directory and not found automatically).

## 8. Running the Experiment

1.  Compile the project as described above to generate an executable `mcs.exe`.
3.  The program will:
    *   Determine the alphabet size `Y`.
    *   Generate (or load if existing) the main text file.
    *   Extract queries.
    *   Generate basic C(N,K) patterns (`tavnits.txt`).
    *   Execute the Naive search (skipped if `naive_matches_details.txt` exists).
    *   For each `M` value specified in `M_VALUES_PARAM`:
        *   Generate and refine MCS_0.
        *   Build the filter map from MCS_0.
        *   Run the "Usual MCS" search.
        *   Generate the sequence of "Positional-Associated MCS" sets.
        *   Run the "Positional MCS" search.
    *   Write all results and intermediate files to the working directory.
    *   Generate `experiment_conclusions.txt` with a summary.

## 9. Output Files and Understanding Results

*   **`main_text_Y<alphabet_size>.txt`**: The randomly generated text used for the experiment.
*   **`tavnits.txt`**: The set of C(N,K) binary patterns starting with '1'.
*   **`Usual_MCS_InitialForms_M<m_val>.txt`**: Filters generated by `chetv_struct_Generation` for MCS_0.
*   **`Usual_MCS_RefinedForms_M<m_val>.txt`**: MCS_0 filters after `Forms1_anal1_Global_Refinement` (if run).
*   **`Usual_Map_M<m_val>.txt`**: The filter map generated from the (refined) MCS_0.
*   **`Pos_MCS_M<m_val>_Filters_Pos<p>.txt`**: The refined MCS set for position `p` for a given `M`.
*   **`Pos_MCS_M<m_val>_ReducedTavnits_Pos<p>.txt`**: (If your older `perform_phase1_pattern_reduction` is still used somewhere or for debugging) The set of patterns remaining after a positional pruning step. With the new logic, this file's role for "flowing patterns" isn't explicitly saved unless you add it. The "problem patterns" used for `refine_mcs_iteratively_phase_b` are not currently saved by default.
*   **`naive_matches_details.txt`**: Detailed match positions for the naive search.
*   **`usual_mcs_matches_M<m_val>.txt`**: Detailed match positions for the Usual MCS search.
*   **`positional_mcs_matches_M<m_val>.txt`**: Detailed match positions for the Positional MCS search.
    *   Each match file typically contains: `QueryID NumMatches MatchPosition_1 MatchPosition_2 ...`
*   **`experiment_conclusions.txt`**:
    *   Summary of input parameters.
    *   Table of Preprocessing Time, Search Time, and Total Matches for each algorithm and M value.
    *   Discussion of theoretical time complexities.

## 10. Key Functions for Understanding the Logic

*   **`Forms::create_mas1` (`forms1.cpp`):** C(N,K) pattern generation.
*   **`Forms::chetv_struct_Generation` (`forms2.cpp`):** Initial MCS filter set generation (core greedy algorithm with substring coverage).
*   **`Forms::Forms1_anal1_Global_Refinement` (`forms2.cpp`):** Optional intensive global refinement for MCS_0.
*   **`build_filter_map` (`make_positions.cpp`):** Creates the hash map of masked text segments.
*   **`perform_positional_pattern_pruning_phase_a` (e.g., `make_positions.cpp`):** Implements the positional prefix-based pattern pruning for the Positional MCS generation.
*   **`Forms::refine_mcs_iteratively_phase_b` (e.g., `forms2.cpp` or `forms_positional_logic.cpp`):** Implements the iterative refinement of an MCS set based on sole substring coverage of a guiding pattern set (either "problem patterns" or all global patterns).
*   Search functions (`search_naive`, `search_usual_mcs`, `search_positional_mcs` in `main.cpp`): Implement the respective search strategies.
