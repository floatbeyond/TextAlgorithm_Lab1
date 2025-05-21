#include "main.h"           // For Forms class, FilterMapCollection
#include "make_positions.h" // For build_filter_map function

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cstdlib>      // For std::rand(), std::srand()
#include <ctime>        // For std::time()
#include <iomanip>      // For std::fixed, std::setprecision, std::setw
#include <cmath>        // For std::ceil, std::pow, std::round
#include <algorithm>    // For std::min, std::max
#include <chrono>       // For timing
#include <random>       // For std::mt19937, std::uniform_int_distribution

// --- Global Parameters (based on X=8, adjust if X changes) ---
const int X_PARAM = 8; // Last digit of ID (example)
const int N_PARAM_VAL = 25 - static_cast<int>(std::round(static_cast<double>(X_PARAM) / 2.0)); // W, Query Word Size (21 for X=8)
const double MATCH_SIMILARITY_THRESHOLD = (60.0 + X_PARAM) / 100.0; // (0.68 for X=8)
const int K_PARAM_VAL = static_cast<int>(std::ceil(N_PARAM_VAL * MATCH_SIMILARITY_THRESHOLD)); // K, ones in base patterns (15 for X=8, W=21, T=0.68)
const std::vector<int> M_VALUES_PARAM = { 3, 4, 5 }; // m values

const long long TEXT_SIZE_PARAM = 10000000; // 10M
// const long long TEXT_SIZE_PARAM = 100000; // 100K for testing
const int NUM_QUERIES_PARAM = 10000;
int ALPHABET_SIZE_Y_PARAM = 4; // Default, will be determined

// For random choices in positional refinement
std::mt19937 global_rng_engine(static_cast<unsigned int>(std::time(0)));


// --- Utility Functions ---
double calculate_similarity(const std::string& s1, const std::string& s2) {
    if (s1.length() != s2.length() || s1.empty()) {
        return 0.0;
    }
    int matches = 0;
    for (size_t i = 0; i < s1.length(); ++i) {
        if (s1[i] == s2[i]) {
            matches++;
        }
    }
    return static_cast<double>(matches) / s1.length();
}

// Helper function to calculate binomial coefficient C(n,k)
double nCr_prob(int n, int r) {
    if (r < 0 || r > n) return 0.0;
    if (r == 0 || r == n) return 1.0;
    if (r > n / 2) r = n - r;
    double res = 1.0;
    for (int i = 1; i <= r; ++i) {
        res = res * (n - i + 1) / i;
    }
    return res;
}

// Probability of exactly k matches in W trials, prob_char_match for one char
double binomial_pmf(int W_len, int k_successes, double prob_char_match) {
    if (k_successes < 0 || k_successes > W_len) return 0.0;
    return nCr_prob(W_len, k_successes) * std::pow(prob_char_match, k_successes) * std::pow(1.0 - prob_char_match, W_len - k_successes);
}

void determine_alphabet_size_Y() {
    std::cout << "Determining appropriate alphabet size Y for W=" << N_PARAM_VAL << ", K_matches_needed=" << K_PARAM_VAL << "..." << std::endl;
    for (int Y_test = 8; Y_test <= 26; ++Y_test) {
        double prob_single_char_match = 1.0 / Y_test;
        double prob_at_least_K_matches = 0.0;
        for (int k = K_PARAM_VAL; k <= N_PARAM_VAL; ++k) {
            prob_at_least_K_matches += binomial_pmf(N_PARAM_VAL, k, prob_single_char_match);
        }

        double expected_random_hits_in_text = (TEXT_SIZE_PARAM - N_PARAM_VAL + 1) * prob_at_least_K_matches;

        std::cout << "  Y = " << std::setw(2) << Y_test
            << ": P(single char match) = " << std::fixed << std::setprecision(4) << prob_single_char_match
            << ", P(>=K matches in W) = " << std::scientific << prob_at_least_K_matches
            << ", Expected random query hits in " << TEXT_SIZE_PARAM << " text: " << std::fixed << std::setprecision(2) << expected_random_hits_in_text
            << std::endl;

        if (expected_random_hits_in_text >= 2.0 && expected_random_hits_in_text <= 20.0) { // Aim for a small number
            ALPHABET_SIZE_Y_PARAM = Y_test;
            std::cout << "  Selected Y = " << ALPHABET_SIZE_Y_PARAM << std::endl;
            return;
        }
    }
    std::cout << "  Could not find Y for 2-20 expected hits. Defaulting to Y = " << ALPHABET_SIZE_Y_PARAM << std::endl;
}

std::string generate_random_text(long long length, int alphabet_size) {
    std::cout << "Generating random text of length " << length << " with alphabet size " << alphabet_size << "..." << std::endl;
    std::string text;
    text.reserve(length);
    std::uniform_int_distribution<> distrib(0, alphabet_size - 1);
    for (long long i = 0; i < length; ++i) {
        text += static_cast<char>('a' + distrib(global_rng_engine));
    }
    std::cout << "Random text generated." << std::endl;

    return text;
}

std::vector<std::string> extract_queries(const std::string& text, int num_queries, int query_length) {
    std::cout << "Extracting " << num_queries << " queries of length " << query_length << "..." << std::endl;
    std::vector<std::string> queries;
    queries.reserve(num_queries);
    for (int i = 0; i < num_queries; ++i) {
        if (static_cast<long long>(i) * query_length + query_length > text.length()) {
            std::cerr << "Warning: Not enough text to extract all " << num_queries << " non-overlapping queries." << std::endl;
            break;
        }
        queries.push_back(text.substr(static_cast<long long>(i) * query_length, query_length));
    }
    std::cout << queries.size() << " queries extracted." << std::endl;
    return queries;
}

// --- Search Algorithm Implementations ---

// Naive Search
std::set<int> search_naive(const std::string& query, const std::string& text) {
    std::set<int> match_positions;
    int query_len = query.length();
    if (query_len == 0 || text.length() < query_len) return match_positions;

    for (long long i = 0; i <= static_cast<long long>(text.length()) - query_len; ++i) {
        std::string text_segment = text.substr(i, query_len);
        if (calculate_similarity(query, text_segment) >= MATCH_SIMILARITY_THRESHOLD) {
            match_positions.insert(static_cast<int>(i));
        }
    }
    return match_positions;
}

// Helper to apply filter from Forms representation
std::string apply_forms_filter_to_string(const std::string& str_to_filter, const int* filter_bits, int filter_span, int N_glob_param) {
    std::string masked_str = "";
    if (filter_span <= 0 || str_to_filter.length() < static_cast<size_t>(filter_span)) {
        // Cannot apply: filter is longer than string or invalid span
        // Return empty or some indicator, or handle upstream
        return "INVALID_FILTER_APPLICATION";
    }
    masked_str.reserve(filter_span);
    for (int k = 0; k < filter_span; ++k) {
        if (filter_bits[k] == 1) {
            masked_str += str_to_filter[k];
        }
        else {
            masked_str += '_';
        }
    }
    return masked_str;
}


// Usual MCS Search
std::set<int> search_usual_mcs(const std::string& query, const Forms& mcs_object,
    const FilterMapCollection& filter_map, const std::string& text) {
    std::set<int> match_positions;
    if (query.length() != N_PARAM_VAL) return match_positions; // Query length must match W

    for (int i = 0; i < mcs_object.Nform1_Glob; ++i) {
        const int* current_filter_bits = mcs_object.forms1Glob[i];
        int filter_span = mcs_object.forms1Glob[i][mcs_object.N_glob]; // Actual span

        if (filter_span <= 0 || filter_span > N_PARAM_VAL) continue;

        std::string masked_query = apply_forms_filter_to_string(query, current_filter_bits, filter_span, mcs_object.N_glob);
        if (masked_query == "INVALID_FILTER_APPLICATION") continue;


        auto it = filter_map.word_to_index_map.find(masked_query);
        if (it != filter_map.word_to_index_map.end()) {
            const std::vector<int>& occurrences = filter_map.entries_vector[it->second].occurrences;
            for (int text_pos : occurrences) {
                if (static_cast<long long>(text_pos) + N_PARAM_VAL > text.length()) continue;
                std::string text_segment = text.substr(text_pos, N_PARAM_VAL);
                if (calculate_similarity(query, text_segment) >= MATCH_SIMILARITY_THRESHOLD) {
                    match_positions.insert(text_pos);
                }
            }
        }
    }
    return match_positions;
}

// Positional Associated MCS Search
std::set<int> search_positional_mcs(const std::string& query,
    const std::vector<Forms>& positional_mcs_sequence, // Sequence of MCS_0, MCS_1, ...
    const FilterMapCollection& filter_map, // Map built from MCS_0
    const std::string& text) {
    std::set<int> aggregated_match_positions;
    if (query.length() != N_PARAM_VAL) return aggregated_match_positions;

    // The PDF implies using MCS_i for a conceptual "position i" in processing.
    // Since our query and filters are fixed length W, we can iterate through the
    // generated MCS sets and apply each one, aggregating unique matches.
    // A more complex interpretation would be needed if queries were much longer than filters.
    int mcs_set_idx = 0;
    for (const Forms& mcs_set : positional_mcs_sequence) {
        // std::cout << "    Searching with Positional MCS_" << mcs_set_idx++ << " (num filters: " << mcs_set.Nform1_Glob << ")" << std::endl;
        for (int i = 0; i < mcs_set.Nform1_Glob; ++i) {
            const int* current_filter_bits = mcs_set.forms1Glob[i];
            int filter_span = mcs_set.forms1Glob[i][mcs_set.N_glob];

            if (filter_span <= 0 || filter_span > N_PARAM_VAL) continue;

            std::string masked_query = apply_forms_filter_to_string(query, current_filter_bits, filter_span, mcs_set.N_glob);
            if (masked_query == "INVALID_FILTER_APPLICATION") continue;


            auto it = filter_map.word_to_index_map.find(masked_query);
            if (it != filter_map.word_to_index_map.end()) {
                const std::vector<int>& occurrences = filter_map.entries_vector[it->second].occurrences;
                for (int text_pos : occurrences) {
                    if (static_cast<long long>(text_pos) + N_PARAM_VAL > text.length()) continue;
                    std::string text_segment = text.substr(text_pos, N_PARAM_VAL);
                    if (calculate_similarity(query, text_segment) >= MATCH_SIMILARITY_THRESHOLD) {
                        aggregated_match_positions.insert(text_pos);
                    }
                }
            }
        }
    }
    return aggregated_match_positions;
}


// --- Main Orchestration ---
int main() {
    std::srand(static_cast<unsigned int>(std::time(0))); // Seed for std::rand() if used by Forms
    // global_rng_engine is already seeded

    std::cout << "--- MCS Algorithm Comparison ---" << std::endl;
    std::cout << "Parameters: W=" << N_PARAM_VAL << " (Query/Window Size)" << std::endl;
    std::cout << "            K=" << K_PARAM_VAL << " (Min 1s in Base Patterns C(N,K))" << std::endl;
    std::cout << "            Match Threshold=" << MATCH_SIMILARITY_THRESHOLD * 100 << "%" << std::endl;
    std::cout << "            M values (1s in Filters): {" << M_VALUES_PARAM[0] << ", "
        << M_VALUES_PARAM[1] << ", " << M_VALUES_PARAM[2] << "}" << std::endl;

    // Step 1: Determine Alphabet Size Y
    // determine_alphabet_size_Y();
	ALPHABET_SIZE_Y_PARAM = 12; // For now, hardcoded for testing
    std::cout << "Using Alphabet Size Y = " << ALPHABET_SIZE_Y_PARAM << std::endl;

    // Step 2: Generate Text and Queries
    std::string main_text = generate_random_text(TEXT_SIZE_PARAM, ALPHABET_SIZE_Y_PARAM);
    std::vector<std::string> queries = extract_queries(main_text, NUM_QUERIES_PARAM, N_PARAM_VAL);

    // --- Generate Full Basic Patterns (Vars_Glob) ONCE ---
    Forms f_global_pattern_holder;
    f_global_pattern_holder.N_glob = N_PARAM_VAL;
    f_global_pattern_holder.N_2_glob = K_PARAM_VAL; // K (ones in base patterns)
    // form1Size will be set internally by methods needing it, or when filters are made

    std::cout << "\nGenerating full basic patterns (Vars_Glob for C(N,K)) ONCE..." << std::endl;
    auto timer_start_base_patterns = std::chrono::high_resolution_clock::now();
    f_global_pattern_holder.create_mas1(); // Creates "tavnits.txt" and populates Vars_Glob
    auto timer_end_base_patterns = std::chrono::high_resolution_clock::now();
    auto duration_base_patterns = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_base_patterns - timer_start_base_patterns);
    if (f_global_pattern_holder.NVars_Glob == 0) {
        std::cerr << "FATAL: No basic patterns generated by create_mas1. Exiting." << std::endl;
        return 1;
    }
    std::cout << "  Full basic patterns (tavnits.txt) generated: " << f_global_pattern_holder.NVars_Glob
        << " (Time: " << duration_base_patterns.count() << " ms)" << std::endl;

    // --- Storage for results for comparison ---
    std::map<std::string, std::chrono::milliseconds> total_times;
    std::map<std::string, int> total_matches_found;

	// --- Storage for results of each algorithm ---
    std::map<std::string, std::vector<std::set<int>>> all_algorithm_results;

    // --- Loop for each M value ---
    for (int m_val : M_VALUES_PARAM) {
        std::string m_str = "M" + std::to_string(m_val);
        std::cout << "\n=====================================================" << std::endl;
        std::cout << "Processing for " << m_str << " (N=" << N_PARAM_VAL << ", K=" << K_PARAM_VAL << ")" << std::endl;
        std::cout << "=====================================================" << std::endl;

        // --- A. "Usual" MCS (MCS_0 for this M) ---
        std::cout << "\nPART A: 'Usual' MCS (MCS_0 for " << m_str << ")" << std::endl;
        Forms mcs_0_object; // Default constructor initializes members
        mcs_0_object.N_glob = N_PARAM_VAL;
        mcs_0_object.N_2_glob = K_PARAM_VAL; // K used for Vars_Glob from which filters are derived
        mcs_0_object.Nsovp1_Glob = m_val;    // M (ones in filter)
        // mcs_0_object.form1Size will be N_PARAM_VAL + 1, set by chetv_struct or when allocating forms1Glob

        // Share the global Vars_Glob (read-only for chetv_struct_Generation)
        // This requires Forms to handle being assigned a Vars_Glob pointer without owning it,
        // OR to make a deep copy if true isolation is needed.
        // For now, let's assume chetv_struct_Generation uses the Vars_Glob already in mcs_0_object
        // (which should be populated by copying from f_global_pattern_holder)
        // This is safer: mcs_0_object gets its own copy of relevant params from f_global_pattern_holder
        //mcs_0_object.Vars_Glob = f_global_pattern_holder.Vars_Glob; // Share pointer
        //mcs_0_object.NVars_Glob = f_global_pattern_holder.NVars_Glob;



        std::string initial_forms_filename = "Usual_MCS_InitialForms_" + m_str + ".txt";
        std::cout << "  1. Generating MCS_0 filters (chetv_struct_Generation)..." << std::endl;
        auto timer_start_mcs0_gen = std::chrono::high_resolution_clock::now();
        mcs_0_object.chetv_struct_Generation(initial_forms_filename, f_global_pattern_holder);
        auto timer_end_mcs0_gen = std::chrono::high_resolution_clock::now();
        auto duration_mcs0_gen = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_mcs0_gen - timer_start_mcs0_gen);
        std::cout << "     Initial filters for " << m_str << ": " << mcs_0_object.Nform1_Glob
            << " (Time: " << duration_mcs0_gen.count() << " ms)" << std::endl;

        if (mcs_0_object.Nform1_Glob == 0) {
            std::cerr << "     WARNING: No MCS_0 filters generated for " << m_str << ". Skipping further steps for this M." << std::endl;
            continue;
        }

        // Optional: Global refinement of MCS_0 using Forms1_anal1_Global_Refinement
        std::string refined_mcs0_filename = "Usual_MCS_RefinedForms_" + m_str + ".txt";
        mcs_0_object.Forms1_anal1_Global_Refinement(
            initial_forms_filename,      // File to READ from (as per prof's code)
            f_global_pattern_holder,     // Source of base patterns
            refined_mcs0_filename        // File to WRITE refined filters to
        );

        std::string mcs0_saved_filename_prefix = "Usual_MCS_" + m_str;
        mcs_0_object.save_forms_to_file(mcs0_saved_filename_prefix + "_Filters.txt", "Usual MCS_0");


        std::cout << "  2. Building Filter Map for MCS_0 (" << m_str << ")..." << std::endl;
        FilterMapCollection usual_mcs_filter_map;
        std::string map_filename = "Usual_Map_" + m_str + ".txt";
        auto timer_start_map_build = std::chrono::high_resolution_clock::now();
        build_filter_map(mcs_0_object, main_text, usual_mcs_filter_map, map_filename, m_val);
        auto timer_end_map_build = std::chrono::high_resolution_clock::now();
        auto duration_map_build = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_map_build - timer_start_map_build);
        std::cout << "     Filter Map for " << m_str << " built. (Time: " << duration_map_build.count() << " ms)" << std::endl;

        std::cout << "  3. Searching with Usual MCS (" << m_str << ")..." << std::endl;
        std::string usual_mcs_algo_name = "Usual_MCS_M" + std::to_string(m_val);
        all_algorithm_results[usual_mcs_algo_name].resize(queries.size()); // Pre-allocate space
        long long current_total_matches = 0;

        auto timer_start_usual_search = std::chrono::high_resolution_clock::now();
        for (size_t q_idx = 0; q_idx < queries.size(); ++q_idx) {
            all_algorithm_results[usual_mcs_algo_name][q_idx] =
                search_usual_mcs(queries[q_idx], mcs_0_object, usual_mcs_filter_map, main_text);
            current_total_matches += all_algorithm_results[usual_mcs_algo_name][q_idx].size();
        }
        auto timer_end_usual_search = std::chrono::high_resolution_clock::now();
        auto duration_usual_search = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_usual_search - timer_start_usual_search);
        total_times[usual_mcs_algo_name] = duration_mcs0_gen + duration_map_build + duration_usual_search; // Assuming these are defined earlier
        total_matches_found[usual_mcs_algo_name] = current_total_matches;
        std::cout << "     Usual MCS Search for " << m_str << " complete. Matches: " << current_total_matches
            << " (Search Time: " << duration_usual_search.count() << " ms)"
            << " (Total Time: " << total_times["Usual_MCS_" + m_str].count() << " ms)" << std::endl;


        // --- B. "Positional-Associated" MCS (for this M value) ---
        std::cout << "\nPART B: 'Positional-Associated' MCS sequence for " << m_str << std::endl;
        std::vector<Forms> positional_mcs_sequence_storage;
        Forms mcs_0_for_pa = mcs_0_object; // Needs DEEP COPY if mcs_0_object is reused/modified
        // Assuming Forms copy constructor handles deep copy
        positional_mcs_sequence_storage.push_back(mcs_0_for_pa); // Add MCS_0
        Forms current_patterns_to_reduce_holder = f_global_pattern_holder; // Deep copy, owns its Vars_Glob for this context

        std::string pa_filename_prefix = "Pos_MCS_" + m_str;
        // MCS_0 for PA is already saved as part of Usual MCS.

        auto timer_start_pa_gen = std::chrono::high_resolution_clock::now();
        const int MAX_POSITIONAL_SETS = N_PARAM_VAL; // Max possible refinement depth
        Forms input_patterns_for_next_phase1 = f_global_pattern_holder; // Deep copy (or careful pointer share)

        for (int pos_idx = 0; pos_idx < MAX_POSITIONAL_SETS - 1; ++pos_idx) { // Generate MCS_1 from MCS_0, etc.
            std::cout << "    Generating Positional MCS Set " << pos_idx + 1 << " for " << m_str << "..." << std::endl;
            const Forms& mcs_i_ref = positional_mcs_sequence_storage.back(); // MCS_i

            Forms mcs_i_plus_1_candidate; // Will become MCS_{i+1}
            // Setup candidate: copy params, copy filters from mcs_i_ref to start with
            mcs_i_plus_1_candidate.N_glob = mcs_i_ref.N_glob;
            mcs_i_plus_1_candidate.N_2_glob = mcs_i_ref.N_2_glob;
            mcs_i_plus_1_candidate.Nsovp1_Glob = mcs_i_ref.Nsovp1_Glob;
            mcs_i_plus_1_candidate.form1Size = mcs_i_ref.form1Size;
            // Deep copy forms1Glob from mcs_i_ref
            if (mcs_i_ref.Nform1_Glob > 0) {
                mcs_i_plus_1_candidate.Nform1_Glob = mcs_i_ref.Nform1_Glob;
                mcs_i_plus_1_candidate.forms1Glob = new int* [mcs_i_ref.Nform1_Glob];
                for (int r = 0; r < mcs_i_ref.Nform1_Glob; ++r) {
                    mcs_i_plus_1_candidate.forms1Glob[r] = new int[mcs_i_ref.form1Size];
                    for (int c = 0; c < mcs_i_ref.form1Size; ++c) {
                        mcs_i_plus_1_candidate.forms1Glob[r][c] = mcs_i_ref.forms1Glob[r][c];
                    }
                }
            }
            else {
                mcs_i_plus_1_candidate.Nform1_Glob = 0;
                mcs_i_plus_1_candidate.forms1Glob = nullptr;
            }

            // Phase 1: Reduce basic patterns
            // The 'input_patterns_for_next_phase1' object holds the patterns to be reduced.
            // For pos_idx = 0, this is f_global_pattern_holder.
            // For pos_idx > 0, this will be the reducedVars_Glob from the previous iteration,
            // packaged into a Forms object.
            std::string reduced_tavnits_filename = pa_filename_prefix + "_ReducedTavnits_Pos" + std::to_string(pos_idx + 1) + ".txt";
            mcs_i_plus_1_candidate.perform_phase1_pattern_reduction(
                mcs_i_ref, 
                f_global_pattern_holder, 
                reduced_tavnits_filename);

            int n_reduced_vars_before_phase2 = mcs_i_plus_1_candidate.NReducedVars_Glob;
            int n_filters_before_phase2 = mcs_i_plus_1_candidate.Nform1_Glob;

            if (n_filters_before_phase2 == 0) { // No filters to refine from MCS_i
                std::cout << "      Candidate MCS_" << pos_idx + 1 << " started with 0 filters. Stopping." << std::endl;
                break;
            }
            if (n_reduced_vars_before_phase2 == 0 && n_filters_before_phase2 > 0) {
                std::cout << "      No basic patterns remained after Phase 1 reduction against MCS_" << pos_idx
                    << ". Stopping positional generation for " << m_str << "." << std::endl;
                break;
            }


            // Phase 2: Iteratively prune filters
            int refinement_prune_steps = 0;
            if (n_filters_before_phase2 > 0) { // Only prune if there are filters
                std::cout << "      Iteratively refining MCS_" << pos_idx + 1 << " (candidate filters: "
                    << mcs_i_plus_1_candidate.Nform1_Glob
                    << ", using " << mcs_i_plus_1_candidate.NReducedVars_Glob << " reduced patterns)" << std::endl;
                while (mcs_i_plus_1_candidate.Nform1_Glob > 0 &&
                    mcs_i_plus_1_candidate.perform_phase2_single_filter_prune_step(global_rng_engine)) {
                    refinement_prune_steps++;
                    // Safety break (optional, but good for debugging)
                    // if (refinement_prune_steps > 2 * n_filters_before_phase2 && n_filters_before_phase2 > 0) {
                    //     std::cerr << "      Warning: Excessive refinement iterations for MCS_" << pos_idx + 1 << ". Breaking." << std::endl;
                    //     break;
                    // }
                }
            }

            std::cout << "      Iteratively refining MCS_" << pos_idx + 1 << " (candidate filters: "
                << mcs_i_plus_1_candidate.Nform1_Glob << ")" << std::endl;
            while (mcs_i_plus_1_candidate.Nform1_Glob > 0 &&
                mcs_i_plus_1_candidate.perform_phase2_single_filter_prune_step(global_rng_engine)) {
                refinement_prune_steps++;
                if (refinement_prune_steps > 2 * mcs_i_ref.Nform1_Glob && mcs_i_ref.Nform1_Glob > 0) { // Safety break
                    std::cerr << "      Warning: Excessive refinement iterations for MCS_" << pos_idx + 1 << ". Breaking." << std::endl;
                    break;
                }
            }
            std::cout << "      Refinement for MCS_" << pos_idx + 1 << " complete (" << refinement_prune_steps << " prunes). Final filters: "
                << mcs_i_plus_1_candidate.Nform1_Glob << std::endl;

            int n_filters_after_phase2 = mcs_i_plus_1_candidate.Nform1_Glob;
            // Stop if:
            // 1. No filters were pruned (n_filters_after_phase2 == n_filters_before_phase2)
            // AND
            // 2. The number of reduced patterns in this step is the same as the number of patterns input to this step's Phase 1
            //    (mcs_i_plus_1_candidate.NReducedVars_Glob == input_patterns_for_next_phase1.NVars_Glob
            //     OR if input_patterns_for_next_phase1 was holding its patterns in reducedVars_Glob, then check that)
            // OR if the candidate MCS became empty.

            bool no_filters_pruned = (refinement_prune_steps == 0);
            // To check if patterns stopped reducing, compare NReducedVars_Glob of mcs_i_plus_1_candidate
            // with the NVars_Glob (or NReducedVars_Glob if that's where they were held) of input_patterns_for_next_phase1.
            // Let's assume input_patterns_for_next_phase1.Vars_Glob was used if it's the global set,
            // and input_patterns_for_next_phase1.reducedVars_Glob if it was a previously reduced set.
            // This gets a bit tricky with how `input_patterns_for_next_phase1` is structured.

            // Let's use your stated conditions:
            // 1. No filters removed in Phase 2 of this step.
            // 2. No patterns removed in Phase 1 of this step (i.e., NReducedVars_Glob == count of patterns fed into Phase 1).

            int patterns_fed_to_phase1_count;

            if (pos_idx == 0) { // First iteration, used global patterns
                patterns_fed_to_phase1_count = f_global_pattern_holder.NVars_Glob;
            }
            else {
                // Subsequent iterations, used the reducedVars_Glob from the *previous* mcs_i_plus_1_candidate,
                // which became input_patterns_for_next_phase1.
                // This assumes input_patterns_for_next_phase1's relevant patterns are in its 'Vars_Glob' if we copy them there.
                // Or, if perform_phase1_pattern_reduction can take patterns from either .Vars_Glob or .reducedVars_Glob of its source.
                // Let's assume 'input_patterns_for_next_phase1' always has its relevant patterns in its '.Vars_Glob' member
                // for simplicity of perform_phase1_pattern_reduction's current signature.
                patterns_fed_to_phase1_count = input_patterns_for_next_phase1.NVars_Glob;
            }
            bool no_patterns_deleted_in_phase1 = (n_reduced_vars_before_phase2 == patterns_fed_to_phase1_count);

            /*std::cout << "      Refinement for MCS_" << pos_idx + 1 << " complete (" << refinement_prune_steps << " prunes). Final filters: "
                << mcs_i_plus_1_candidate.Nform1_Glob << std::endl;*/

            if (mcs_i_plus_1_candidate.Nform1_Glob == 0) {
                std::cout << "      MCS_" << pos_idx + 1 << " became empty after refinement. Stopping." << std::endl;
                // Destructor of mcs_i_plus_1_candidate will clean its memory
                break;
            }

            positional_mcs_sequence_storage.push_back(mcs_i_plus_1_candidate); // Adds a COPY
            mcs_i_plus_1_candidate.save_forms_to_file(pa_filename_prefix + "_Filters_Pos" + std::to_string(pos_idx + 1) + ".txt", "Positional MCS");
            // --- Prepare for the NEXT iteration's Phase 1 input ---
            // The reducedVars_Glob from the *current* mcs_i_plus_1_candidate (now stored in the sequence)
            // becomes the set of patterns to be further reduced in the *next* iteration.
            // We need to transfer these into a Forms object that `perform_phase1_pattern_reduction` can read as a source.
            // Let's make `input_patterns_for_next_phase1` take these.
            // This requires `input_patterns_for_next_phase1` to clear its old patterns and take new ones.

            // Clear previous content of input_patterns_for_next_phase1 (if it was a deep copy holder)
            // (Its destructor will handle its own members if it goes out of scope and is re-created,
            // or we manage it explicitly)

            // Create a new Forms object to hold the reduced patterns that will be the input for the next iteration.
            Forms next_input_patterns;
            next_input_patterns.N_glob = positional_mcs_sequence_storage.back().N_glob; // N of the patterns
            next_input_patterns.N_2_glob = positional_mcs_sequence_storage.back().N_2_glob; // K context
            next_input_patterns.NVars_Glob = positional_mcs_sequence_storage.back().NReducedVars_Glob; // Count of patterns

            if (next_input_patterns.NVars_Glob > 0) {
                next_input_patterns.Vars_Glob = new bool* [next_input_patterns.NVars_Glob];
                for (int k = 0; k < next_input_patterns.NVars_Glob; ++k) {
                    next_input_patterns.Vars_Glob[k] = new bool[next_input_patterns.N_glob];
                    for (int l = 0; l < next_input_patterns.N_glob; ++l) {
                        // Copy from the reducedVars_Glob of the *just added* MCS set
                        next_input_patterns.Vars_Glob[k][l] = positional_mcs_sequence_storage.back().reducedVars_Glob[k][l];
                    }
                }
            }
            else {
                next_input_patterns.Vars_Glob = nullptr;
            }
            input_patterns_for_next_phase1 = next_input_patterns; // Assignment operator (deep copy)


            // New stopping condition check
            if (pos_idx > 0) { // Don't check for stopping on the very first generated MCS_1
                if (no_filters_pruned && no_patterns_deleted_in_phase1) {
                    std::cout << "      Stabilized: No filters pruned in this step AND no patterns deleted in Phase 1 of this step." << std::endl;
                    std::cout << "      Stopping positional generation for " << m_str << "." << std::endl;
                    break;
                }
            }
        } // End of pos_idx loop
        
        auto timer_end_pa_gen = std::chrono::high_resolution_clock::now();
        auto duration_pa_gen = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_pa_gen - timer_start_pa_gen);
        std::cout << "    Positional MCS sequence generation for " << m_str << " complete. Generated "
            << positional_mcs_sequence_storage.size() << " sets (including MCS_0)."
            << " (Time: " << duration_pa_gen.count() << " ms)" << std::endl;


        std::cout << "  3. Searching with Positional MCS (" << m_str << ")..." << std::endl;
        long long total_pos_mcs_matches_count = 0;
        // The filter map for PA search is the one built from MCS_0 (usual_mcs_filter_map)
        auto timer_start_pa_search = std::chrono::high_resolution_clock::now();
        for (const auto& q_str : queries) {
            total_pos_mcs_matches_count += search_positional_mcs(q_str, positional_mcs_sequence_storage, usual_mcs_filter_map, main_text).size();
        }
        auto timer_end_pa_search = std::chrono::high_resolution_clock::now();
        auto duration_pa_search = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_pa_search - timer_start_pa_search);
        // Total time for PA includes MCS_0 gen + PA sequence gen + map build (same as usual) + PA search
        total_times["Pos_MCS_" + m_str] = duration_mcs0_gen + duration_pa_gen + duration_map_build + duration_pa_search;
        total_matches_found["Pos_MCS_" + m_str] = total_pos_mcs_matches_count;
        std::cout << "     Positional MCS Search for " << m_str << " complete. Matches: " << total_pos_mcs_matches_count
            << " (Search Time: " << duration_pa_search.count() << " ms)"
            << " (Total Time: " << total_times["Pos_MCS_" + m_str].count() << " ms)" << std::endl;

        // mcs_0_object destructor will be called here
        // positional_mcs_sequence_storage destructor will call destructor for each Forms object
    } // End of M-value loop


    // --- Naive Search (Run Once) ---
    std::cout << "\nRunning Naive Search Algorithm..." << std::endl;
    long long total_naive_matches_count = 0;
    auto timer_start_naive_search = std::chrono::high_resolution_clock::now();
    for (const auto& q_str : queries) {
        total_naive_matches_count += search_naive(q_str, main_text).size();
    }
    auto timer_end_naive_search = std::chrono::high_resolution_clock::now();
    auto duration_naive_search = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_naive_search - timer_start_naive_search);
    total_times["Naive"] = duration_naive_search;
    total_matches_found["Naive"] = total_naive_matches_count;
    std::cout << "  Naive Search complete. Matches: " << total_naive_matches_count
        << " (Time: " << duration_naive_search.count() << " ms)" << std::endl;


    // --- Task 4 & 5: Compare Results & Theoretical Time ---
    std::cout << "\n--- Final Results Summary ---" << std::endl;
    std::cout << std::setw(25) << std::left << "Algorithm"
        << std::setw(15) << std::right << "Total Matches"
        << std::setw(15) << std::right << "Total Time (ms)" << std::endl;
    std::cout << std::string(55, '-') << std::endl;

    for (const auto& m_val : M_VALUES_PARAM) {
        std::string m_str = "M" + std::to_string(m_val);
        std::string usual_key = "Usual_MCS_" + m_str;
        std::string pos_key = "Pos_MCS_" + m_str;
        if (total_times.count(usual_key)) {
            std::cout << std::setw(25) << std::left << usual_key
                << std::setw(15) << std::right << total_matches_found[usual_key]
                << std::setw(15) << std::right << total_times[usual_key].count() << std::endl;
        }
        if (total_times.count(pos_key)) {
            std::cout << std::setw(25) << std::left << pos_key
                << std::setw(15) << std::right << total_matches_found[pos_key]
                << std::setw(15) << std::right << total_times[pos_key].count() << std::endl;
        }
    }
    std::cout << std::setw(25) << std::left << "Naive"
        << std::setw(15) << std::right << total_matches_found["Naive"]
        << std::setw(15) << std::right << total_times["Naive"].count() << std::endl;

    std::cout << "\nTheoretical Time Discussion:" << std::endl;
    std::cout << "Naive: O(num_queries * text_length * W) = O("
        << NUM_QUERIES_PARAM << " * " << TEXT_SIZE_PARAM << " * " << N_PARAM_VAL << ")" << std::endl;
    std::cout << "MCS-based methods are more complex:" << std::endl;
    std::cout << "  - Base Pattern Gen (create_mas1): Depends on C(N,K) or 2^N iteration." << std::endl;
    std::cout << "  - MCS_0 Gen (chetv_struct): O(NVars_Glob * N_glob * avg_mcs_size_check_complexity)." << std::endl;
    std::cout << "  - Map Build: O(TEXT_SIZE_PARAM * num_filters_in_mcs0 * filter_span_avg)." << std::endl;
    std::cout << "  - Search: O(num_queries * num_filters_in_mcs_set * filter_application + num_candidates * W)." << std::endl;
    std::cout << "  - Positional MCS Gen: Adds O(MAX_POS_SETS * (NVars_Glob_reduction + Nform1_Glob_pruning_loop * NReducedVars * Nform1_Glob_candidate))." << std::endl;
    std::cout << "Actual performance depends heavily on data, number of filters, and candidate hits." << std::endl;


    // f_global_pattern_holder destructor will be called, cleaning Vars_Glob
    std::cout << "\nAll processing complete." << std::endl;
    return 0;
}