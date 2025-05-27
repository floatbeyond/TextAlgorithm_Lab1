// main.cpp
#include "main.h"           
#include "make_positions.h" 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cstdlib>      
#include <ctime>        
#include <iomanip>      
#include <cmath>        
#include <algorithm>    
#include <chrono>       
#include <random>       

// --- Global Parameters (based on X=8) ---
const int X_PARAM = 8;
const int N_PARAM_VAL = 25 - static_cast<int>(std::round(static_cast<double>(X_PARAM) / 2.0));
const double MATCH_SIMILARITY_THRESHOLD = (60.0 + X_PARAM) / 100.0;
const int K_PARAM_VAL = static_cast<int>(std::round(N_PARAM_VAL * MATCH_SIMILARITY_THRESHOLD)); 
const std::vector<int> M_VALUES_PARAM = { 3, 4, 5 };

const long long TEXT_SIZE_PARAM = 10000000;
const int NUM_QUERIES_PARAM = 10000;
int ALPHABET_SIZE_Y_PARAM = 4; // Default, will be determined or overridden

std::mt19937 global_rng_engine; // Declare global, seed in main


// --- Utility Functions ---
double nCr_prob_util(int n, int r) { // Renamed to avoid conflict if Forms has one
    if (r < 0 || r > n) return 0.0;
    if (r == 0 || r == n) return 1.0;
    if (r > n / 2) r = n - r;
    double res = 1.0;
    for (int i = 1; i <= r; ++i) res = res * (n - i + 1) / i;
    return res;
}

double binomial_pmf_util(int W_len, int k_successes, double prob_char_match) {
    if (k_successes < 0 || k_successes > W_len) return 0.0;
    if (prob_char_match == 0.0) return (k_successes == 0) ? 1.0 : 0.0;
    if (prob_char_match == 1.0) return (k_successes == W_len) ? 1.0 : 0.0;

    // Use logs for stability if numbers are extreme, direct for typical cases
    double log_nCr = std::log(nCr_prob_util(W_len, k_successes));
    double log_p_k = k_successes * std::log(prob_char_match);
    double log_1_minus_p_nk = (W_len - k_successes) * std::log(1.0 - prob_char_match);

    if (std::isinf(log_nCr) || std::isinf(log_p_k) || std::isinf(log_1_minus_p_nk)) { // check for -inf
        if (log_nCr == -std::numeric_limits<double>::infinity() && k_successes > 0 && k_successes < W_len) return 0; // C(n,k) was 0
        // other underflow/overflow checks might be needed if direct method is used
    }
    try {
        return std::exp(log_nCr + log_p_k + log_1_minus_p_nk);
    }
    catch (const std::overflow_error& e) {
        // Fallback to direct if exp overflows, though logs should help prevent this
        return nCr_prob_util(W_len, k_successes) * std::pow(prob_char_match, k_successes) * std::pow(1.0 - prob_char_match, W_len - k_successes);
    }
}

void determine_alphabet_size_Y() {
    std::cout << "Determining appropriate alphabet size Y for W=" << N_PARAM_VAL << ", K_matches_needed=" << K_PARAM_VAL << "..." << std::endl;
    bool found_y = false;
    for (int Y_test = 2; Y_test <= 26; ++Y_test) { // Changed from 12 to 2, as per original prompt idea
        double prob_single_char_match = 1.0 / Y_test;
        double prob_at_least_K_matches = 0.0;
        for (int k = K_PARAM_VAL; k <= N_PARAM_VAL; ++k) {
            prob_at_least_K_matches += binomial_pmf_util(N_PARAM_VAL, k, prob_single_char_match);
        }
        double expected_random_hits_in_text = (TEXT_SIZE_PARAM - N_PARAM_VAL + 1) * prob_at_least_K_matches;
        std::cout << "  Y = " << std::setw(2) << Y_test
            << ": P(single char match) = " << std::fixed << std::setprecision(4) << prob_single_char_match
            << ", P(>=K matches in W) = " << std::scientific << std::setprecision(4) << prob_at_least_K_matches
            << ", Expected random query hits in " << TEXT_SIZE_PARAM << " text: " << std::fixed << std::setprecision(2) << expected_random_hits_in_text
            << std::endl;
        if (!found_y && expected_random_hits_in_text >= 2.0 && expected_random_hits_in_text <= 20.0) {
            ALPHABET_SIZE_Y_PARAM = Y_test;
            std::cout << "  Selected Y = " << ALPHABET_SIZE_Y_PARAM << " (based on 2-20 expected hits criteria)" << std::endl;
            found_y = true; // Take the first one in range
        }
    }
    if (!found_y) {
        std::cout << "  Could not find Y for 2-20 expected hits in range 2-26. Using hardcoded Y." << std::endl;
        ALPHABET_SIZE_Y_PARAM = 12; // Fallback to user's test value if criteria not met
        std::cout << "  Using fallback Y = " << ALPHABET_SIZE_Y_PARAM << std::endl;
    }
}

std::string generate_random_text(long long length, int alphabet_size, const std::string& output_filename, std::mt19937& rng_engine_ref) {
    std::cout << "Generating random text of length " << length << " with alphabet size " << alphabet_size << "..." << std::endl;
    std::string text;
    if (length <= 0) { std::cout << "Text length 0, no text generated." << std::endl; return text; }
    text.reserve(length);
    std::uniform_int_distribution<> distrib(0, alphabet_size - 1);
    for (long long i = 0; i < length; ++i) {
        text += static_cast<char>('a' + distrib(rng_engine_ref));
    }
    std::cout << "Random text generation complete." << std::endl;
    if (!output_filename.empty()) {
        std::ofstream outfile(output_filename);
        if (outfile.is_open()) {
            outfile << text; outfile.close();
            std::cout << "  Generated text saved to " << output_filename << std::endl;
        }
        else {
            std::cerr << "  ERROR: Could not open " << output_filename << " to save text." << std::endl;
        }
    }
    return text;
}

std::vector<std::string> extract_queries(const std::string& text, int num_queries, int query_length) {
    std::cout << "Extracting " << num_queries << " queries of length " << query_length << "..." << std::endl;
    std::vector<std::string> queries;
    if (query_length <= 0) { std::cerr << "Query length 0, no queries extracted." << std::endl; return queries; }
    queries.reserve(num_queries);
    for (int i = 0; i < num_queries; ++i) {
        if (static_cast<long long>(i) * query_length + query_length > text.length()) {
            std::cerr << "Warning: Not enough text to extract all " << num_queries << " non-overlapping queries." << std::endl;
            break;
        }
        queries.push_back(text.substr(static_cast<size_t>(i) * query_length, query_length)); // Use size_t for substr
    }
    std::cout << queries.size() << " queries extracted." << std::endl;
    return queries;
}

// Place this utility function somewhere accessible in main.cpp or a utility section
void save_match_positions(const std::string& filename,
    const std::string& query_id, // Or the query string itself
    const std::set<int>& match_positions) {
    std::ofstream outfile(filename, std::ios_base::app); // Open in append mode
    if (!outfile.is_open()) {
        std::cerr << "  ERROR: Could not open file " << filename << " to save match positions." << std::endl;
        return;
    }
    // outfile << "Query: " << query_id << "\tMatches (" << match_positions.size() << "):\t"; // Old more verbose
    outfile << query_id << "\t" << match_positions.size(); // Query ID and count
    for (int pos : match_positions) {
        outfile << "\t" << pos; // Tab-separated positions
    }
    outfile << std::endl;
    // No need to close explicitly here if only appending one line, destructor will close.
    // But if appending many times in a loop, periodic close/reopen or flush might be good.
    // For this use case (once per query), it's fine.
}

double calculate_similarity(const std::string& s1, const std::string& s2) {
    if (s1.length() != s2.length() || s1.empty()) return 0.0;
    int matches = 0;
    for (size_t i = 0; i < s1.length(); ++i) if (s1[i] == s2[i]) matches++;
    return static_cast<double>(matches) / s1.length();
}

std::string apply_forms_filter_to_string(const std::string& str_to_filter, const int* filter_bits, int filter_span, int N_glob_param_filter_storage_length) {
    std::string masked_str = "";
    if (filter_span <= 0 || str_to_filter.length() < static_cast<size_t>(filter_span) || filter_span > N_glob_param_filter_storage_length) {
        return "INVALID_FILTER_APPLICATION";
    }
    masked_str.reserve(filter_span);
    for (int k = 0; k < filter_span; ++k) {
        masked_str += (filter_bits[k] == 1 ? str_to_filter[k] : '_');
    }
    return masked_str;
}

// --- Search Algorithm Implementations ---
// Performs a naive, brute-force search for a query in a text.
// Slides a window of size W across the text and calculates similarity at each position.
std::set<int> search_naive(const std::string& query,         // The query string to search for
    const std::string& text,          // The text to search within
    int W,                           // The expected length of the query (and text segments)
    double threshold)                 // The similarity threshold for a match
{
    std::set<int> match_positions; // Stores starting positions of confirmed matches

    // Basic validation: query length must match W, and text must be long enough.
    if (W == 0 || query.length() != static_cast<size_t>(W) || text.length() < static_cast<size_t>(W)) {
        // std::cerr << "Warning (search_naive): Invalid parameters or query/text length." << std::endl;
        return match_positions; // Return empty set if params are invalid
    }

    // Iterate through all possible starting positions for a segment of length W in the text.
    // The loop condition ensures that `text.substr(i, W)` is always valid.
    for (long long i = 0; i <= static_cast<long long>(text.length()) - W; ++i) {
        // Extract the current text segment of length W.
        std::string text_segment = text.substr(i, W);

        // Calculate the similarity between the query and the current text segment.
        if (calculate_similarity(query, text_segment) >= threshold) {
            // If similarity meets or exceeds the threshold, it's a match.
            match_positions.insert(static_cast<int>(i)); // Store the starting position of the match.
        }
    }
    return match_positions; // Return the set of all found match positions.
}

// --- MCS Search Algorithm ---
// Performs a search using the "Usual" MCS (Multiple Correlated Substrings) approach.
// It uses a precomputed filter map to quickly identify candidate positions,
// then verifies these candidates using full similarity calculation.
std::set<int> search_usual_mcs(
    const std::string& query,         // The query string (length W)
    const Forms& mcs_object,        // MCS filters
    const FilterMapCollection& filter_map, // Precomputed map from build_filter_map
    const std::string& text,          // The text to search within
    int W,                           // Query length (should be query.length())
    double threshold)                 // Similarity threshold
{
    std::set<int> match_positions; // Stores starting positions of confirmed matches

    // Validate query length against W.
    if (query.length() != static_cast<size_t>(W)) {
        // std::cerr << "Warning (search_usual_mcs): Query length mismatch." << std::endl;
        return match_positions;
    }

    // Loop 1: For each filter in the MCS object
    for (int filter_idx = 0; filter_idx < mcs_object.Nform1_Glob; ++filter_idx) {
        const int* current_filter_bits = mcs_object.forms1Glob[filter_idx];
        int filter_span = mcs_object.forms1Glob[filter_idx][mcs_object.N_glob]; // Actual span of this filter

        // Basic filter validation
        if (filter_span <= 0 || filter_span > W) { // Filter cannot be longer than the query
            continue;
        }

        // Loop 2: Slide the current filter across the query string
        // query_slide_pos is the starting position of the filter on the query
        for (int query_slide_pos = 0; query_slide_pos <= W - filter_span; ++query_slide_pos) {

            // Extract the segment of the query that the filter currently covers
            std::string query_segment_to_mask = query.substr(query_slide_pos, filter_span);

            // Apply the filter to this specific query segment
            // The last argument to apply_forms_filter_to_string is the max storage length for filter_bits,
            // actual masking uses filter_span.
            std::string masked_query_segment = apply_forms_filter_to_string(
                query_segment_to_mask, // Pass the query *segment*
                current_filter_bits,
                filter_span,
                mcs_object.N_glob // N_glob from the Forms object storing the filter
            );

            if (masked_query_segment == "INVALID_FILTER_APPLICATION") {
                continue;
            }

            // Look up this masked_query_segment in the precomputed filter_map
            auto it = filter_map.word_to_index_map.find(masked_query_segment);

            if (it != filter_map.word_to_index_map.end()) {
                // If the masked_query_segment (e.g., "B_DE" from query) is found in the map...
                // ...it means this masked pattern also appeared in the original TEXT.
                // 'occurrences_in_text' are the starting positions in TEXT where this masked pattern appeared.
                const std::vector<int>& occurrences_in_text = filter_map.entries_vector[it->second].occurrences;

                // STAGE 2: Candidate Verification
                for (int text_start_of_masked_segment : occurrences_in_text) {
                    // `text_start_of_masked_segment` is where the *masked segment* started in the text.
                    // We need to figure out where the *full W-length query* would align in the text
                    // if its `query_segment_to_mask` (which started at `query_slide_pos` within the query)
                    // matches the text at `text_start_of_masked_segment`.
                    //
                    // The start of the full W-length query in the text would be:
                    // `text_start_of_masked_segment - query_slide_pos`
                    long long potential_full_query_start_in_text_ll =
                        static_cast<long long>(text_start_of_masked_segment) - query_slide_pos;

                    // Boundary checks for the W-length alignment in the text
                    if (potential_full_query_start_in_text_ll < 0 ||
                        potential_full_query_start_in_text_ll + W > text.length()) {
                        continue; // This alignment would place the full query outside text bounds
                    }
                    int potential_full_query_start_in_text = static_cast<int>(potential_full_query_start_in_text_ll);

                    // Now, extract the W-length text segment for full, unmasked comparison
                    std::string text_segment_for_verification = text.substr(potential_full_query_start_in_text, W);

                    // Compare the ORIGINAL full query with the ORIGINAL W-length text segment
                    if (calculate_similarity(query, text_segment_for_verification) >= threshold) {
                        match_positions.insert(potential_full_query_start_in_text); // Confirmed match
                    }
                } // End loop over occurrences_in_text
            } // End if masked_query_segment found in map
        } // End of Loop 2 (sliding filter across query)
    } // End of Loop 1 (iterating through filters)
    return match_positions;
}

// --- Positional MCS Search Algorithm ---
// Performs a search using the "Positional-Associated" MCS approach.
// It iterates through a sequence of MCS sets (MCS_0, MCS_1, ..., MCS_p).
// For each filter in each MCS set, it applies the filter to the query,
// looks up the masked query in a filter map (typically built from MCS_0),
// and then verifies candidate positions. Matches from all sets are aggregated.
std::set<int> search_positional_mcs(const std::string& query,                     // The query string
    const std::vector<Forms>& positional_mcs_sequence, // Sequence of MCS sets
    const FilterMapCollection& map_from_mcs0,       // Filter map, usually built from MCS_0
    const std::string& text,                      // The text to search within
    int W,                                       // Expected query/window length
    double threshold)                             // Similarity threshold
{
    std::set<int> aggregated_match_positions; // Stores unique match positions from all MCS sets

    // Validate query length.
    if (query.length() != static_cast<size_t>(W)) {
        // std::cerr << "Warning (search_positional_mcs): Query length mismatch." << std::endl;
        return aggregated_match_positions;
    }

    // Iterate through each MCS set in the positional sequence (MCS_0, MCS_1, etc.).
    for (const Forms& mcs_set : positional_mcs_sequence) {
        // --- STAGE 1 (per MCS set): Candidate Identification ---
        // For each filter within the current mcs_set...
        for (int i = 0; i < mcs_set.Nform1_Glob; ++i) {
            const int* current_filter_bits = mcs_set.forms1Glob[i];
            int filter_span = mcs_set.forms1Glob[i][mcs_set.N_glob]; // N_glob of this mcs_set

            // Validate filter span against query length W.
            if (filter_span <= 0 || filter_span > W) {
                continue;
            }

            // Apply the current filter to the query string.
            // mcs_set.N_glob is passed as the filter storage length.
            std::string masked_query = apply_forms_filter_to_string(query, current_filter_bits, filter_span, mcs_set.N_glob);

            if (masked_query == "INVALID_FILTER_APPLICATION") {
                continue;
            }

            // Look up the masked_query in the filter map (which was typically built from MCS_0).
            // The assumption is that filters in subsequent positional sets (MCS_1, MCS_2...)
            // are still meaningful when their masked versions are looked up in MCS_0's map.
            auto it = map_from_mcs0.word_to_index_map.find(masked_query);

            if (it != map_from_mcs0.word_to_index_map.end()) {
                const std::vector<int>& occurrences = map_from_mcs0.entries_vector[it->second].occurrences;

                // --- STAGE 2 (per MCS set): Candidate Verification ---
                for (int text_pos : occurrences) {
                    // Bounds check.
                    if (static_cast<long long>(text_pos) + W > text.length()) {
                        continue;
                    }

                    // Extract original text segment and verify against original query.
                    std::string text_segment = text.substr(text_pos, W);
                    if (calculate_similarity(query, text_segment) >= threshold) {
                        aggregated_match_positions.insert(text_pos); // Add to set (handles uniqueness).
                    }
                }
            }
        }
    }
    return aggregated_match_positions; // Return all unique matches found across all positional MCS sets.
}

// --- NEW Helper Functions for file operations ---
bool file_exists(const std::string& filename) {
    std::ifstream infile(filename);
    return infile.good();
}

std::string read_text_from_file(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "  ERROR (read_text_from_file): Could not open file " << filename << " for reading." << std::endl;
        return ""; // Return empty string on failure
    }
    // Read the whole file into a string
    std::string content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
    infile.close();
    std::cout << "  Text successfully read from " << filename << " (Length: " << content.length() << ")" << std::endl;
    return content;
}

// --- Main Orchestration ---
int main() {
    unsigned int seed = static_cast<unsigned int>(std::time(0));
    global_rng_engine.seed(seed);
    std::srand(seed);

    std::cout << std::fixed << std::setprecision(4); // Default precision for cout

    std::cout << "--- MCS Algorithm Comparison ---" << std::endl;
    std::cout << "Parameters: W(N_PARAM_VAL)=" << N_PARAM_VAL
        << ", K(K_PARAM_VAL)=" << K_PARAM_VAL
        << ", Match Threshold=" << MATCH_SIMILARITY_THRESHOLD * 100 << "%" << std::endl;
    std::cout << "            M values: {" << M_VALUES_PARAM[0] << ", "
        << M_VALUES_PARAM[1] << ", " << M_VALUES_PARAM[2] << "}" << std::endl;

    determine_alphabet_size_Y(); // Sets ALPHABET_SIZE_Y_PARAM
	// ALPHABET_SIZE_Y_PARAM = 6; // Hardcoded for testing purposes
    std::cout << "Using final Alphabet Size Y = " << ALPHABET_SIZE_Y_PARAM << std::endl;

	// Old text generation/loading logic commented out
    //std::string main_text_file = "main_text_Y" + std::to_string(ALPHABET_SIZE_Y_PARAM) + ".txt";
    //std::string main_text = generate_random_text(TEXT_SIZE_PARAM, ALPHABET_SIZE_Y_PARAM, main_text_file, global_rng_engine);
    //if (main_text.empty() && TEXT_SIZE_PARAM > 0) {
    //    std::cerr << "FATAL: Main text generation failed. Exiting." << std::endl; return 1;
    //}
    //std::vector<std::string> queries = extract_queries(main_text, NUM_QUERIES_PARAM, N_PARAM_VAL);
    //if (queries.empty() && NUM_QUERIES_PARAM > 0) {
    //    std::cerr << "FATAL: Query extraction failed. Exiting." << std::endl; return 1;
    //}

    // --- Conditional Text Generation/Loading ---
    std::string main_text_file = "main_text_Y" + std::to_string(ALPHABET_SIZE_Y_PARAM) + ".txt";
    std::string main_text;

    if (file_exists(main_text_file)) {
        std::cout << "\nText file " << main_text_file << " found. Loading text..." << std::endl;
        main_text = read_text_from_file(main_text_file);
        if (main_text.empty() && TEXT_SIZE_PARAM > 0) { // Check if read failed or file was empty but expected content
            std::cerr << "FATAL: Failed to read text from existing file or file is empty. Exiting." << std::endl; return 1;
        }
    }
    else {
        std::cout << "\nText file " << main_text_file << " not found. Generating new text..." << std::endl;
        main_text = generate_random_text(TEXT_SIZE_PARAM, ALPHABET_SIZE_Y_PARAM, main_text_file, global_rng_engine);
        if (main_text.empty() && TEXT_SIZE_PARAM > 0) {
            std::cerr << "FATAL: Main text generation failed. Exiting." << std::endl; return 1;
        }
    }

    std::vector<std::string> queries = extract_queries(main_text, NUM_QUERIES_PARAM, N_PARAM_VAL);
    if (queries.empty() && NUM_QUERIES_PARAM > 0) {
        std::cerr << "FATAL: Query extraction failed. Exiting." << std::endl; return 1;
    }


    Forms f_global_pattern_holder;
    f_global_pattern_holder.N_glob = N_PARAM_VAL;
    f_global_pattern_holder.N_2_glob = K_PARAM_VAL;
    std::cout << "\nGenerating full basic patterns (Vars_Glob for C(N,K)) ONCE..." << std::endl;
    auto timer_start_base_patterns = std::chrono::high_resolution_clock::now();
    f_global_pattern_holder.create_mas1();
    auto timer_end_base_patterns = std::chrono::high_resolution_clock::now();
    auto duration_base_patterns = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_base_patterns - timer_start_base_patterns);
    if (f_global_pattern_holder.NVars_Glob == 0 && K_PARAM_VAL > 0) { // If K=0, NVars_Glob can be 0 (e.g. pattern "000...")
        std::cerr << "FATAL: No basic patterns generated by create_mas1. Exiting." << std::endl; return 1;
    }
    std::cout << "  Full basic patterns (tavnits.txt) generated: " << f_global_pattern_holder.NVars_Glob
        << " (Time: " << duration_base_patterns.count() << " ms)" << std::endl;

    std::map<std::string, std::chrono::milliseconds> total_times_preprocess;
    std::map<std::string, std::chrono::milliseconds> total_times_search;
    std::map<std::string, long long> total_matches_found;


    for (int m_val : M_VALUES_PARAM) {
        std::string m_str = "M" + std::to_string(m_val);
        std::cout << "\n================ PROCESSING FOR " << m_str << " ================" << std::endl;

        // A. "Usual" MCS (MCS_0)
        std::cout << "\nPART A: 'Usual' MCS (MCS_0 for " << m_str << ")" << std::endl;
        Forms mcs_0_object;
        mcs_0_object.N_glob = N_PARAM_VAL;
        mcs_0_object.N_2_glob = K_PARAM_VAL;
        mcs_0_object.Nsovp1_Glob = m_val;

        std::string initial_forms_filename = "Usual_MCS_InitialForms_" + m_str + ".txt";
        std::cout << "  1. Generating MCS_0 filters (chetv_struct_Generation)..." << std::endl;
        auto timer_start_mcs0_gen = std::chrono::high_resolution_clock::now();
        mcs_0_object.chetv_struct_Generation(initial_forms_filename, f_global_pattern_holder);
        auto timer_end_mcs0_gen = std::chrono::high_resolution_clock::now();
        auto duration_mcs0_gen = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_mcs0_gen - timer_start_mcs0_gen);
        std::cout << "     Initial filters for " << m_str << ": " << mcs_0_object.Nform1_Glob
            << " (Time: " << duration_mcs0_gen.count() << " ms)" << std::endl;

        bool run_global_refinement = false; // Set to true to run Forms1_anal1
        if (run_global_refinement && mcs_0_object.Nform1_Glob > 0) {
            std::cout << "  (Optional) Globally refining MCS_0 for " << m_str << "..." << std::endl;
            std::string refined_mcs0_filename = "Usual_MCS_RefinedForms_" + m_str + ".txt";
            auto timer_start_mcs0_refine = std::chrono::high_resolution_clock::now();
            mcs_0_object.Forms1_anal1_Global_Refinement(initial_forms_filename, f_global_pattern_holder, refined_mcs0_filename);
            auto timer_end_mcs0_refine = std::chrono::high_resolution_clock::now();
            duration_mcs0_gen += std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_mcs0_refine - timer_start_mcs0_refine); // Add to gen time
            std::cout << "     MCS_0 for " << m_str << " refined. Final count: " << mcs_0_object.Nform1_Glob
                << " (Refinement Time: " << (timer_end_mcs0_refine - timer_start_mcs0_refine).count() << " ms)" << std::endl;
        }
        else if (mcs_0_object.Nform1_Glob == 0) {
            std::cout << "     Skipping global refinement as MCS_0 is empty for " << m_str << "." << std::endl;
        }


        if (mcs_0_object.Nform1_Glob == 0) {
            std::cerr << "     WARNING: No MCS_0 filters for " << m_str << ". Skipping further steps for this M." << std::endl;
            continue;
        }
        // save_forms_to_file is called inside Forms1_anal1_Global_Refinement if run, or by chetv_struct if not.
        // If neither saves with final name, add a save here:
        // mcs_0_object.save_forms_to_file("Usual_MCS_Final_" + m_str + ".txt", "Final Usual MCS_0");


        std::cout << "  2. Building Filter Map for MCS_0 (" << m_str << ")..." << std::endl;
        FilterMapCollection usual_mcs_filter_map;
        std::string map_filename = "Usual_Map_" + m_str + ".txt";
        auto timer_start_map_build = std::chrono::high_resolution_clock::now();
        build_filter_map(mcs_0_object, main_text, usual_mcs_filter_map, map_filename, m_val);
        auto timer_end_map_build = std::chrono::high_resolution_clock::now();
        auto duration_map_build = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_map_build - timer_start_map_build);
        std::cout << "     Filter Map for " << m_str << " built. (Time: " << duration_map_build.count() << " ms)" << std::endl;
        total_times_preprocess["Usual_MCS_" + m_str] = duration_base_patterns + duration_mcs0_gen + duration_map_build;


        std::cout << "  3. Searching with Usual MCS (" << m_str << ")..." << std::endl;
        long long current_usual_mcs_matches_count_total = 0; // Renamed
        std::string usual_mcs_matches_filename = "usual_mcs_matches_" + m_str + ".txt";
        // Clear file and write header
        std::ofstream usual_outfile(usual_mcs_matches_filename, std::ios_base::trunc);
        if (usual_outfile.is_open()) {
            usual_outfile << "QueryID\tNumMatches\tMatchPosition_1\tMatchPosition_2\t..." << std::endl;
            usual_outfile.close();
        }
        else {
            std::cerr << "ERROR: Could not open " << usual_mcs_matches_filename << " to write header for M=" << m_val << std::endl;
        }

        auto timer_start_usual_search = std::chrono::high_resolution_clock::now();
        int query_idx_usual = 0;
        for (const auto& q_str : queries) {
            std::set<int> current_query_matches = search_usual_mcs(q_str, mcs_0_object, usual_mcs_filter_map, main_text, N_PARAM_VAL, MATCH_SIMILARITY_THRESHOLD);
            current_usual_mcs_matches_count_total += current_query_matches.size();

            save_match_positions(usual_mcs_matches_filename, "Query_" + std::to_string(query_idx_usual), current_query_matches);
            query_idx_usual++;
        }
        auto timer_end_usual_search = std::chrono::high_resolution_clock::now();
        total_times_search["Usual_MCS_" + m_str] = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_usual_search - timer_start_usual_search);
        total_matches_found["Usual_MCS_" + m_str] = current_usual_mcs_matches_count_total; // Update map
        std::cout << "     Usual MCS Search for " << m_str << " complete. Matches: " << current_usual_mcs_matches_count_total
            << " (Search Time: " << total_times_search["Usual_MCS_" + m_str].count() << " ms)" << std::endl;

        // B. "Positional-Associated" MCS
        std::cout << "\nPART B: 'Positional-Associated' MCS sequence for " << m_str << std::endl;
        std::vector<Forms> positional_mcs_sequence_storage;
        positional_mcs_sequence_storage.push_back(mcs_0_object); // MCS_0 is the first element (deep copy by vector)

        std::string pa_filename_prefix = "Pos_MCS_" + m_str;
        std::string current_input_patterns_file = "tavnits.txt"; // Initial patterns from global holder
        int patterns_fed_to_current_phase1 = f_global_pattern_holder.NVars_Glob;

        auto timer_start_pa_gen_seq = std::chrono::high_resolution_clock::now();
        for (int pos_idx = 0; pos_idx < N_PARAM_VAL - 1; ++pos_idx) { // Max N-1 refinements
            std::cout << "    Generating Positional MCS Set " << pos_idx + 1 << " for " << m_str << "..." << std::endl;
            const Forms& mcs_i_ref = positional_mcs_sequence_storage.back();

            Forms mcs_i_plus_1_candidate;
            mcs_i_plus_1_candidate.N_glob = N_PARAM_VAL; // Will be updated by Read_Patterns_FromFile inside perform_phase1
            mcs_i_plus_1_candidate.N_2_glob = K_PARAM_VAL; // Context
            mcs_i_plus_1_candidate.Nsovp1_Glob = mcs_i_ref.Nsovp1_Glob;
            // mcs_i_plus_1_candidate.form1Size will be N_glob+1, N_glob from read patterns

            // Deep copy filters from mcs_i_ref
            if (mcs_i_ref.Nform1_Glob > 0) {
                mcs_i_plus_1_candidate.form1Size = mcs_i_ref.form1Size; // Important for allocation
                mcs_i_plus_1_candidate.Nform1_Glob = mcs_i_ref.Nform1_Glob;
                mcs_i_plus_1_candidate.forms1Glob = new int* [mcs_i_ref.Nform1_Glob];
                for (int r = 0; r < mcs_i_ref.Nform1_Glob; ++r) {
                    mcs_i_plus_1_candidate.forms1Glob[r] = new int[mcs_i_ref.form1Size];
                    for (int c = 0; c < mcs_i_ref.form1Size; ++c) {
                        mcs_i_plus_1_candidate.forms1Glob[r][c] = mcs_i_ref.forms1Glob[r][c];
                    }
                }
            } // else Nform1_Glob = 0, forms1Glob = nullptr

            std::string output_reduced_patterns_filename = pa_filename_prefix + "_ReducedTavnits_Pos" + std::to_string(pos_idx + 1) + ".txt";

            if (patterns_fed_to_current_phase1 == 0 && mcs_i_plus_1_candidate.Nform1_Glob > 0) {
                std::cout << "      Phase 1: No input patterns from file '" << current_input_patterns_file << "' to reduce, but filters exist. Stopping." << std::endl;
                mcs_i_plus_1_candidate.clear_reducedVars_Glob(); // Ensure NReducedVars_Glob = 0
                // Save empty file for consistency
                std::ofstream empty_fout(output_reduced_patterns_filename);
                if (empty_fout.is_open()) { empty_fout << "NReducedVars_Glob = 0" << std::endl; empty_fout.close(); }
            }
            else if (patterns_fed_to_current_phase1 > 0 || mcs_i_plus_1_candidate.Nform1_Glob == 0) { // Proceed if patterns exist or no filters to refine anyway
                mcs_i_plus_1_candidate.perform_phase1_pattern_reduction(mcs_i_ref, current_input_patterns_file, output_reduced_patterns_filename);
            }


            int n_reduced_vars_this_step = mcs_i_plus_1_candidate.NReducedVars_Glob;
            int n_filters_at_start_of_phase2 = mcs_i_plus_1_candidate.Nform1_Glob;

            if (n_filters_at_start_of_phase2 == 0) { std::cout << "      MCS_" << pos_idx + 1 << " candidate has 0 filters before Phase 2. Stopping." << std::endl; break; }
            if (n_reduced_vars_this_step == 0 && n_filters_at_start_of_phase2 > 0) {
                std::cout << "      No patterns remained after Phase 1 for MCS_" << pos_idx + 1 << ". Stopping." << std::endl; break;
            }

            int prunes_in_phase2 = 0;
            if (n_filters_at_start_of_phase2 > 0 && n_reduced_vars_this_step > 0) {
                std::cout << "      Iteratively refining MCS_" << pos_idx + 1 << " (filters: " << n_filters_at_start_of_phase2
                    << ", guiding patterns: " << n_reduced_vars_this_step << ")" << std::endl;
                while (mcs_i_plus_1_candidate.Nform1_Glob > 0 && mcs_i_plus_1_candidate.perform_phase2_single_filter_prune_step(global_rng_engine)) {
                    prunes_in_phase2++;
                }
                std::cout << "      Refinement for MCS_" << pos_idx + 1 << " complete (" << prunes_in_phase2 << " prunes). Final filters: " << mcs_i_plus_1_candidate.Nform1_Glob << std::endl;
            }
            else if (n_filters_at_start_of_phase2 > 0) {
                std::cout << "      Skipping Phase 2 for MCS_" << pos_idx + 1 << " (no guiding patterns)." << std::endl;
            }

            if (mcs_i_plus_1_candidate.Nform1_Glob == 0) { std::cout << "      MCS_" << pos_idx + 1 << " became empty. Stopping." << std::endl; break; }

            positional_mcs_sequence_storage.push_back(mcs_i_plus_1_candidate);
            mcs_i_plus_1_candidate.save_forms_to_file(pa_filename_prefix + "_Filters_Pos" + std::to_string(pos_idx + 1) + ".txt", "Positional MCS");

            bool no_filters_pruned_this_step = (prunes_in_phase2 == 0);
            bool no_patterns_deleted_this_phase1 = (n_reduced_vars_this_step == patterns_fed_to_current_phase1);

            current_input_patterns_file = output_reduced_patterns_filename; // Output of this phase 1 is input for next
            patterns_fed_to_current_phase1 = n_reduced_vars_this_step;

            if (no_filters_pruned_this_step && no_patterns_deleted_this_phase1 && patterns_fed_to_current_phase1 > 0) { // check patterns_fed > 0
                std::cout << "      Stabilized at Positional Set " << pos_idx + 1 << ". Stopping." << std::endl; break;
            }
            if (patterns_fed_to_current_phase1 == 0 && positional_mcs_sequence_storage.back().Nform1_Glob > 0) { // Check last added MCS
                std::cout << "      Stopping: No patterns to guide further refinement, but filters exist in MCS_" << pos_idx + 1 << "." << std::endl; break;
            }
        }
        auto timer_end_pa_gen_seq = std::chrono::high_resolution_clock::now();
        auto duration_pa_gen_seq = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_pa_gen_seq - timer_start_pa_gen_seq);
        std::cout << "    Positional MCS sequence generation for " << m_str << " complete. Generated "
            << positional_mcs_sequence_storage.size() << " sets (including MCS_0)."
            << " (Time: " << duration_pa_gen_seq.count() << " ms)" << std::endl;
        total_times_preprocess["Pos_MCS_" + m_str] = duration_base_patterns + duration_mcs0_gen + duration_map_build + duration_pa_gen_seq;

        std::cout << "  3. Searching with Positional MCS (" << m_str << ")..." << std::endl;
        long long current_pos_mcs_matches_count_total = 0; // Renamed
        std::string pos_mcs_matches_filename = "positional_mcs_matches_" + m_str + ".txt";
        // Clear file and write header
        std::ofstream pos_outfile(pos_mcs_matches_filename, std::ios_base::trunc);
        if (pos_outfile.is_open()) {
            pos_outfile << "QueryID\tNumMatches\tMatchPosition_1\tMatchPosition_2\t..." << std::endl;
            pos_outfile.close();
        }
        else {
            std::cerr << "ERROR: Could not open " << pos_mcs_matches_filename << " to write header for M=" << m_val << std::endl;
        }

        auto timer_start_pa_search = std::chrono::high_resolution_clock::now();
        int query_idx_pos = 0;
        for (const auto& q_str : queries) {
            std::set<int> current_query_matches = search_positional_mcs(q_str, positional_mcs_sequence_storage, usual_mcs_filter_map, main_text, N_PARAM_VAL, MATCH_SIMILARITY_THRESHOLD);
            current_pos_mcs_matches_count_total += current_query_matches.size();

            save_match_positions(pos_mcs_matches_filename, "Query_" + std::to_string(query_idx_pos), current_query_matches);
            query_idx_pos++;
        }
        auto timer_end_pa_search = std::chrono::high_resolution_clock::now();
        total_times_search["Pos_MCS_" + m_str] = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_pa_search - timer_start_pa_search);
        total_matches_found["Pos_MCS_" + m_str] = current_pos_mcs_matches_count_total; // Update map
        std::cout << "     Positional MCS Search for " << m_str << " complete. Matches: " << current_pos_mcs_matches_count_total
            << " (Search Time: " << total_times_search["Pos_MCS_" + m_str].count() << " ms)" << std::endl;
    } // End of M-value loop

    // Naive Search
    std::cout << "\n--- PART C: Naive Search Algorithm ---" << std::endl;
    std::string naive_matches_filename = "naive_matches_details.txt";
    bool naive_search_performed = false;

    if (file_exists(naive_matches_filename)) {
        std::cout << "  Naive search output file (" << naive_matches_filename << ") found. Skipping naive search." << std::endl;
        // Optionally, you could try to parse total_matches_found["Naive"] from the file if needed for comparison,
        // but for now, we'll just mark it as skipped for timing.
        total_matches_found["Naive"] = -2; // Indicate skipped due to existing file
    }
    else {
        naive_search_performed = true;
        long long current_naive_matches_count_total = 0;
        std::ofstream naive_outfile_header(naive_matches_filename, std::ios_base::trunc); // Create/Truncate file for header
        if (naive_outfile_header.is_open()) {
            naive_outfile_header << "QueryID\tNumMatches\tMatchPosition_1\tMatchPosition_2\t..." << std::endl;
            naive_outfile_header.close();
        }
        else {
            std::cerr << "ERROR: Could not open " << naive_matches_filename << " to write header." << std::endl;
        }

        auto timer_start_naive_search = std::chrono::high_resolution_clock::now();
        int query_idx_naive = 0;
        const int total_queries_for_naive = queries.size();
        const int naive_progress_interval = std::max(1, total_queries_for_naive / 100);

        for (const auto& q_str : queries) {
            std::set<int> current_query_matches = search_naive(q_str, main_text, N_PARAM_VAL, MATCH_SIMILARITY_THRESHOLD);
            current_naive_matches_count_total += current_query_matches.size();

            // Save matches for this query
            save_match_positions(naive_matches_filename, "Query_" + std::to_string(query_idx_naive), current_query_matches);

            query_idx_naive++;
            if (query_idx_naive % naive_progress_interval == 0 || query_idx_naive == total_queries_for_naive) {
                auto current_time_naive = std::chrono::high_resolution_clock::now();
                auto elapsed_ms_naive = std::chrono::duration_cast<std::chrono::milliseconds>(current_time_naive - timer_start_naive_search).count();
                double percent_done = (static_cast<double>(query_idx_naive) / total_queries_for_naive) * 100.0;

                long long estimated_total_ms = 0;
                if (query_idx_naive > 0 && percent_done > 0.1) { // Avoid division by zero and unstable early estimates
                    estimated_total_ms = static_cast<long long>((elapsed_ms_naive / (percent_done / 100.0)));
                }

                std::cout << "\r  Naive Search: Processed query " << std::setw(5) << query_idx_naive << "/" << total_queries_for_naive
                    << " (" << std::fixed << std::setprecision(1) << percent_done << "%)"
                    << ". Elapsed: " << elapsed_ms_naive / 1000.0 << "s."
                    << (estimated_total_ms > 0 ? " Est. Total: " + std::to_string(estimated_total_ms / 1000.0) + "s." : "")
                    << std::flush;
            }
        }
        std::cout << std::endl; // Newline after the progress indicator finishes

        auto timer_end_naive_search = std::chrono::high_resolution_clock::now();
        total_times_search["Naive"] = std::chrono::duration_cast<std::chrono::milliseconds>(timer_end_naive_search - timer_start_naive_search);
        total_matches_found["Naive"] = current_naive_matches_count_total; // Update map
        std::cout << "  Naive Search complete. Total Matches: " << current_naive_matches_count_total
            << " (Total Search Time: " << total_times_search["Naive"].count() << " ms)" << std::endl;
    }

    std::string conclusions_filename = "experiment_conclusions.txt";
    std::ofstream conclusion_fout(conclusions_filename);
    if (!conclusion_fout.is_open()) {
        std::cerr << "ERROR: Could not open " << conclusions_filename << " to write conclusions." << std::endl;
    }
    else {
        std::cout << "\nWriting conclusions to " << conclusions_filename << "..." << std::endl;
        conclusion_fout << "--- MCS Algorithm Comparison Experiment Conclusions ---" << std::endl;
        conclusion_fout << "Parameters:" << std::endl;
        conclusion_fout << "  X_PARAM: " << X_PARAM << std::endl;
        conclusion_fout << "  N_PARAM_VAL (W - Query/Window Size): " << N_PARAM_VAL << std::endl;
        conclusion_fout << "  K_PARAM_VAL (Min 1s in Base Patterns C(N,K)): " << K_PARAM_VAL << std::endl;
        conclusion_fout << "  MATCH_SIMILARITY_THRESHOLD: " << MATCH_SIMILARITY_THRESHOLD * 100 << "%" << std::endl;
        conclusion_fout << "  M_VALUES_PARAM (1s in Filters): {";
        for (size_t i = 0; i < M_VALUES_PARAM.size(); ++i) {
            conclusion_fout << M_VALUES_PARAM[i] << (i == M_VALUES_PARAM.size() - 1 ? "" : ", ");
        }
        conclusion_fout << "}" << std::endl;
        conclusion_fout << "  TEXT_SIZE_PARAM: " << TEXT_SIZE_PARAM << std::endl;
        conclusion_fout << "  NUM_QUERIES_PARAM: " << NUM_QUERIES_PARAM << std::endl;
        conclusion_fout << "  ALPHABET_SIZE_Y_PARAM: " << ALPHABET_SIZE_Y_PARAM << std::endl;
        conclusion_fout << "  Base patterns C(" << N_PARAM_VAL - 1 << "," << K_PARAM_VAL - 1 << ") generated: "
            << f_global_pattern_holder.NVars_Glob << std::endl;


        conclusion_fout << "\n--- Performance Summary ---" << std::endl;
        conclusion_fout << std::setw(25) << std::left << "Algorithm"
            << std::setw(20) << std::right << "Preprocessing(ms)"
            << std::setw(15) << std::right << "Search(ms)"
            << std::setw(20) << std::right << "Total Matches" << std::endl;
        conclusion_fout << std::string(80, '-') << std::endl;

        for (const auto& m_val_loop : M_VALUES_PARAM) { // Use a different loop var name
            std::string m_str_loop = "M" + std::to_string(m_val_loop);
            std::string usual_key_loop = "Usual_MCS_" + m_str_loop;
            std::string pos_key_loop = "Pos_MCS_" + m_str_loop;

            if (total_times_search.count(usual_key_loop)) {
                conclusion_fout << std::setw(25) << std::left << usual_key_loop
                    << std::setw(20) << std::right << (total_times_preprocess.count(usual_key_loop) ? total_times_preprocess[usual_key_loop].count() : 0)
                    << std::setw(15) << std::right << total_times_search[usual_key_loop].count()
                    << std::setw(20) << std::right << total_matches_found[usual_key_loop] << std::endl;
            }
            if (total_times_search.count(pos_key_loop)) {
                conclusion_fout << std::setw(25) << std::left << pos_key_loop
                    << std::setw(20) << std::right << (total_times_preprocess.count(pos_key_loop) ? total_times_preprocess[pos_key_loop].count() : 0)
                    << std::setw(15) << std::right << total_times_search[pos_key_loop].count()
                    << std::setw(20) << std::right << total_matches_found[pos_key_loop] << std::endl;
            }
        }
        conclusion_fout << std::setw(25) << std::left << "Naive"
            << std::setw(20) << std::right << 0
            << std::setw(15) << std::right << total_times_search["Naive"].count()
            << std::setw(20) << std::right << total_matches_found["Naive"] << std::endl;

        conclusion_fout << "\n--- Theoretical Time Complexity Discussion ---" << std::endl;
        long long n_q = NUM_QUERIES_PARAM;
        long long t_len = TEXT_SIZE_PARAM;
        long long w_len = N_PARAM_VAL;
        long long n_vars_glob = f_global_pattern_holder.NVars_Glob;

        conclusion_fout << "Key: Q=NumQueries (" << n_q << "), T=TextLength (" << t_len
            << "), W=QueryLength (" << w_len << "), N_base=NumBasePatterns (" << n_vars_glob << ")" << std::endl;
        conclusion_fout << "     N_filters_M=NumFiltersForM, Span_avg=AvgFilterSpan, Hits_cand=CandidateHitsFromMap" << std::endl;

        conclusion_fout << "\nNaive Search:" << std::endl;
        conclusion_fout << "  Theory: O(Q * T * W)" << std::endl;
        conclusion_fout << "  Values: O(" << n_q << " * " << t_len << " * " << w_len << ") = O("
            << n_q * t_len * w_len << ")" << std::endl; // This product can be huge, might overflow std::to_string or be unreadable.
        // Better to keep as factors or use scientific notation if very large.
        // For now, let's print the factors.
        conclusion_fout << "  Calculated: O(" << n_q << " * " << t_len << " * " << w_len << ")" << std::endl;


        conclusion_fout << "\nMCS-based Methods (General Components):" << std::endl;
        conclusion_fout << "  Base Pattern Gen (create_mas1): Using prev_permutation for C(W-1, K-1)" << std::endl;
        conclusion_fout << "    Theory: Roughly O( (W-1) * C(W-1, K-1) ) if permutation is dominant." << std::endl;
        // C(20,14) = 38760. 20 * 38760 = 775200. This is feasible.
        conclusion_fout << "    Values: C(" << w_len - 1 << ", " << K_PARAM_VAL - 1 << ") = " << n_vars_glob
            << ". Approx O(" << (w_len - 1) * n_vars_glob << ")" << std::endl;

        conclusion_fout << "  MCS_0 Gen (chetv_struct_Generation) per M:" << std::endl;
        conclusion_fout << "    Theory: O(N_base * W * AvgFilterStrSetInsert)" << std::endl;
        conclusion_fout << "    Values: O(" << n_vars_glob << " * " << w_len << " * log(N_filters_M_avg))" << std::endl;

        conclusion_fout << "  Forms1_anal1_Global_Refinement (if run) per M:" << std::endl;
        conclusion_fout << "    Theory (Professor's version): Complex. Pass1 is O(N_base * N_filters_M * W * Span_avg). Pass2 is O(N_filters_M * N_base * N_filters_M_active * W * Span_avg)." << std::endl;
        conclusion_fout << "    This is very intensive. N_filters_M_active also changes. Dominated by N_base * (N_filters_M)^2 * W * Span_avg roughly for Pass2." << std::endl;


        conclusion_fout << "  Filter Map Building per M:" << std::endl;
        conclusion_fout << "    Theory: O(T * N_filters_M * Span_avg * AvgMapInsert)" << std::endl;
        conclusion_fout << "    Values: O(" << t_len << " * N_filters_M_avg * Span_avg_M * log(UniqueMaskedWords_M))" << std::endl;

        conclusion_fout << "  Search (Usual MCS) per M:" << std::endl;
        conclusion_fout << "    Theory: O(Q * (N_filters_M * W_filter_apply + Hits_cand_avg_M * W_similarity_check))" << std::endl;

        conclusion_fout << "  Positional MCS Sequence Generation per M:" << std::endl;
        conclusion_fout << "    Phase 1 (perform_phase1_pattern_reduction) per Positional Step:" << std::endl;
        conclusion_fout << "      Theory: O(N_patterns_input_to_step * N_filters_MCS_i * Span_avg_MCS_i)" << std::endl;
        conclusion_fout << "    Phase 2 (perform_phase2_single_filter_prune_step loop) per Positional Step:" << std::endl;
        conclusion_fout << "      Theory: O(N_filters_candidate * N_reduced_patterns_this_step * N_glob_pattern * N_filters_candidate * AvgFilterSpan)"
            << " (Rough, as it's iterative with random choices)" << std::endl;

        conclusion_fout << "  Search (Positional MCS) per M:" << std::endl;
        conclusion_fout << "    Theory: O(Q * Sum_over_Pos_Sets(N_filters_in_set * W_filter_apply) + Hits_cand_avg_M0 * W_similarity_check)" << std::endl;

        conclusion_fout << "\nNote: N_filters_M, Span_avg, UniqueMaskedWords, Hits_cand are data-dependent and vary per M." << std::endl;
        conclusion_fout << "The 'Values' section for MCS provides the symbolic formula with parameters; actual numbers would require per-M stats (e.g., actual N_filters for M3, M4, M5)." << std::endl;

        conclusion_fout.close();
        std::cout << "Conclusions saved to " << conclusions_filename << std::endl;
    }
    std::cout << "\nAll processing complete." << std::endl;
    return 0;
}