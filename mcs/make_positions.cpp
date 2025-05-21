#include "make_positions.h" // Include its own header
#include <iostream>
#include <fstream>
#include <iomanip> // For std::setw

// Note: The Forms class definition is included via "main.h" inside "make_positions.h"

void build_filter_map(
    const Forms& f_mcs,                       // Input MCS object containing the filters
    const std::string& text_to_map,           // The main text to scan
    FilterMapCollection& out_filter_map_data, // Output: The populated filter map
    const std::string& output_map_filename,   // File to save the map (can be empty if no save desired)
    int m_value_for_logging)                  // M-value, used for logging purposes
{
    std::cout << "  build_filter_map: Building filter map for M=" << m_value_for_logging
        << " using " << f_mcs.Nform1_Glob << " filters..." << std::endl;

    // Clear any previous data in the output map
    out_filter_map_data.word_to_index_map.clear();
    out_filter_map_data.entries_vector.clear();

    if (f_mcs.Nform1_Glob <= 0 || f_mcs.forms1Glob == nullptr) { // Check Nform1_Glob <= 0 for robustness
        std::cerr << "    WARNING (build_filter_map): No filters in MCS to build map for M="
            << m_value_for_logging << ". Map will be empty." << std::endl;
        if (!output_map_filename.empty()) {
            std::ofstream outfile_map(output_map_filename);
            if (outfile_map.is_open()) { // Check if file opened successfully
                outfile_map << "# No filters in MCS for M=" << m_value_for_logging << std::endl;
                outfile_map.close();
            }
            else {
                std::cerr << "    ERROR (build_filter_map): Could not open map file for writing empty map note: " << output_map_filename << std::endl;
            }
        }
        return;
    }

    long long text_length = text_to_map.length();

    // Your original DEBUGGING PRINT (useful, can be conditionally compiled)
    // std::cout << "  DEBUG build_filter_map: Text length = " << text_length
    //     << ", Num filters (Nform1_Glob) = " << f_mcs.Nform1_Glob
    //     << ", N_glob (filter array storage length) = " << f_mcs.N_glob << std::endl;

    for (int filter_idx = 0; filter_idx < f_mcs.Nform1_Glob; ++filter_idx) {
        if (f_mcs.forms1Glob[filter_idx] == nullptr) {
            std::cerr << "    Warning (build_filter_map): Filter " << filter_idx << " (M=" << m_value_for_logging
                << ") data is null. Skipping." << std::endl;
            continue;
        }

        // f_mcs.N_glob is the N value used for filter storage (max pattern length for filters)
        // The actual span of this filter is stored at forms1Glob[filter_idx][f_mcs.N_glob]
        int current_filter_span = f_mcs.forms1Glob[filter_idx][f_mcs.N_glob];

        if (current_filter_span <= 0 || current_filter_span > f_mcs.N_glob) {
            std::cerr << "    Warning (build_filter_map): Filter " << filter_idx << " (M=" << m_value_for_logging
                << ") has invalid span " << current_filter_span
                << ". Max expected: " << f_mcs.N_glob << ". Skipping this filter." << std::endl;
            continue;
        }

        // Your original DEBUGGING PRINT for each filter (can be verbose)
        // std::cout << "    DEBUG Filter " << filter_idx << ": Span = " << current_filter_span << " Bits: ";
        // for (int k = 0; k < current_filter_span; ++k) std::cout << f_mcs.forms1Glob[filter_idx][k];
        // std::cout << std::endl;

        // Slide this filter over the text
        // The window in the text will have length current_filter_span
        for (long long text_pos = 0; text_pos <= text_length - current_filter_span; ++text_pos) {
            std::string masked_word_str = "";
            masked_word_str.reserve(current_filter_span); // Pre-allocate for efficiency

            for (int k = 0; k < current_filter_span; ++k) { // k is index within the filter's span
                if (f_mcs.forms1Glob[filter_idx][k] == 1) { // If filter bit is 1 (keep char)
                    masked_word_str += text_to_map[text_pos + k];
                }
                else { // Filter bit is 0 (mask char)
                    masked_word_str += '_'; // Use underscore for masked positions
                }
            }

            // Add the masked word and its position to the map
            auto it_map = out_filter_map_data.word_to_index_map.find(masked_word_str);
            if (it_map == out_filter_map_data.word_to_index_map.end()) { // New masked word encountered
                MapEntryData new_entry;
                new_entry.masked_word = masked_word_str;
                new_entry.occurrences.push_back(static_cast<int>(text_pos));
                out_filter_map_data.entries_vector.push_back(new_entry);
                // Store the index of the new entry in the map
                out_filter_map_data.word_to_index_map[masked_word_str] = out_filter_map_data.entries_vector.size() - 1;
            }
            else { // Masked word already exists, add current position to its occurrences
                out_filter_map_data.entries_vector[it_map->second].occurrences.push_back(static_cast<int>(text_pos));
            }
        } // End of text sliding loop (text_pos)

        // Progress logging
        if ((filter_idx + 1) % 100 == 0 || filter_idx == f_mcs.Nform1_Glob - 1) {
            std::cout << "\r    build_filter_map: Processed filter " << std::setw(5) << filter_idx + 1 << "/" << f_mcs.Nform1_Glob
                << " for M=" << m_value_for_logging << ". Current unique masked words: "
                << out_filter_map_data.word_to_index_map.size() << std::flush; // Use std::flush for \r
        }
    } // End of filter loop (filter_idx)
    std::cout << std::endl; // Ensure newline after progress indicator

    // Your original DEBUGGING PRINTS for final map stats
    std::cout << "  DEBUG build_filter_map FINAL: Total unique masked words = "
        << out_filter_map_data.word_to_index_map.size() << std::endl;
    std::cout << "  DEBUG build_filter_map FINAL: Total entries in vector = "
        << out_filter_map_data.entries_vector.size() << std::endl;

    long long total_occurrences_stored = 0;
    for (const auto& entry : out_filter_map_data.entries_vector) {
        total_occurrences_stored += entry.occurrences.size();
    }
    std::cout << "  DEBUG build_filter_map FINAL: Total integer occurrences stored = "
        << total_occurrences_stored << std::endl;
    std::cout << "  DEBUG build_filter_map FINAL: Estimated memory for occurrences (ints only) = "
        << (total_occurrences_stored * sizeof(int) / (1024.0 * 1024.0)) << " MB" << std::endl;

    // Save the map to file if a filename is provided
    if (!output_map_filename.empty()) {
        std::ofstream outfile_map(output_map_filename);
        if (!outfile_map.is_open()) {
            std::cerr << "    ERROR (build_filter_map): Could not open map file for writing: " << output_map_filename << std::endl;
        }
        else {
            std::cout << "  build_filter_map: Writing filter map for M=" << m_value_for_logging
                << " to " << output_map_filename << "..." << std::endl;
            outfile_map << "# Filter Map for M=" << m_value_for_logging
                << ", N_glob_filters=" << f_mcs.N_glob << ", K_filters=" << f_mcs.N_2_glob << std::endl;
            outfile_map << "# Unique masked words: " << out_filter_map_data.word_to_index_map.size() << std::endl;

            for (size_t i = 0; i < out_filter_map_data.entries_vector.size(); ++i) {
                outfile_map << "\n" << i << " " << out_filter_map_data.entries_vector[i].masked_word << ":";
                for (size_t j = 0; j < out_filter_map_data.entries_vector[i].occurrences.size(); ++j) {
                    outfile_map << " " << out_filter_map_data.entries_vector[i].occurrences[j];
                }
            }
            outfile_map.close();
            std::cout << "  build_filter_map: Filter map for M=" << m_value_for_logging << " written." << std::endl;
        }
    }
    std::cout << "  build_filter_map: Finished. Map contains "
        << out_filter_map_data.word_to_index_map.size() << " unique masked words." << std::endl;
}