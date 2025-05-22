// make_positions.cpp
#include "make_positions.h" 
#include <iostream>
#include <fstream>
#include <iomanip> // For std::setw, std::flush

void build_filter_map(
    const Forms& f_mcs,
    const std::string& text_to_map,
    FilterMapCollection& out_filter_map_data,
    const std::string& output_map_filename,
    int m_value_for_logging)
{
    std::cout << "  build_filter_map: Building filter map for M=" << m_value_for_logging
        << " using " << f_mcs.Nform1_Glob << " filters..." << std::endl;

    out_filter_map_data.word_to_index_map.clear();
    out_filter_map_data.entries_vector.clear();

    if (f_mcs.Nform1_Glob <= 0 || f_mcs.forms1Glob == nullptr) {
        std::cerr << "    WARNING (build_filter_map): No filters in MCS to build map for M="
            << m_value_for_logging << ". Map will be empty." << std::endl;
        if (!output_map_filename.empty()) {
            std::ofstream outfile_map(output_map_filename);
            if (outfile_map.is_open()) {
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

    for (int filter_idx = 0; filter_idx < f_mcs.Nform1_Glob; ++filter_idx) {
        if (f_mcs.forms1Glob[filter_idx] == nullptr) {
            std::cerr << "    Warning (build_filter_map): Filter " << filter_idx << " (M=" << m_value_for_logging
                << ") data is null. Skipping." << std::endl;
            continue;
        }
        int current_filter_span = f_mcs.forms1Glob[filter_idx][f_mcs.N_glob];

        if (current_filter_span <= 0 || current_filter_span > f_mcs.N_glob) {
            std::cerr << "    Warning (build_filter_map): Filter " << filter_idx << " (M=" << m_value_for_logging
                << ") has invalid span " << current_filter_span
                << ". Max expected: " << f_mcs.N_glob << ". Skipping this filter." << std::endl;
            continue;
        }

        for (long long text_pos = 0; text_pos <= text_length - current_filter_span; ++text_pos) {
            std::string masked_word_str = "";
            masked_word_str.reserve(current_filter_span);
            for (int k = 0; k < current_filter_span; ++k) {
                if (f_mcs.forms1Glob[filter_idx][k] == 1) {
                    masked_word_str += text_to_map[text_pos + k];
                }
                else {
                    masked_word_str += '_';
                }
            }
            auto it_map = out_filter_map_data.word_to_index_map.find(masked_word_str);
            if (it_map == out_filter_map_data.word_to_index_map.end()) {
                MapEntryData new_entry;
                new_entry.masked_word = masked_word_str;
                new_entry.occurrences.push_back(static_cast<int>(text_pos));
                out_filter_map_data.entries_vector.push_back(new_entry);
                out_filter_map_data.word_to_index_map[masked_word_str] = out_filter_map_data.entries_vector.size() - 1;
            }
            else {
                out_filter_map_data.entries_vector[it_map->second].occurrences.push_back(static_cast<int>(text_pos));
            }
        }
        if ((filter_idx + 1) % 100 == 0 || filter_idx == f_mcs.Nform1_Glob - 1) {
            std::cout << "\r    build_filter_map: Processed filter " << std::setw(5) << filter_idx + 1 << "/" << f_mcs.Nform1_Glob
                << " for M=" << m_value_for_logging << ". Current unique masked words: "
                << out_filter_map_data.word_to_index_map.size() << std::flush;
        }
    }
    std::cout << std::endl;

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