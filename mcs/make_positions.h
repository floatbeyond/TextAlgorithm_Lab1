#ifndef MAKE_POSITIONS_H
#define MAKE_POSITIONS_H

#include <string>
#include <vector>
#include <map>
#include "main.h" // For Forms class definition

// Structures for the map (can be kept here or moved to a common utility header if used elsewhere)
struct MapEntryData // Renamed from www for clarity in this context
{
	std::string masked_word;
	std::vector<int> occurrences;
};

struct FilterMapCollection // Renamed from MP
{
	std::map<std::string, int> word_to_index_map; // Maps masked_word to index in 'entries_vector'
	std::vector<MapEntryData> entries_vector;      // Stores the actual MapEntryData
};

// Function to build and save the filter map
void build_filter_map(
	const Forms& f_mcs,                          // Contains the generated MCS (forms1Glob, Nform1_Glob, N_glob)
	const std::string& text_to_map,              // The large text to process
	FilterMapCollection& out_filter_map_data,    // Output: The populated filter map
	const std::string& output_map_filename,      // Filename to save the resulting map
	int m_value_for_logging                      // For logging which M this map is for
);

#endif // MAKE_POSITIONS_H