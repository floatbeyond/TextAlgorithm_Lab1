// make_positions.h
#ifndef MAKE_POSITIONS_H
#define MAKE_POSITIONS_H

#include <string>
// main.h (included below) will bring in vector, map, Forms, FilterMapCollection, MapEntryData

#include "main.h" // For Forms class definition and FilterMapCollection structs

void build_filter_map(
    const Forms& f_mcs,
    const std::string& text_to_map,
    FilterMapCollection& out_filter_map_data,
    const std::string& output_map_filename,
    int m_value_for_logging
);

#endif // MAKE_POSITIONS_H