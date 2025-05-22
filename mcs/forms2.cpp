// forms2.cpp
#include "main.h"
#include <vector>
#include <string>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>   // For std::ostringstream in Forms1_anal1
#include <iomanip>   // For std::setw, std::flush

// --- chetv_struct_Generation ---
void Forms::chetv_struct_Generation(const std::string& output_initial_forms_filename, const Forms& pattern_source) {
    std::cout << "chetv_struct_Generation: Generating initial MCS (M=" << this->Nsovp1_Glob
        << ") from " << pattern_source.NVars_Glob << " basic patterns..." << std::endl;

    if (pattern_source.Vars_Glob == nullptr || pattern_source.NVars_Glob == 0) {
        std::cerr << "  ERROR (chetv_struct): No basic patterns (pattern_source.Vars_Glob) available." << std::endl;
        this->Nform1_Glob = 0; this->forms1Glob = nullptr; return;
    }
    if (this->Nsovp1_Glob <= 0) {
        std::cerr << "  ERROR (chetv_struct): Invalid M value (Nsovp1_Glob=" << this->Nsovp1_Glob << ")." << std::endl;
        this->Nform1_Glob = 0; this->forms1Glob = nullptr; return;
    }

    clear_forms1Glob();
    this->form1Size = this->N_glob + 1; // N_glob of 'this' object, for filter storage

    std::vector<int*> collected_filters_list;
    std::set<std::string> unique_candidate_filter_strings;
    int M = this->Nsovp1_Glob;

    for (int i = 0; i < pattern_source.NVars_Glob; ++i) {
        const bool* current_basic_pattern = pattern_source.Vars_Glob[i];
        std::string candidate_filter_str_derived = "";
        int ones_in_candidate = 0;
        int candidate_span = 0;

        for (int bit_idx = 0; bit_idx < pattern_source.N_glob; ++bit_idx) { // Use pattern_source.N_glob for length
            candidate_span++;
            if (current_basic_pattern[bit_idx]) {
                candidate_filter_str_derived += '1'; ones_in_candidate++;
            }
            else {
                candidate_filter_str_derived += '0';
            }
            if (ones_in_candidate == M) break;
        }

        if (ones_in_candidate < M) continue;

        if (unique_candidate_filter_strings.find(candidate_filter_str_derived) == unique_candidate_filter_strings.end()) {
            unique_candidate_filter_strings.insert(candidate_filter_str_derived);
            int* new_filter_data = new int[this->form1Size]();
            for (int bit_idx = 0; bit_idx < candidate_span; ++bit_idx) {
                new_filter_data[bit_idx] = (candidate_filter_str_derived[bit_idx] == '1' ? 1 : 0);
            }
            new_filter_data[this->N_glob] = candidate_span; // Span stored using this->N_glob
            collected_filters_list.push_back(new_filter_data);
        }
        if (i > 0 && i % 10000 == 0) {
            std::cout << "\r    chetv_struct: Processed " << i << "/" << pattern_source.NVars_Glob
                << " basic patterns. Current unique filters: " << collected_filters_list.size() << std::flush;
        }
    }
    std::cout << "\r    chetv_struct: Processed " << pattern_source.NVars_Glob << "/" << pattern_source.NVars_Glob
        << " basic patterns. Total unique filters generated: " << collected_filters_list.size() << std::endl;

    this->Nform1_Glob = collected_filters_list.size();
    if (this->Nform1_Glob > 0) {
        this->forms1Glob = new int* [this->Nform1_Glob];
        for (size_t j = 0; j < collected_filters_list.size(); ++j) {
            this->forms1Glob[j] = collected_filters_list[j];
        }
    }
    else {
        this->forms1Glob = nullptr;
    }
    this->save_forms_to_file(output_initial_forms_filename, "Initial MCS_0 (chetv_struct)");
    std::cout << "chetv_struct_Generation: Finished." << std::endl;
}

// --- save_forms_to_file ---
void Forms::save_forms_to_file(const std::string& filename, const std::string& comment_prefix) const {
    if (this->Nform1_Glob <= 0 || this->forms1Glob == nullptr) {
        std::cout << "    Note: MCS for " << comment_prefix << " (" << filename << ") is empty or not allocated, not saving." << std::endl;
        std::ofstream fout_empty(filename);
        if (fout_empty.is_open()) {
            fout_empty << "# MCS data is empty for " << comment_prefix << std::endl;
            fout_empty.close();
        }
        return;
    }
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        std::cerr << "    ERROR: Could not open file to save MCS: " << filename << std::endl;
        return;
    }
    fout << "# " << comment_prefix << std::endl;
    fout << "Nform1 = " << this->Nform1_Glob
        << "\tform1Size = " << this->form1Size
        << "\tNsovp1 = " << this->Nsovp1_Glob
        << "\tN_2 = " << this->N_2_glob
        << "\tN_glob = " << this->N_glob << std::endl;
    for (int i = 0; i < this->Nform1_Glob; ++i) {
        fout << "\n>" << i << "\t";
        for (int j = 0; j < this->form1Size; ++j) {
            fout << this->forms1Glob[i][j] << " ";
        }
    }
    fout.close();
    std::cout << "    " << comment_prefix << " (" << this->Nform1_Glob << " filters) saved to " << filename << std::endl;
}

// --- Read_Forms1 ---
void Forms::Read_Forms1(std::string fname) {
    std::ifstream fin(fname);
    if (!fin.is_open()) {
        std::cerr << "  ERROR (Read_Forms1): Could not open file: " << fname << std::endl;
        this->Nform1_Glob = 0; this->forms1Glob = nullptr; return;
    }
    std::cout << "Read_Forms1: Reading filters from " << fname << "..." << std::endl;
    clear_forms1Glob(); // Clear existing before loading

    if (!FindWord(fin, "Nform1 =")) { std::cerr << "Read_Forms1: 'Nform1 =' not found." << std::endl; fin.close(); return; }
    fin >> this->Nform1_Glob;
    if (!FindWord(fin, "form1Size =")) { std::cerr << "Read_Forms1: 'form1Size =' not found." << std::endl; fin.close(); return; }
    fin >> this->form1Size;
    if (!FindWord(fin, "Nsovp1 =")) { std::cerr << "Read_Forms1: 'Nsovp1 =' not found." << std::endl; fin.close(); return; }
    fin >> this->Nsovp1_Glob;
    if (!FindWord(fin, "N_2 =")) { /* Optional, might not always be in all files */ FindWord(fin, "N_glob ="); } // Try to find N_glob if N_2 is missing
    fin >> this->N_2_glob; // This might read N_glob value if N_2 was missing and N_glob followed. Be careful.
    // A more robust header parsing is needed if format varies. For now, assume structure.
    // It's better if FindWord also consumes the value or if we parse line by line.
    // Let's assume N_glob is also in header or can be derived (form1Size -1)

    // Try to read N_glob if present
    std::string temp_line;
    std::getline(fin, temp_line); // Read rest of the header line
    size_t n_glob_pos = temp_line.find("N_glob =");
    if (n_glob_pos != std::string::npos) {
        std::string n_glob_str = temp_line.substr(n_glob_pos + 8); // Length of "N_glob = "
        try {
            this->N_glob = std::stoi(n_glob_str);
        }
        catch (const std::exception& e) {
            std::cerr << "  Warning (Read_Forms1): Could not parse N_glob from header. Relying on form1Size." << std::endl;
            if (this->form1Size > 0) this->N_glob = this->form1Size - 1; else this->N_glob = 0;
        }
    }
    else {
        if (this->form1Size > 0) this->N_glob = this->form1Size - 1; else this->N_glob = 0;
    }


    std::cout << "  Read Header: Nform1=" << Nform1_Glob << ", form1Size=" << form1Size
        << ", Nsovp1(M)=" << Nsovp1_Glob << ", N_2(K)=" << N_2_glob << ", N_glob(N)=" << N_glob << std::endl;

    if (Nform1_Glob <= 0 || form1Size <= 0) {
        std::cerr << "  Read_Forms1: Invalid header values. Nform1 or form1Size is zero or negative." << std::endl;
        Nform1_Glob = 0; fin.close(); return;
    }

    this->forms1Glob = new int* [Nform1_Glob];
    for (int i = 0; i < Nform1_Glob; ++i) {
        this->forms1Glob[i] = new int[form1Size]();
    }
    char dummy_char; int dummy_idx;
    for (int n = 0; n < Nform1_Glob; ++n) {
        if (!(fin >> dummy_char >> dummy_idx)) {
            std::cerr << "  Read_Forms1: Error reading filter line prefix for filter " << n << std::endl;
            Nform1_Glob = n; fin.close(); return;
        }
        if (dummy_char != '>') {
            std::cerr << "  Read_Forms1: Expected '>' at start of filter line " << n << ", found '" << dummy_char << "'" << std::endl;
            Nform1_Glob = n; fin.close(); return;
        }
        for (int i = 0; i < form1Size; ++i) {
            if (!(fin >> this->forms1Glob[n][i])) {
                std::cerr << "  Read_Forms1: Error reading data for filter " << n << ", item " << i << std::endl;
                Nform1_Glob = n; fin.close(); return;
            }
        }
    }
    fin.close();
    std::cout << "Read_Forms1: Successfully read " << Nform1_Glob << " filters." << std::endl;
}


// --- Forms1_anal1_Global_Refinement ---
// Performs a global refinement of the filters currently in this->forms1Glob.
// It uses the full set of base patterns from global_pattern_source.Vars_Glob.
// The logic is based on a two-pass system to identify and keep "essential" filters.
// Pass 1: Marks filters as essential if they are the *sole* cover for any base pattern.
// Pass 2: Iteratively and tentatively removes non-essential filters, and if this removal
//         causes another filter to become the sole cover for a base pattern, that other
//         filter is then promoted to essential.
// Finally, only filters marked essential (mas4[i] == 1) are kept.
// Assuming Forms class definition is available

// Assuming Forms class definition is available

void Forms::Forms1_anal1_Global_Refinement(
	const std::string& input_initial_forms_filename, // Original filename, used for Read_Forms1
	const Forms& global_pattern_source,
	const std::string& output_refined_forms_filename) // Target filename for results
{
	// Use class members for K (N_2_glob) and M (Nsovp1_Glob)
	int current_K = this->N_2_glob;     // K from the object whose filters are being refined
	int current_M = this->Nsovp1_Glob;  // M from the object whose filters are being refined

	// These should ideally match the K and N of the global_pattern_source
	// but we use global_pattern_source members directly for clarity with patterns.
	int pattern_N = global_pattern_source.N_glob; // Length of patterns in Vars_Glob
	int pattern_K = global_pattern_source.N_2_glob; // Number of 1s in patterns in Vars_Glob

	std::cerr << "\nForms1_anal1_Global_Refinement: Starting for M=" << current_M
		<< " (filters) and K=" << pattern_K << ", N=" << pattern_N << " (base patterns)." << std::endl;

	// 1. Read the initial set of forms that this object is supposed to refine.
	this->Read_Forms1(input_initial_forms_filename);
	if (this->Nform1_Glob == 0) {
		std::cerr << "  ERROR (Forms1_anal1): No forms loaded from " << input_initial_forms_filename << " to refine." << std::endl;
		return;
	}
	if (global_pattern_source.Vars_Glob == nullptr || global_pattern_source.NVars_Glob == 0) {
		std::cerr << "  ERROR (Forms1_anal1): No global base patterns provided for refinement." << std::endl;
		return;
	}

	std::cerr << "  Read " << this->Nform1_Glob << " initial filters. Refining against "
		<< global_pattern_source.NVars_Glob << " base patterns." << std::endl;

	std::ostringstream temp_fname_ss;
	temp_fname_ss << "dbg_forms_c_K" << pattern_K << "_M" << current_M << ".txt";
	this->save_forms_to_file(temp_fname_ss.str(), "Initial forms for refinement (debug copy)");
	temp_fname_ss.str("");

	int filter_storage_N = this->N_glob;

	int* mas3 = new int[this->Nform1_Glob]();
	int* mas4 = new int[this->Nform1_Glob]();

	// --- First Pass - Identify Initially Essential Filters ---
	std::cerr << "  Forms1_anal1 Pass 1 (M=" << current_M << "): Identifying essential filters..." << std::endl;
	for (int i = 0; i < global_pattern_source.NVars_Glob; ++i) {
		if (i % 1000 == 0 && i > 0) {
			std::cerr << "\r    Pass 1: Processing base pattern " << i << "/" << global_pattern_source.NVars_Glob << std::flush;
		}
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) mas3[i1] = 0;

		int covering_filters_count_for_pattern_i = 0;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			int filter_span = this->forms1Glob[i1][filter_storage_N];
			if (filter_span <= 0) continue;
			for (int j1 = 0; j1 <= pattern_N - filter_span; ++j1) {
				if (global_pattern_source.Vars_Glob[i][j1] == false && this->forms1Glob[i1][0] == 1) continue;
				bool current_alignment_matches = true;
				for (int j2 = 0; j2 < filter_span; ++j2) {
					if (this->forms1Glob[i1][j2] == 1) {
						if (global_pattern_source.Vars_Glob[i][j1 + j2] == false) {
							current_alignment_matches = false; break;
						}
					}
				}
				if (current_alignment_matches) {
					mas3[i1] = 1;
					covering_filters_count_for_pattern_i++; // Not strictly needed here, but was in original
					break;
				}
			}
		}
		int k_sum_mas3 = 0;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) k_sum_mas3 += mas3[i1];
		if (k_sum_mas3 == 1) {
			for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
				if (mas3[i1] == 1) {
					mas4[i1] = 1;
					// DEBUG: Optional print for Pass 1 essential
					// std::cerr << "\n    DEBUG Pass 1: Filter " << i1 << " is essential for pattern " << i << std::endl;
				}
			}
		}
	}
	std::cerr << "\n  Forms1_anal1 Pass 1 (M=" << current_M << ") complete." << std::endl;

	temp_fname_ss << "dbg_forms_an1_K" << pattern_K << "_M" << current_M << ".txt";
	// ... (file saving logic for an1) ... (omitted for brevity, same as before)
	std::ofstream fout_an1_dbg(temp_fname_ss.str()); // Renamed to avoid conflict if original fout_an1 is used
	if (fout_an1_dbg.is_open()) {
		fout_an1_dbg << "Nform1 = " << this->Nform1_Glob << "\tform1Size = " << this->form1Size
			<< "\tNsovp1 = " << this->Nsovp1_Glob << "\tN_2(K_pattern) = " << pattern_K;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			fout_an1_dbg << "\n>" << i1 << "\t";
			for (int j1 = 0; j1 < this->form1Size; ++j1) fout_an1_dbg << this->forms1Glob[i1][j1] << " ";
			fout_an1_dbg << "\t" << mas4[i1]; // Essentiality marker
		}
		fout_an1_dbg.close();
		std::cerr << "  Debug info (an1) for M=" << current_M << " written to " << temp_fname_ss.str() << std::endl;
	}
	temp_fname_ss.str("");


	// --- Second Pass - Iterative Refinement ---
	std::cerr << "  Forms1_anal1 Pass 2 (M=" << current_M << ") starting iterative refinement..." << std::endl;
	for (int ii1 = 0; ii1 < this->Nform1_Glob; ++ii1) {
		// Regular progress update
		// Replaced the \r version with a new line for each main filter ii1, to not interfere with detailed debug prints
		std::cerr << "\n  --- Pass 2: Iteration for filter ii1 = " << ii1 << " / " << (this->Nform1_Glob - 1)
			<< " (Status: mas4[" << ii1 << "]=" << mas4[ii1] << ") ---" << std::endl;

		if (mas4[ii1] == 1) {
			std::cerr << "    DEBUG[ii1=" << ii1 << "]: Already essential from Pass 1. Skipping." << std::endl;
			continue;
		}

		int previous_mas4_ii1_val = mas4[ii1]; // Should be 0 if not essential from Pass 1
		mas4[ii1] = -1;
		std::cerr << "    DEBUG[ii1=" << ii1 << "]: Tentatively marked mas4[" << ii1 << "] = -1. (Was " << previous_mas4_ii1_val << ")" << std::endl;

		bool made_change_promoting_to_essential_this_iter = false;
		// Scan 1: Check if removing ii1 promotes other filters
		std::cerr << "      DEBUG[ii1=" << ii1 << "]: Starting Scan 1 (check for promotion of OTHERS due to ii1 removal)." << std::endl;
		for (int i = 0; i < global_pattern_source.NVars_Glob; ++i) {
			if (global_pattern_source.NVars_Glob > 1000 && i % (global_pattern_source.NVars_Glob / 100) == 0 && i > 0) { // Progress for long pattern scans
				std::cerr << "\r        DEBUG[ii1=" << ii1 << "]: Scan 1, pattern " << i << "/" << global_pattern_source.NVars_Glob << std::flush;
			}
			for (int r = 0; r < this->Nform1_Glob; ++r) mas3[r] = 0;

			int active_covering_filters_count = 0;
			for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
				if (mas4[i1] == -1) continue;

				int filter_span = this->forms1Glob[i1][filter_storage_N];
				if (filter_span <= 0) continue;

				for (int j1 = 0; j1 <= pattern_N - filter_span; ++j1) {
					if (global_pattern_source.Vars_Glob[i][j1] == false && this->forms1Glob[i1][0] == 1) continue;

					bool current_alignment_matches = true;
					for (int j2 = 0; j2 < filter_span; ++j2) {
						if (this->forms1Glob[i1][j2] == 1) {
							if (global_pattern_source.Vars_Glob[i][j1 + j2] == false) {
								current_alignment_matches = false; break;
							}
						}
					}
					if (current_alignment_matches) {
						mas3[i1] = 1;
						active_covering_filters_count++;
						break;
					}
				}
			}

			int k_sum_mas3_active = 0;
			int newly_essential_idx = -1;
			for (int r = 0; r < this->Nform1_Glob; ++r) {
				if (mas4[r] != -1 && mas3[r] == 1) {
					k_sum_mas3_active++;
					// If it becomes 1, this 'r' is the one. If more than one, newly_essential_idx will be the last one.
					// This is okay because we only care if k_sum_mas3_active == 1.
					newly_essential_idx = r;
				}
			}

			if (k_sum_mas3_active == 1) {
				if (mas4[newly_essential_idx] == 0) {
					mas4[newly_essential_idx] = 1;
					made_change_promoting_to_essential_this_iter = true;
					std::cerr << "\n        DEBUG[ii1=" << ii1 << "]: Filter " << newly_essential_idx
						<< " (was mas4=0) PROMOTED to essential (mas4=1) because it uniquely covers pattern " << i
						<< " after ii1 was tentatively removed." << std::endl;
				}
				// Optional: else if mas4[newly_essential_idx] == 1, it was already essential and still is.
				// else if mas4[newly_essential_idx] == -1, this case should not happen due to the check `mas4[r] != -1`.
			}
		}
		if (global_pattern_source.NVars_Glob > 1000) std::cerr << "\r                                                                                           \r"; // Clear line
		std::cerr << "      DEBUG[ii1=" << ii1 << "]: Scan 1 finished. made_change_promoting_to_essential_this_iter = " << (made_change_promoting_to_essential_this_iter ? "true" : "false") << std::endl;


		bool ii1_must_be_restored = false;
		if (!made_change_promoting_to_essential_this_iter) {
			std::cerr << "      DEBUG[ii1=" << ii1 << "]: Starting Scan 2 (check if ii1 itself is essential for coverage)." << std::endl;
			for (int i = 0; i < global_pattern_source.NVars_Glob; ++i) {
				if (global_pattern_source.NVars_Glob > 1000 && i % (global_pattern_source.NVars_Glob / 100) == 0 && i > 0) { // Progress for long pattern scans
					std::cerr << "\r        DEBUG[ii1=" << ii1 << "]: Scan 2, pattern " << i << "/" << global_pattern_source.NVars_Glob << std::flush;
				}
				bool pattern_i_is_covered = false;
				for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
					if (mas4[i1] == -1) continue;

					int filter_span = this->forms1Glob[i1][filter_storage_N];
					if (filter_span <= 0) continue;
					for (int j1 = 0; j1 <= pattern_N - filter_span; ++j1) {
						if (global_pattern_source.Vars_Glob[i][j1] == false && this->forms1Glob[i1][0] == 1) continue;
						bool matches = true;
						for (int j2 = 0; j2 < filter_span; ++j2) {
							if (this->forms1Glob[i1][j2] == 1 && global_pattern_source.Vars_Glob[i][j1 + j2] == false) {
								matches = false; break;
							}
						}
						if (matches) { pattern_i_is_covered = true; break; }
					}
					if (pattern_i_is_covered) break;
				}
				if (!pattern_i_is_covered) {
					ii1_must_be_restored = true;
					std::cerr << "\n        DEBUG[ii1=" << ii1 << "]: Filter " << ii1
						<< " (tentatively removed) IS NOW ESSENTIAL as its removal uncovers pattern " << i << "." << std::endl;
					break;
				}
			}
			if (global_pattern_source.NVars_Glob > 1000) std::cerr << "\r                                                                                           \r"; // Clear line
			std::cerr << "      DEBUG[ii1=" << ii1 << "]: Scan 2 finished. ii1_must_be_restored = " << (ii1_must_be_restored ? "true" : "false") << std::endl;
		}
		else {
			std::cerr << "      DEBUG[ii1=" << ii1 << "]: Skipping Scan 2 because another filter was promoted due to ii1's removal." << std::endl;
		}


		// Decision for ii1
		if (ii1_must_be_restored) { // This implies !made_change_promoting_to_essential_this_iter was true
			mas4[ii1] = 1;
			std::cerr << "    DEBUG[ii1=" << ii1 << "]: Filter " << ii1 << " RESTORED AND PROMOTED to essential. Final mas4[" << ii1 << "]=1." << std::endl;
		}
		else {
			// If mas4[ii1] is still -1 here, it means either:
			// 1. made_change_promoting_to_essential_this_iter was TRUE (so ii1's removal was beneficial for promoting another)
			// OR
			// 2. made_change_promoting_to_essential_this_iter was FALSE AND ii1_must_be_restored was FALSE (ii1 is redundant)
			if (mas4[ii1] == -1) { // Check explicitly to only print if it's still -1
				if (made_change_promoting_to_essential_this_iter) {
					std::cerr << "    DEBUG[ii1=" << ii1 << "]: Filter " << ii1 << " REMAINS REMOVED (mas4[" << ii1 << "]=" << mas4[ii1]
						<< ") because its removal promoted another filter." << std::endl;
				}
				else { // This means !made_change and !ii1_must_be_restored
					std::cerr << "    DEBUG[ii1=" << ii1 << "]: Filter " << ii1 << " deemed REDUNDANT and REMAINS REMOVED (mas4[" << ii1 << "]=" << mas4[ii1]
						<< "). Its removal did not promote others and did not uncover any pattern." << std::endl;
				}
			}
			// If mas4[ii1] was somehow changed from -1 by another mechanism (should not happen here), this part is skipped.
		}
	} // End of Pass 2 loop (ii1)
	std::cerr << "\n  Forms1_anal1 Pass 2 (M=" << current_M << ") complete." << std::endl;

	temp_fname_ss << "dbg_forms_an2_K" << pattern_K << "_M" << current_M << ".txt";
	// ... (file saving logic for an2) ... (omitted for brevity, same as before)
	std::ofstream fout_an2_dbg(temp_fname_ss.str()); // Renamed
	if (fout_an2_dbg.is_open()) {
		fout_an2_dbg << "Nform1 = " << this->Nform1_Glob << "\tform1Size = " << this->form1Size
			<< "\tNsovp1 = " << this->Nsovp1_Glob << "\tN_2(K_pattern) = " << pattern_K;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			fout_an2_dbg << "\n>" << i1 << "\t";
			for (int j1 = 0; j1 < this->form1Size; ++j1) fout_an2_dbg << this->forms1Glob[i1][j1] << " ";
			fout_an2_dbg << "\t" << mas4[i1]; // Essentiality marker after Pass 2
		}
		fout_an2_dbg.close();
		std::cerr << "  Debug info (an2) for M=" << current_M << " written to " << temp_fname_ss.str() << std::endl;
	}
	temp_fname_ss.str("");

	// --- Final Selection and Compaction ---
	std::vector<int*> refined_forms_list;
	std::vector<int> original_indices_kept;

	for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
		if (mas4[i1] == 1) {
			refined_forms_list.push_back(this->forms1Glob[i1]);
			original_indices_kept.push_back(i1);
		}
		else {
			delete[] this->forms1Glob[i1];
			this->forms1Glob[i1] = nullptr;
		}
	}

	if (this->forms1Glob != nullptr) {
		delete[] this->forms1Glob;
		this->forms1Glob = nullptr;
	}

	this->Nform1_Glob = refined_forms_list.size();
	if (this->Nform1_Glob > 0) {
		this->forms1Glob = new int* [this->Nform1_Glob];
		for (size_t i = 0; i < this->Nform1_Glob; ++i) {
			this->forms1Glob[i] = refined_forms_list[i];
		}
	}
	else {
		this->forms1Glob = nullptr;
	}

	std::cerr << "  Forms1_anal1: Refinement resulted in " << this->Nform1_Glob << " essential filters for M=" << current_M << "." << std::endl;
	std::cerr << "    Original indices of filters kept: ";
	for (size_t i = 0; i < original_indices_kept.size(); ++i) {
		std::cerr << original_indices_kept[i] << (i == original_indices_kept.size() - 1 ? "" : ", ");
	}
	std::cerr << std::endl;


	this->save_forms_to_file(output_refined_forms_filename, "Globally Refined MCS_0 (Forms1_anal1)");

	std::string ch_output_filename = output_refined_forms_filename;
	size_t dot_pos = ch_output_filename.rfind('.');
	if (dot_pos != std::string::npos) {
		ch_output_filename.insert(dot_pos, "_ch");
	}
	else {
		ch_output_filename += "_ch";
	}
	this->save_forms_to_file(ch_output_filename, "Globally Refined MCS_0 (Forms1_anal1_ch)");

	delete[] mas3;
	delete[] mas4;

	std::cerr << "Forms1_anal1_Global_Refinement for M=" << current_M << " finished.\n";
}