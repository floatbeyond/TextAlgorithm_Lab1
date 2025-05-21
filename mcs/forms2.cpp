#ifndef _INC_forms2
#define _INC_forms2

#include "main.h"
#include <vector>
#include <algorithm> // For std::find_if or similar if needed for choosing which to remove
#include <atlstr.h> 
#include <fstream> 

void Forms::save_forms_to_file(const std::string& filename, const std::string& comment_prefix) const {
	if (this->Nform1_Glob == 0 && this->forms1Glob == nullptr) { // Corrected to &&, or check Nform1_Glob <=0
		std::cout << "    Note: MCS for " << comment_prefix << " (" << filename << ") is empty or not allocated, not saving." << std::endl;
		std::ofstream fout_empty(filename);
		if (fout_empty.is_open()) {
			fout_empty << "# MCS data is empty for " << comment_prefix << std::endl;
			fout_empty.close();
		}
		return;
	}

	std::ofstream fout(filename.c_str()); // .c_str() is for older C++ fstream, modern fstream takes std::string directly
	if (!fout.is_open()) {
		std::cerr << "    ERROR: Could not open file to save MCS: " << filename << std::endl;
		return;
	}

	fout << "# " << comment_prefix << std::endl;
	fout << "Nform1 = " << this->Nform1_Glob
		<< "\tform1Size = " << this->form1Size
		<< "\tNsovp1 = " << this->Nsovp1_Glob // M
		<< "\tN_2 = " << this->N_2_glob    // K
		<< "\tN_glob = " << this->N_glob << std::endl; // N

	for (int i = 0; i < this->Nform1_Glob; ++i) {
		fout << "\n>" << i << "\t";
		// int filter_actual_span = this->forms1Glob[i][this->N_glob]; // Not used here but good for context
		// Print all N_glob bits + the span at the end (as per your original file format)
		for (int j = 0; j < this->form1Size; ++j) { // form1Size is N_glob + 1
			fout << this->forms1Glob[i][j] << " ";
		}
	}
	fout.close();
	std::cout << "    " << comment_prefix << " (" << this->Nform1_Glob << " filters) saved to " << filename << std::endl;
}

void Forms::chetv_struct_Generation(const std::string& output_initial_forms_filename, const Forms& pattern_source)
{
	int i, i1, j1, j2, k, n1, Nform, flag1;
	std::ofstream fout;

	// N (N) - pattern length
	// N_2 (K) - number of 1s in pattern
	// Nsovp (M) - number of 1s to look for within the pattern

	int N = pattern_source.N_glob;
	int N_2 = pattern_source.N_2_glob;
	int Nsovp = pattern_source.Nsovp1_Glob;

	CString qqq1, qqq2;

	// "array of forms"
	int NformMax = 60000; // max forms to generate
	int** forms; // "array of forms; indices from 0 to N-1 can be 1 or 0; N - length of the form"

	forms = new int* [NformMax];
	for (i = 0; i < NformMax; i++)
		forms[i] = new int[N + 1];

	for (i = 0; i < NformMax; i++)
		for (j1 = 0; j1 < N; j1++)
			forms[i][j1] = 0;

	for (i = 0; i < NformMax; i++)
	{
		forms[i][N] = -1;
		forms[i][0] = 1;
	}

	k = 0;
	Nform = 0;
	for (i = 0; i < pattern_source.NVars_Glob; i++)
	{
		if (i % 100 == 0)
			std::cerr << "\n " << i;

		// comparison with existing forms
		flag1 = 1;
		for (i1 = 0; i1 < Nform; i1++)
		{
			for (j1 = 0; j1 <= N - forms[i1][N]; j1++)
			{
				if (pattern_source.Vars_Glob[i][j1] == 0) continue;

				flag1 = 0;
				for (j2 = 1; j2 < forms[i1][N]; j2++)
				{
					if (forms[i1][j2] == 0) continue;

					if (pattern_source.Vars_Glob[i][j1 + j2] != 1)
					{
						flag1 = 1;
						break;
					}
				}

				if (flag1 == 0) break;
			}

			if (flag1 == 0) break;
		}

		if (flag1 == 0) continue;

		// if there is no existing form - add a new form (can be optimized)
		for (j1 = 0; j1 < N; j1++)
		{
			if (pattern_source.Vars_Glob[i][j1] == 0) continue;

			n1 = 1;
			for (j2 = 1; j2 < N; j2++)
			{
				forms[Nform][j2] = pattern_source.Vars_Glob[i][j2 + j1];

				if (forms[Nform][j2] == 1)
					n1++;

				if (n1 == Nsovp)
				{
					forms[Nform][N] = j2 + 1; // zapis' razmera
					Nform++;
					break;
				}
			}

			break;
		}

		if (Nform == NformMax)
			break;
	}

	fout.open(output_initial_forms_filename.c_str(), std::ios::out);
	if (!fout.is_open()) {
		std::cerr << "\nERROR: Could not open forms file for writing: " << output_initial_forms_filename << std::endl;
		return;
	}
	fout << "Nform1 = " << Nform << "\t" << "form1Size = " << N + 1 << " Nsovp1 = " << Nsovp << " N_2 = " << N_2;

	for (i1 = 0; i1 < Nform; i1++)
	{
		fout << "\n>" << i1 << "\t";
		for (j1 = 0; j1 < N + 1; j1++)
			fout << forms[i1][j1] << " ";
	}

	fout.close();
	std::cerr << "\n   Initial MCS for M=" << pattern_source.Nsovp1_Glob << " written to " << output_initial_forms_filename;
	std::cerr << "\nk = " << Nform;

	// Copy to forms_Glob
	this->Nform1_Glob = Nform;
	this->form1Size = N + 1;

	this->forms1Glob = new int* [this->Nform1_Glob];
	for (i = 0; i < this->Nform1_Glob; i++)
		this->forms1Glob[i] = new int[this->form1Size];

	for (i1 = 0; i1 < this->Nform1_Glob; i1++)
		for (j1 = 0; j1 < this->form1Size; j1++)
			this->forms1Glob[i1][j1] = forms[i1][j1];

	// Clean up
	for (i = 0; i < NformMax; i++)
		delete[] forms[i];
	delete[] forms;
}



void  Forms::Read_Forms1(std::string fname)
{
	int i,n,k;
	char a,buf[1000];
	
	std::ofstream fout;
	std::ifstream fin;
	
	fin.open(fname.c_str(), std::ios::in);
	if (!fin.good()) return;

	std::cerr << "\nT00";

	FindWord(fin, "Nform1 =");

	std::cerr << "\nT01";
	fin>>Nform1_Glob;

	

	FindWord(fin, "form1Size =");
		fin>>form1Size;

	FindWord(fin, "Nsovp1 =");
		fin>>Nsovp1_Glob;

	FindWord(fin, "N_2 =");
		fin>>N_2_glob;

	std::cerr << "\nT01";

	forms1Glob = new int *[Nform1_Glob];
	for (i = 0; i<Nform1_Glob; i++)
		forms1Glob[i] = new int[form1Size];

	std::cerr<<"\n"<<Nform1_Glob<<" "<<form1Size;
//		cin>>k;

	for (n = 0; n<Nform1_Glob; n++)
		for (i = 0; i<form1Size; i++)
			forms1Glob[n][i] = 0;
	
	for (n = 0; n<Nform1_Glob; n++)
	{
		fin.getline(buf,sizeof(buf));
		
		fin>>a;
		fin>>k;
		
		for (i = 0; i<form1Size; i++)
			fin>>forms1Glob[n][i];

	}

	fin.close();

	std::cerr << "\nT03";
}


#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream> // For std::ostringstream
#include <iomanip> // For std::setw (if needed for cerr)

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
	//    The professor's code constructs the filename based on M.
	//    Let's use input_initial_forms_filename which should be like "Usual_MCS_InitialForms_M3.txt"
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

	// Save a copy of the initial forms (like forms_c_...txt)
	std::ostringstream temp_fname_ss;
	temp_fname_ss << "dbg_forms_c_K" << pattern_K << "_M" << current_M << ".txt";
	this->save_forms_to_file(temp_fname_ss.str(), "Initial forms for refinement (debug copy)");
	temp_fname_ss.str(""); // Clear stream

	// N for filter related indexing (span is at forms1Glob[i1][this->N_glob])
	int filter_storage_N = this->N_glob; // This should be N_PARAM_VAL

	int* mas3 = new int[this->Nform1_Glob](); // Zero-initialize
	int* mas4 = new int[this->Nform1_Glob](); // Zero-initialize (0 means not essential yet)

	// --- First Pass - Identify Initially Essential Filters ---
	std::cerr << "  Forms1_anal1 Pass 1 (M=" << current_M << "): Identifying essential filters..." << std::endl;
	for (int i = 0; i < global_pattern_source.NVars_Glob; ++i) { // Iterate over base patterns
		if (i % 1000 == 0 && i > 0) {
			std::cerr << "\r    Pass 1: Processing base pattern " << i << "/" << global_pattern_source.NVars_Glob << std::flush;
		}

		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			mas3[i1] = 0; // Reset for this pattern
		}

		int covering_filters_count_for_pattern_i = 0;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) { // Iterate over current filters
			int filter_span = this->forms1Glob[i1][filter_storage_N];
			if (filter_span <= 0) continue;

			// Sliding window match of filter i1 within pattern i
			for (int j1 = 0; j1 <= pattern_N - filter_span; ++j1) { // j1 is start pos in pattern
				// Professor's check: pattern must have '1' at the start of this potential match window
				// if the filter implicitly starts with '1' (which they do).
				if (global_pattern_source.Vars_Glob[i][j1] == false && this->forms1Glob[i1][0] == 1) {
					continue;
				}

				bool current_alignment_matches = true;
				for (int j2 = 0; j2 < filter_span; ++j2) { // j2 is index within filter
					if (this->forms1Glob[i1][j2] == 1) { // If filter bit is '1'
						if (global_pattern_source.Vars_Glob[i][j1 + j2] == false) { // Pattern must have '1'
							current_alignment_matches = false;
							break;
						}
					}
					// If filter bit is '0', it's a wildcard - matches.
				}

				if (current_alignment_matches) {
					mas3[i1] = 1; // Filter i1 covers pattern i (at this alignment)
					covering_filters_count_for_pattern_i++;
					break; // Filter i1 covers pattern i, no need to check other alignments for this filter
				}
			}
		}

		// If exactly one filter (from the original set) covers this pattern, mark it essential
		int k_sum_mas3 = 0;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			k_sum_mas3 += mas3[i1];
		}

		if (k_sum_mas3 == 1) {
			for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
				if (mas3[i1] == 1) {
					mas4[i1] = 1; // Mark as essential
				}
			}
		}
	}
	std::cerr << "\n  Forms1_anal1 Pass 1 (M=" << current_M << ") complete." << std::endl;

	// Save intermediate results (like forms_an1_...txt)
	temp_fname_ss << "dbg_forms_an1_K" << pattern_K << "_M" << current_M << ".txt";
	std::ofstream fout_an1(temp_fname_ss.str());
	if (fout_an1.is_open()) {
		fout_an1 << "Nform1 = " << this->Nform1_Glob << "\tform1Size = " << this->form1Size
			<< "\tNsovp1 = " << this->Nsovp1_Glob << "\tN_2(K_pattern) = " << pattern_K;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			fout_an1 << "\n>" << i1 << "\t";
			for (int j1 = 0; j1 < this->form1Size; ++j1) fout_an1 << this->forms1Glob[i1][j1] << " ";
			fout_an1 << "\t" << mas4[i1]; // Essentiality marker
		}
		fout_an1.close();
		std::cerr << "  Debug info (an1) for M=" << current_M << " written to " << temp_fname_ss.str() << std::endl;
	}
	temp_fname_ss.str("");


	// --- Second Pass - Iterative Refinement ---
	std::cerr << "  Forms1_anal1 Pass 2 (M=" << current_M << ") starting iterative refinement..." << std::endl;
	for (int ii1 = 0; ii1 < this->Nform1_Glob; ++ii1) { // Iterate over each filter to potentially remove
		if (ii1 > 0 && ii1 % 100 == 0 && this->Nform1_Glob > 100) {
			std::cerr << "\r    Pass 2: Refinement iteration " << ii1 << "/" << this->Nform1_Glob << std::flush;
		}

		if (mas4[ii1] == 1) continue; // Skip if already marked essential from Pass 1

		mas4[ii1] = -1; // Tentatively mark for removal

		bool made_change_promoting_to_essential_this_iter = false;
		for (int i = 0; i < global_pattern_source.NVars_Glob; ++i) { // Re-evaluate all base patterns
			for (int r = 0; r < this->Nform1_Glob; ++r) mas3[r] = 0; // Reset mas3 for this pattern

			int active_covering_filters_count = 0;
			for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) { // Iterate over filters
				if (mas4[i1] == -1) continue; // Skip if tentatively removed (this includes current ii1 if not restored)

				int filter_span = this->forms1Glob[i1][filter_storage_N];
				if (filter_span <= 0) continue;

				for (int j1 = 0; j1 <= pattern_N - filter_span; ++j1) { // Sliding window
					if (global_pattern_source.Vars_Glob[i][j1] == false && this->forms1Glob[i1][0] == 1) {
						continue;
					}
					bool current_alignment_matches = true;
					for (int j2 = 0; j2 < filter_span; ++j2) {
						if (this->forms1Glob[i1][j2] == 1) {
							if (global_pattern_source.Vars_Glob[i][j1 + j2] == false) {
								current_alignment_matches = false;
								break;
							}
						}
					}
					if (current_alignment_matches) {
						mas3[i1] = 1; // Active filter i1 covers pattern i
						active_covering_filters_count++;
						break;
					}
				}
			}

			// If after tentatively removing ii1 (and any others marked -1),
			// a pattern is now covered by exactly one *active* filter, that filter becomes essential.
			int k_sum_mas3_active = 0;
			int newly_essential_idx = -1;
			for (int r = 0; r < this->Nform1_Glob; ++r) {
				if (mas4[r] != -1 && mas3[r] == 1) { // if active and covers
					k_sum_mas3_active++;
					newly_essential_idx = r; // If k_sum_mas3_active ends up 1, this is the one.
				}
			}

			if (k_sum_mas3_active == 1) {
				if (mas4[newly_essential_idx] == 0) { // If it wasn't already essential (mas4=1) or removed (mas4=-1)
					mas4[newly_essential_idx] = 1; // Promote it
					made_change_promoting_to_essential_this_iter = true;
				}
				if (newly_essential_idx == ii1) { // If the filter we tentatively removed became essential
					// This means removing it would make some pattern uncovered by any other single filter.
					// So, we must keep ii1.
					// However, the logic of professor's code is more about promoting others.
					// If ii1 itself becomes essential this way, it gets restored from -1 to 1.
				}
			}
		}
		// If filter ii1 (marked -1) was found to be "needed" because its removal caused
		// another filter to become uniquely essential for some pattern, then ii1 might still be removed.
		// The professor's code doesn't seem to explicitly "restore" ii1 if its removal doesn't lead to uncoverage.
		// It relies on the fact that if removing ii1 leaves a pattern uncovered by *any* filter,
		// then the `k_sum_mas3_active` for that pattern would be 0, and nothing gets promoted.
		// The crucial part is that `mas4[ii1]` remains -1 unless it itself gets promoted to 1.
		// Let's refine this: if removing `ii1` makes *any* pattern completely uncovered by *active* filters,
		// `ii1` must be restored.

		bool ii1_must_be_restored = false;
		if (!made_change_promoting_to_essential_this_iter) { // Only check if no other filter was promoted due to ii1's removal
			for (int i = 0; i < global_pattern_source.NVars_Glob; ++i) {
				// Recalculate coverage for pattern i, assuming ii1 is removed (mas4[ii1] == -1)
				bool pattern_i_is_covered = false;
				for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
					if (mas4[i1] == -1) continue; // Skip removed filters

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
					ii1_must_be_restored = true; // Removing ii1 uncovers this pattern
					break;
				}
			}
		}
		if (ii1_must_be_restored) {
			mas4[ii1] = 0; // Restore it (or to 1 if it was uniquely covering)
			// For simplicity, restore to 0, it might become 1 again if it's sole cover.
			// Or, if it leads to another becoming 1, it means ii1 was a key part of that.
			// The original logic is subtle. Let's stick to: if removing it uncovers something, restore it.
			// More accurately, if mas4[ii1] was 0 and removing it uncovers something, it becomes 1.
			// This state restoration is the tricky part.
			// The provided code doesn't explicitly restore `mas4[ii1]` from -1 in this loop.
			// It seems to rely on the idea that essential filters are found.
			// If `mas4[ii1]` is -1, and some pattern is ONLY covered by it, then `k_sum_mas3_active` would be 0
			// for that pattern, and `ii1` doesn't get promoted. This implies `ii1` stays removed.
			// This implies that if a filter is not marked `mas4[i]=1` after pass 1 and pass 2, it's gone.
		}

	} // End of Pass 2 loop (ii1)
	std::cerr << "\n  Forms1_anal1 Pass 2 (M=" << current_M << ") complete." << std::endl;

	// Save intermediate results (like forms_an2_...txt)
	temp_fname_ss << "dbg_forms_an2_K" << pattern_K << "_M" << current_M << ".txt";
	std::ofstream fout_an2(temp_fname_ss.str());
	if (fout_an2.is_open()) {
		fout_an2 << "Nform1 = " << this->Nform1_Glob << "\tform1Size = " << this->form1Size
			<< "\tNsovp1 = " << this->Nsovp1_Glob << "\tN_2(K_pattern) = " << pattern_K;
		for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
			fout_an2 << "\n>" << i1 << "\t";
			for (int j1 = 0; j1 < this->form1Size; ++j1) fout_an2 << this->forms1Glob[i1][j1] << " ";
			fout_an2 << "\t" << mas4[i1]; // Essentiality marker after Pass 2
		}
		fout_an2.close();
		std::cerr << "  Debug info (an2) for M=" << current_M << " written to " << temp_fname_ss.str() << std::endl;
	}
	temp_fname_ss.str("");

	// --- Final Selection and Compaction ---
	std::vector<int*> refined_forms_list;
	std::vector<int> original_indices_kept; // For debugging if needed

	for (int i1 = 0; i1 < this->Nform1_Glob; ++i1) {
		if (mas4[i1] == 1) { // Only keep filters marked as essential
			refined_forms_list.push_back(this->forms1Glob[i1]);
			original_indices_kept.push_back(i1);
		}
		else {
			delete[] this->forms1Glob[i1]; // Delete data of filters not kept
			this->forms1Glob[i1] = nullptr; // Avoid dangling pointer in original array
		}
	}

	delete[] this->forms1Glob; // Delete the old array of pointers

	this->Nform1_Glob = refined_forms_list.size();
	if (this->Nform1_Glob > 0) {
		this->forms1Glob = new int* [this->Nform1_Glob];
		for (size_t i = 0; i < this->Nform1_Glob; ++i) {
			this->forms1Glob[i] = refined_forms_list[i]; // Transfer pointer
		}
	}
	else {
		this->forms1Glob = nullptr;
	}

	std::cerr << "  Forms1_anal1: Refinement resulted in " << this->Nform1_Glob << " essential filters for M=" << current_M << "." << std::endl;

	// Save the final refined forms
	// The professor's code saves two files, forms_Result_... and forms_Result_..._ch. They seem identical.
	this->save_forms_to_file(output_refined_forms_filename, "Globally Refined MCS_0 (Forms1_anal1)");

	// Construct the "_ch" filename variant
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
	// delete[] nums; // nums array was for compaction, handled differently here with vector

	std::cerr << "Forms1_anal1_Global_Refinement for M=" << current_M << " finished.\n";
}


#endif  // _INC_forms2