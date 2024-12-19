#include <iostream>
#include <string>
#include <fstream>

#include "./deps/CLI11.hpp"
#include "./deps/structs.h"
#include "./deps/addons.h"


int main(int argc, char** argv) {
    // CLI Interface starts----------------------------//
    
    CLI::App app{"Another FASTA parser"};
    argv = app.ensure_utf8(argv);

    string path;
    app.add_option("-i,--fasta-input", path, "Path to fasta file. Format: .fna, .fasta and .fa or .gz versions")
                    ->required();

    bool header;
    app.add_flag("-H,--print-headers", header, "Print Headers of a FASTA sequence.");

    string reference;
    app.add_option("-r,--reference-id", reference, "ID of the reference sequence found in FASTA file.");


    string subsequence;
    app.add_option("-s,--subsequence-pattern", subsequence, "A nucleotide sequence that you want to file in the fasta sequences of fasta file. Start inclusive, End non-inclusive.");


    CLI11_PARSE(app, argc, argv);
    // CLI Interface ends----------------------------//


    BLAST_info my_fasta;
    map<string, pair<string, string>> FastaStruct; // This variable is commonly shared. Structured to only create one representation of fasta im memory.
    
    
    // Prints headers if flag is passed.
    if( header == true){
        my_fasta.New(path);
        my_fasta.FastaHeaders();
        
    }
    
    // Carrys out alignment with all sequences in the file if one of the sequences ID id given as a reference.
    if (reference.empty() ==false){
        FastaStruct = my_fasta.New(path); // Struct: ID (metadata,sequence)

        auto reference_obj = FastaStruct.find(reference);
        if (reference_obj != FastaStruct.end()) {
            
            cout << "Reference"<<endl;
            cout << FastaStruct[reference].second<<endl;
            cout << "Aligned"<<endl;
            for (const auto& fasta_entry : FastaStruct){ 
                
                if (fasta_entry.first != reference){
                    
                    string seq1=FastaStruct[fasta_entry.first].second;
                    string seq2 = FastaStruct[fasta_entry.first].second;

                    Align rad_seqs;
                    
                    vector<string> my_alignment = rad_seqs.Global_NW(seq1,seq2);

                    cout << "#################################################" << endl;
                    cout << "Ref: " << reference << " Target: " << fasta_entry.first<< endl; 
                    cout << my_alignment[0] << endl;
                    cout << my_alignment[1] << endl;
                }


            }


        } else {
            cout << "Reference: '" << reference << "' does not exist in the file!" << endl;
            exit (2);
    }

    }

    // Finds first occurrences of subsequence in all sequences of fasta file.
    if (subsequence.empty() ==false){
        if (reference.empty() ==false){ // Uses already loaded file to find sequences.
            
            cout << "ID\tStart\tEnd"<<endl;

            for (const auto& fasta_entry : FastaStruct){
                
                vector<pair<int,int>> match = FirstMatch(FastaStruct[fasta_entry.first].second,subsequence);

                if (match[0].first == -1 && match[0].second){

                    cout << fasta_entry.first << "\t" << "NULL" << "\t" << "NULL" << endl;

                } else { 

                    cout << fasta_entry.first << "\t" << match[0].first << "\t" << match[0].second << endl;
                }

            } 
        } else if (reference.empty() ==true){ // Loads in file and Finds sequences.

            FastaStruct = my_fasta.New(path);
            cout << "ID\tStart\tEnd"<< endl;
            for (const auto& fasta_entry : FastaStruct){
                
                vector<pair<int,int>> match = FirstMatch(FastaStruct[fasta_entry.first].second,subsequence);

                if (match[0].first == -1 && match[0].second == -1){

                    cout << fasta_entry.first << "\t" << "NULL" << "\t" << "NULL" << endl;

                } else { 
                    cout << fasta_entry.first << "\t" << match[0].first << "\t" << match[0].second << endl;
                }

            }

        }

    }

    exit(0); // Exit program when finished with no errors 

}

