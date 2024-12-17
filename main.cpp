#include <iostream>
#include <string>
#include <fstream>

#include "./deps/CLI11.hpp"
#include "./deps/structs.h"



int main(int argc, char** argv) {
    // CLI Interface starts----------------------------//
    
    CLI::App app{"Another FASTA parser"};
    argv = app.ensure_utf8(argv);

    string path;
    app.add_option("-i,--fasta-input", path, "Path to fasta file")
                    ->required();

    bool header;
    app.add_flag("-H,--print-headers", header, "Print Headers of a FASTA sequence.");

    string reference;
    app.add_option("-r,--reference-id", reference, "ID of the reference sequence found in FASTA file.");

    CLI11_PARSE(app, argc, argv);
    // CLI Interface ends----------------------------//


    BLAST_info my_fasta;
    
    
    // Prints headers if flag is passed.
    if( header == true){
        my_fasta.New(path);
        my_fasta.FastaHeaders();
        
    }
    
    // Carrysout alignment with all sequences in the file give one of the sequences ID as a reference.
    if (reference.empty() ==false){
        map<string, pair<string, string>> FastaStruct = my_fasta.New(path); // Struct: ID (metadata,sequence)

        auto it = FastaStruct.find(reference);
        if (it != FastaStruct.end()) {
            
            cout << "Reference"<<endl;
            cout << FastaStruct[reference].second<<endl;
            cout << "Aligned"<<endl;
            for (const auto& fasta_entry : FastaStruct){ 
                
                if (fasta_entry.first != reference){
                    string ID1=reference;
                    string ID2=fasta_entry.first;

                    string seq1=FastaStruct[ID1].second;
                    string seq2 = FastaStruct[fasta_entry.first].second;

                    Align rad_seqs;
                    
                    vector<string> my_alignment = rad_seqs.Global_NW(seq1,seq2);

                    cout << "#################################################" << endl;
                    cout << "Ref: " << ID1 << " Target: " << ID2<< endl; 
                    cout << my_alignment[0] << endl;
                    cout << my_alignment[1] << endl;
                }


            }


        } else {
            cout << "Reference: '" << reference << "' does not exist in the file!" << endl;
            exit (2);
    }

    }


    return 0;

}

