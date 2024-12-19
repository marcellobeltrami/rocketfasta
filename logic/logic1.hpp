#include <iostream>
#include <string>
#include <fstream>

#include "../deps/CLI11.hpp"
#include "../deps/structs.h"
#include "../deps/addons.h"


void alignment(string path, string reference, map<string, pair<string, string>> FastaStruct, BLAST_info my_fasta, Align rad_seqs){
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

                    rad_seqs;
                    
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