#include <iostream>
#include <string>
#include <fstream>

#include "../deps/structs.h"

using namespace std;
// Carrys out alignment if argument is provided.
void alignment(string path_to_fasta,string reference_seq, map<string, pair<string, string>> FastaFileStructure, BLAST_info blast_info){
 
 if (reference_seq.empty() ==false){
        FastaFileStructure = blast_info.New(path_to_fasta); // Struct: ID (metadata,sequence)

        auto reference_seq_obj = FastaFileStructure.find(reference_seq);
        if (reference_seq_obj != FastaFileStructure.end()) {
            
            cout << "reference_seq"<<endl;
            cout << FastaFileStructure[reference_seq].second<<endl;
            cout << "Aligned"<<endl;
            for (const auto& fasta_entry : FastaFileStructure){ 
                
                if (fasta_entry.first != reference_seq){
                    
                    string seq1=FastaFileStructure[reference_seq].second;
                    string seq2 = FastaFileStructure[fasta_entry.first].second;
                    
                    Align rad_seqs;

                    vector<string> my_alignment = rad_seqs.Global_NW(seq1,seq2);

                    cout << "#################################################" << endl;
                    cout << "Ref: " << reference_seq << " Target: " << fasta_entry.first<< endl; 
                    cout << my_alignment[0] << endl;
                    cout << my_alignment[1] << endl;
                    }
            } 
        } else {
                cout << "reference_seq: '" << reference_seq << "' does not exist in the file!" << endl;
                exit (2);
            }
    } 
}






