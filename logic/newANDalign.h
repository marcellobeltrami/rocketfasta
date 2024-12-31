#include <iostream>
#include <string>
#include <fstream>

#include "../deps/structs.h"

using namespace std;


void global_alignment( map<string, pair<string, string>> FastaFileStructure, string reference_seq){
    for (const auto& fasta_entry : FastaFileStructure){ 
                
                if (fasta_entry.first != reference_seq){ // Prevents alignent with reference.
                    
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

}

void local_alignment(map<string, pair<string, string>> FastaFileStructure, string reference_seq){

    for (const auto& fasta_entry : FastaFileStructure){ 
                
                if (fasta_entry.first != reference_seq){ // Prevents alignent with reference.
                    
                    string seq1=FastaFileStructure[reference_seq].second;
                    string seq2 = FastaFileStructure[fasta_entry.first].second;
                    
                    Align rad_seqs;

                    vector<string> my_alignment = rad_seqs.Local_NW(seq1,seq2);

                    cout << "#################################################" << endl;
                    cout << "Ref: " << reference_seq << " Target: " << fasta_entry.first<< endl; 
                    cout << my_alignment[0] << endl;
                    cout << my_alignment[1] << endl;
                    }
            }

}


// Carrys out alignment if argument is provided.
void alignment(string path_to_fasta,string reference_seq, map<string, pair<string, string>> FastaFileStructure, BLAST_info blast_info, string al_type ){
 
 if (reference_seq.empty() ==false){
        FastaFileStructure = blast_info.New(path_to_fasta); // Struct: ID (metadata,sequence)

        auto reference_seq_obj = FastaFileStructure.find(reference_seq);
        if (reference_seq_obj != FastaFileStructure.end()) { 
            
            cout << "reference_seq"<<endl;
            cout << FastaFileStructure[reference_seq].second<<endl;
            cout << "Aligned"<<endl;
            
            if (al_type == "g"){
                global_alignment(FastaFileStructure, reference_seq);
            } else if (al_type == "l"){
                local_alignment(FastaFileStructure, reference_seq);
            } else {
                cout << "Invalid alignment type. Use 'g' for global and 'l' for local alignment." << endl;
                exit (2);
            }  
        
        } else {
                cout << "reference_seq: '" << reference_seq << "' does not exist in the file!" << endl;
                exit (2);
            }
    } 
}






