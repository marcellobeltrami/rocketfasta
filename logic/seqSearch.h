#include <iostream>
#include <string>
#include <fstream>

#include "../deps/addons.h"

using namespace std;

void subseqsearch(string subsequence, string path, string reference, map<string, pair<string, string>> FastaStruct, BLAST_info my_fasta){

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
}