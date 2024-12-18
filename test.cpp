#include <iostream>
#include <string>
#include <fstream>

#include "./deps/addons.h"


int main(){
    
    string sequence = "ATCGAGCATGCATGAAAAAAAJAAAACTAGCTAGCTAGC";
    string subseq = "AAAAAAAAAA";
    
    
    vector<pair<int,int>> positions = FirstMatch(sequence,subseq);
    
    
    if (positions[0].first == -1 &&positions[0].second == -1 ){
        cout << "No sequences found" << endl;
        exit(404);
    }
    cout << "Start: " << positions[0].first
         <<" End: " << positions[0].second << endl;


    return 0;

}