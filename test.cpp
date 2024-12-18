#include <iostream>
#include <string>
#include <fstream>

#include "./deps/addons.h"


int main(){
    
    string sequence = "ATCGAGCATGCATGAAAAAAAAAAACTAGCTAGCTAGC";
    string subseq = "AAAAAAAAAA";

    cout << "Location" << PartialMatch(sequence,subseq)<< endl;


    return 0;

}