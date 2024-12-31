#include <iostream>
#include <string>
#include <fstream>

#include "./deps/CLI11.hpp"
#include "./logic/newANDalign.h"
#include "./logic/seqSearch.h"

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

    string al_type; // Function to be implemented. (logic is there, need to connect to CLI)
    app.add_option("-a,--alignment-type", al_type, "Specifies a global alignment of all sequences in the fasta file against reference. g = global, l = local");
    
    string subsequence;
    app.add_option("-s,--subsequence-pattern", subsequence, "A nucleotide sequence that you want to file in the fasta sequences of fasta file. Start inclusive, End non-inclusive.");


    CLI11_PARSE(app, argc, argv);
    // CLI Interface ends----------------------------//

    if (al_type == ""){ // Default alignment is global.
        al_type = "g"; 
    }


    BLAST_info my_fasta;
    map<string, pair<string, string>> FastaStruct; // This variable is commonly shared. Structured to only create one representation of fasta im memory.
    
    
    // Prints headers if flag is passed.
    if( header == true){
        my_fasta.New(path);
        my_fasta.FastaHeaders();
        
    }
    
    // Carrys out alignment when reference is provided.
    alignment(path,reference,FastaStruct, my_fasta, al_type); // add alignement type to function.
    
    // Finds first occurrences of subsequence in all sequences of fasta file.
    subseqsearch(subsequence, path, reference, FastaStruct, my_fasta);

    exit(0); // Exit program when finished with no errors 

}

