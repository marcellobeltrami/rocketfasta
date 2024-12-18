#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include "./zlib-1.3.1/zlib.h"
#include <array>
#include <sstream> // For string trimming

using namespace std;

class BLAST_info  
{   
    private:
        
    string FastaFile_path;
    map<string, pair<string, string>> fastaObj;

    // Helper function to trim whitespace when defining ID
    string trim(const string& str) {
        size_t first = str.find_first_not_of(" \t\n\r");
        size_t last = str.find_last_not_of(" \t\n\r");
        return (first == string::npos || last == string::npos) ? "" : str.substr(first, (last - first + 1));
    }

    // Function to check if a file path ends with ".gz"
    bool endsWithGZ(const string& path) {
        return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
    }

   
    public: 

   map<string, pair<string, string>> New(string& path) {
    auto FastaFile = static_cast<void*>(nullptr); // defines a pointer to potentially hold the address of an object in the future
    string fastatxt;
    string ID;
    string fastaseq;
    string metadata;
    string temp_id;
    map<string, pair<string, string>> output;

    // Open the file based on its type (gzipped or plain)
    bool isGzipped = endsWithGZ(path);
    if (isGzipped) {
        FastaFile = gzopen(path.c_str(), "rb");
    } else {
        FastaFile = fopen(path.c_str(), "r");
    }

    // Check if the file was opened successfully
    if (!FastaFile) {
        cerr << "Error: File does not exist or cannot be opened: " << path << endl;
        exit(1);
    }

    cout << "Reading FASTA..." << endl;

    if (isGzipped) {
        // For gzipped files
        char buffer[8192];
        while (gzgets(static_cast<gzFile>(FastaFile), buffer, sizeof(buffer))) {
            fastatxt = string(buffer);
            fastatxt = trim(fastatxt); // Trim whitespace and newlines

            if (!fastatxt.empty() && fastatxt[0] == '>') {
                // Process header
                size_t spacePos = fastatxt.find(' ');
                temp_id = fastatxt.substr(1);  // Remove the '>' from the ID
                ID = trim(temp_id.substr(0, spacePos));
                metadata = spacePos != string::npos ? trim(temp_id.substr(spacePos )) : "";

                // If there is a previous sequence, store it in the map
                if (!ID.empty() && !fastaseq.empty()) {
                    output[ID] = make_pair(metadata, fastaseq);
                    fastaseq.clear();
                }
            } else if (!fastatxt.empty()) {
                // Append sequence line
                fastaseq += fastatxt;
            }
        }
        // Add the last record
        if (!ID.empty() && !fastaseq.empty()) {
            output[ID] = make_pair(metadata, fastaseq);
        }
        gzclose(static_cast<gzFile>(FastaFile));
    } else {
        // For plain text files
        char buffer[8192];
        while (fgets(buffer, sizeof(buffer), static_cast<FILE*>(FastaFile))) {
            fastatxt = string(buffer);
            fastatxt = trim(fastatxt); // Trim whitespace and newlines

            if (!fastatxt.empty() && fastatxt[0] == '>') {
                // Process header
                size_t spacePos = fastatxt.find(' ');
                temp_id = fastatxt.substr(1);  // Remove the '>' from the ID
                ID = trim(temp_id.substr(0, spacePos));
                metadata = spacePos != string::npos ? trim(temp_id.substr(spacePos)) : "";

                // If there is a previous sequence, store it in the map
                if (!ID.empty() && !fastaseq.empty()) {
                    output[ID] = make_pair(metadata, fastaseq);
                    fastaseq.clear();
                }
            } else if (!fastatxt.empty()) {
                // Append sequence line
                fastaseq += fastatxt;
            }
        }
        // Add the last record
        if (!ID.empty() && !fastaseq.empty()) {
            output[ID] = make_pair(metadata, fastaseq);
        }
        fclose(static_cast<FILE*>(FastaFile));
    }

    cout << "Fasta Object created." << endl;
    fastaObj = output;
    return output;

   }
    
    // Print Fasta Headers
    void FastaHeaders(){
        for (const auto& entry : fastaObj) {
            cout << entry.first <<  " " << entry.second.first << endl;
        }

    }
   
};



class Align{

    private: 

        // ################ Global_NW |Needleman-Wunsch| functions ################      
        int match = 1;   // Match score
        int mismatch = -1; // Mismatch penalty
        int gap = -2;     // Gap penalty

        // Function to initialize the score matrix
        void initializeMatrix_NW(vector<vector<int>>& matrix, const string& seq1, const string& seq2) {
            int m = seq1.size();
            int n = seq2.size();
            
            for (int i = 0; i <= m; ++i) {
                matrix[i][0] = i * gap;  // Gap penalty for column
            }
            
            for (int j = 0; j <= n; ++j) {
                matrix[0][j] = j * gap;  // Gap penalty for row
            }
        }

        // Function to fill the score matrix using dynamic programming
        void fillMatrix_NW(vector<vector<int>>& matrix, const string& seq1, const string& seq2) {
            int m = seq1.size();
            int n = seq2.size();
            
            for (int i = 1; i <= m; ++i) {
                for (int j = 1; j <= n; ++j) {
                    int score_diagonal = matrix[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch);
                    int score_up = matrix[i-1][j] + gap;
                    int score_left = matrix[i][j-1] + gap;
                    
                    matrix[i][j] = max({score_diagonal, score_up, score_left});
                }
            }
        }

        // Function to perform the traceback and get the alignment
        vector<string> traceback_NW(const vector<vector<int>>& matrix, const string& seq1, const string& seq2) {
            int i = seq1.size();
            int j = seq2.size();
            string alignedSeq1 = "";
            string alignedSeq2 = "";
            
            vector<string> output_strings;

            while (i > 0 && j > 0) {
                if (matrix[i][j] == matrix[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)) {
                    alignedSeq1 = seq1[i-1] + alignedSeq1;
                    alignedSeq2 = seq2[j-1] + alignedSeq2;
                    --i; --j;
                }
                else if (matrix[i][j] == matrix[i-1][j] + gap) {
                    alignedSeq1 = seq1[i-1] + alignedSeq1;
                    alignedSeq2 = "-" + alignedSeq2;
                    --i;
                }
                else {
                    alignedSeq1 = "-" + alignedSeq1;
                    alignedSeq2 = seq2[j-1] + alignedSeq2;
                    --j;
                }
            }

            // Handle remaining gaps if any sequence ends first
            while (i > 0) {
                alignedSeq1 = seq1[i-1] + alignedSeq1;
                alignedSeq2 = "-" + alignedSeq2;
                --i;
            }
            while (j > 0) {
                alignedSeq1 = "-" + alignedSeq1;
                alignedSeq2 = seq2[j-1] + alignedSeq2;
                --j;
            }

            output_strings.push_back(alignedSeq1);
            output_strings.push_back(alignedSeq2);

            return output_strings;
        }
        



        // ################ Local_SW |Smith-Waterman| functions ###################



    public: 

        // Returns a vector with the 2 aligned sequences.
        vector<string> Global_NW(string seq1,string seq2){
            
            vector<vector<int>> matrix(seq1.size() + 1, vector<int>(seq2.size() + 1, 0));
            
            initializeMatrix_NW(matrix,seq1,seq2);
            fillMatrix_NW(matrix,seq1,seq2);
            vector<string> alignment = traceback_NW(matrix,seq1,seq2);

            return alignment;
            
        }


        //vector<string> Local_SW(string seq1,string seq2) {} // IMPLEMENT


};


