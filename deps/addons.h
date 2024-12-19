#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <optional>

using namespace std;


// Finds partially matching subsequences in a sequence. TODO: returns starting position, ensure it also returns the ending position 
// Implementation of Boyer-Moore Algorithm. Finds only the first subseq in sequence.
vector<pair<int,int>>  FirstMatch(string& sequence, string& subseq){

    int n = sequence.length();
    int m = subseq.length();

    map<char,int> last_occurrence;

    for (int i =0; i < m; i++){
        
        last_occurrence[subseq[i]] = i;
    }

    int i = m - 1;
    int j = m - 1;


    while (i < n){

        if (subseq[j]==sequence[i]){

            if (j==0){
                vector<pair<int,int>> seq_range;  
                seq_range.push_back({i, i + m+1}); // BED format. Start inclusive, End non-inclusive.
                return seq_range;
            }

            i -= 1;
            j -= 1;
        } else {
            
            // Check if the current character in `sequence[i]` exists in the `last_occurrence` map.
            // If it exists, retrieve its value (the last occurrence index) from the map.
            // Otherwise, use -1 as the default value.
            // This ensures we have a valid value to work with, even if the character is not in the map.

            int lastIndex = (last_occurrence.find(sequence[i]) != last_occurrence.end())? last_occurrence[sequence[i]] : -1;

            // Compute the minimum value
            i+= m - min(j, 1 + lastIndex);

            j = m - 1;
        }
    }
    
     vector<pair<int,int>> no_found;  
                no_found.push_back({-1, -1});

    return no_found;  //Pattern not found

}

