# rocketfasta 🚀

## Overview

 The **FASTA Parser** is a command-line tool designed to process FASTA files, retrieve sequence headers, and extract sequences based on a specified reference ID. 

## Command-line Arguments

### 1. **Input File** (`-i`, `--fasta-input`)
- **Description**: Specifies the path to the input FASTA file.
- **Usage**: This argument is **required** for the program to run.
- **Syntax**: 
  ```bash
  -i <path_to_file> 
  --fasta-input <path_to_file>
  ```
- **Example**:
  ```bash
  ./fasta_parser -i input.fasta
  ```

---

### 2. **Print Headers** (`-H`, `--print-headers`)
- **Description**: A flag to print the headers of the FASTA sequences. When this option is specified, the program will display the headers from the input FASTA file.
- **Usage**: This argument is optional and acts as a flag (no additional value required).
- **Syntax**:
  ```bash
  -H 
  --print-headers
  ```
- **Example**:
  ```bash
  ./fasta_parser -i input.fasta -H
  ```

---

### 3. **Reference ID** (`-r`, `--reference-id`)
- **Description**: Specifies the ID of the reference sequence to be located in the FASTA file. The program will search for and operate on the sequence with the given reference ID.
- **Usage**: This argument is optional but requires a string value corresponding to the reference ID.
- **Syntax**:
  ```bash
  -r <reference_id>
  --reference-id <reference_id>
  ```
- **Example**:
  ```bash
  ./fasta_parser -i input.fasta -r ref123
  ```

---

### 4. **Subseqeunce** (`-s`, `--subsequence-pattern`)
- **Description**: "A nucleotide sequence that you want to file in the fasta sequences of fasta file. Return BED style Start and End coordinates.
- **Usage**: This argument is optional but requires a string value corresponding to the wanted subsequence to be searched.
- **Syntax**:
  ```bash
  -s <nucleotide_sequence> 
  --subsequence-pattern <nucleotide_sequence>
  ```
- **Example**:
  ```bash
  ./fasta_parser -i input.fasta -s ATGT
  ```

---

### 5. **Alignemnt type** (`-a`, `--alignment-type`)
- **Description**: Allows to specify alignment type between global (g) or local (l).
- **Usage**: If this argument is not provided, alignment defaults to global.
- **Syntax**:
  ```bash
  -a <g> or <l>
  --alignment-type <g> or <l>
  ```
- **Example**:
  ```bash
  ./fasta_parser -i input.fasta -r ID 7 -a l // for local alignment
  ./fasta_parser -i input.fasta -r ID 7 -a g // for global alignment
  ```

---


## Exit Codes

- **0**: Program executed successfully.
- **1**: Missing required arguments or invalid input.
- **2**: Invalid argument for alignment type / Alignment reference doesnt exist 
---

## Notes
- The input file must be in valid FASTA format.
- The program ensures UTF-8 compliance for all input arguments.

---

## Credits

> [CLI11](https://cliutils.github.io/CLI11/) library.
> [zlib-1.3.1](https://github.com/madler/zlib) library