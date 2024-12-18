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

## Usage Examples

### Example 1: Parse a FASTA file and print all headers
```bash
./fasta_parser -i input.fasta -H
```

### Example 2: Parse a FASTA file and find a sequence with a specific reference ID
```bash
./fasta_parser -i input.fasta -r ref123
```

### Example 3: Parse a FASTA file without additional options
```bash
./fasta_parser -i input.fasta
```

---

## Exit Codes

- **0**: Program executed successfully.
- **1**: Missing required arguments or invalid input.

---

## Notes
- The input file must be in valid FASTA format.
- The program ensures UTF-8 compliance for all input arguments.

---

## Credits

Developed using the [CLI11](https://cliutils.github.io/CLI11/) library.
