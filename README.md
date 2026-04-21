# **MutSeq**

MutSeq is a Python CLI tool for smart DNA sequence mutagenesis. It automates the process of mutating DNA sequences based on protein-level requirements while preventing common synthesis issues like DNA slippage.

## **Features**

* **Smart Codon Selection**: Uses synonymous codons to avoid repetitive sequences (anti-slippage) and unwanted restriction sites.  
* **Flexible Referencing**: Mutate directly on a target sequence or align a target DNA sequence to a reference (Local FASTA or UniProt ID).  
* **Validation**: Automatically verifies residue identity and handles gaps in alignment.  
* **Batch Processing**: Process multiple mutations at once and receive a detailed CSV audit log.

## **Installation**

Bash

git clone https://github.com/ignaciojci/mutseq.git  
cd mutseq  
pip install .

## **Quick Start**

### **1\. Direct Mutation**

Mutate the DNA by providing coordinates relative to the translation of the input file.

Bash

mutseq target.fasta \--muts M1A R24G \--outpre my\_mutation

### **2\. Reference Alignment (UniProt)**

Align your DNA to a known protein sequence and mutate based on the protein's numbering.

Bash

mutseq target.fasta \--uniprot P04637 \--muts R24A \--avoid GAATTC

## **CLI Arguments**

| Argument | Description |
| :---- | :---- |
| target | Path to the target DNA FASTA file. |
| \--muts | List of mutations (e.g., R24A Y102F). **Required.** |
| \--uniprot | UniProt Accession ID for reference alignment. |
| \--ref | Path to a local reference protein FASTA. |
| \--outpre | Prefix for output files (default: mutated\_output). |
| \--avoid | List of DNA sequences/restriction sites to avoid. |
| \--debug | Enable detailed debug logging and tracebacks. |

## **Outputs**

1. **\[outpre\].fasta**: The mutated DNA sequence.  
2. **\[outpre\]\_log.csv**: A manifest detailing the success or reason for skipping every requested mutation.