import argparse
import sys
import os
import csv
import requests
from io import StringIO
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Data import CodonTable
import logging
import traceback

# Setup logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger("mutseq")

class DNASequenceMutator:
    def __init__(self, target_dna_path, avoid_sites=None):
        self.target_record = SeqIO.read(target_dna_path, "fasta")
        self.target_dna = str(self.target_record.seq).upper()
        self.target_aa = str(self.target_record.seq.translate(to_stop=True))
        self.avoid_sites = avoid_sites or []
        self.mutation_log = []

    def fetch_uniprot_ref(self, accession):
        """Fetches AA sequence from UniProt."""
        print(f"[*] Fetching Reference from UniProt: {accession}...")
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        response = requests.get(url)
        response.raise_for_status()
        ref_record = list(SeqIO.parse(StringIO(response.text), "fasta"))[0]
        return str(ref_record.seq)

    def align_and_map(self, ref_aa):
        """Aligns Ref AA to Target AA and maps indices."""
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        alignments = aligner.align(ref_aa, self.target_aa)
        if not alignments:
            raise ValueError("Alignment failed: No significant similarity found between Ref and Target.")
        
        best_alignment = alignments[0]
        # Check if indices were actually generated
        if not hasattr(best_alignment, 'aligned') or len(best_alignment.aligned) == 0:
            raise IndexError("Alignment mapping is empty. Check your input sequences.")
        
        # Mapping: Ref_Index -> Target_Index
        # best_alignment.indices provides the mapping for the alignment
        ref_indices, target_indices = best_alignment.aligned
        mapping = {}
        
        # Convert alignment segments into a flat dictionary
        for r_range, t_range in zip(ref_indices, target_indices):
            for i in range(r_range[1] - r_range[0]):
                mapping[r_range[0] + i] = t_range[0] + i
        return mapping

    def get_anti_slippage_codon(self, amino_acid, index, current_dna):
        """Selects a codon to avoid repeats and restriction sites."""
        table = CodonTable.unambiguous_dna_by_id[1]
        possible_codons = [c for c, aa in table.forward_table.items() if aa == amino_acid]
        
        # Get preceding 6 bases to check for repeats
        context_start = max(0, index * 3 - 6)
        prev_sequence = current_dna[context_start : index * 3]

        for codon in possible_codons:
            test_seq = prev_sequence + codon
            
            # 1. Anti-slippage: Avoid making a 3-unit homopolymer of codons (e.g., GCA GCA GCA)
            if len(prev_sequence) >= 6 and prev_sequence[-3:] == codon:
                continue
                
            # 2. Restriction site check
            site_found = any(site in test_seq for site in self.avoid_sites)
            if site_found:
                continue
            
            return codon
        
        return possible_codons[0] # Fallback to first available

    def process_mutations(self, ref_aa, mutations):
        mapping = self.align_and_map(ref_aa)
        mutated_dna_list = list(self.target_dna)
        
        for mut in mutations:
            # Parse: R24A -> ('R', 24, 'A')
            orig_res = mut[0]
            new_res = mut[-1]
            
            try:
                pos = int(mut[1:-1]) - 1 # 0-indexed
            except ValueError:
                continue
            
            # Before accessing the index, check length
            if pos < 0 or pos >= len(ref_aa):
                logger.warning(f"Mutation {mut} is invalid: Position {pos+1} is out of bounds for Ref (len: {len(ref_aa)})")
                continue

            # 1. Identity Check
            if pos >= len(ref_aa) or ref_aa[pos] != orig_res:
                actual = ref_aa[pos] if pos < len(ref_aa) else "OUT_OF_BOUNDS"
                self.log_mutation(mut, "SKIPPED", f"Identity Mismatch: Ref has {actual}")
                continue

            # 2. Gap Check
            target_res_idx = mapping.get(pos)
            if target_res_idx is None:
                self.log_mutation(mut, "SKIPPED", "Gap: Position aligns to gap in target")
                continue

            # 3. Apply Mutation
            new_codon = self.get_anti_slippage_codon(new_res, target_res_idx, "".join(mutated_dna_list))
            dna_pos = target_res_idx * 3
            mutated_dna_list[dna_pos:dna_pos+3] = list(new_codon)
            
            self.log_mutation(mut, "SUCCESS", f"Mutated to {new_codon}")

        return "".join(mutated_dna_list)

    def log_mutation(self, mut, status, note):
        self.mutation_log.append({"Mutation": mut, "Status": status, "Note": note})

    def export_results(self, final_dna, fasta_out, csv_out, annotation=None):
        # Save FASTA
        with open(fasta_out, "w") as f:
            if annotation:
                f.write(f">{annotation}\n{final_dna}\n")
            else:
                f.write(f">mutated_{self.target_record.id}\n{final_dna}\n")
        
        # Save CSV Log
        keys = self.mutation_log[0].keys()
        with open(csv_out, "w", newline="") as f:
            dict_writer = csv.DictWriter(f, fieldnames=keys)
            dict_writer.writeheader()
            dict_writer.writerows(self.mutation_log)
        print(f"[+] Process complete. DNA saved to {fasta_out}, Log saved to {csv_out}")

def run_cli():
    parser = argparse.ArgumentParser(description="MutSeq: Smart DNA Mutagenesis")
    parser.add_argument("target", help="Path to target DNA FASTA")
    parser.add_argument("--uniprot", help="UniProt Accession ID for reference")
    parser.add_argument("--ref", help="Path to local reference AA FASTA")
    parser.add_argument("--muts", nargs="+", required=True, help="Mutations (e.g. R24A)")
    parser.add_argument("--outpre", default="mutated_output", help="Output file prefix")
    parser.add_argument("--annotation", default="", help="Custom FASTA annotation")
    parser.add_argument("--avoid", nargs="+", help="Restriction sites to avoid")
    parser.add_argument("--debug", action='store_true', help="Enable debug logging")
    
    args = parser.parse_args()
    
    if args.debug:
        logger.setLevel(logging.DEBUG)

    try:
        # Step 1: Initialize Mutator
        logger.debug(f"Loading target DNA from {args.target}")
        mutator = DNASequenceMutator(args.target, avoid_sites=args.avoid)
        
        # Step 2: Reference Fetching Logic (with Direct Mode Fallback)
        if args.uniprot:
            ref_seq = mutator.fetch_uniprot_ref(args.uniprot)
            logger.debug(f"Fetched UniProt Ref: {len(ref_seq)} residues")
        elif args.ref:
            ref_seq = str(SeqIO.read(args.ref, "fasta").seq)
            logger.debug(f"Loaded Local Ref: {len(ref_seq)} residues")
        else:
            # FALLBACK: Use the target sequence itself as the reference
            logger.info("No reference provided. Using translated target DNA as reference (Direct Mode).")
            ref_seq = mutator.target_aa

        # Step 3: Mutation Processing
        logger.info(f"Starting mutation processing for {len(args.muts)} sites...")
        final_seq = mutator.process_mutations(ref_seq, args.muts)
        
        # Step 4: Export with sanitized filenames
        fasta_out = f"{args.outpre}.fasta"
        csv_out = f"{args.outpre}_log.csv"
        
        mutator.export_results(final_seq, fasta_out, csv_out, args.annotation)
        
    except FileNotFoundError:
        logger.error(f"File not found: Ensure '{args.target}' or your reference path is correct.")
    except IndexError:
        logger.error("Index Error: A mutation position is outside the sequence range.")
        if args.debug: logger.error(traceback.format_exc())
    except Exception as e:
        logger.critical(f"Unexpected Failure: {e}")
        if args.debug: logger.error(traceback.format_exc())

if __name__ == "__main__":
    run_cli()

# # --- CLI EXECUTION LOGIC ---
# if __name__ == "__main__":
#     # Configuration (In a real CLI, these come from argparse)
#     TARGET_FILE = "target.fasta"   # Path to your DNA file
#     REF_UNIPROT = "P04637"         # Example: p53
#     MUTATIONS = ["R24A", "Y102F", "E201H"]
#     AVOID = ["GAATTC"]             # EcoRI

#     if os.path.exists(TARGET_FILE):
#         mutator = DNASequenceMutator(TARGET_FILE, avoid_sites=AVOID)
#         ref_sequence = mutator.fetch_uniprot_ref(REF_UNIPROT)
        
#         final_seq = mutator.process_mutations(ref_sequence, MUTATIONS)
#         mutator.export_results(final_seq, "mutated.fasta", "mutation_log.csv")
#     else:
#         print("Please provide a valid target.fasta file.")

