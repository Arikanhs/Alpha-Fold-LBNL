#!/usr/bin/env python
# coding: utf-8

# In[3]:


# This file includes the codes in pdb_manipulation.ipynb file as functions
# Functions take ppdb_codes such as 2C0K as input 

import sys
from biopandas.pdb import PandasPdb
import numpy as np
import os
import requests

from Bio import SeqIO
from collections import OrderedDict
import re

# Output destination folder location
destination = "../Outputs/data_manipulation/"

# Dataset destination folder location
destination2 = "../Dataset/test_pdbs/"


# In[4]:


# get_sequence function extracts the sequence from the ATOM section in pdb data file
def get_sequence(ppdb):
    amino_acid_sequence = ppdb.amino3to1().residue_name
    amino_acid_sequence = ''.join(amino_acid_sequence)
    return amino_acid_sequence


# In[5]:


# get_seqres_sequence function extracts amino acids from the SEQRES section, converts them from 3 letter code to 1 letter code and returns the sequence
def get_seqres_sequence(ppdb):
    df_others = ppdb.df["OTHERS"]
    seqres_entry = df_others[df_others['record_name'] == 'SEQRES']['entry'].values

    # Initialize an empty dictionary
    result_dict = {}

    # Process each string in the array
    for item in seqres_entry:
        # Split the string into words
        words = item.split()
        
        # Extract the second letter to get the chain identifier "A"
        second_letter = words[1]
        
        # Extract the substring after the '46', which starts from the 4th word
        substring = ' '.join(words[3:])
        
        # Append the substring to the corresponding key in the dictionary
        if second_letter in result_dict:
            result_dict[second_letter].append(substring)
        else:
            result_dict[second_letter] = [substring]

    # print(result_dict)

    # Initialize a new dictionary to store the merged strings
    merged_dict = {}

    # Iterate through each key in the dictionary
    for key, strings in result_dict.items():
        # Merge all strings in the list into one string
        merged_string = ' '.join(strings)
        # Add the merged string to the new dictionary
        merged_dict[key] = merged_string

    # print(merged_dict)

    # Mapping dictionary from three-letter to one-letter amino acid codes
    three_to_one = {
        # 21 Amino-acids
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U",
        
        "ACE": "X", "DA": "", "DT": "", "DG": "", "DC": "", "DU": ""
    }

    # Function to convert three-letter codes to one-letter codes
    def convert_to_one_letter(amino_acid_sequence):
        one_letter_sequence = ''.join(three_to_one[aa] for aa in amino_acid_sequence.split())
        return one_letter_sequence

    # Convert the sequences in the dictionary
    converted_amino_acids = {key: convert_to_one_letter(seq) for key, seq in merged_dict.items()}

    # Print the converted sequences
    for key, seq in converted_amino_acids.items():
        print(f"{key}: {seq}")

    # Store the converted sequences in strings
    # sequence_A = converted_amino_acids['A']
    # sequence_B = converted_amino_acids['B']

    merged_sequence = "" 
    seen_chains = set()

    for chain in converted_amino_acids:
        if converted_amino_acids[chain] not in seen_chains:
            merged_sequence += converted_amino_acids[chain]
            # Add the seen chains to set
            seen_chains.add(converted_amino_acids[chain])

    return merged_sequence



# In[6]:

# Extracts each chain and writes their PDB record to individual files. Creates a folder for each protein
def extract_chains_to_pdb(ppdb, ppdb_code):

    sequence = ppdb.amino3to1()

    for chain_id in sequence['chain_id'].unique():
        print('Chain ID: %s' % chain_id)
        print(''.join(sequence.loc[sequence['chain_id'] == chain_id, 'residue_name']))

    print("")    

    output_folder = f"{destination}{ppdb_code}"
    os.makedirs(output_folder, exist_ok=True)


    # Get the unique chain IDs
    chain_ids = ppdb.df['ATOM']['chain_id'].unique()

    # Iterate over each chain ID
    for chain_id in chain_ids:
        # Filter the dataframe for the current chain
        chain_df = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == chain_id]

        # print(chain_df)

        # Create a new PandasPdb object for the chain
        chain_ppdb = PandasPdb()
        chain_ppdb.df['ATOM'] = chain_df

        # Write the chain to a new PDB file
        # output_file = f"../Outputs/data_manipulation/{ppdb_code}_chain_{chain_id}.pdb"
        output_file = os.path.join(output_folder, f"{ppdb_code}_chain_{chain_id}.pdb")
        chain_ppdb.to_pdb(path=output_file, records=None, gz=False, append_newline=True)

# Downloads the Fatsa file and stores it in the destination2 folder 
def download_file(ppdb_code):

    # Example usage
    # https://www.rcsb.org/fasta/entry/1K8K
    url = f"https://www.rcsb.org/fasta/entry/{ppdb_code}"  # URL to the PDB file

    output_file = os.path.join(destination2, f"{ppdb_code}.fasta")  # Path to save the downloaded file
    response = requests.get(url)

    if response.status_code == 200:
        with open(output_file, 'wb') as f:
            f.write(response.content)
        print("File downloaded successfully.")
    else:
        print(f"Failed to download file. Status code: {response.status_code}")

    return output_file



def read_fatsa(fasta_file_path):

    # Read the FASTA file and get the final sequence

    # -------------------------   EDGE CASES  ------------------------
    # Two or more chains with same sequence, Chains A,B = AANAAA:    Code chooses the chain that comes first in alphabetical order (A)
    # Chain names with one or more letter, Chains B, BA, BB, BC:     Code chooses the chain that comes first in alphabetical order (B) 
    # Chains with multiple names, Chains A[auth P], B[auth Q]:       Code removes anything between the brackets 
    # DNA Chains:                                                    Code identifies the "DNA" word and skips the sequence
    # RNA Chains:                                                    Code identifies the "RNA" word and skips the sequence
    # Long chain sequences with multiple lines:                      Multiple lines of sequences are automatically extracted as one piece by the record.seq function


    sequences_dict = OrderedDict()

    for record in SeqIO.parse(fasta_file_path, "fasta"):
        desc = record.description # Gets each line starts with > symbol aka description part
        words = desc.split("|")   # Split the string into words using the "|" delimiter

        # Check if "DNA" or "RNA" exists in the third index of the line
        if len(words) > 2 and ("DNA" in words[2] or "RNA" in words[2]):
            continue  # Skip this line and move to the next one
        
        # Process the line if "DNA" is not found
        
        # print(words)
        chains = words[1]                       # Get the chain info
        chains = chains.replace(",", "")        # Remove the comma between chain IDs

        # print(chains)
        cleaned_chains = re.sub(r'\[.*?\]', '', chains)         # Use regular expression to remove everything inside square brackets and the brackets
        cleaned_chains = cleaned_chains.strip()                 # Strip any extra spaces
        splitted_chains = cleaned_chains.split(" ")             # Split everything by white spaces

        # print(splitted_chains)
        
        splitted_chains = splitted_chains[1:]

        sorted_chains = sorted(splitted_chains)

        sequence = str(record.seq)

        for chain_id in splitted_chains:
            sequences_dict[chain_id] = sequence

    for chain_id in sequences_dict:
        print('Chain ID: %s' % chain_id)
        print(sequences_dict[chain_id])

    # all_sequence = ""
    # seen_chains = set()
    # short_sequence = ""

    # # sorted_sequences_dict = sequences_dict(sorted(sequences_dict.items()))
    # sorted_sequences_dict = OrderedDict(sorted(sequences_dict.items()))

    # for chain_id in sorted_sequences_dict:
    #     if sorted_sequences_dict[chain_id] not in seen_chains:
    #         chain_seq = sorted_sequences_dict[chain_id]

    #         print(f"Chain: {chain_id}")
    #         print(f"Sequence: {chain_seq}\n")

    #         all_sequence += chain_seq
    #         short_sequence += chain_seq[0] + chain_seq[1] + chain_seq[2] + "-"

    #         # Add the seen chains to set
    #         seen_chains.add(chain_seq)

    # print(f"Final sequence: {all_sequence}")
    # print(f"Shortened final sequence (first 3 aminoacids of each chain): {short_sequence}")        


# In[ ]:


if len(sys.argv) != 2:
    print("Usage: python script.py <PDB_code>")
    sys.exit(1)

# Get the PDB code from command-line argument
ppdb_code = sys.argv[1]

# Call the functions with the provided PDB code
print(f"PDB code: {ppdb_code}\n")

ppdb = PandasPdb().fetch_pdb(ppdb_code)
fatsa_loc = download_file(ppdb_code)

sequence = get_sequence(ppdb)
print(f"Sequence from ATOM section: \n{sequence}\n")

print("Sequence from SEQRES section: ")
seqres_sequence = get_seqres_sequence(ppdb)
print(f"Sequence: \n{seqres_sequence}\n")

print(f"Get each chain from ATOM section and write them to individual PDBs:")
extract_chains_to_pdb(ppdb, ppdb_code)

print(f"Get each chain from FATSA file:")
read_fatsa(fatsa_loc)





