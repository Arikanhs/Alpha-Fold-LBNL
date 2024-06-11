#!/usr/bin/env python
# coding: utf-8

import unittest
import difflib
from biopandas.pdb import PandasPdb
import numpy as np
import os

# get_sequence function extracts the sequence from the ATOM section in pdb data file
def get_sequence(ppdb_code):
    ppdb = PandasPdb().fetch_pdb(ppdb_code)
    amino_acid_sequence = ppdb.amino3to1().residue_name
    amino_acid_sequence = ''.join(amino_acid_sequence)
    return amino_acid_sequence

# get_seqres_sequence function extracts amino acids from the SEQRES section, converts them from 3 letter code to 1 letter code and returns the sequence
def get_seqres_sequence(ppdb_code):
    ppdb = PandasPdb().fetch_pdb(ppdb_code)
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


    merged_sequence = "" 
    seen_chains = set()

    for chain in converted_amino_acids:
        if converted_amino_acids[chain] not in seen_chains:
            merged_sequence += converted_amino_acids[chain]
            # Add the seen chains to set
            seen_chains.add(converted_amino_acids[chain])

    return merged_sequence

class TestGetSequence(unittest.TestCase):

    def test_get_sequence(self):

        #   -----------------------        EDGE CASES        -----------------------------
        #   
        #   Missing Amino acid:           Amino acid ACE replaced with X (or nothing depending on need)
        #   Same chains different names:  A,B = VLSPAD C,D = VHLTPE or A = KETAA B = KETAA (if the chain was seen before, it won't be included in the final sequence)
        #   Nucleotide bases (A-T-G-C):   They are not included in the SEQRES section
        
        protein_codes = ["1L2Y", "1UBQ", "2GB1", "1CTF", "4HHB", "1HRC", "3KUY", "1YCR", "2W5I", "2LZM", "2BG9", ]
        expected_sequences = [
            "NLYIQWLKDGGPSSGRPPPS",
            "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
            "MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE",
            "AAEEKTEFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK",
            "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
            "XGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
            "PEPAKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAKARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVALFEDTNLCAIHAKRVTIMPKDIQLARRIRGERSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGGSGRGKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKK",
            "SQIPASEQETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVVVNQQESSDSGTSVSENSQETFSDLWKLLPEN",
            "KETAAAKFERQHMDSSTSAASSSNYCNQMMKSRNLTKDRCKPVNTFVHESLADVQAVCSQKNVACKNGQTNCYQSYSTMSITDCRETGSSKYPNCAYKTTQANKHIIVACEGNPYVPVHFDASV",
            "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNCNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRCALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL",
            "SEHETRLVANLLENYNKVIRPVEHHTHFVDITVGLQLIQLINVDEVNQIVETNVRLRQQWIDVRLRWNPADYGGIKKIRLPSDDVWLPDLVLYNNADGDFAIVHMTKLLLDYTGKIMWTPPAIFKSYCEIIVTHFPFDQQNCTMKLGIWTYDGTKVSISPESDRPDLSTFMESGEWVMKDYRGWKHWVYYTCCPDTPYLDITYHFIMQRIPLYFVVNVIIPCLLFSFLTVLVFYLPTDSGEKMTLSISVLLSLTVFLLVIVELIPSTSSAVPLIGKYMLFTMIFVISSIIVTVVVINTHHRSPSTHSAIEGVKYIAEHMKSDEESSNAAEEWKYVAMVIDHILLCVFMLICIIGTVSVFAGRLIELSQEGSVMEDTLLSVLFENYNPKVRPSQTVGDKVTVRVGLTLTSLLILNEKNEEMTTSVFLNLAWTDYRLQWDPAAYEGIKDLSIPSDDVWQPDIVLMNNNDGSFEITLHVNVLVQHTGAVSWHPSAIYRSSCTIKVMYFPFDWQNCTMVFKSYTYDTSEVILQHALDAMINQDAFTENGQWSIEHKPSRKNWRSDDPSYEDVTFYLIIQRKPLFYIVYTIVPCILISILAILVFYLPPDAGEKMSLSISALLALTVFLLLLADKVPETSLSVPIIISYLMFIMILVAFSVILSVVVLNLHHRSPNTHEAVEAIKYIAEQLESASEFDDLKKDWQYVAMVADRLFLYIFITMCSIGTFSIFLDASHNVPPDNPFAVNEEERLINDLLIVNKYNKHVRPVKHNNEVVNIALSLTLSNLISLKETDETLTTNVWMDHAWYDHRLTWNASEYSDISILRLRPELIWIPDIVLQNNNDGQYNVAYFCNVLVRPNGYVTWLPPAIFRSSCPINVLYFPFDWQNCSLKFTALNYNANEISMDLIIDPEAFTENGEWEIIHKPAKKNIYGDKFPNGTNYQDVTFYLIIRRKPLFYVINFITPCVLISFLAALAFYLPAESGEKMSTAICVLLAQAVFLLLTSQRLPETALAVPLIGKYLMFIMSLVTGVVVNCGIVLNFHFRTPSTHSGIDSTNYIVKQIKEKNAYDEEVGNWNLVGQTIDRLSMFIITPVMVLGTIFIFVMGNFNRPPAKNEEGRLIEKLLGDYDKRIKPAKTLDHVIDVTLKLTLTNLISLNEKEEALTTNVWIEIQWNDYRLSWNTSEYEGIDLVRIPSELLWLPDVVLENNVDGQFEVAYYANVLVYNDGSMYWLPPAIYRSTCPIAVTYFPFDWQNCSLVFRSQTYNAHEVNLQLSAEEGIDPEDFTENGEWTIRHRPAKKNYNWQLTKDDIDFQEIIFFLIIQRKPLFYIINIIAPCVLISSLVVLVYFLPAQAGGQKCTLSISVLLAQTIFLFLIAQKVPETSLNVPLIGKYLIFVMFVSLVIVTNCVIVLNVSLRTPNTHSCVEACNFIAKSTKEQNDSGSENENWVLIGKVIDKACFWIALLLFSLGTLAIFLTGHLNQVPE",
        

        ]
        
        threshold = 1.0  # Define a threshold for similarity

        for i, code in enumerate(protein_codes):
            with self.subTest(code=code):
                # result = get_sequence(code)
                result = get_seqres_sequence(code)
                similarity = difflib.SequenceMatcher(None, result, expected_sequences[i]).ratio()
                print(f"\nSequence similarity for {code}:\t\t\t {similarity}")
                 
                if similarity < threshold:
                    print(f"Sequence:\t\t\t\t\t {result}")
                    print(f"Expected sequence taken from Fatsa file:\t {expected_sequences[i]}\n")

                # self.assertGreaterEqual(similarity, threshold, f"Sequence similarity for {code} below threshold. Similarity: {similarity}")

if __name__ == '__main__':
    unittest.main()

