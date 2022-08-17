import Bio
from Bio.PDB import PDBList
import os
import subprocess
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

class eval:
    def __init__(self, pdb_a, pdb_b, lddt_path = "/home/tanghan/hirah_tools/protein/lddt-linux/lddt") -> None:
        
        self.pdb_a = pdb_a
        self.pdb_b = pdb_b
        self.lddt_path = lddt_path
        self.pdb_a_name = self.pdb_a.split('/')[-1].split('.')[0]
        self.pdb_b_name = self.pdb_b.split('/')[-1].split('.')[0]


    def lddt(self):
        """
        Return the lddt score of the predicted pdb to the target pdb
        Output: lddt score, lddt report
        """
        execute_command = f"{self.lddt_path} -c {self.pdb_a} {self.pdb_b}"
        report = subprocess.check_output(execute_command, shell = True).decode("utf-8")
        prefix = "Global LDDT score: "
        for line in report.split("\n"):
            if line.startswith(prefix):
                score = float(line[len(prefix) :])
        return score, report



    def seq_len(self):
        """
        Return the pdb length of the pdbs, separatly, for the cases they are different
        Output: seq_a(predicted pdb), seq_a_len, seq_b(reference pdb), seq_b_len
        """

        parser = PDBParser(QUIET=True)
        
        structure_a = parser.get_structure('{}_pred'.format(self.pdb_a_name), self.pdb_a)
        structure_b = parser.get_structure('{}_ref'.format(self.pdb_b_name), self.pdb_b)

        seq_a, seq_a_len = self.seq_from_pdb(structure_a)
        seq_b, seq_b_len = self.seq_from_pdb(structure_b)

        return seq_a, seq_a_len, seq_b, seq_b_len
        
                #print('>some_header\n',''.join(seq))
    
    def get_pdb_chain(self):
        
        """
        From the nomenclature of PDB Entry ID,Entity ID,Chain Asym ID,Chain Auth Asym ID
        Output: PDB Entry ID, Chain ID
        """
        
        pdb_info = self.pdb_a.split('/')[-1].split('.')
        return pdb_info[0], pdb_info[2]

    def seq_from_pdb(structure):
        
        """
        Get sequence in string and len(seq) from PDB input (parsed by Bio.PDBParser)
        """

        d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        for model in structure:
            for chain in model:
                seq = []
                for residue in chain:
                    #print(residue)
                    try:
                        seq.append(d3to1[residue.resname])
                    except:
                        continue
                seq_out = ''.join(seq)

        return seq_out, len(seq_out)

        
