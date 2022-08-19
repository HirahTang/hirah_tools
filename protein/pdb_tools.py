import difflib
from Bio.PDB import PDBParser
from .pdb_eval import eval
from Bio.PDB.PDBIO import PDBIO



class alignment:
    def __init__(self, retrieve_pdb_path, orphan_ref_path, pdb_pred_path) -> None:
        self.retrieve_pdb_path = retrieve_pdb_path
        self.orphan_ref_path = orphan_ref_path
        self.pdb_pred_path = pdb_pred_path


    def desire_chain(self, pdb_id, chain_id):
        """
        Check whether specific chain id exists in the PDB
        Input: PDB_id, Chain_id
        Output: Save the single chain PDB to [orphan_ref_path]/[pdb_id]_[chain_id].pdb
        """
        parser = PDBParser()
        name = "pdb{}.ent".format(pdb_id.lower())
        data = parser.get_structure(pdb_id, self.retrieve_pdb_path+'/'+name)
        io = PDBIO()
        for chain in data[0]:
            print(chain.get_id())
        io.set_structure(data[0][chain_id])
        io.save(f'{self.orphan_ref_path}/{pdb_id}_{chain_id}.pdb')

    def get_overlap(s1, s2):
        
        """
        Get the maximum overlap part between input sequence s1 and input sequence s2
        """

        s = difflib.SequenceMatcher(None, s1, s2)
        pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
        return s1[pos_a:pos_a+size], pos_a, pos_b, size

    def reindex(self, pdb_id, chain_id, pos_a=0):

        """
        Reindex the chain to start from pos_a, and continuously index the rest sequence.
        Input: pdb_id, chain_id, pos_a=0, un-reindexed pdb at [orphan_ref_path/pdb_id_chain_id.pdb]
        Output: Reindexed pdb at [orphan_ref_path/pdb_id_chain_id.pdb]
        """

        parser = PDBParser()
        name = "{}_{}.pdb".format(pdb_id, chain_id)
        structure = parser.get_structure(pdb_id, f'{self.orphan_ref_path}/{name}')
        
        for res in structure[0][chain_id]:
            res.id = (' ', pos_a, ' ')
            pos_a += 1
        io = PDBIO()
        io.set_structure(structure[0][chain_id])
        io.save(f'{self.orphan_ref_path}/{pdb_id}_{chain_id}.pdb')


    def align(self, pdb_id, chain_id, record):

        pdb_ref = f'{self.orphan_ref_path}/{pdb_id}_{chain_id}.pdb'
        pdb_pred = self.pdb_pred_path + '/' + record.id + '.pdb'
        
        eval_tool = eval(pdb_pred, pdb_ref)
        seq_a, seq_a_len, seq_b, seq_b_len = eval_tool.seq_len()

        if seq_a_len == seq_b_len:
            self.reindex(pdb_id, chain_id, 0)
        else:
            overlap, pos_a, pos_b, size = alignment.get_overlap(seq_a, seq_b)
            print(pos_a)
            self.reindex(pdb_id, chain_id, pos_a)

        

    
