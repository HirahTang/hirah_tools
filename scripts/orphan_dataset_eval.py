from hirah_tools.protein.pdb_eval import eval
from hirah_tools.protein.pdb_tools import alignment
import os
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
import pandas as pd
from tqdm import tqdm
seq_path = '/home/tanghan/OmegaFold/orphan_sequence'
pdb_pred_path = "/home/tanghan/OmegaFold/orphan_pred"
retrieve_pdb_path = "/home/tanghan/PDB"
orphan_ref_path = "/home/tanghan/OmegaFold/orphan_ref"

pdbl = PDBList()
seq_list = os.listdir(seq_path)
print(seq_list)
#protein(pd)
parser = PDBParser()
alignment_tool = alignment(retrieve_pdb_path, orphan_ref_path, pdb_pred_path)



result_df = {"seq_cat":[],
             "record":[],
             "pdb_id":[],
             "chain_id":[],
             "pdb_pred_omegafold":[],
             "pdb_ref":[],
             "seq_pred":[],
             "seq_pred_len":[],
             "seq_ref":[],
             "seq_ref_len":[],
             "lddt_omegafold":[],
             "tm_align_ref_omegafold":[],
             "align_length_omegafold":[],
             "rmsd_omegafold":[],
             }




def desire_chain(pdb_id, chain_id):
    name = "pdb{}.ent".format(pdb_id.lower())
    data = parser.get_structure(pdb_id, retrieve_pdb_path+'/'+name)
    io = PDBIO()
    for chain in data[0]:
        print(chain.get_id())
    io.set_structure(data[0][chain_id])
    io.save(f'{orphan_ref_path}/{pdb_id}_{chain_id}.pdb')
    #return data[0][chain_id]

for orig in seq_list[1:]:

    for record in tqdm(SeqIO.parse(seq_path+'/'+orig, "fasta")):
        #print(record.seq, record.id)
        pdb_pred = pdb_pred_path + '/' + record.id + '.pdb'
        pdb_info = record.id.split('_')
        pdb_id, chain_id = pdb_info[0].strip(','), pdb_info[2].strip(',')
        
        #---
        # Retrieve from PDB
        #--

        try:
            pdbl.retrieve_pdb_file(pdb_id, pdir = 'PDB', file_format = 'pdb')
        except:
            print("no available PDB", pdb_id)

            result_df['seq_cat'].append(orig)
            result_df['record'].append(record.id)
            result_df['pdb_id'].append(pdb_id)
            result_df['chain_id'].append(chain_id)
            result_df['pdb_pred_omegafold'].append(pdb_pred)
            result_df['pdb_ref'].append('NA')
            result_df['seq_pred'].append(record.seq)
            result_df['seq_pred_len'].append(len(record.seq))
            result_df['seq_ref'].append('NA')
            result_df['seq_ref_len'].append('NA')
            result_df['lddt_omegafold'].append('NA')
            result_df['tm_align_ref_omegafold'].append('NA')
            result_df['align_length_omegafold'].append('NA')
            result_df['rmsd_omegafold'].append('NA')

            continue
        
        #--
        # Process to leave the only chain
        #--
        
        try:
            alignment_tool.desire_chain(pdb_id, chain_id)

            alignment_tool.align(pdb_id, chain_id, record)

        
        except:
            print("There is no desired chain or problems occur during saving the chain", pdb_id, chain_id)

            result_df['seq_cat'].append(orig)
            result_df['record'].append(record.id)
            result_df['pdb_id'].append(pdb_id)
            result_df['chain_id'].append(chain_id)
            result_df['pdb_pred_omegafold'].append(pdb_pred)
            result_df['pdb_ref'].append('NA')
            result_df['seq_pred'].append(record.seq)
            result_df['seq_pred_len'].append(len(record.seq))
            result_df['seq_ref'].append('NA')
            result_df['seq_ref_len'].append('NA')
            result_df['lddt_omegafold'].append('NA')
            result_df['tm_align_ref_omegafold'].append('NA')
            result_df['align_length_omegafold'].append('NA')
            result_df['rmsd_omegafold'].append('NA')

            continue
        
        
        pdb_ref = f'{orphan_ref_path}/{pdb_id}_{chain_id}.pdb'
        eval_machine = eval(pdb_pred, pdb_ref, lddt_path="/home/tanghan/psp-pipeline/pipeline/tool/lddt-linux/lddt")


        result_df['seq_cat'].append(orig)
        result_df['record'].append(record.id)
        result_df['pdb_id'].append(pdb_id)
        result_df['chain_id'].append(chain_id)
        result_df['pdb_pred_omegafold'].append(pdb_pred)
        result_df['pdb_ref'].append(pdb_ref)
        result_df['seq_pred'].append(record.seq)
        result_df['seq_pred_len'].append(len(record.seq))
        try:
            eq_pred, seq_pred_len, seq_ref, seq_ref_len = eval_machine.seq_len()
            result_df['seq_ref'].append(seq_ref)
            result_df['seq_ref_len'].append(seq_ref_len)
        except:
            print("Could not retrieve ref sequence & length")
            result_df['seq_ref'].append('NA')
            result_df['seq_ref_len'].append('NA')
            result_df['lddt_omegafold'].append('NA')
            result_df['tm_align_ref_omegafold'].append('NA')
            result_df['align_length_omegafold'].append('NA')
            result_df['rmsd_omegafold'].append('NA')
            continue

        print(eval_machine.seq_len())
        try:
            lddt_score, lddt_report = eval_machine.lddt()
            result_df['lddt_omegafold'].append(lddt_score)
        except:
            print("Could not calculate lddt")
            result_df['lddt_omegafold'].append('NA')
            
        try:

            tm_align_res = eval_machine.tmalign()
            result_df['tm_align_ref_omegafold'].append(tm_align_res['tm_score'][1])
            result_df['align_length_omegafold'].append(tm_align_res['align_length'])
            result_df['rmsd_omegafold'].append(tm_align_res['rmsd'])
            
            print("A new line added")
        except:
            print('Could not calculate TM align')
            result_df['tm_align_ref_omegafold'].append("NA")
            result_df['align_length_omegafold'].append("NA")
            result_df['rmsd_omegafold'].append("NA")





print(result_df)
df = pd.DataFrame(data = result_df)
df.to_csv('/home/tanghan/OmegaFold/orphan_eval.csv')
