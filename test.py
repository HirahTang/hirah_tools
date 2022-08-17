
from protein.pdb_eval import eval
pdb_pred = "/home/tanghan/OmegaFold/orphan_pred/6BYJ_2_G,.pdb"
pdb_target = "/home/tanghan/OmegaFold/6BYJ_G.pdb"

eval_machine = eval(pdb_pred, pdb_target)

print(eval_machine.lddt())
print(eval_machine.seq_len())
print(eval_machine.get_pdb_chain())
print(eval_machine.tmalign())