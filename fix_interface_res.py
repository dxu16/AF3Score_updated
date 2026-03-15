from Bio.PDB import PDBParser, MMCIFParser
from collections import defaultdict
import numpy as np

def get_interface_res_from_pdb(pdb_file, chain1="A", chain2="B", dist_cutoff=10):
    if str(pdb_file).endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    model = structure[0]
    
    chain_coords = defaultdict(dict)
    for chain in model:
        cid = chain.id
        if cid not in (chain1, chain2):
            continue
        for residue in chain:
            if "CA" in residue:
                chain_coords[cid][residue.id[1]] = residue["CA"].get_coord()
                
    chain_1_res = sorted(chain_coords[chain1].keys())
    chain_2_res = sorted(chain_coords[chain2].keys())
    if not chain_1_res or not chain_2_res:
        return [], []
    
    chain_1_coords = np.array([chain_coords[chain1][res] for res in chain_1_res])
    chain_2_coords = np.array([chain_coords[chain2][res] for res in chain_2_res])
    
    dist = np.sqrt(np.sum((chain_1_coords[:, None, :] - chain_2_coords[None, :, :]) ** 2, axis=2))
    interface_residues = np.where(dist < dist_cutoff)
    
    interface_1 = sorted(set(chain_1_res[i] for i in interface_residues[0]))
    interface_2 = sorted(set(chain_2_res[i] for i in interface_residues[1]))
    return interface_1, interface_2
