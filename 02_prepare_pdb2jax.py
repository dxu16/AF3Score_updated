import os
import sys
import time
from typing import Dict, Tuple, List
import numpy as np
import jax
import jax.numpy as jnp
from Bio import PDB
from Bio.PDB.Structure import Structure
import h5py
import multiprocessing as mp
from tqdm import tqdm
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import argparse
import io

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

sys.path.insert(0, '/home/dxu39/scr4_jgray21/dxu39/software/AF3Score_updated/src')
from alphafold3.structure import parsing as structure_parsing
from alphafold3.constants import chemical_components
from alphafold3.model import features
from alphafold3.model.pipeline import structure_cleaning
from alphafold3.model.atom_layout import atom_layout

# Initialize CCD (lazy load overhead once)
ccd = None

def get_ccd():
    global ccd
    if ccd is None:
        ccd = chemical_components.cached_ccd()
    return ccd

def find_bucket_size(seq_length: int, buckets: List[int]) -> int:
    for bucket in buckets:
        if bucket >= seq_length:
            return bucket
    return buckets[-1]

def save_traced_array(
    traced_array: jax.Array,
    seq_length: int,
    save_path: str,
    metadata: Dict = None,
) -> None:
    if not save_path.endswith(".h5"):
        save_path = save_path + ".h5"

    numpy_array = np.array(traced_array)

    with h5py.File(save_path, "w") as f:
        f.create_dataset("coordinates", data=numpy_array)
        f.create_dataset("seq_length", data=seq_length)
        f.create_dataset("shape", data=numpy_array.shape)

        if metadata:
            metadata_grp = f.create_group("metadata")
            for key, value in metadata.items():
                metadata_grp.attrs[key] = value

def load_structure_as_cif_string(file_path: str) -> str:
    """Loads a PDB/CIF and returns CIF string content, normalizing PDBs."""
    if file_path.endswith(".cif"):
        with open(file_path, "r") as f:
            return f.read()
    else:
        # Convert PDB to CIF string using Biopython
        parser = PDB.PDBParser(QUIET=True)
        structure_id = os.path.basename(file_path).split(".")[0]
        struct = parser.get_structure(structure_id, file_path)
        io_writer = PDB.MMCIFIO()
        io_writer.set_structure(struct)
        s_io = io.StringIO()
        io_writer.save(s_io)
        return s_io.getvalue()

def pdb_to_traced_array(
    pdb_path: str,
    chain_ids: List[str] = None,
    num_copies: int = 1,
    save_path: str = None,
    max_length: int = 3072,
) -> Tuple[jnp.ndarray, int]:
    
    
    cif_str = load_structure_as_cif_string(pdb_path)
    
    struct = structure_parsing.from_mmcif(cif_str)
    c = get_ccd()
    
    cleaned_struc, _ = structure_cleaning.clean_structure(
        struct,
        ccd=c,
        drop_non_standard_atoms=True,
        drop_missing_sequence=True,
        filter_clashes=False,
        filter_crystal_aids=False,
        filter_waters=True,
        filter_hydrogens=True,
        filter_leaving_atoms=False,
        only_glycan_ligands_for_leaving_atoms=True,
        covalent_bonds_only=True,
        remove_polymer_polymer_bonds=True,
        remove_bad_bonds=True,
        remove_nonsymmetric_bonds=False,
    )

    flat_output_layout = atom_layout.atom_layout_from_structure(cleaned_struc)
    _, layout, _ = features.tokenizer(
        flat_output_layout=flat_output_layout,
        ccd=c,
        max_atoms_per_token=24,
        flatten_non_standard_residues=True,
        logging_name="eval_target",
    )

    num_tokens = layout.shape[0]
    seq_length = num_tokens

    if seq_length > 3072:
        print(f"Warning: {os.path.basename(pdb_path)} - Sequence length {seq_length} exceeds maximum allowed length 3072")
        return None, seq_length

    target_length = find_bucket_size(seq_length, [3072])
    
    # Build fast mapper for coordinates
    mapper = {}
    c_x = np.array(cleaned_struc.atom_x)
    c_y = np.array(cleaned_struc.atom_y)
    c_z = np.array(cleaned_struc.atom_z)
    c_aname = np.array(cleaned_struc.atom_name)
    c_rname = np.array(cleaned_struc.res_name)
    c_rid = np.array(cleaned_struc.res_id)
    c_cid = np.array(cleaned_struc.chain_id)

    for i in range(len(c_aname)):
        cn = c_cid[i].decode('utf-8') if isinstance(c_cid[i], bytes) else c_cid[i]
        rn = c_rname[i].decode('utf-8') if isinstance(c_rname[i], bytes) else c_rname[i]
        an = c_aname[i].decode('utf-8') if isinstance(c_aname[i], bytes) else c_aname[i]
        key = (cn, int(c_rid[i]), rn, an)
        mapper[key] = (c_x[i], c_y[i], c_z[i])

    coords = np.zeros((num_tokens, 24, 3), dtype=np.float32)

    for t in range(num_tokens):
        for a in range(24):
            an = layout.atom_name[t, a]
            if not an: continue
            cn = layout.chain_id[t, a]
            rn = layout.res_name[t, a]
            
            an = an.decode('utf-8') if isinstance(an, bytes) else an
            cn = cn.decode('utf-8') if isinstance(cn, bytes) else cn
            rn = rn.decode('utf-8') if isinstance(rn, bytes) else rn
            
            key = (cn, int(layout.res_id[t, a]), rn, an)
            if key in mapper:
                coords[t, a] = mapper[key]

    # Pad with zeros to meet bucket size
    if seq_length < target_length:
        padding = np.zeros((target_length - seq_length, 24, 3), dtype=np.float32)
        coords = np.concatenate([coords, padding], axis=0)

    # Tile the array to create multiple copies
    coords_repeated = np.stack([coords] * num_copies)
    traced_array = coords_repeated

    
    
        

    

    if save_path:
        metadata = {
            "pdb_file": os.path.basename(pdb_path),
            "chain_ids": ",".join(chain_ids) if chain_ids else "all",
            "num_copies": num_copies,
            "original_length": seq_length,
            "padded_length": target_length,
        }
        save_traced_array(traced_array, seq_length, save_path, metadata)
    
    return traced_array, seq_length

def process_single_file(args):
    input_path, output_path, chain_ids, num_copies = args
    try:
        result = pdb_to_traced_array(
            pdb_path=input_path,
            chain_ids=chain_ids,
            num_copies=num_copies,
            save_path=output_path,
            max_length=3072,
        )
        if result[0] is None:
            return (False, input_path)
        return (True, input_path)
    except Exception as e:
        print(f"Error processing {os.path.basename(input_path)}: {str(e)}")
        return (False, input_path)

def process_pdb_folder(
    pdb_folder: str,
    output_folder: str,
    chain_ids: List[str] = None,
    num_copies: int = 1,
    num_workers: int = None,
) -> None:
    os.makedirs(output_folder, exist_ok=True)

    processing_args = []
    for filename in os.listdir(pdb_folder):
        if filename.endswith(".pdb") or filename.endswith(".cif"):
            input_path = os.path.join(pdb_folder, filename)
            output_path = os.path.join(
                output_folder, f"{os.path.splitext(filename)[0]}.h5"
            )
            processing_args.append(
                (input_path, output_path, chain_ids, num_copies)
            )

    if not processing_args:
        print("No valid files to process.")
        return

    print(f"🔵 Preparing to process {len(processing_args)} files...")
    print(f"🔵 Using {num_workers} processes for parallel execution.")

    # Pass sequential to test safely.
    # multiprocessing with JAX can be an issue if JAX grabs all GPU memory, 
    # but the old script used map anyway.
    
    # We'll use mp pool:
    if True:
        results = []
        for args in tqdm(processing_args, desc="Processing files"):
            results.append(process_single_file(args))
        
    failed = [path for success, path in results if not success]
    print(f"🔴 Failed: {len(failed)} files")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_folder", type=str, required=True)
    parser.add_argument("--output_folder", type=str, required=True)
    parser.add_argument("--num_workers", type=int, default=4)
    args = parser.parse_args()

    try:
        bucket = int(os.path.basename(args.pdb_folder).split("_")[-1])
    except (ValueError, IndexError):
        bucket = 3072

    
    BUCKETS = [bucket]
    
    process_pdb_folder(args.pdb_folder, args.output_folder, num_workers=args.num_workers)

if __name__ == "__main__":
    main()
