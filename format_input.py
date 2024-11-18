from typing import Union
from pathlib import Path
import yaml
import dnaio
import numpy as np
from collections import defaultdict, OrderedDict
import click
import sys


def get_protein_yaml(id:Union[str, list], sequence: str, msa: Union[Path,str]):
    # if type(id) == list:
    #     id = ",".join(id)
    #     id = f"[id]

    output = {"protein":
              {"id":id,
              "sequence":sequence,
              "msa":str(msa)}
            }
    return output
    
def get_ligand_yaml(id:Union[str, list], smiles: str):
    if type(id) == list:
        id = ",".join(id)
        id = f"[{id}]"
    output = {"ligand":
              {"id":id,
              "smiles":sequence}
            }
    return output

def make_yaml(protein_sequences: Union[list,dict] = None,ligand_sequences: Union[list,dict] = None):
    document = {"version":1}
    document["sequences"]=[]
    if protein_sequences:
        if type(protein_sequences) == dict:
            document["sequences"].append(protein_sequences)
        elif type(protein_sequences) == list:
            for ps in protein_sequences:
                document["sequences"].append(ps)
                
    if ligand_sequences:
        if type(ligand_sequences) == dict:
            document["sequences"].append(ligand_sequences)
        elif type(ligand_sequences) == list:
            for ls in ligand_sequences:
                document["sequences"].append(ps)
    return document

def concat_list(list_in):
    list_out = []
    for x in list_in:
        if type(x) == list:
            list_out+=concat_list(x)
        else:
            list_out.append(x)
    return list_out


def get_from_folder(folder:Union[Path,str]):
    if type(folder) == str:
        folder = Path(folder)
    msa = sorted(folder.rglob("*.a3m"))
    msa = msa[np.argmax([x.stat().st_size for x in msa])].resolve()
    fasta = Path(f"{folder}.fasta")
    chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    seq_chains = defaultdict(list)
    sequence = [x.sequence.split(":") for x in dnaio.FastaReader(str(fasta))]
    sequence = concat_list(sequence)
    
    for xi,x in enumerate(sequence):
        c = chains[xi]
        seq_chains[str(x)].append(c)
    
    return [get_protein_yaml(id,sequence,msa=msa) for sequence, id in seq_chains.items()]



@click.command()
@click.argument('folder', required=False)  # Positional argument
def main(folder):
    # print("input: ",folder)

    document = make_yaml(get_from_folder(folder))
    output_yaml = Path(folder) / f"run.yaml"
    # print("output: ",output_yaml)
    with open(output_yaml,"w+") as h:
        yaml.safe_dump(document,h)
    sys.stdout.write(str(output_yaml.resolve()))
    
    # documents = make_yaml(get_from_folder(folder))
    # output_yaml = Path(folder) / f"run_{np.random.random():.3f}.yaml"
    # print(documents)
    # with open(str(output_yaml),"w+") as h:
    #     yaml.safe_dump(documents)
        

if __name__ == '__main__':
    main()