import gzip
import pickle
import json
from bioservices.apps.fasta import FASTA
from rdkit import Chem

from Bio.PDB import *

from BioIO.utility import Read


def read_pickle(pickle_file):
    p = pickle.load(open(pickle_file, 'rb'))
    return p


def write_pickle(pickle_file, object):
    pickle.dump(object, open(pickle_file, 'wb'))


def write_zip(file_name, object_name):
    with gzip.open(file_name, 'wb') as f:
        f.write(pickle.dumps(object_name, protocol=4))
        f.close()


def read_zip(file_name):
    with gzip.open(file_name, 'rb') as f:
        output_object = pickle.load(f)
    return output_object


def open_file(file_name):
    f = open(file_name, 'rb')
    lines = f.readlines()
    return lines


def write_file(file_name, object_name):
    f = open(file_name, 'wb')
    f.writelines(object_name)
    f.close()
    
def pd_read_table(file_name):
    t = pd.read_table(file_name)
    return t

def read_json(file_name):
    f = open(file_name).read()
    data = json.loads(f)
    return data


def read_fasta(file_name):
    '''
    :param file_name:
    :return: the fasta raw data
    '''
    f = FASTA.read_fasta(file_name)
    return f._fasta


# Leos: the output can be different from the input one because of removing errors
def read_sdf(file_name):
    '''
    :param file_name:
    :return: the output is the mol object list obtained from rdkit
    '''
    r = Read(input=file_name, format='sdf')
    r.read()
    return r.mols


def read_mol2(file_name):
    '''
    :param file_name:
    :return: the output is the mol object list obtained from rdkit
    Note: the format is useless here
    '''
    r = Read(input=file_name, format='sdf')
    r.read()
    return r.mols


def read_mol(file_name):
    '''
    :param file_name:
    :return: the output is the mol object list obtained from rdkit
    Note: the format is useless here
    '''
    r = Read(input=file_name, format='sdf')
    r.read()
    return r.mols


def read_pandas_smiles(file_name):
    '''
    :param file_name:
    :return: the output is the mol object list obtained from rdkit
    '''
    r = Read(input=file_name, format='smiles')
    r.read()
    return r.mols


def read_pdb_file(pdb_file):
    '''
    :param pdb_file:
    :return: the structure object from biopython
    '''
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_file)
    return structure


def write_pdb():
    pass


def write_sdf(mols, file_name):
    """
    :param mols: the rdkit format of mol or list of mols
    :param file_name: str, the name of the file to write
    :return: no return
    """
    w = Chem.SDWriter(file_name)
    if type(mols) == list:
        for m in mols:
            w.write(m)
    else:
        w.write(mols)


def read_gro_file():
    pass


def read_cms_file():
    pass
