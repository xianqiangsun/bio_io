from .fileIO import read_zip,read_pickle
from rdkit import Chem
import os
from rdkit.Chem.SaltRemover import SaltRemover
import pandas as pd

smile_data = ['%', 't', 'e', '9', 'i', "7", "6", "o", "]", "3", "s", \
              "(", "-", "S", "/", "B", "4", "[", ")", "#", "I", "l", \
              "O", "H", "c", "1", "@", "=", "n", "P", "8", "C", "2", \
              "F", "5", "r", "N", "+", "\\", " "]


def get_suffix(file_name):
    file_name, file_extension = os.path.splitext(file_name)
    return file_extension


def read_sdf(sdf_file, start, slice, stride):
    m = [i for i in Chem.ForwardSDMolSupplier(sdf_file)][start:slice:stride]
    return m


def read_smiles(smile_file, start, slice, stride):
    # m = []
    # for i in Chem.SmilesMolSupplier(smile_file):
    #    i=CanonicalRankAtoms(i)
    #    m.append(i)
    m = [i for i in Chem.SmilesMolSupplier(smile_file)][start:slice:stride]
    return m


def read_csv_smiles(pd_project, start, slice, stride):
    all_s = pd_project['smiles'][start:slice:stride]
    m = []
    for each_mol in all_s:
        each_mol = each_mol.rstrip()
        em = Chem.MolFromSmiles(each_mol)
        m.append(em)
    return m


def check_smiles(smiles):
    output = True
    for each_latter in smiles:
        if each_latter not in smile_data:
            output = False
            print(smiles, "does not pass the check because it is not in our atom database")
            break
    return output


def remove_salt(mol):
    remover = SaltRemover(defnData="[te+,K+,Tc-3,Br,Br-,I-,Cl,Cl-,Li+,Na+]")
    res = remover(mol)
    remover = SaltRemover(defnData="CS(=O)(=O)O")
    res = remover(res)
    remover = SaltRemover(defnData="O=S(=O)(O)O")
    res = remover(res)
    remover = SaltRemover(defnData="CS(=O)(=O)O")
    res = remover(res)
    remover = SaltRemover(defnData="O=C(O)C(=O)O")
    res = remover(res)
    remover = SaltRemover(defnData="O=C(O)C=CC(=O)O")
    res = remover(res)
    remover = SaltRemover(defnData="O=C(O)C(F)(F)F")
    res = remover(res)
    remover = SaltRemover(defnData="CCN(CC)CC")
    res = remover(res)
    remover = SaltRemover(defnData="CC(=O)O")
    res = remover(res)
    remover = SaltRemover(defnData="O=C(O)CCC(=O)O")
    res = remover(res)
    remover = SaltRemover(defnData="NCCO")
    res = remover(res)
    return res


class Read(object):
    def __init__(self, input, format, **kwargs):
        '''
        :param input: the input files in smiles or sdf, mol, mol2
        :param format: {smiles, sdf}
        '''
        if "stride" in kwargs:
            self.stride = kwargs["stride"]
        else:
            self.stride = 1
        if "slice" in kwargs:
            self.slice = kwargs["slice"]
        else:
            self.slice = None
        if "start" in kwargs:
            self.start = kwargs["start"]
        else:
            self.start = 0

        file_extension = get_suffix(input)
        if file_extension in ['.gz', '.gzip']:
            self.input = read_zip(input)
        elif file_extension in ['.pkl', '.pickle']:
            self.input = read_pickle(input)
        elif file_extension in ['.csv']:
            self.input = pd.read_csv(input)
        else:
            self.input = input

        self.format = format
        self.mols = None
        if 'cutoff' in kwargs:
            self.cutoff = kwargs['cutoff']
        else:
            self.cutoff = None
        self.can_mol = None
        self.can_smi = None
        if 'max_no' in kwargs:
            self.max_no = kwargs['max_no']
        else:
            self.max_no = None

    def read(self):
        print('reading the molecules from file...')
        if self.format == 'smiles':
            self.mols = read_smiles(self.input, self.start, self.slice, self.stride)
        elif self.format == 'sdf':
            self.mols = read_sdf(self.input, self.start, self.slice, self.stride)
        elif self.format == 'csv':
            self.mols = read_csv_smiles(self.input, self.start, self.slice, self.stride)

    def remove_error(self):
        print('preforming remove errors of the molecules...')
        can_mol = []
        can_smi = []
        for i in self.mols:
            try:
                res = remove_salt(i)
                atoms_numbers = res.GetNumAtoms()
                each_smi = Chem.MolToSmiles(res)
                if check_smiles(each_smi) and atoms_numbers != 0:
                    if self.max_no:
                        if all([len(each_smi) <= self.max_no, '[' not in each_smi]):
                            # print (atoms_numbers)
                            can_smi.append(each_smi)
                            each_mol = Chem.MolFromSmiles(each_smi)
                            can_mol.append(each_mol)
                    else:
                        can_smi.append(each_smi)
                        each_mol = Chem.MolFromSmiles(each_smi)
                        can_mol.append(each_mol)

            except:
                pass
        print("Total:", len(can_smi), 'molecules are loaded...')
        self.can_mol = can_mol
        self.can_smi = can_smi
