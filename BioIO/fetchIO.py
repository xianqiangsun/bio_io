from bioservices import UniProt
from bioservices import PDB
from bioservices.apps.fasta import FASTA
from bioservices import ChEMBL

# Leos from uniprot id to the protein name
f = FASTA()
u = UniProt(cache=True, verbose=False)
s = PDB()
c = ChEMBL(verbose=False)


def uniprot_id_2_protein_name(uni_id):
    '''
    :param uni_id:
    :return: a list of protein names to the corresponding uniprot id
    '''
    all_names = u.search(str(uni_id), columns='protein names').split('\n')[1:-1]
    print(all_names[0])
    protein_names = []
    for each_name in all_names:
        each_name = each_name.rstrip()
        each_name = each_name.split(' (')[0].upper()
        if each_name not in protein_names:
            protein_names.append(each_name)
    print(uni_id, protein_names)
    return protein_names

def uniprot_id_2_protein_family(uni_id):
    '''
    :param uni_id:
    :return: a list of protein family to the corresponding uniprot id
    check all coloumns  bioservices.uniprot.UniProt._valid_columns
    '''
    all_families = u.search(str(uni_id), columns='families').split('\n')[1:-1]
    protein_families = []
    for each_name in all_families:
        if each_name!= '':
            protein_families.append(each_name)
            print (uni_id,each_name)
    return protein_families

def uniprot_id_2_pdb_id(uni_id):
    """
    :param uni_id:
    :return: a list of pdb structures corresponding to the pdb
    """
    output = u.mapping("ID", "PDB_ID", uni_id)
    output_pdb_id = list(output.values())[0]
    return output_pdb_id


def uniprot_id_2_gene_name(uni_id):
    '''
    :param uni_id:
    :return: the gene names for the corresponding uni_id
    '''
    gene_names = u.search(str(uni_id), columns='genes').split('\n')[1].split()
    return gene_names


def uniprot_id_2_sequence(uni_id):
    '''
    :param uni_id:
    :return: the sequence of the uni_id
    '''
    sequence = u.retrieve(str(uni_id), "fasta")
    return sequence


def uniprot_id_2_fasta(uni_id):
    '''
    :param uni_id:
    :return: the sequence of the uni_id
    '''
    fasta_object = f.load(uni_id)
    return fasta_object


def pdb_id_2_pdb(pdb_id):
    '''
    :param pdb_id:
    :return: pdb structure of the pdb_id the output is in xml format
    '''
    pdb_file = s.get_file(pdb_id, "pdb")
    return pdb_file


def chembl_id_structures(chembl_id):
    '''
    Note: refer to the function in bioservises
    :param chembl_id:
    :return: json format of the chemblid
    '''
    resjson = c.get_compounds_by_chemblId(str(chembl_id))
    # Note json format of the return
    return resjson


def chembl_id_targets(chembl_id):
    '''
    Note: refer to the function in bioservises
    :param chembl_id:
    :return:
    '''
    resjson = c.get_target_by_chemblId(chembl_id)
    return resjson


def chembl_id_activities(chembl_id):
    '''
    Note: refer to the function in bioservises
    :param chembl_id:
    :return:
    '''
    resjson = c.get_compounds_activities(chembl_id)
    return resjson
