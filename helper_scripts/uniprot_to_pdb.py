import numpy as np
import requests
import pypdb
import os 
import gzip
import json

pdbs_folder = '<name>'

def create_pdb_list(pdb_folder):

  pdb_files = sorted(os.listdir(pdb_folder))
  filenames = sorted([i for i in pdb_files if i[0]!='.'])
  pdb_names = [i.split('.')[0] for i in pdb_files if i[0]!='.']
  rcsb_accession_codes = []

  for index in range(len(filenames)):
    with gzip.open(pdb_folder+'/'+filenames[index],'r') as fin:
      lines_pdb = fin.readlines()
      ATOM_only = [line for line in lines_pdb if line.startswith(b'ATOM')]
      chain_name = ATOM_only[0][21:22].decode('utf-8')
      rcsb_accession_code = pdb_names[index]+'.'+chain_name
      rcsb_accession_codes.append(rcsb_accession_code)
  return pdb_names, rcsb_accession_codes

pdb_names, rcsb_codes = create_pdb_list(pdbs_folder)

def map_pdb_uniprot(query_id):
  dict_mapping = {}
  query = """
  {{
    annotations(reference:PDB_INSTANCE sources:[UNIPROT] queryId:"{query_id}"){{target_id features {{feature_id description name provenance_source type feature_positions {{beg_ori_id beg_seq_id end_ori_id end_seq_id}}}}}}
  }}
  """.format(query_id=query_id)
  accesions_list = []
  url = 'https://1d-coordinates.rcsb.org/graphql?{}'
  response = requests.get(url, params={'query':query}).json()
  success = len(response['data']['annotations'])
  if success > 0:
    accession_uniprot = response['data']['annotations'][0]['target_id']
    if query_id not in dict_mapping.keys():
      dict_mapping[query_id] = accession_uniprot
    else: 
      dict_mapping[query_id].append(accession_uniprot)
  else: 
    print('--- failed to find UniProt accession ---')
    print(response)
    dict_mapping[query_id] = None
  return dict_mapping

results_mapping = []
for query_id in rcsb_codes:
  x = map_pdb_uniprot(query_id)
  results_mapping.append(x)

# with open('loc/<name>, 'w') as pdb_uniprot_mapping:
#    json.dump(results_mapping, pdb_uniprot_mapping)
