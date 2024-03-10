import numpy as np
import requests
import pypdb
import os 
import gzip
import json
import shutil
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

pdbs_folder = '/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/pdbs'
af2_folder = '/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/AF2_predictions'

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

def map_pdb_uniprot(rcsb_codes):
	results_mapping = []
	for query_id in rcsb_codes:
		result = get_pdb_uniprot(query_id)
		results_mapping.append(result)
	def get_pdb_uniprot(query_id):
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

# with open('/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/test_mapping.txt', 'w') as pdb_uniprot_mapping:
#   json.dump(results_mapping, pdb_uniprot_mapping)

with open('/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/test_mapping.txt', 'r') as json_data:
	test = json.load(json_data)

selected_100 = test[:100]

def download_af_db(accession_codes):
	folder_pdbs_out = '/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/AF2_predictions'
	missing_predictions = []
	for entry in accession_codes:
		version='v4'
		try:
			entry_code = entry.split('\n')[0]
		except AttributeError:
			print('NoneType error! Perhaps the UniProt accession is missing')
		alphafold_ID = 'AF-{}-F1-model_{}'.format(str(entry_code), str(version))
		model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_ID}.pdb'
		pdb_out = os.path.join(folder_pdbs_out, alphafold_ID+'.pdb')
		os.system(f'curl --max-time 900 {model_url} -o {alphafold_ID}'+'.pdb')
		shutil.move('/Users/ncyx'+'/'+str(alphafold_ID)+'.pdb', pdb_out)
		file_stats = os.stat(pdb_out)
		file_size = file_stats.st_size
		if file_size == 127:
			missing_predictions.append(entry)
	def write_missing_predictions(missing_predictions):
		folder_out = '/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/missing_predictions.txt'
		with open(folder_out, 'w') as missing_predictions_save:
			for index in range(len(missing_predictions)):
				if index == len(missing_predictions)-1:
					missing_predictions_save.write(missing_predictions[index])
				else:
					missing_predictions_save.write(missing_predictions[index]+'\n')
	write_missing_predictions(missing_predictions)
	return missing_predictions

def find_missing_af2pred(af2_folder):
	missing_predictions = []
	filenames = os.listdir(af2_folder)
	for filename in filenames:
		if filename[0] != '.':
			file_stats = os.stat(af2_folder+'/'+filename)
			file_size = file_stats.st_size
			if file_size == 127:
				missing_predictions.append(filename.split('-')[1])
	return missing_predictions
missing_pred = find_missing_af2pred(af2_folder)

def gen_uniprot_for_pdbs(json_data, missing_af2_pred):
	uniprot_codes = [list(i.values())[0] for i in json_data]
	upd_uniprot_names = []
	for index in range(len(json_data)):
		for key, value in json_data[index].items():
				if value not in upd_uniprot_names:
					upd_uniprot_names.append(value)
	entries_present = []
	trunc_json = []
	for i in json_data:
		if list(i.values())[0] in upd_uniprot_names:
				if list(i.values())[0] in entries_present:
					continue
				else:
					if list(i.values())[0] not in missing_pred:
						if list(i.values())[0] != None:
							entries_present.append(list(i.values())[0])
							trunc_json.append(i)
	return trunc_json

dict_trunc_json = gen_uniprot_for_pdbs(selected_100, missing_pred)

def find_af2_for_pdbs(dict_pdbs_uniprot, af2_folder):
	list_af2_locs = []
	uniprot_pdb_accession = [list(i.values())[0] for i in dict_pdbs_uniprot]
	# print(uniprot_pdb_accession)
	for filename in os.listdir(af2_folder):
		if filename[0] != '.':
			if filename.split('-')[1] in uniprot_pdb_accession:
				list_af2_locs.append(af2_folder+'/'+filename)
	return list_af2_locs

def select_mapped_pdbs(dict_pdbs_uniprot, pdb_folder):
	filenames = [i for i in os.listdir(pdb_folder) if i[0]!='.']
	pdb_accessions = [list(i.keys())[0].split('.')[0] for i in dict_pdbs_uniprot]
	if filenames[0].split('.')[-1] == 'gz':
		list_pdbs_subset = [pdb_folder+'/'+i+'.pdb.gz' for i in pdb_accessions]
	else:
		list_pdbs_subset = [pdb_folder+'/'+i+'.pdb' for i in pdb_accessions]
	return list_pdbs_subset

list_mapped_af2 = find_af2_for_pdbs(dict_trunc_json, af2_folder)
list_mapped_pdbs = select_mapped_pdbs(dict_trunc_json, pdbs_folder)

def open_pdb(location_pdb):
	'''
	location_pdb: absolute path to the file
	'''
	if location_pdb.split('.')[-1] == 'gz':
		with gzip.open(location_pdb,'r') as fin:
			lines_pdb = fin.readlines()
		ATOM_only = [line.decode('utf-8') for line in lines_pdb if line.startswith(b'ATOM')]
		return ATOM_only
	else:
		with open(location_pdb, 'r') as pdb:
			lines_pdb = pdb.readlines()
		ATOM_only = [line for line in lines_pdb if line.startswith('ATOM')]
		return(ATOM_only)

model_af2 = list_mapped_af2[1]
model_pdb = list_mapped_pdbs[1]
test_pdb = open_pdb(model_pdb)
test_af2 = open_pdb(model_af2)

def find_res_exp_af2(exp_pdb, af2_pdb):
	resnumbs_exp = []
	resnumbs_af2 = []
	resnumb_zero = int(''.join(exp_pdb[0][23:31].split(' ')))
	for residue in exp_pdb:
		if resnumb_zero == 0:
			resnumb_int = int(''.join(residue[23:31].split(' ')))+1 #residue definition
			if resnumb_int not in resnumbs_exp:
				resnumbs_exp.append(resnumb_int)
		else:
			for residue in exp_pdb:
				resnumb_int = int(''.join(residue[23:31].split(' '))) #residue definition
				if resnumb_int not in resnumbs_exp:
					resnumbs_exp.append(resnumb_int)
	for residue in af2_pdb:
		resnumb_int = int(''.join(residue[23:31].split(' ')))
		if resnumb_int not in resnumbs_af2:
			resnumbs_af2.append(resnumb_int) 
	return resnumbs_exp, resnumbs_af2

find_res_exp_af2(test_pdb, test_af2)

def missing_termini_pdb(exp_pdb, af2_pdb):
    '''
    compares residue number of af2 prediction and experimental structure and
    returns 3 lists with missing residues from N- and C- termini, breaks and all exp-missing residues
    '''
    exp_missing_res = []
    termini_missing_res = []
    breaks_res_missing = []
    resnumbs_exp, resnumbs_af2 = find_res_exp_af2(exp_pdb, af2_pdb)

    if len(resnumbs_af2) > len(resnumbs_exp):
        for index in resnumbs_af2:
            if index not in resnumbs_exp:
                exp_missing_res.append(index)
                if int(index) < int(resnumbs_exp[0]) or int(index) > int(resnumbs_exp[-1]): 
                    termini_missing_res.append(index)
                else: 
                    breaks_res_missing.append(index)

    return termini_missing_res, breaks_res_missing, exp_missing_res

t, b, all_miss = missing_termini_pdb(test_pdb, test_af2)

def find_plddt_resolved(exp_pdb, af2_pdb):
    residues_exp, residues_af2 = find_res_exp_af2(exp_pdb, af2_pdb)
    dict_plddts = {}
    for index in af2_pdb:
        res_numb = ''.join(index[23:31].split(' '))
        plddt = ''.join(index[60:67].split(' '))
        if int(res_numb) in residues_exp:
            if res_numb not in dict_plddts.keys():
                dict_plddts[res_numb] = float(plddt)
    return dict_plddts

def find_plddt_missing(af2_pdb, missing):
    dict_plddts_missing = {}
    for index in af2_pdb:
        res_numb = ''.join(index[23:31].split(' '))
        plddt = ''.join(index[60:67].split(' '))
        if int(res_numb) in missing:
            if res_numb not in dict_plddts_missing.keys():
                dict_plddts_missing[res_numb] = float(plddt)
    return dict_plddts_missing

def find_Bfac(exp_pdb):
	dict_bfacs = {}
	resnumb_zero = int(''.join(exp_pdb[0][23:31].split(' ')))
	if resnumb_zero == 0:
		for index in exp_pdb:
			res_numb = ''.join(index[23:31].split(' '))
			bfac = ''.join(index[60:67].split(' '))
			if bfac not in dict_bfacs.keys():
				dict_bfacs[int(res_numb)+1] = float(bfac)
		return dict_bfacs
	else: 
		for index in exp_pdb:
			res_numb = ''.join(index[23:31].split(' '))
			bfac = ''.join(index[60:67].split(' '))
			if bfac not in dict_bfacs.keys():
				dict_bfacs[res_numb] = float(bfac)
		return dict_bfacs

bfacs = find_Bfac(test_pdb)

plddt_resolved = find_plddt_resolved(test_pdb, test_af2)
pldtt_termini = find_plddt_missing(test_af2, t)
pldtt_breaks = find_plddt_missing(test_af2, b)

def plot_plddts(plddt_rd, plddt_ms, plddt_bk):

    merged_plddts = {**plddt_rd, **plddt_ms, **plddt_bk}
    merged_plddts_sorted = sorted(merged_plddts.items(), key=lambda x: int(x[0]))
    sorted_plddts = [i[1] for i in merged_plddts_sorted]

    annot = []
    for i in range(len(sorted_plddts)):
        res = i+1
        if str(res) in plddt_rd.keys():
            annot.append('resolved')
        if str(res) in plddt_ms.keys():
            annot.append('missing_termini')
        if str(res) in plddt_bk.keys(): 
            annot.append('breaks')

    plot_df = pd.DataFrame(sorted_plddts, columns=['plddt'])
    plot_df['annot'] = annot
    order = ['missing_termini', 'breaks', 'resolved']
    sns.barplot(data=plot_df, x='annot', y='plddt', order=order, capsize=.1, linewidth=1.5, 
                palette='crest')
    sns.despine()
    plt.xlabel('')
    plt.ylabel('plddt', fontsize=14, labelpad=10)
    plt.ylim(0,100)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()


def plot_plddt_R_corr(bfacs, af2_plddts):
	xvals = np.array([value for value in list(af2_plddts.values())])
	yvals = np.array([value for value in list(bfacs.values())])
	pearson = stats.pearsonr(xvals, yvals)
	print(pearson)
	plot_df = pd.DataFrame([xvals, yvals], index=['plddt', 'bfac'])
	plot_df_T = plot_df.T
	# sns.kdeplot(data=plot_df_T, x='plddt', y='bfac', fill=True, palette="crest")
	plt.scatter(xvals, yvals, s=10)
	plt.ylabel('bfac', fontsize=14, labelpad=10)
	plt.xlabel('plddt', fontsize=14, labelpad=10)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.ylim(min(yvals)-5, max(yvals)+5)
	plt.xlim(90, 100)
	# plt.text(91, 90, '$p=${}'.format(round(pearson[0],3)))
	# plt.axline(xy1=(100,min(yvals)-5), xy2=(70,max(yvals)+5), linestyle='--', color='black')
	plt.show()

plot_plddt_R_corr(bfacs, plddt_resolved)
plot_plddts(plddt_resolved, pldtt_termini, pldtt_breaks)
