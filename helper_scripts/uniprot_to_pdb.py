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
import re

pdbs_folder = '/home/viliugvd/Desktop/project/pdb_dataset/pdbs'
af2_folder = '/home/viliugvd/Desktop/project/pdb_dataset/af2'
pdbs_selected_folder = '/Volumes/Seagate/HITS_project/pdbs_3.5_monomeric/pdbs_selected'
pdbs_mapped = '/home/viliugvd/Desktop/project/pdb_dataset/pdbs_mapped'

def create_pdb_list(pdb_folder):
	pdb_files = sorted(os.listdir(pdb_folder))
	filenames = sorted([i for i in pdb_files if i[0]!='.'])
	pdb_names = [i.split('.')[0] for i in pdb_files if i[0]!='.']
	rcsb_accession_codes = []

	for index in range(len(filenames)):
		print(index)
		with gzip.open(pdb_folder+'/'+filenames[index],'r') as fin:
			lines_pdb = fin.readlines()
			ATOM_only = [line for line in lines_pdb if line.startswith(b'ATOM')]
			chain_name = ATOM_only[0][21:22].decode('utf-8')
			rcsb_accession_code = pdb_names[index]+'.'+chain_name
			rcsb_accession_codes.append(rcsb_accession_code)
	return pdb_names, rcsb_accession_codes

# pdb_names, rcsb_codes = create_pdb_list(pdbs_folder)

def map_pdb_uniprot(rcsb_codes):
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
	
	results_mapping = []
	for query_id in rcsb_codes:
		print(query_id)
		result = get_pdb_uniprot(query_id)
		results_mapping.append(result)
	return results_mapping
	
def write_pdb_uniprot_mapping(rcsb_codes):
	pdb_uniprot_mapping = map_pdb_uniprot(rcsb_codes)
	with open('/home/viliugvd/Desktop/project/pdb_dataset/01_10k_CATH_3.5_mapping.txt', 'w') as res:
  		json.dump(pdb_uniprot_mapping, res)

with open('/home/viliugvd/Desktop/project/pdb_dataset/01_10k_CATH_3.5_mapping.txt', 'r') as json_data:
	test = json.load(json_data)

selected_100 = test[:100]

def find_mapped_pdbs(json_mapped):
	'''
	the problem here is that some PDB files have chains beginning not with "A" e.g. 
	first there's a B chain (ATOMS 0-n) and then A chain with (ATOMS n-...k)
	'''
	exclude_indices = []
	for index, item in enumerate(selected_100):
		value = list(item.values())[0]
		if value == None:
			exclude_indices.append(index)
	upd_list_filenames = [json_mapped[i] for i in range(len(json_mapped)) if i not in exclude_indices] 
	return upd_list_filenames
mapped_100_pdbs = find_mapped_pdbs(selected_100)
mapped_all = find_mapped_pdbs(test)

def download_af_db(mapped_pdbs_dict):
	uniprot_codes = []
	for item in mapped_pdbs_dict:
		value = list(item.values())[0] 
		if value not in uniprot_codes:
			uniprot_codes.append(value) 
	folder_pdbs_out = '/home/viliugvd/Desktop/project/pdb_dataset/af2'
	missing_predictions = []
	for entry in uniprot_codes:
		version='v4'
		try:
			entry_code = entry.split('\n')[0]
		except AttributeError:
			print('NoneType error! Perhaps the UniProt accession is missing')
		alphafold_ID = 'AF-{}-F1-model_{}'.format(str(entry_code), str(version))
		model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_ID}.pdb'
		pdb_out = os.path.join(folder_pdbs_out, alphafold_ID+'.pdb')
		os.system(f'curl --max-time 900 {model_url} -o {alphafold_ID}'+'.pdb')
		shutil.move('/home/viliugvd'+'/'+str(alphafold_ID)+'.pdb', pdb_out)
# download_af_db(mapped_all)

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

# test_multiple_pdbs_one_uniprot = mapped_100_pdbs[:5]

def gen_uniprot_for_pdbs(json_data, missing_af2_pred, remove_dub=False):
	'''
	eliminates dublicates in pdb -> uniprot mapping; only individual pdb structures are compared
	this is bad bc # of atoms is different in different conformations -> i want to compare everything
	i find to the ONE af2 prediction -> think of a way of indexing certain regions
	'''
	uniprot_codes = [list(i.values())[0] for i in json_data]
	upd_uniprot_names = []
	for index in range(len(json_data)):
		for key, value in json_data[index].items():
			if value not in upd_uniprot_names:
				upd_uniprot_names.append(value)
	entries_present = []
	trunc_json = []
	if remove_dub == True:
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
	else: 
		for i in json_data:
			if list(i.values())[0] in upd_uniprot_names:
				if list(i.values())[0] not in missing_pred:
					if list(i.values())[0] != None:
						entries_present.append(list(i.values())[0])
						trunc_json.append(i)
		return trunc_json
	
mapped_all_upd = gen_uniprot_for_pdbs(mapped_all, missing_pred)
mapped_all_upd_nodub = gen_uniprot_for_pdbs(mapped_all, missing_pred, remove_dub = True)

def find_af2_for_pdbs(dict_pdbs_uniprot, af2_folder):
	list_af2_locs = []
	uniprot_pdb_accession = [list(i.values())[0] for i in dict_pdbs_uniprot]
	filenames = [i for i in os.listdir(af2_folder) if i[0]!='.']
	# print(uniprot_pdb_accession)
	for filename in filenames:
		list_af2_locs.append(af2_folder+'/'+filename)
			# if filename.split('-')[1] in uniprot_pdb_accession: #eliminates dublicates!
			# 	list_af2_locs.append(af2_folder+'/'+filename)
	return list_af2_locs

def select_mapped_pdbs(dict_pdbs_uniprot, pdb_folder):
	filenames = [i for i in os.listdir(pdb_folder) if i[0]!='.']
	pdb_accessions = [list(i.keys())[0].split('.')[0] for i in dict_pdbs_uniprot]
	if filenames[0].split('.')[-1] == 'gz':
		list_pdbs_subset = [pdb_folder+'/'+i+'.pdb.gz' for i in pdb_accessions]
	else:
		list_pdbs_subset = [pdb_folder+'/'+i+'.pdb' for i in pdb_accessions]
	return list_pdbs_subset

list_mapped_af2 = find_af2_for_pdbs(mapped_all_upd, af2_folder)
list_mapped_pdbs = select_mapped_pdbs(mapped_all_upd, pdbs_folder)

list_mapped_af2_nodub = find_af2_for_pdbs(mapped_all_upd_nodub, af2_folder)
list_mapped_pdbs_nodub = select_mapped_pdbs(mapped_all_upd_nodub, pdbs_folder)

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

def find_res_exp_af2(exp_pdb, af2_pdb):
	resnumbs_exp = []
	resnumbs_af2 = []
	try:
		resnumb_zero = int(''.join(exp_pdb[0][23:31].split(' ')))
		for residue in exp_pdb:
			if resnumb_zero == 0:
				try:	
					resnumb_int = int(''.join(residue[22:26].split(' ')))+1 #residue definition
				except ValueError:
					print('value error encountered, perhaps residue numb contains a letter')
					break
				if resnumb_int not in resnumbs_exp:
					resnumbs_exp.append(resnumb_int)
			else:
				try:
					resnumb_int = int(''.join(residue[22:26].split(' ')))
				except ValueError:
					print('value error encountered, perhaps residue numb contains a letter')
					break
				if resnumb_int not in resnumbs_exp:
					resnumbs_exp.append(resnumb_int)
	except ValueError:
		print('value error encountered, perhaps residue numb contains a letter')
	if len(resnumbs_exp) != 0:
		for residue in af2_pdb:
			resnumb_int = int(''.join(residue[22:26].split(' ')))
			if resnumb_int not in resnumbs_af2:
				resnumbs_af2.append(resnumb_int) 
	return resnumbs_exp, resnumbs_af2, len(resnumbs_exp), len(resnumbs_af2)

def missing_termini_pdb(exp_pdb, af2_pdb):
	'''
	compares residue number of af2 prediction and experimental structure and
	returns 3 lists with missing residues from N- and C- termini, breaks and all exp-missing residues
	'''
	exp_missing_res = []
	termini_missing_res = []
	breaks_res_missing = []
	resnumbs_exp, resnumbs_af2, length_exp, length_af2 = find_res_exp_af2(exp_pdb, af2_pdb)

	if len(resnumbs_af2) > len(resnumbs_exp):
		for index in resnumbs_af2:
			if index not in resnumbs_exp:
				exp_missing_res.append(index)
				if int(index) < int(resnumbs_exp[0]) or int(index) > int(resnumbs_exp[-1]): 
					termini_missing_res.append(index)
				else: 
					breaks_res_missing.append(index)

	return termini_missing_res, breaks_res_missing, exp_missing_res

# t, b, all_miss = missing_termini_pdb(test_pdb, test_af2)

def find_plddt_resolved(exp_pdb, af2_pdb):
	residues_exp, residues_af2, length_exp, length_af2 = find_res_exp_af2(exp_pdb, af2_pdb)
	dict_plddts = {}
	for index in af2_pdb:
		res_numb = ''.join(index[22:26].split(' '))
		plddt = ''.join(index[60:67].split(' '))
		if int(res_numb) in residues_exp:
			if res_numb not in dict_plddts.keys():
				dict_plddts[res_numb] = float(plddt)
	return dict_plddts

def find_plddt_missing(af2_pdb, missing):
	dict_plddts_missing = {}
	for index in af2_pdb:
		res_numb = ''.join(index[22:26].split(' '))
		plddt = ''.join(index[60:67].split(' '))
		if int(res_numb) in missing:
			if res_numb not in dict_plddts_missing.keys():
				dict_plddts_missing[res_numb] = float(plddt)
	return dict_plddts_missing

def find_Bfac(exp_pdb):
	dict_bfacs = {}
	try:
		resnumb_zero = int(''.join(exp_pdb[0][22:26].split(' ')))
		if resnumb_zero == 0:
			for index in exp_pdb:
				try:
					resnumb_int = int(''.join(index[22:26].split(' ')))
					bfac = ''.join(index[60:67].split(' '))
					if bfac not in dict_bfacs.keys():
						dict_bfacs[int(resnumb_int)+1] = float(bfac)
				except ValueError:
					print('b-factor is missing or smth is wrong with it')
					break
			return dict_bfacs
	except ValueError: 
		print('could not convert residue number to integer..')
	else: 
		for index in exp_pdb:
			try:
				res_numb = ''.join(index[22:26].split(' '))
				bfac = ''.join(index[60:67].split(' '))
				if bfac not in dict_bfacs.keys():
					dict_bfacs[res_numb] = float(bfac)
			except ValueError: 
				print('b-factor is missing or something is wrong with it')
		return dict_bfacs

test_file = list_mapped_pdbs[0]
test_pdb = open_pdb(test_file)
bfacs = find_Bfac(test_pdb)

# plddt_resolved = find_plddt_resolved(test_pdb, test_af2)
# pldtt_termini = find_plddt_missing(test_af2, t)
# pldtt_breaks = find_plddt_missing(test_af2, b)

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
	plt.scatter(xvals, yvals, s=30, color='#3E7D83')
	plt.ylabel('bfac', fontsize=14, labelpad=10)
	plt.xlabel('plddt', fontsize=14, labelpad=10)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.ylim(min(yvals)-5, max(yvals)+5)
	plt.xlim(60, 100)
	plt.title('P62149 / 1AHR.A', fontsize=14, pad=10)
	# plt.text(91, 90, '$p=${}'.format(round(pearson[0],3)))
	# plt.axline(xy1=(100,min(yvals)-5), xy2=(70,max(yvals)+5), linestyle='--', color='black')
	plt.show()

# plot_plddt_R_corr(bfacs, plddt_resolved)
# plot_plddts(plddt_resolved, pldtt_termini, pldtt_breaks)

# problems: 
# indexing: not 0 or 1-based; maybe I should discard them from analyses and from subsequent training then?
# AF2 predicts the whole sequence which doesn't make too much sense: too big difference bw/ exp and predictions -> discard them as well? This will ruin the analysis for plddt and bfac
# solution: write an algo that accounts for the difference in numb of missing residues: if too high -> compare only resolved stretch to the AF2 prediction!
# problem: if you throw away the dublicates, you cannot cover the whole exp solved protein sometimes -> you need more structures and only their focused regions for AF2 comparison 
# observation:
# for smaller proteins where differences bw/ exp and AF2 not too large, plddt are significantly lower for missing termini and breaks; for bfac sometimes the correlation is not bad either
# example: model 1A29 pdb (uniprot P62157) or 1AHR (P62149) dict_trunc_json[7] dict_trunc_json[33]

def compute_plddt_bfac_batch(dict_mapped, list_mapped_pdbs, list_mapped_af2, len_filter = True):
	list_lens_models = []
	dict_all_plddts = {}
	dict_all_bfacs = {}
	count=0
	for index in range(len(dict_mapped)):
		pdb_name = list(dict_mapped[index].keys())[0].split('.')[0]
		af2_name = list(dict_mapped[index].values())[0]
		search_name_pdb = re.compile('.*'+str(pdb_name))
		search_name_af2 = re.compile('.*'+str(af2_name))
		file_pdb = list(filter(search_name_pdb.match, list_mapped_pdbs))[0]
		file_af2 = list(filter(search_name_af2.match, list_mapped_af2))[0]
		model_pdb = open_pdb(file_pdb)
		model_af2 = open_pdb(file_af2)
		try:
			res_numb_exp, res_numb_af2, len_exp, len_af2 = find_res_exp_af2(model_pdb, model_af2)
		except ValueError:
			print('encountered ValueError at', index)
			print(file_pdb)
			print(file_af2)
		if len_filter == True: 
			if (len_af2 - len_exp) >= 100:
				count+=1
				continue
			else:
				list_lens_models.append([len_exp, len_af2])
				print('exp:',len_exp,'af2:',len_af2) #discard everything which exceeds the length of diff is higher than threshold
				t, b, all_miss = missing_termini_pdb(model_pdb, model_af2)
				plddt_resolved = find_plddt_resolved(model_pdb, model_af2)
				pldtt_termini = find_plddt_missing(model_af2, t)
				pldtt_breaks = find_plddt_missing(model_af2, b)
				bfacs = find_Bfac(model_pdb)
				dict_all_plddts[index] = [plddt_resolved, pldtt_termini, pldtt_breaks]
				dict_all_bfacs[index] = bfacs
		else: 
			list_lens_models.append([len_exp, len_af2])
			print('exp:',len_exp,'af2:',len_af2) #discard everything which exceeds the length of diff is higher than threshold
			t, b, all_miss = missing_termini_pdb(model_pdb, model_af2)
			plddt_resolved = find_plddt_resolved(model_pdb, model_af2)
			pldtt_termini = find_plddt_missing(model_af2, t)
			pldtt_breaks = find_plddt_missing(model_af2, b)
			bfacs = find_Bfac(model_pdb)
			dict_all_plddts[index] = [plddt_resolved, pldtt_termini, pldtt_breaks]
			dict_all_bfacs[index] = bfacs
	return dict_all_plddts, dict_all_bfacs
dict_all_plddts, dict_all_bfacs = compute_plddt_bfac_batch(mapped_all_upd_nodub, list_mapped_pdbs_nodub, list_mapped_af2_nodub,len_filter=True)

def plot_plddt_batch(dict_all_plddts):

	keys = list(dict_all_plddts.keys())
	plddt_resolved_all = []
	plddt_termini_all = []
	plddt_breaks_all = []

	for i in keys:
		for j in list(dict_all_plddts[i][0].values()):
			plddt_resolved_all.append(j)
		for k in list(dict_all_plddts[i][1].values()):
			plddt_termini_all.append(k)
		for v in list(dict_all_plddts[i][2].values()):
			plddt_breaks_all.append(v)

	plddts = plddt_resolved_all+plddt_termini_all+plddt_breaks_all
	annot_res = ['resolved' for i in plddt_resolved_all]
	annot_termini = ['missing_termini' for i in plddt_termini_all]
	annot_breaks = ['breaks' for i in plddt_breaks_all]
	annots = annot_res+annot_termini+annot_breaks

	df_plot = pd.DataFrame([plddts, annots], index=['plddt', 'annot'])
	df_plot_T = df_plot.T
	n_resolved = len(df_plot_T[df_plot_T['annot'] == 'resolved'])
	n_breaks = len(df_plot_T[df_plot_T['annot'] == 'breaks'])
	n_missing_term = len(df_plot_T[df_plot_T['annot'] == 'missing_termini'])
	print('num missing_termini:', n_missing_term)
	print('num breaks:', n_breaks)
	print('num resolved:', n_resolved)
	order = ['missing_termini', 'breaks', 'resolved']
	sns.barplot(data=df_plot_T, x='annot', y='plddt', order=order, capsize=.1, linewidth=1.5, 
					palette='crest')
	sns.despine()
	plt.ylim(0,100)
	plt.xlabel('')
	plt.ylabel('plddt', fontsize=14, labelpad=10)
	plt.ylim(0,100)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	# plt.title('With dublicates no filtering (n=8708)', fontsize=14, pad=10)
	# plt.savefig('/home/viliugvd/Desktop/project/pdb_dataset/plddt_batch_CATH_3.5_monomeric/plddt_batch_10k_dub_nofiltering.png', 
	# 		bbox_inches = 'tight', dpi = 300, transparent=True)
	plt.show()
	return df_plot_T
df_plot = plot_plddt_batch(dict_all_plddts)

# with open('/home/viliugvd/Desktop/project/pdb_dataset/plddt_batch_CATH_3.5_monomeric/plddts_dub_diff_100_01_10k.json', 'w') as file:
# 	json.dump(dict_all_plddts, file)

print(len(dict_all_plddts))

df_plot[df_plot['annot'] == 'resolved']
np.mean(df_plot[df_plot['annot'] == 'breaks']['plddt'])
np.mean(df_plot[df_plot['annot'] == 'missing_termini']['plddt'])

#### testing here to plot bfacs to check whether i filtered everything out

def bfacs_kdeplot(dict_bfacs):
	pal = sns.color_palette('crest')
	palette = pal.as_hex()
	all_bfacs = [j for i in dict_all_bfacs.keys() for j in dict_all_bfacs[i].values()]
	sns.kdeplot(data=all_bfacs, fill=True, color = palette[1], linewidth=0, alpha=0.6)
	plt.xlabel('B-factor $[\AA]$', fontsize = 14)
	plt.ylabel('Density', fontsize = 14)
	plt.xticks(fontsize=12)
	plt.yticks(fontsize=12)
	plt.savefig('/home/viliugvd/Desktop/project/pdb_dataset/plddt_batch_CATH_3.5_monomeric/bfacs_nodub_filtering_01_10k.png', 
			bbox_inches = 'tight', dpi = 300, transparent=True)
	plt.show()
bfacs_kdeplot(dict_all_bfacs)