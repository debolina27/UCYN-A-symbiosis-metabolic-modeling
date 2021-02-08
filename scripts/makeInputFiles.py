'''this code makes GAMS input files from .sbml files'''

import cobra 

def __writelist2file(inlist, fname):
	f1 = open(fname,'w')
	for item in inlist:
		f1.write('\''+item+'\'\n')
	f1.close()

def __writestoic2file(indict, fname):
	f1 = open(fname,'w')
	for item in indict.keys():
		for item2 in indict[item].keys():
			f1.write('\''+item2+'\'.\''+item+'\' '+indict[item][item2]+'\n')
	f1.close()

def __writerxnbnds2file(indict, fname):
	f1 = open(fname,'w')
	for rxn in indict.keys():
		f1.write('LB(\''+rxn+'\') = '+str(indict[rxn]['LB'])+';\n')
		f1.write('UB(\''+rxn+'\') = '+str(indict[rxn]['UB'])+';\n')
	f1.close()

def __gettransfermets(infile):
	outmets = dict()  ###key is no=-comp

	f1 = open(infile,'r')
	data = f1.readlines()
	f1.close()

	for i in range(0,len(data)):
		line = data[i]
		if (line.strip())=='':
			pass
		elif (line.strip()).split()[0] == '<species':
			met = line[line.find(' id="'):line.find(' metaid=')]
			met = met.strip()
			met = met.replace('id=','')
			met = met.replace('"','')
			if 'value="True"' in data[i+3]:
				outmets[met[:met.find('_')]] = met
				if 'cpd00293' in met:
						print (met+'\t'+outmets[met[:met.find('_')]])
	return outmets


def __maketransferrxns(infile1,org1,infile2,org2):
	print (org1)
	mets_org1 = __gettransfermets(infile1)
	print (org2)
	mets_org2 = __gettransfermets(infile2)


	####make rxn list and sij 
	rxnlist_org1 = []
	rxnlist_org2 = []
	rxn2stoic = dict()

	fout = open('./../data/OrgTransRxns_BinaryCons.txt','w')
	fout2 = open('./../data/OrgTransRxns_BinaryConsList.txt','w')
	consnum = 0
	for met in list(mets_org1.keys()):
		if met in (list(mets_org2.keys())):

			rxn1 = 'Trans_'+met+'_'+org1+'2'+org2
			rxn2stoic[rxn1] = {}
			rxn2stoic[rxn1][mets_org1[met]+'['+org1+']'] = '-1'
			rxn2stoic[rxn1][mets_org2[met]+'['+org2+']'] = '1'
			rxnlist_org1.append(rxn1)

			rxn2 = 'Trans_'+met+'_'+org2+'2'+org1
			rxn2stoic[rxn2] = {}
			rxn2stoic[rxn2][mets_org2[met]+'['+org2+']'] = '-1'
			rxn2stoic[rxn2][mets_org1[met]+'['+org1+']'] = '1'
			rxnlist_org2.append(rxn2)

			###now the binary list 
			consnum = consnum + 1
			fout2.write('\nLB1_'+str(consnum)+'\n')
			fout.write('LB1_'+str(consnum)+'..\tv(\''+rxn1+'\') =g= LB(\''+rxn1+'\')*ytrans(\''+rxn1+'\');\n')

			fout2.write('UB1_'+str(consnum)+'\n')
			fout.write('UB1_'+str(consnum)+'..\tv(\''+rxn1+'\') =l= UB(\''+rxn1+'\')*ytrans(\''+rxn1+'\');\n')

			fout2.write('\nLB2_'+str(consnum)+'\n')
			fout.write('LB2_'+str(consnum)+'..\tv(\''+rxn2+'\') =g= LB(\''+rxn2+'\')*ytrans(\''+rxn2+'\');\n')

			fout2.write('UB2_'+str(consnum)+'\n')
			fout.write('UB2_'+str(consnum)+'..\tv(\''+rxn2+'\') =l= UB(\''+rxn2+'\')*ytrans(\''+rxn2+'\');\n')

			fout2.write('OrgTransCons'+str(consnum)+'\n')
			fout.write('OrgTransCons'+str(consnum)+'..\tytrans(\''+rxn1+'\') + ytrans(\''+rxn2+'\') =l= 1;\n')

	fout.close()
	fout2.close()

	return rxnlist_org1, rxnlist_org2, rxn2stoic


def __parsemodel(infile,org):
	'''parses .sbml model file to give rxnlist, metablist, rxnbnds, mets2transfer'''
	metlist = []
	rxnlist = []
	rxn2stoic = dict()
	rxn2bnd = dict()

	model = cobra.io.read_sbml_model(infile)
	
	for met in model.metabolites:
		metlist.append(met.id+'['+org+']')

	for rxn in model.reactions:
		#print (dir(rxn))
		rxnlist.append(rxn.id+'['+org+']')
		rxn2bnd[rxn.id+'['+org+']'] = {}
		rxn2bnd[rxn.id+'['+org+']']['LB'] = rxn.lower_bound
		rxn2bnd[rxn.id+'['+org+']']['UB'] = rxn.upper_bound

		###make the stoichiometry
		if '<=>' in rxn.reaction:
			stoic = str(rxn.reaction).split('<=>')
		elif '-->' in rxn.reaction:
			stoic = str(rxn.reaction).split('-->')
		elif '<--' in rxn.reaction:
			stoic = str(rxn.reaction).split('<--')
		lhs = stoic[0]
		rhs = stoic[-1]

		rxn2stoic[rxn.id+'['+org+']'] = {}
		if lhs.strip()!='':
			lhs = lhs.split('+')
			for met in lhs:
				met = met.strip()
				#check if stoic exists
				if len(met.split())==1:
					rxn2stoic[rxn.id+'['+org+']'][met+'['+org+']'] = '-1'
				else:
					rxn2stoic[rxn.id+'['+org+']'][met.split()[-1]+'['+org+']'] = str(-1.0*float(met.split()[0]))

		if rhs.strip()!='':
			rhs = rhs.split('+')
			for met in rhs:
				met = met.strip()
				#check if stoic exists
				if len(met.split())==1:
					rxn2stoic[rxn.id+'['+org+']'][met+'['+org+']'] = '1'
				else:
					rxn2stoic[rxn.id+'['+org+']'][met.split()[-1]+'['+org+']'] = met.split()[0]

	return metlist, rxnlist, rxn2stoic,rxn2bnd



##UCYN model
mets, rxns, stoic_dict, bnd_dict  = __parsemodel('./../GSMs/iUCYN401.sbml', 'ucyn')
##write stuff to file 
__writelist2file(mets,'./../data/UCYN_metabolites.txt')
__writelist2file(rxns,'./../data/UCYN_reactions.txt')
__writestoic2file(stoic_dict,'./../data/UCYN_sij.txt')
__writerxnbnds2file(bnd_dict, './../data/UCYN_reaction_bounds.txt')

##Host model
mets, rxns, stoic_dict, bnd_dict  = __parsemodel('./../GSMs/iCtobin714.sbml', 'ctobin')
##write stuff to file 
__writelist2file(mets,'./../data/Ctobin_metabolites.txt')
__writelist2file(rxns,'./../data/Ctobin_reactions.txt')
__writestoic2file(stoic_dict,'./../data/Ctobin_sij.txt')
__writerxnbnds2file(bnd_dict, './../data/Ctobin_reaction_bounds.txt')

###make transfer reactions 
rxns_ucyn2host, rxns_host2ucyn, stoic_dict = __maketransferrxns('./../GSMs/iUCYN401.sbml','UCYN','./../GSMs/iCtobin714.sbml','Ctobin')
__writelist2file(rxns_ucyn2host,'./../data/OrgTransRxns_UCYN2Ctobin.txt')
__writelist2file(rxns_host2ucyn, './../data/OrgTransRxns_Ctobin2UCYN.txt')
__writestoic2file(stoic_dict, './../data/OrgTransRxns_sij.txt')


