$INLINECOM /*  */

/*
this code finds the least number of metabolites that must be transferred between the two organisms to produce ALL biomass precursors for UCYN-A

*/


OPTION decimals = 8
       sysout = off
       solprint = on
       reslim = 100000
       iterlim = 10000000
       domlim = 1000
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 0.0
       work = 10000000
       nlp = minos5
       mip = cplex; 
       

$set data_path ./../data/
$set out_path ./../results/

SETS
i     Name of the metabolites in the model
/
$include %data_path%UCYN_metabolites.txt
$include %data_path%Ctobin_metabolites.txt
/

j     Name of the rxns in the two models and inter organism transfer
/
$include %data_path%UCYN_reactions.txt
$include %data_path%Ctobin_reactions.txt

***Exchange and transfer rxns 
$include %data_path%OrgTransRxns_UCYN2Ctobin.txt
$include %data_path%OrgTransRxns_Ctobin2UCYN.txt

***the biomass exchange rxns 
$include %data_path%BiomassMetabolitesCutoff_rxns_simultaneousBios.txt
/

UCYNbioEx(j)  list of UCYNA biomass precursors
/
$include %data_path%BiomassMetabolitesCutoff_rxns_simultaneousBios.txt
/

	
***carbon Ctobin -> UCYN
ucyn_carbon_exchanges(j)
/
$include %data_path%OrgTransRxns_Ctobin2UCYN.txt
/

****nitrogen UCYN -> Ctobin
tobin_nitrogen_exchanges(j)
/
$include %data_path%OrgTransRxns_UCYN2Ctobin.txt
/

org_trans(j)  all inter organism transfer reactions
/
$include %data_path%OrgTransRxns_UCYN2Ctobin.txt
$include %data_path%OrgTransRxns_Ctobin2UCYN.txt
/

;


********************* Define parameters *******************
PARAMETERS
UB(j)  Upper bound on fluxes

LB(j)  Lower bound on fluxes

S(i,j)   Stoichiomteric matrix 
/
$include %data_path%UCYN_sij.txt
$include %data_path%Ctobin_sij.txt

***individual exchange rxns for every precursor, to find blocked ones 
$include %data_path%BiomassMetabolitesCutoff_sij.txt

****the organism transfer 
$include %data_path%OrgTransRxns_sij.txt

***UCYN ATPM
*'cpd00002[c][ucyn]'.'R_ATPM_NGAM[ucyn]' -1
*'cpd00001[c][ucyn]'.'R_ATPM_NGAM[ucyn]' -1
*'cpd00008[c][ucyn]'.'R_ATPM_NGAM[ucyn]' 1
*'cpd00009[c][ucyn]'.'R_ATPM_NGAM[ucyn]' 1
*'cpd00067[c][ucyn]'.'R_ATPM_NGAM[ucyn]' 1
/

carbon_in_ex_rxns(j)
nitrogen_in_ex_rxns(j)

c(j)


epsilon  A small value

n_min  The minimum number of added rxns

counter

done  Represents when we are done with finding the rxns to resol
;

epsilon = 0.001;

********************* Define variables *******************
VARIABLES
	v(j)      Fluxes
	u(j)      positive fluxes
	Z         objective	
	a(j) 		for decomposing the lumped product v(j)y_neg(j)
	fluxsum
	max_org_trans
;

biNARY VARIABLE
	ytrans(j)     1 if a metabolite transfer equation is activated 
;

**********************define bounds on reactions ***********************************
scalar vmax /10000/;

$include %data_path%Ctobin_reaction_bounds.txt
$include %data_path%UCYN_reaction_bounds.txt


LB(j)$(ucynBioEx(j)) = 0.0;
UB(j)$(ucynBioEx(j)) = 1000.0;

UB(j)$(org_trans(j)) = 1000.0;
LB(j)$(org_trans(j)) = 0.0;


***make water be in excess
LB('Exchange_cpd00001_ctobin[ctobin]') = -100000.0;
LB('Exchange_cpd00001_ucyn[ucyn]') = -100000.0;


*****define the C and N composition of exchanges 
loop(j,
carbon_in_ex_rxns(j) = 0;
);
$include %data_path%AddedAllExRxns2_NumCarbon_OrgTransRxns.txt

loop(j,
nitrogen_in_ex_rxns(j) = 0;
);
$include %data_path%AddedAllExRxns2_NumNitrogen_OrgTransRxns.txt

***CO2 limited growth 
LB('Exchange_cpd00011_ctobin[ctobin]') = -100.0;

$ontext
**N2 supply
LB('Exchange_cpd00528_ucyn[ucyn]') = -10000;
***photons
LB('Exchange_hvphoton1_ucyn[ucyn]') = -100.0;
LB('Exchange_hvphoton1_ctobin[ctobin]') = -100.0;
LB('Exchange_hvphoton2_ctobin[ctobin]') = -100.0;
****add minimal oxygen to test o2-scavenging systems
LB('Exchange_cpd00007_ucyn[ucyn]') = -0.1;
$offtext

v.lo(j) = LB(j);
v.up(j) = UB(j);



********************* Define equations *******************
EQUATIONS
	obj
	massbalance_all(i)
	symb_C_ex
	symb_N_ex

	ucyn_bio_cutoff
	tobin_bio_cutoff
	n2fix_cutoff
	testcutoff
	pfba_ub(j)
	pfba_lb(j)
	pfba_cons
	obj_pfba
	max_org_trans_cons
	
$include %data_path%OrgTransRxns_BinaryConsList.txt
	
;


* Minimize the number of reactions added from the database
obj ..              z =e= sum(j$(org_trans(j)), ytrans(j));

max_org_trans_cons..		sum(j$(org_trans(j)), ytrans(j)) =l= max_org_trans ;

pfba_ub(j)..		u(j) =g= v(j);
pfba_lb(j)..		u(j) =g= -v(j);
pfba_cons..			sum(j, u(j)) =l= fluxsum;
obj_pfba..			z =e= sum(j, u(j));

* Mass balance constraints on metabolites
massbalance_all(i)..					   sum(j,S(i,j)*v(j)) =e= 0;

*****the C N exchanges 
symb_C_ex..         sum(j$(ucyn_carbon_exchanges(j)), carbon_in_ex_rxns(j)*v(j)) =l= -0.17*v('Exchange_cpd00011_ctobin[ctobin]');
symb_N_ex..         sum(j$(tobin_nitrogen_exchanges(j)), nitrogen_in_ex_rxns(j)*v(j)) =g= 0.95*2*v('rxn06874[ucyn]');
$include %data_path%OrgTransRxns_BinaryCons.txt



********** Define and solve the model. Store the results in output files  *************

Model minOrgTransBins
/
obj
massbalance_all

symb_C_ex
symb_N_ex

$include %data_path%OrgTransRxns_BinaryConsList.txt
/;
minOrgTransBins.optfile=1;

Model pfba
/
obj_pfba
pfba_ub
pfba_lb
massbalance_all

symb_C_ex
symb_N_ex

$include %data_path%OrgTransRxns_BinaryConsList.txt
max_org_trans_cons
/;
pfba.optfile=1;

FILE resultfile /"%out_path%FindHost2UCYNtransfer_UCYNsimultaneousBioProd.txt"/;

resultfile.pw=500;
PUT resultfile;

****************

alias(j,j1);

***set a 0.001 LB for all UCYN-A biomass precursors
loop(j1$(UCYNbioEx(j1)),
v.lo(j1) = epsilon;
v.up(j1) = 1000.0;
);

put //;

*****uptake of oxygen 
v.lo('Exchange_cpd00007_ucyn[ucyn]') = -10.0;
v.up('Exchange_cpd00007_ucyn[ucyn]') = 1000.0;

put "Minimizing inter-organism metabolite transfer"//;

Solve minOrgTransBins using mip minimizing z;
max_org_trans.fx = z.l;
put "total number of metabolites transferred: ",z.l:0:10,"  ",minOrgTransBins.modelstat//;
*loop(j$(org_trans(j) and (ytrans.l(j) eq 1.0)),
*put j.tl:0:20,"  ",v.l(j):10:5," ",y.l(j):10:5/;
*);

put "Finding pFBA flux distribution at minimal metabolite transfer"//;

Solve pfba using  mip minimizing z;
put "pFBA objective: ", z.l:0:10,"   ",pfba.modelstat/;
*loop(j$(org_trans(j) and (y.l(j) eq 1.0)),
*put j.tl:0:20,"  ",v.l(j):10:5," ",y.l(j):10:5/;);

loop(j$(v.l(j) ne 0.0),
put j.tl:0:20,"   ",v.l(j):0:10,"    ",v.lo(j):0:10,"   ",v.up(j):0:10/);



