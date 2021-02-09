# UCYN-A symbiosis metabolic modeling

## This is a public repo for 'Elucidation of trophic interactions in an unusual single-cell nitrogen-fixing symbiosis using metabolic modeling'

### Getting Started

- GSMs: stores genome-scale metabolic models for UCYN-A and the host (SBML and JSON formats)
- data: stores input data for GAMS scripts
- results: stores GAMS outputs
- scripts: 
  - makeInputFiles.py: makes GAMS-readable input files from the 
  - MinimalMetTransfers_UCYNbiomass.gms: finds the least number of metabolites that must be transferred between the two organisms to produce all UCYN-A biomass precursors
  - MinimalMetTransfers_UCYNbiomass_Icuts.gms: uses integer cuts to iteratively determine alternate inter-organism metabolite transfer sets to produce all UCYN-A biomass precursors
  - MinimalMetTransfers_maxN2fix.gms: finds the least number of metabolites that must be transferred between the two organisms for maximal UCYN-A nitrogen fixation flux
  - MinimalMetTransfers_maxHostBiomass.gms: finds the least number of metabolites that must be transferred between the two organisms for maximal host biomass production flux
