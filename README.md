# 2hdm
BDT and DNN jupyter notebooks for signal vs background classification for different signal mass ranges of the Pseudoscalar Higgs $A$:
* Low masses: $A = 400$ $GeV$
* Medium masses: $A = 500$, $600$ $GeV$
* High masses: $A = 700+$ $GeV$
---
Samples can be located at the group EOS: `/eos/atlas/atlascerngroupdisk/phys-hdbs/hbsm/ANA-HDBS-2023-05/` and `/eos/user/r/rjoshi/ntuples`
* Signal Information: https://hplus-wminus-to-bbww-analysis-docs.docs.cern.ch/samples/
* Background Information: [in progress]
---
### Features Used
Baseline features:
* "bjet1_pt_NOSYS", "bjet1_phi_NOSYS", "bjet1_eta_NOSYS", "bjet1_mass_NOSYS" 
* "bjet2_pt_NOSYS", "bjet2_phi_NOSYS", "bjet2_eta_NOSYS", "bjet2_mass_NOSYS"> 
* "A_pt_fitted_NOSYS", "A_phi_fitted_NOSYS", "A_eta_fitted_NOSYS", "A_mass_fitted_NOSYS"
* "Hp_pt_fitted_NOSYS", "Hp_phi_fitted_NOSYS", "Hp_eta_fitted_NOSYS", "Hp_mass_fitted_NOSYS"

Extra 14 features for Low Mass Model:
* "MET_NOSYS"
* "deltaR_bjet1_bjet2"
* "Wb_nonTop_mass_fitted_NOSYS"
* "lepton_pt_NOSYS"
* "deltaR_ljet1_ljet2"
* "deltaR_ljet2_ljet3"
* "deltaR_ljet2_lepton"
* "deltaR_ljet4_lepton"
* "lepton_eta_NOSYS"
* "deltaR_ljet3_ljet4"
* "deltaR_bjet2_ljet1"
* "deltaR_bjet1_lepton"
* "deltaR_top_ljet2"
* "deltaR_Wb_nonTop_ljet1"
 
