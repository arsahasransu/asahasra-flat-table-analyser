
for i in DYToLL_M50_PU0  DYToLL_M50_PU200  QCD_Pt120To170  QCD_Pt20To30  QCD_Pt30To50  QCD_Pt50To80  QCD_Pt80To120
do
  scp -r mercury13:/opt/ppd/scratch/asahasra/condor_scratch/Phase2EGPuppiIso_142Xv0_Ntuples/$i ./puppiIsoData/
done
