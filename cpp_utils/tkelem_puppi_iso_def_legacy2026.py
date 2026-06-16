"""
The isolation as computed below was verified to be equivalent to the
puppi isolation for egamma as defined in the firmware simulation
before summer 2026. This will be replaced by the modern definition
at a later date.

TODO - Note the meeting number and the other details of the changed
definiton below.
"""

import ROOT


def tkelem_puppi_iso_def_legacy2026():

    STRCPPFUNC_calciso_legacy26_single = """
        std::tuple< float, float> calciso_legacy26_single(float sig_pt,
                                                        float sig_eta,
                                                        float sig_calo_eta,
                                                        float sig_phi,
                                                        float sig_calo_phi,
                                                        float sig_vz,
                                                        ROOT::VecOps::RVec<float> &bkg_pt,
                                                        ROOT::VecOps::RVec<float> &bkg_eta,
                                                        ROOT::VecOps::RVec<float> &bkg_phi,
                                                        ROOT::VecOps::RVec<float> &bkg_pid,
                                                        ROOT::VecOps::RVec<float> &bkg_z0,
                                                        float dRmin, float dRmax, float ptmin,
                                                        float dzmax) {
                                                        
            float isotot = 0.0;
            for(int i=0; i<bkg_pt.size(); i++) {

                if(bkg_pt[i] < ptmin) continue;
                if(fabs(bkg_pid[i]) == 11  || fabs(bkg_pid[i]) == 13  || fabs(bkg_pid[i]) == 211 ) {
                    if(abs(bkg_z0[i]) >= dzmax) continue;
                }
                // else {
                //     if(abs(sig_pt - bkg_pt[i]) < 0.5) continue;
                // }

                float calodR = getdR(sig_calo_eta, bkg_eta[i], sig_calo_phi, bkg_phi[i]);
                // float nonCalodR = getdR(sig_eta, bkg_eta[i], sig_phi, bkg_phi[i]);
                // float dR = (fabs(bkg_pid[i]) == 22) || (fabs(bkg_pid[i]) == 130) ? calodR : nonCalodR;
                float dR = calodR;
                // float deta = sig_calo_eta - bkg_eta[i];
                // float dphi = sig_calo_phi - bkg_phi[i];
                // float dR = sqrt(deta*deta + dphi*dphi);

                if (dR >= dRmin && dR < dRmax) {
                        isotot += bkg_pt[i];
                    }
            }

            return std::make_tuple(isotot, isotot/sig_pt);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_calciso_legacy26_single)

    STRCPPFUNC_calciso_legacy26 = """
        std::tuple< ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float> > calciso_legacy26(ROOT::VecOps::RVec<float> &sig_pt,
                                                                 ROOT::VecOps::RVec<float> &sig_eta,
                                                                 ROOT::VecOps::RVec<float> &sig_calo_eta,
                                                                 ROOT::VecOps::RVec<float> &sig_phi,
                                                                 ROOT::VecOps::RVec<float> &sig_calo_phi,
                                                                 ROOT::VecOps::RVec<float> &sig_vz,
                                                                 ROOT::VecOps::RVec<float> &bkg_pt,
                                                                 ROOT::VecOps::RVec<float> &bkg_eta,
                                                                 ROOT::VecOps::RVec<float> &bkg_phi,
                                                                 ROOT::VecOps::RVec<float> &bkg_pid,
                                                                 ROOT::VecOps::RVec<float> &bkg_z0,
                                                                 float dRmin, float dRmax, float ptmin,
                                                                 float dzmax) {
            ROOT::VecOps::RVec<float> isotot(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> relisotot(sig_pt.size(), -1);

            for (int i = 0; i < sig_pt.size(); i++) {
                std::tuple isotuple = calciso_legacy26_single( sig_pt[i], sig_eta[i], sig_calo_eta[i],
                                                               sig_phi[i], sig_calo_phi[i], sig_vz[i],
                                                               bkg_pt, bkg_eta, bkg_phi, bkg_pid, bkg_z0,
                                                               dRmin, dRmax, ptmin, dzmax );
                isotot[i] = std::get<0>(isotuple);
                relisotot[i] = std::get<1>(isotuple);
            }
            return std::make_tuple(isotot, relisotot);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_calciso_legacy26)

