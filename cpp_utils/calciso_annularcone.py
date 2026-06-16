import ROOT

def calciso_annularcone():
    STRCPPFUNC_calcisoanncone_singleobj = """
        std::tuple< float, float, float,  float, float, float, float, float, float, float,  float, float, float,
            float, float, float, float > calcisoanncone_singleobj(float sig_pt,
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
            float isotot = 0.0, iso11 = 0.0, iso13 = 0.0, iso22 = 0.0, iso130 = 0.0, iso211 = 0.0, isooth = 0.0;
              //  std::cout<<"Calculate Iso ------"<<std::endl;
              
            for(int i=0; i<bkg_pt.size(); i++) {

                if(bkg_pt[i] < ptmin) continue;
                if(fabs(bkg_pid[i]) == 11  || fabs(bkg_pid[i]) == 13  || fabs(bkg_pid[i]) == 211 ) {
                    if(abs(bkg_z0[i]) >= dzmax) continue;
                }

                float calodR = getdR(sig_calo_eta, bkg_eta[i], sig_calo_phi, bkg_phi[i]);
                float nonCalodR = getdR(sig_eta, bkg_eta[i], sig_phi, bkg_phi[i]);
                float dR = (fabs(bkg_pid[i]) == 22) || (fabs(bkg_pid[i]) == 130) ? calodR : nonCalodR;
                // float dR = calodR;

                if (dR > dRmin && dR < dRmax) {

                 //                  std::cout<<"L1PuppiCand: "<<bkg_pt[i]<<"\t"<<
                 // bkg_eta[i]<<"\t"<<bkg_phi[i]<<"\t"<<bkg_pid[i]<<"\t"<<dR<<std::endl;

                        isotot += bkg_pt[i];
                        if( fabs(bkg_pid[i]) == 11 ) iso11 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 13 ) iso13 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 22 ) iso22 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 130 ) iso130 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 211 ) iso211 += bkg_pt[i];
                        else isooth += bkg_pt[i];
                    }
            }

            // std::cout<<sig_calo_eta<<"\t"<<sig_calo_phi<<"\t"<<sig_pt<<"\t"<<isotot<<"\t"<<isotot/sig_pt<<"\t"<<quantise_relisoval(isotot/sig_pt)<<std::endl;
            //if(isotot != 0) {
            //    for(int i=0; i<bkg_pt.size(); i++) {
            //       std::cout<<"L1PuppiCand: "<<bkg_pt[i]<<"\t"<<
            //      bkg_eta[i]<<"\t"<<bkg_phi[i]<<"\t"<<bkg_pid[i]<<"\t"<<bkg_z0[i]<<std::endl;
             //   }
            //}
            return std::make_tuple(isotot, iso11, iso13, iso22, iso130, iso211, isooth,
                                   (iso11+iso13+iso211), (iso22+iso130),
                                   isotot/sig_pt, iso11/sig_pt, iso13/sig_pt,
                                   iso22/sig_pt, iso130/sig_pt, iso211/sig_pt,
                                   (iso11+iso13+iso211)/sig_pt, (iso22+iso130)/sig_pt);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_calcisoanncone_singleobj)

    STRCPPFUNC_calcisoannularcone = """
        std::tuple< ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float> >calcisoannularcone(ROOT::VecOps::RVec<float> &sig_pt,
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
            ROOT::VecOps::RVec<float> iso11(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> iso13(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> iso22(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> iso130(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> iso211(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> isooth(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> isochg(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> isonut(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> relisotot(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> reliso11(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> reliso13(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> reliso22(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> reliso130(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> reliso211(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> relisochg(sig_pt.size(), -1);
            ROOT::VecOps::RVec<float> relisonut(sig_pt.size(), -1);

            for (int i = 0; i < sig_pt.size(); i++) {
                // std::cout<<"Calc iso. for: "<<sig_pt[i]<<"\t"<<sig_calo_eta[i]<<"\t"<<sig_calo_phi[i]<<std::endl;
                std::tuple isotuple = calcisoanncone_singleobj( sig_pt[i], sig_eta[i], sig_calo_eta[i],
                                                                sig_phi[i], sig_calo_phi[i], sig_vz[i],
                                                                bkg_pt, bkg_eta, bkg_phi, bkg_pid, bkg_z0,
                                                                dRmin, dRmax, ptmin, dzmax );
                isotot[i] = std::get<0>(isotuple);
                iso11[i] = std::get<1>(isotuple);
                iso13[i] = std::get<2>(isotuple);
                iso22[i] = std::get<3>(isotuple);
                iso130[i] = std::get<4>(isotuple);
                iso211[i] = std::get<5>(isotuple);
                isooth[i] = std::get<6>(isotuple);
                isochg[i] = std::get<7>(isotuple);
                isonut[i] = std::get<8>(isotuple);
                relisotot[i] = std::get<9>(isotuple);
                reliso11[i] = std::get<10>(isotuple);
                reliso13[i] = std::get<11>(isotuple);
                reliso22[i] = std::get<12>(isotuple);
                reliso130[i] = std::get<13>(isotuple);
                reliso211[i] = std::get<14>(isotuple);
                relisochg[i] = std::get<15>(isotuple);
                relisonut[i] = std::get<16>(isotuple);
            }
            return std::make_tuple(isotot, iso11, iso13, iso22, iso130, iso211, isooth, isochg, isonut,
                                   relisotot, reliso11, reliso13, reliso22, reliso130, reliso211,
                                   relisochg, relisonut);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_calcisoannularcone)
