import ROOT


def define_cpp_utils():
    STRCPPFUNC_getdR = """
        float getdR(float eta1, float eta2, float phi1, float phi2) {
            ROOT::Math::PtEtaPhiMVector p1(1.0, eta1, phi1, 1.0);
            ROOT::Math::PtEtaPhiMVector p2(1.0, eta2, phi2, 1.0);
            return ROOT::Math::VectorUtil::DeltaR(p1, p2);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_getdR)

    STRCPPFUNC_getminangs = """
        std::unordered_map< std::string,
            ROOT::VecOps::RVec<double> > getminangs(ROOT::VecOps::RVec<float> &geta, // the output is of the size of geta
                                                    ROOT::VecOps::RVec<float> &gphi,
                                                    ROOT::VecOps::RVec<float> &eta,
                                                    ROOT::VecOps::RVec<float> &phi) {

            ROOT::VecOps::RVec<double> mindeta;
            ROOT::VecOps::RVec<double> mindphi;
            ROOT::VecOps::RVec<double> mindR;

            for (int i = 0; i < geta.size(); i++) {
                double min_dR = 99999;
                double min_deta = 99999;
                double min_dphi = 99999;
                if(eta.size() == 0) continue;

                ROOT::Math::PtEtaPhiMVector pi(1.0, geta[i], gphi[i], 1.0);
                ROOT::Math::PtEtaPhiMVector pj0(1.0, eta[0], phi[0], 1.0);
                min_dR = ROOT::Math::VectorUtil::DeltaR(pi, pj0);
                min_deta = geta[i] - eta[0];
                min_dphi = ROOT::Math::VectorUtil::DeltaPhi(pi, pj0);

                for (int j = 1; j < eta.size(); j++) {
                    ROOT::Math::PtEtaPhiMVector pj(1.0, eta[j], phi[j], 1.0);
                    double deta = geta[i] - eta[j];
                    double dphi = ROOT::Math::VectorUtil::DeltaPhi(pi, pj);
                    double dR = ROOT::Math::VectorUtil::DeltaR(pi, pj);
                    if (dR < min_dR) min_dR = dR;
                    if (fabs(deta) < fabs(min_deta)) min_deta = deta;
                    if (fabs(dphi) < fabs(min_dphi)) min_dphi = dphi;
                }
                mindeta.push_back(min_deta);
                mindphi.push_back(min_dphi);
                mindR.push_back(min_dR);
            }

            std::unordered_map< std::string, ROOT::VecOps::RVec<double> > angs;
            angs["deta"] = mindeta;
            angs["dphi"] = mindphi;
            angs["dR"] = mindR;
            return angs;
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_getminangs)

    STRCPPFUNC_getmatchedidxs = """
        std::tuple< ROOT::VecOps::RVec<int>,
            ROOT::VecOps::RVec<int> > getmatchedidxs(ROOT::VecOps::RVec<float> &geta,
                                                    ROOT::VecOps::RVec<float> &gphi,
                                                    ROOT::VecOps::RVec<float> &eta,
                                                    ROOT::VecOps::RVec<float> &phi,
                                                    float dRcut) {
            ROOT::VecOps::RVec<int> idxs(geta.size(), -1); // matched reco index. The size of gen collection
            ROOT::VecOps::RVec<int> gidxs(eta.size(), -1); // matched gen index. The size of reco collection

            for (int i = 0; i < geta.size(); i++) {
                double min_dR = 99999;
                for (int j = 0; j < eta.size(); j++) {
                    double dR = getdR(geta[i], eta[j], gphi[i], phi[j]);
                    if (dR < min_dR && dR < dRcut) {
                        min_dR = dR;
                        gidxs[j] = i;
                        idxs[i] = j;
                    }
                }
            }

            return std::make_tuple(idxs, gidxs);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_getmatchedidxs)

    STRCPPFUNC_calcisoannularcone = """
        ROOT::VecOps::RVec<float> calcisoannularcone(ROOT::VecOps::RVec<float> &sig_pt,
                                                     ROOT::VecOps::RVec<float> &sig_eta,
                                                     ROOT::VecOps::RVec<float> &sig_phi,
                                                     ROOT::VecOps::RVec<float> &bkg_pt,
                                                     ROOT::VecOps::RVec<float> &bkg_eta,
                                                     ROOT::VecOps::RVec<float> &bkg_phi,
                                                     float dRmin, float dRmax) {
            ROOT::VecOps::RVec<float> iso(sig_pt.size(), -1);

            for (int i = 0; i < sig_pt.size(); i++) {
                float sum = 0;
                for (int j = 0; j < bkg_pt.size(); j++) {
                    float dR = getdR(sig_eta[i], bkg_eta[j], sig_phi[i], bkg_phi[j]);
                    if (dR > dRmin && dR < dRmax) {
                        sum += bkg_pt[j];
                    }
                }
                iso[i] = sum / sig_pt[i];
            }
            return iso;
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_calcisoannularcone)

    STRCPPFUNC_calcisoanncone_singleobj = """
        std::tuple< float, float, float,
            float, float, float, float > calcisoanncone_singleobj(float sig_pt,
                                                                  float sig_eta,
                                                                  float sig_calo_eta,
                                                                  float sig_phi,
                                                                  float sig_calo_phi,
                                                                  ROOT::VecOps::RVec<float> &bkg_pt,
                                                                  ROOT::VecOps::RVec<float> &bkg_eta,
                                                                  ROOT::VecOps::RVec<float> &bkg_phi,
                                                                  ROOT::VecOps::RVec<float> &bkg_pid,
                                                                  float dRmin, float dRmax) {
            float isotot = 0.0, iso11 = 0.0, iso13 = 0.0, iso22 = 0.0, iso130 = 0.0, iso211 = 0.0, isooth = 0.0;
            for(int i=0; i<bkg_pt.size(); i++) {
                float calodR = getdR(sig_calo_eta, bkg_eta[i], sig_calo_phi, bkg_phi[i]);
                float nonCalodR = getdR(sig_eta, bkg_eta[i], sig_phi, bkg_phi[i]);
                float dR = (fabs(bkg_pid[i]) == 22) || (fabs(bkg_pid[i]) == 130) ? calodR : nonCalodR;
                if (dR > dRmin && dR < dRmax) {
                        isotot += bkg_pt[i];
                        if( fabs(bkg_pid[i]) == 11 ) iso11 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 13 ) iso13 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 22 ) iso22 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 130 ) iso130 += bkg_pt[i];
                        else if( fabs(bkg_pid[i]) == 211 ) iso211 += bkg_pt[i];
                        else isooth += bkg_pt[i];
                    }
            }
            return std::make_tuple(isotot/sig_pt, iso11/sig_pt, iso13/sig_pt,
                                   iso22/sig_pt, iso130/sig_pt, iso211/sig_pt,
                                   isooth/sig_pt);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_calcisoanncone_singleobj)

    STRCPPFUNC_checksorting = """
        template <typename T>
        unsigned int checksorting(ROOT::VecOps::RVec<T> &collection, bool descending=true) {
            bool is_sorted = true;
            unsigned int collection_size = collection.size();
            if(collection_size <= 1) return 2;
            for (unsigned int i=1; i<collection_size; i++) {
                bool check_sort_order = descending ? collection[i] <= collection[i-1] :
                                                     collection[i] >= collection[i-1];
                if(!check_sort_order) {
                    is_sorted = false;
                    break;
                }
            }
            return is_sorted ? 1 : 0;
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_checksorting)

    STRCPPFUNC_getmask_annulardR = """
        ROOT::VecOps::RVec<bool> getpuppimask_annulardR(ROOT::VecOps::RVec<float> &alleta,
                                                        ROOT::VecOps::RVec<float> &allphi,
                                                        ROOT::VecOps::RVec<float> &allpid,
                                                        ROOT::VecOps::RVec<float> &refeta,
                                                        ROOT::VecOps::RVec<float> &refcaloeta,
                                                        ROOT::VecOps::RVec<float> &refphi,
                                                        ROOT::VecOps::RVec<float> &refcalophi,
                                                        float dRmin=0.1, float dRmax=0.4) {
            ROOT::VecOps::RVec<bool> mask(alleta.size(), false);

            for(unsigned int mcnt=0; mcnt<alleta.size(); mcnt++) {
                for(unsigned int rcnt=0; rcnt<refeta.size(); rcnt++) {
                    float nonCalodR = getdR(alleta[mcnt], refeta[rcnt], allphi[mcnt], refphi[rcnt]);
                    float calodR = getdR(alleta[mcnt], refcaloeta[rcnt], allphi[mcnt], refcalophi[rcnt]);
                    float dR = (fabs(allpid[mcnt]) == 22) || (fabs(allpid[mcnt]) == 130) ? calodR : nonCalodR;
                    if(dR > dRmin && dR < dRmax) {
                        mask[mcnt] = true;
                        break;
                    }
                }
            }
            
            return mask;
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_getmask_annulardR)

    STRCPPFUNC_getmindR_PuppiCandSize_TkElRef = """
        ROOT::VecOps::RVec<double> getmindR_PuppiCandSize_TkElRef(ROOT::VecOps::RVec<float> &peta, // the output is of the size of peta
                                                                  ROOT::VecOps::RVec<float> &pphi,
                                                                  ROOT::VecOps::RVec<float> &ppdg,
                                                                  ROOT::VecOps::RVec<float> &eta,
                                                                  ROOT::VecOps::RVec<float> &phi,
                                                                  ROOT::VecOps::RVec<float> &caloeta,
                                                                  ROOT::VecOps::RVec<float> &calophi) {

            ROOT::VecOps::RVec<double> mindr_vals(peta.size(), -0.01);

            // std::cout<<"Start of a new event"<<std::endl;
            // for(int i=0; i<peta.size(); i++) {
            //     std::cout<<ppdg[i]<<"\t"<<peta[i]<<"\t"<<pphi[i]<<std::endl;
            // }

            for(int tki=0; tki<eta.size(); tki++) {
                int psize = peta.size();
                if(psize > 0) {
                    int mindr_idx = 0;
                    double mindr = (ppdg[0] == 22 or ppdg[0] == 130) ? getdR(peta[0], caloeta[tki], pphi[0], calophi[tki])
                                                                     : getdR(peta[0], eta[tki], pphi[0], phi[tki]);

                    for(int pi=1; pi<psize; pi++) {
                        double dr = (ppdg[pi] == 22 or ppdg[pi] == 130) ? getdR(peta[pi], caloeta[tki], pphi[pi], calophi[tki])
                                                                        : getdR(peta[pi], eta[tki], pphi[pi], phi[tki]);
                        if(abs(dr) < abs(mindr)) {
                            mindr = dr;
                            mindr_idx = pi;
                        }
                    }
                    mindr_vals[mindr_idx] = mindr;
                }
            }

            return mindr_vals;

        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_getmindR_PuppiCandSize_TkElRef)
