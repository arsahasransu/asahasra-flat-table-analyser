import ROOT

from cpp_utils.cmssw_cpp_utils import cmssw_cpp_utils
from cpp_utils.calciso_annularcone import calciso_annularcone
from cpp_utils.tkelem_puppi_iso_def_legacy2026 import tkelem_puppi_iso_def_legacy2026


def define_cpp_utils():

    # Order is important here due to inter-dependent functions
    cmssw_cpp_utils()
    calciso_annularcone()

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

    # Import and define the other specialised CPP utils
    tkelem_puppi_iso_def_legacy2026()