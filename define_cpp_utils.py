import ROOT


def define_cpp_utils():

    STRCPPFUNC_getminangs = """
        std::tuple< ROOT::VecOps::RVec<double>, ROOT::VecOps::RVec<double>,
            ROOT::VecOps::RVec<double> > getminangs(ROOT::VecOps::RVec<float> &geta,
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
                for (int j = 0; j < eta.size(); j++) {
                    double deta = geta[i] - eta[j];
                    double dphi = gphi[i] - phi[j];
                    dphi = abs(dphi) > M_PI ? dphi - 2 * M_PI : dphi;
                    double dR = sqrt(deta * deta + dphi * dphi);
                    if (dR < min_dR) min_dR = dR;
                    if (fabs(deta) < fabs(min_deta)) min_deta = deta;
                    if (fabs(dphi) < fabs(min_dphi)) min_dphi = dphi;
                }
                mindeta.push_back(min_deta);
                mindphi.push_back(min_dphi);
                mindR.push_back(min_dR);
            }
            return std::make_tuple(mindeta, mindphi, mindR);
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
                    double deta = geta[i] - eta[j];
                    double dphi = gphi[i] - phi[j];
                    dphi = dphi > M_PI ? dphi - 2 * M_PI : dphi;
                    double dR = sqrt(deta * deta + dphi * dphi);
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
                    float deta = sig_eta[i] - bkg_eta[j];
                    float dphi = sig_phi[i] - bkg_phi[j];
                    dphi = abs(dphi) > M_PI ? dphi - 2 * M_PI : dphi;
                    float dR = sqrt(deta * deta + dphi * dphi);
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
                float deta = sig_eta - bkg_eta[i];
                float dcaloeta = sig_calo_eta - bkg_eta[i];
                float dphi = sig_phi - bkg_phi[i];
                float dcalophi = sig_calo_phi - bkg_phi[i];
                dphi = abs(dphi) > M_PI ? dphi - 2 * M_PI : dphi;
                dcalophi = abs(dcalophi) > M_PI ? dcalophi - 2 * M_PI : dcalophi;
                float dR = (fabs(bkg_pid[i]) == 22) || (fabs(bkg_pid[i]) == 130) ?
                           sqrt(dcaloeta * dcaloeta + dcalophi * dcalophi) : sqrt(deta * deta + dphi * dphi);
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
