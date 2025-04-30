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
                    if (dphi < min_dphi) min_dphi = dphi;
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