import ROOT

def define_cpp_utils():

    calculate_angdiff = """
        #include <tuple>

        std::tuple< ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>,
                    ROOT::VecOps::RVec<float> > calculate_angdiff(ROOT::VecOps::RVec<float>& p1_pt,
                                                                   ROOT::VecOps::RVec<float>& p1_eta,
                                                                   ROOT::VecOps::RVec<float>& p1_phi,
                                                                   ROOT::VecOps::RVec<int>& p1_q,
                                                                   ROOT::VecOps::RVec<float>& p2_pt,
                                                                   ROOT::VecOps::RVec<float>& p2_eta,
                                                                   ROOT::VecOps::RVec<float>& p2_phi,
                                                                   ROOT::VecOps::RVec<int>& p2_q) {

            ROOT::VecOps::RVec<float> deta;
            ROOT::VecOps::RVec<float> dphi;
            ROOT::VecOps::RVec<float> p1qdphi;
            ROOT::VecOps::RVec<float> p2qdphi;
            ROOT::VecOps::RVec<float> dR;

            for(unsigned int i=0; i<p1_pt.size(); i++) {

                TVector3 p1;
                p1.SetPtEtaPhi(p1_pt[i], p1_eta[i], p1_phi[i]);

                float min_dR = 999999.0;
                unsigned int min_dR_pos = 0;

                for(unsigned int j=0; j<p2_pt.size(); j++) {

                    TVector3 p2;
                    p2.SetPtEtaPhi(p2_pt[j], p2_eta[j], p2_phi[j]);

                    if(j==0) {
                        min_dR = p1.DeltaR(p2);
                    }

                    if( p1.DeltaR(p2) < min_dR ){
                        min_dR = p1.DeltaR(p2);
                        min_dR_pos = j;
                    }
                }

                TVector3 p2_mindr;
                p2_mindr.SetPtEtaPhi(p2_pt[min_dR_pos], p2_eta[min_dR_pos], p2_phi[min_dR_pos]);
                deta.push_back( p1.Eta() - p2_mindr.Eta() );
                dphi.push_back( p1.DeltaPhi(p2_mindr) );
                p1qdphi.push_back( p1_q[i]*(p1.DeltaPhi(p2_mindr)) );
                p2qdphi.push_back( p2_q[min_dR_pos]*(p1.DeltaPhi(p2_mindr)) );
                dR.push_back( p1.DeltaR(p2_mindr) );
            }

            return std::make_tuple(deta, dphi, p1qdphi, p2qdphi, dR);
        }
    """

    ROOT.gInterpreter.Declare(calculate_angdiff)

    perform_genmatch = """
        #include <tuple>

        std::tuple< ROOT::VecOps::RVec<int>, 
                    ROOT::VecOps::RVec<int> > perform_genmatch(ROOT::VecOps::RVec<float>& gen_pt,
                                                                ROOT::VecOps::RVec<float>& gen_eta,
                                                                ROOT::VecOps::RVec<float>& gen_phi,
                                                                ROOT::VecOps::RVec<float>& pup_pt,
                                                                ROOT::VecOps::RVec<float>& pup_eta,
                                                                ROOT::VecOps::RVec<float>& pup_phi,
                                                                float dR_cut) {

            ROOT::VecOps::RVec<int> pupidx;
            ROOT::VecOps::RVec<int> genidx(pup_pt.size(), -1);

            for(unsigned int i=0; i<gen_pt.size(); i++) {

                if(pup_pt.size() == 0) {
                    genidx.push_back(-1);
                    pupidx.push_back(-1);
                    continue;
                }

                TVector3 gen, pup0;
                gen.SetPtEtaPhi(gen_pt[i], gen_eta[i], gen_phi[i]);
                pup0.SetPtEtaPhi(pup_pt[0], pup_eta[0], pup_phi[0]);

                float min_dR = gen.DeltaR(pup0);
                int pupmatch_idx = min_dR < dR_cut ? 0 : -1;

                for(unsigned int j=1; j<pup_pt.size(); j++) {

                    TVector3 pup;
                    pup.SetPtEtaPhi(pup_pt[j], pup_eta[j], pup_phi[j]);

                    if( gen.DeltaR(pup) < dR_cut && gen.DeltaR(pup) < min_dR ) {
                        min_dR = gen.DeltaR(pup);
                        pupmatch_idx = j;
                    }

                }

                pupidx.push_back(pupmatch_idx);
                genidx[pupmatch_idx] = i;
            }
            
            // Check for duplicates
            for(unsigned int i=0; i<pupidx.size(); i++) {
                for(unsigned int j=i+1; j<pupidx.size(); j++) {
                    if(pupidx[i] == pupidx[j] && pupidx[i] != -1) {
                        pupidx[j] = -1;
                        genidx[pupidx[j]] = -1;
                    }
                }
            }

            return std::make_tuple(pupidx, genidx);
        }
    """

    ROOT.gInterpreter.Declare(perform_genmatch)