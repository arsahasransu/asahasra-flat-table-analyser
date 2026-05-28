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


    STRCPPFUNC_quantise_relisoval = """
        // Helper: reinterpret uint32_t bits as float
        inline float bits_to_float(std::uint32_t bits) {
            static_assert(sizeof(float) == sizeof(std::uint32_t), "Size mismatch");
            float result;
            std::memcpy(&result, &bits, sizeof(float)); // safe copy
            return result;
        }

        // Helper: extract float components
        inline void unpack_float(float x, std::uint32_t& sign, std::uint32_t& exp, std::uint32_t& frac) {
            std::uint32_t bits;
            std::memcpy(&bits, &x, sizeof(float));
            sign = (bits >> 31) & 1;
            exp  = (bits >> 23) & 0xFF;
            frac = bits & 0x7FFFFF;
        }

        float quantise_relisoval(float iso) {
            // Handle NaN and Inf gracefully: return as-is or clamp
            if (std::isnan(iso) || std::isinf(iso)) {
                return iso;
            }

            // Handle zero
            if (iso == 0.0f) return 0.0f;

            std::uint32_t sign, orig_exp, orig_frac;
            unpack_float(iso, sign, orig_exp, orig_frac);

            // Extract unbiased exponent
            std::int32_t unbiased_exp = static_cast<std::int32_t>(orig_exp) - 127;

            // Clamp exponent to representable range for 1–8–8 format:
            // Valid exponents: -126 (smallest normalized) to +127 (largest finite)
            // With 8-bit exponent (bias=127), max biased exp = 254 → unbiased 127
            // min normalized: biased=1 → unbiased=-126
            if (unbiased_exp < -126) {
                // Underflow: clamp to smallest normalized
                unbiased_exp = -126;
                orig_exp = 1;
                orig_frac = 0;
            } else if (unbiased_exp > 127) {
                // Overflow: clamp to largest finite (biased = 254)
                unbiased_exp = 127;
                orig_exp = 254;
                orig_frac = 0x7FFFFF; // max mantissa
            }

            // For 1–8–8 format:
            // - exponent is 8-bit (bias=127)
            // - mantissa: take top 8 bits of original 23-bit mantissa (truncation)
            //   or round: std::round(orig_frac * 2^(8-23)) / 2^(8-23)
            // Let's do rounding to nearest:

            // Shift and scale original fraction to 8 bits:
            // orig_frac has 23 bits; we want to round to 8 bits.
            // Rounding: consider bit 8 (i.e., bit 23−8 = 15 in orig_frac)
            const int bits_to_keep = 8;
            const int shift = 23 - bits_to_keep;  // 15

            std::uint32_t mantissa_8bit;
            if (shift > 0) {
                // Get top 8 bits, with rounding
                std::uint32_t guard_bit = (orig_frac >> (shift - 1)) & 1;  // bit after MSB of 8-bit window
                std::uint32_t round_bit = (orig_frac >> (shift - 2)) & 1;  // next bit for round-to-even? Let's use simple round-half-up
                std::uint32_t remainder_mask = (1u << shift) - 1u;
                std::uint32_t remainder = orig_frac & remainder_mask;

                // Round: if remainder >= half, round up
                std::uint32_t half = 1u << (shift - 1);
                if (remainder >= half) {
                    mantissa_8bit = (orig_frac >> shift) + 1u;
                } else {
                    mantissa_8bit = orig_frac >> shift;
                }

                // Handle overflow: if mantissa overflows (e.g., 0xFF → need 9 bits)
                if (mantissa_8bit >= (1u << bits_to_keep)) {
                    mantissa_8bit = (1u << bits_to_keep) - 1u;  // clamp to max (0xFF)
                    // Note: in real 1–8–8, this would increment exponent — but we skip normalization for simplicity
                    // or handle it by re-clamping exponent if desired
                }
            } else {
                mantissa_8bit = orig_frac & 0xFFu;
            }

            // If original was subnormal or zero, handle specially (here we assume normalized)
            // For normalized numbers: leading '1' is implicit → stored mantissa = 8 bits of fractional part
            // So we store mantissa_8bit as-is (fraction only)

            // Build 1–8–8 representation (17 bits total):
            // sign (1) | exponent (8) | mantissa (8)
            std::uint32_t packed = 0;
            packed |= (sign << 16);             // sign at bit 16 (top of 17-bit)
            packed |= (orig_exp << 8);          // exponent at bits 8–15
            packed |= (mantissa_8bit);          // mantissa at bits 0–7

            // Now reinterpret this 17-bit value as 32-bit float
            // But 17-bit value must be zero-extended to 32 bits, and placed in the *low 17 bits*
            // of a float's bit pattern — that's unconventional, but matches "emulating" the format.

            // Alternative: map to standard float's 32-bit layout:
            //   sign = 1 bit, exponent = 8 bits, mantissa = 23 bits.
            // To emulate 8-bit mantissa: put our 8-bit mantissa in the top 8 bits of the 23-bit fraction.
            std::uint32_t float_bits = 0;
            float_bits |= (sign << 31);
            float_bits |= (orig_exp << 23);                     // 8-bit exponent in top 8 bits of 8-bit exponent field
            float_bits |= (mantissa_8bit << (23 - 8));          // place 8-bit mantissa in top 8 bits of fraction (bits 22–15)

            return bits_to_float(float_bits);
        }
    """

    ROOT.gInterpreter.Declare(STRCPPFUNC_quantise_relisoval)

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
                    // if(abs(bkg_z0[i]-sig_vz) > dzmax) continue; INCORRECT
                    if(abs(bkg_z0[i]) >= dzmax) continue;
                }

                float calodR = getdR(sig_calo_eta, bkg_eta[i], sig_calo_phi, bkg_phi[i]);
                float nonCalodR = getdR(sig_eta, bkg_eta[i], sig_phi, bkg_phi[i]);
                // float dR = (fabs(bkg_pid[i]) == 22) || (fabs(bkg_pid[i]) == 130) ? calodR : nonCalodR;
                float dR = calodR;
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
