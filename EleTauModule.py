import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from TreeProducerEleTau import *
from CorrectionTools.ElectronSFs import *
from CorrectionTools.PileupWeightTool import *
from CorrectionTools.LeptonTauFakeSFs import *


class declareVariables(TreeProducerEleTau):
    
    def __init__(self, name):
        
        super(declareVariables, self).__init__(name)


class EleTauProducer(Module):
    
    def __init__(self, name, dataType, **kwargs):
        
        year = kwargs.get('year',  2017 )
        tes  = kwargs.get('tes',   1.0  )
        
        self.name = name
        self.out = declareVariables(name)
        
        if dataType=='data':
            self.isData = True
        else:
            self.isData = False
        
        setYear(year)
        self.eleSFs = ElectronSFs(year=year)
        self.puTool = PileupWeightTool(year=year)
        self.ltfSFs = LeptonTauFakeSFs('loose','tight',year=year)
        
        self.Nocut = 0
        self.Trigger = 1
        self.GoodElectrons = 2
        self.GoodTaus = 3
        self.GoodDiLepton = 4
        self.TotalWeighted = 15
        self.TotalWeighted_no0PU = 16
        
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.Nocut,               "no cut"                 )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.Trigger,             "trigger"                )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodElectrons,       "electron object"        )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodTaus,            "tau object"             )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodDiLepton,        "eletau pair"            )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.TotalWeighted,       "no cut, weighted"       )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.TotalWeighted_no0PU, "no cut, weighted, PU>0" )
        self.out.cutflow.GetXaxis().SetLabelSize(0.041)
    
    def beginJob(self):
        pass

    def endJob(self):
        self.out.outputfile.Write()
        self.out.outputfile.Close()

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):        
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        #electrons = Collection(event, "Electron")
        
        
        #####################################
        self.out.cutflow.Fill(self.Nocut)
        if self.isData:
          self.out.cutflow.Fill(self.TotalWeighted, 1.)
          if event.PV_npvs>0:
            self.out.cutflow.Fill(self.TotalWeighted_no0PU, 1.)
          else:
            return False
        else:
          self.out.cutflow.Fill(self.TotalWeighted, event.genWeight)
          if event.Pileup_nTrueInt>0:
            self.out.cutflow.Fill(self.TotalWeighted_no0PU, event.genWeight)
          else:
            return False
        #####################################
        
        
        # HLT_Ele32_WPTight_Gsf_L1DoubleEG
        # HLT_Ele32_WPTight_Gsf
        if not (event.HLT_Ele35_WPTight_Gsf or event.HLT_Ele32_WPTight_Gsf):
            return False
        
        
        #####################################
        self.out.cutflow.Fill(self.Trigger)
        #####################################
        
        
        idx_goodelectrons = [ ]
        for ielectron in range(event.nElectron):
            if event.Electron_pt[ielectron] < 36: continue
            if abs(event.Electron_eta[ielectron]) > 2.1: continue
            if abs(event.Electron_dz[ielectron]) > 0.2: continue
            if abs(event.Electron_dxy[ielectron]) > 0.045: continue
            if event.Electron_convVeto[ielectron] !=1: continue
            if ord(event.Electron_lostHits[ielectron]) > 1: continue
            if event.Electron_mvaFall17Iso_WP80[ielectron] < 0.5: continue
            idx_goodelectrons.append(ielectron)
        
        if len(idx_goodelectrons)==0:
            return False
        
        
        #####################################
        self.out.cutflow.Fill(self.GoodElectrons)
        #####################################
        
        
        idx_goodtaus = [ ]
        for itau in range(event.nTau):
            if event.Tau_pt[itau] < 20: continue
            if abs(event.Tau_eta[itau]) > 2.3: continue
            if abs(event.Tau_dz[itau]) > 0.2: continue
            if event.Tau_decayMode[itau] not in [0,1,10]: continue
            if abs(event.Tau_charge[itau])!=1: continue
            #if ord(event.Tau_idAntiEle[itau])<1: continue
            #if ord(event.Tau_idAntiMu[itau])<1: continue
            idx_goodtaus.append(itau)
        
        if len(idx_goodtaus)==0:
            return False
        
        
        #####################################
        self.out.cutflow.Fill(self.GoodTaus)
        #####################################
        
        
        # to check dR matching
        electrons = Collection(event, 'Electron')
        taus = Collection(event, 'Tau')
        ltaus = [ ]
        for idx1 in idx_goodelectrons:
            for idx2 in idx_goodtaus:
                dR = taus[idx2].p4().DeltaR(electrons[idx1].p4())
                if dR < 0.5: continue
                electron_reliso = event.Electron_pfRelIso03_all[idx1]
                ltau = LeptonTauPair(idx1, event.Electron_pt[idx1], electron_reliso, idx2, event.Tau_pt[idx2], event.Tau_rawMVAoldDM2017v2[idx2])
                ltaus.append(ltau)
        
        if len(ltaus)==0:
            return False
        
        ltau = bestDiLepton(ltaus)
        #print 'chosen tau1 (idx, pt) = ', ltau.id1, ltau.tau1_pt, 'check', taus[ltau.id1].p4().Pt()
        #print 'chosen tau2 (idx, pt) = ', ltau.id2, ltau.tau2_pt, 'check', taus[ltau.id2].p4().Pt()
        
        
        #####################################
        self.out.cutflow.Fill(self.GoodDiLepton)
        #####################################
        
        
        jetIds = [ ]
        jets = Collection(event, 'Jet')
        #jets = filter(self.jetSel,jets):
        nfjets = 0
        ncjets = 0
        nbtag = 0
        for ijet in range(event.nJet):  
        #for j in filter(self.jetSel,jets):
            if event.Jet_pt[ijet] < 30: continue
            if abs(event.Jet_eta[ijet]) > 4.7: continue
            dR = electrons[ltau.id1].p4().DeltaR(jets[ijet].p4())
            if dR < 0.5: continue
            dR = taus[ltau.id2].p4().DeltaR(jets[ijet].p4())
            if dR < 0.5: continue
            
            #print '#', ijet, 'pt = ', jets[ijet].p4().Pt(), event.Jet_pt[ijet]
            jetIds.append(ijet)
            
            if abs(event.Jet_eta[ijet]) > 2.4:
                nfjets += 1
            else:
                ncjets += 1
            
            if event.Jet_btagCSVV2[ijet] > 0.8838:
                nbtag += 1
        
        #eventSum = ROOT.TLorentzVector()
        #
        #for lep in electrons :
        #    eventSum += lep.p4()
        #for lep in electrons :
        #    eventSum += lep.p4()
        #for j in filter(self.jetSel,jets):
        #    eventSum += j.p4()
        
        
        # ELECTRON
        self.out.pt_1[0]                       = event.Electron_pt[ltau.id1]
        self.out.eta_1[0]                      = event.Electron_eta[ltau.id1]
        self.out.phi_1[0]                      = event.Electron_phi[ltau.id1]
        self.out.mass_1[0]                     = event.Electron_mass[ltau.id1]
        self.out.dxy_1[0]                      = event.Electron_dxy[ltau.id1]
        self.out.dz_1[0]                       = event.Electron_dz[ltau.id1]         
        self.out.q_1[0]                        = event.Electron_charge[ltau.id1]
        self.out.pfRelIso03_all_1[0]           = event.Electron_pfRelIso03_all[ltau.id1]
        self.out.cutBased_1[0]                 = event.Electron_cutBased[ltau.id1]
        self.out.mvaFall17Iso_1[0]             = event.Electron_mvaFall17Iso[ltau.id1]
        self.out.mvaFall17Iso_WP80_1[0]        = event.Electron_mvaFall17Iso_WP80[ltau.id1]
        self.out.mvaFall17Iso_WP90_1[0]        = event.Electron_mvaFall17Iso_WP90[ltau.id1]
        self.out.mvaFall17Iso_WPL_1[0]         = event.Electron_mvaFall17Iso_WPL[ltau.id1]
        
        
        # TAU
        self.out.pt_2[0]                       = event.Tau_pt[ltau.id2]
        self.out.eta_2[0]                      = event.Tau_eta[ltau.id2]
        self.out.phi_2[0]                      = event.Tau_phi[ltau.id2]
        self.out.mass_2[0]                     = event.Tau_mass[ltau.id2]
        self.out.dxy_2[0]                      = event.Tau_dxy[ltau.id2]
        self.out.dz_2[0]                       = event.Tau_dz[ltau.id2]         
        self.out.leadTkPtOverTauPt_2[0]        = event.Tau_leadTkPtOverTauPt[ltau.id2]
        self.out.chargedIso_2[0]               = event.Tau_chargedIso[ltau.id2]
        self.out.neutralIso_2[0]               = event.Tau_neutralIso[ltau.id2]
        self.out.photonsOutsideSignalCone_2[0] = event.Tau_photonsOutsideSignalCone[ltau.id2]
        self.out.puCorr_2[0]                   = event.Tau_puCorr[ltau.id2]
        self.out.rawAntiEle_2[0]               = event.Tau_rawAntiEle[ltau.id2]
        self.out.rawIso_2[0]                   = event.Tau_rawIso[ltau.id2]
        self.out.rawMVAnewDM2017v2_2[0]        = event.Tau_rawMVAnewDM2017v2[ltau.id2]
        self.out.rawMVAoldDM_2[0]              = event.Tau_rawMVAoldDM[ltau.id2]
        self.out.rawMVAoldDM2017v1_2[0]        = event.Tau_rawMVAoldDM2017v1[ltau.id2]
        self.out.rawMVAoldDM2017v2_2[0]        = event.Tau_rawMVAoldDM2017v2[ltau.id2]
        self.out.q_2[0]                        = event.Tau_charge[ltau.id2]
        self.out.decayMode_2[0]                = event.Tau_decayMode[ltau.id2]
        self.out.rawAntiEleCat_2[0]            = event.Tau_rawAntiEleCat[ltau.id2]
        self.out.idAntiEle_2[0]                = ord(event.Tau_idAntiEle[ltau.id2])
        self.out.idAntiMu_2[0]                 = ord(event.Tau_idAntiMu[ltau.id2])
        self.out.idDecayMode_2[0]              = event.Tau_idDecayMode[ltau.id2]
        self.out.idDecayModeNewDMs_2[0]        = event.Tau_idDecayModeNewDMs[ltau.id2]
        self.out.idMVAnewDM2017v2_2[0]         = ord(event.Tau_idMVAnewDM2017v2[ltau.id2])
        self.out.idMVAoldDM_2[0]               = ord(event.Tau_idMVAoldDM[ltau.id2])
        self.out.idMVAoldDM2017v1_2[0]         = ord(event.Tau_idMVAoldDM2017v1[ltau.id2])
        self.out.idMVAoldDM2017v2_2[0]         = ord(event.Tau_idMVAoldDM2017v2[ltau.id2])
        
        
        # GENERATOR
        #print type(event.Tau_genPartFlav[ltau.id2])        
        if not self.isData:
            self.out.genPartFlav_1[0]          = ord(event.Electron_genPartFlav[ltau.id1])
            self.out.genPartFlav_2[0]          = ord(event.Tau_genPartFlav[ltau.id2])
            
            genvistau = Collection(event, 'GenVisTau')
            dRmax = 1000
            gendm = -1
            genpt = -1
            geneta = -1
            genphi = -1
            for igvt in range(event.nGenVisTau):
                dR = genvistau[igvt].p4().DeltaR(taus[ltau.id2].p4())
                if dR < 0.5 and dR < dRmax:
                    dRmax = dR
                    gendm = event.GenVisTau_status[igvt]
                    genpt = event.GenVisTau_pt[igvt]
                    geneta = event.GenVisTau_eta[igvt]
                    genphi = event.GenVisTau_phi[igvt]
            
            self.out.gendecayMode_2[0]         = gendm
            self.out.genvistaupt_2[0]          = genpt
            self.out.genvistaueta_2[0]         = geneta
            self.out.genvistauphi_2[0]         = genphi
        
        
        # EVENT
        self.out.isData[0]                     = self.isData
        self.out.run[0]                        = event.run
        self.out.luminosityBlock[0]            = event.luminosityBlock
        self.out.event[0]                      = event.event & 0xffffffffffffffff
        self.out.MET_pt[0]                     = event.MET_pt
        self.out.MET_phi[0]                    = event.MET_phi
        self.out.PuppiMET_pt[0]                = event.PuppiMET_pt
        self.out.PuppiMET_phi[0]               = event.PuppiMET_phi
        self.out.MET_significance[0]           = event.MET_significance
        self.out.MET_covXX[0]                  = event.MET_covXX
        self.out.MET_covXY[0]                  = event.MET_covXY
        self.out.MET_covYY[0]                  = event.MET_covYY
        self.out.fixedGridRhoFastjetAll[0]     = event.fixedGridRhoFastjetAll
        self.out.npvs[0]                       = event.PV_npvs
        self.out.npvsGood[0]                   = event.PV_npvsGood
        
        if not self.isData:
          self.out.GenMET_pt[0]                = event.GenMET_pt
          self.out.GenMET_phi[0]               = event.GenMET_phi
          self.out.nPU[0]                      = event.Pileup_nPU
          self.out.nTrueInt[0]                 = event.Pileup_nTrueInt
          self.out.genWeight[0]                = event.genWeight
          #self.out.LHE_Njets[0]                = event.LHE_Njets
          self.out.LHE_Njets[0]                = getattr(event,'LHE_Njets',-1)
        
        
        # JETS
        if len(jetIds)>0:
          self.out.jpt_1[0]                    = event.Jet_pt[jetIds[0]]
          self.out.jeta_1[0]                   = event.Jet_eta[jetIds[0]]
          self.out.jphi_1[0]                   = event.Jet_phi[jetIds[0]]
          self.out.jcsvv2_1[0]                 = event.Jet_btagCSVV2[jetIds[0]]
          self.out.jdeepb_1[0]                 = event.Jet_btagDeepB[jetIds[0]]
        else:
          self.out.jpt_1[0]                    = -9.
          self.out.jeta_1[0]                   = -9.
          self.out.jphi_1[0]                   = -9.
          self.out.jcsvv2_1[0]                 = -9.
          self.out.jdeepb_1[0]                 = -9.
        
        if len(jetIds)>1:
          self.out.jpt_2[0]                    = event.Jet_pt[jetIds[1]]
          self.out.jeta_2[0]                   = event.Jet_eta[jetIds[1]]
          self.out.jphi_2[0]                   = event.Jet_phi[jetIds[1]]
          self.out.jcsvv2_2[0]                 = event.Jet_btagCSVV2[jetIds[1]]
          self.out.jdeepb_2[0]                 = event.Jet_btagDeepB[jetIds[1]]
        else:
          self.out.jpt_2[0]                    = -9.
          self.out.jeta_2[0]                   = -9.
          self.out.jphi_2[0]                   = -9.
          self.out.jcsvv2_2[0]                 = -9.
          self.out.jdeepb_2[0]                 = -9.
        
        self.out.njets[0]                      = len(jetIds)
        self.out.njets50[0]                    = len([j for j in jetIds if event.Jet_pt[j]>50])
        self.out.nfjets[0]                     = nfjets
        self.out.ncjets[0]                     = ncjets
        self.out.nbtag[0]                      = nbtag
        
        self.out.pfmt_1[0]                     = math.sqrt( 2 * self.out.pt_1[0] * self.out.MET_pt[0] * ( 1 - math.cos(deltaPhi(self.out.phi_1[0], self.out.MET_phi[0])) ) );
        self.out.pfmt_2[0]                     = math.sqrt( 2 * self.out.pt_2[0] * self.out.MET_pt[0] * ( 1 - math.cos(deltaPhi(self.out.phi_2[0], self.out.MET_phi[0])) ) );
        
        self.out.m_vis[0]                      = (electrons[ltau.id1].p4() + taus[ltau.id2].p4()).M()
        self.out.pt_tt[0]                      = (electrons[ltau.id1].p4() + taus[ltau.id2].p4()).Pt()
        
        self.out.dR_ll[0]                      = electrons[ltau.id1].p4().DeltaR(taus[ltau.id2].p4())
        self.out.dphi_ll[0]                    = deltaPhi(self.out.phi_1[0], self.out.phi_2[0])
        
        
        # PZETA
        leg1 = ROOT.TVector3(electrons[ltau.id1].p4().Px(), electrons[ltau.id1].p4().Py(), 0.)
        leg2 = ROOT.TVector3(taus[ltau.id2].p4().Px(), taus[ltau.id2].p4().Py(), 0.)
        met_tlv = ROOT.TLorentzVector()
        met_tlv.SetPxPyPzE(self.out.MET_pt[0]*math.cos(self.out.MET_phi[0]), 
                           self.out.MET_pt[0]*math.sin(self.out.MET_phi[0]),
                           0, 
                           self.out.MET_pt[0])
        metleg = met_tlv.Vect()
        zetaAxis = ROOT.TVector3(leg1.Unit() + leg2.Unit()).Unit()
        pZetaVis_ = leg1*zetaAxis + leg2*zetaAxis
        pZetaMET_ = metleg*zetaAxis
        self.out.pzetamiss[0]  = pZetaMET_
        self.out.pzetavis[0]   = pZetaVis_
        self.out.pzeta_disc[0] = pZetaMET_ - 0.5*pZetaVis_
        
        
        # VETOS
        self.out.extramuon_veto[0], self.out.extraelec_veto[0], self.out.dilepton_veto[0]  = extraLeptonVetos(event, [-1], [ltau.id1], self.name)        
        
        
        # WEIGHTS
        if not self.isData:
          self.out.puweight[0]      = self.puTool.getWeight(event.Pileup_nTrueInt)
          self.out.trigweight[0]    = self.eleSFs.getTriggerSF(self.out.pt_1[0], self.out.eta_1[0])
          self.out.idisoweight_1[0] = self.eleSFs.getIdIsoSF(self.out.pt_1[0],self.out.eta_1[0])
          self.out.idisoweight_2[0] = self.ltfSFs.getSF(self.out.genPartFlav_2[0],self.out.eta_2[0])
          self.out.weight[0]        = self.out.trigweight[0]*self.out.puweight[0]
        
        
        self.out.tree.Fill()
        return True
