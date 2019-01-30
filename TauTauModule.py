import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from TreeProducerTauTau import *
from CorrectionTools.TauTauSFs import TauTauSFs
from CorrectionTools.PileupWeightTool import PileupWeightTool
from CorrectionTools.LeptonTauFakeSFs import LeptonTauFakeSFs
from CorrectionTools.BTaggingTool import BTagWeightTool, BTagWPs
from CorrectionTools.RecoilCorrectionTool import RecoilCorrectionTool, getZPTMass


class TauTauProducer(Module):

    def __init__(self, name, dataType, **kwargs):
        
        year           = kwargs.get('year',  2017 )
        tes            = kwargs.get('tes',   1.0  )
        channel        = 'tautau'
        
        self.name      = name
        self.year      = year
        self.tes       = tes
        self.out       = TreeProducerTauTau(name)
        self.isData    = dataType=='data'
        self.doZpt     = 'DY' in self.name
        
        setYear(year)
        self.vlooseIso = getVLooseTauIso(year)
        if year==2017:
          self.trigger = lambda e: e.HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg or e.HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg or e.HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg
        else:
          if self.isData:
            self.trigger = lambda e: e.HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg or e.HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg or e.HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg \
                                       if e.run<317509 else e.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg
          else:
            self.trigger = lambda e: e.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg
        self.tauCutPt = 40
        
        if not self.isData:
          self.tauSFs   = TauTauSFs('tight',year=year)
          self.tauSFsVT = TauTauSFs('vtight',year=year)
          self.ltfSFs   = LeptonTauFakeSFs('loose','vloose',year=year)
          self.puTool   = PileupWeightTool(year=year)
          self.btagTool      = BTagWeightTool('CSVv2','medium',channel=channel,year=year)
          self.btagTool_deep = BTagWeightTool('DeepCSV','medium',channel=channel,year=year)
          if self.doZpt:
            self.recoilTool  = RecoilCorrectionTool(year=year)
        self.csvv2_wp   = BTagWPs('CSVv2',year=year)
        self.deepcsv_wp = BTagWPs('DeepCSV',year=year)
        
        self.Nocut = 0
        self.Trigger = 1
        self.GoodTaus = 2
        self.GoodDiTau = 3
        
        self.Nocut_GT = 20
        self.Trigger_GT = 21
        self.GoodTaus_GT = 22
        self.GoodDiTau_GT = 23
        
        self.TotalWeighted = 15
        self.TotalWeighted_no0PU = 16
        
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.Nocut,               "no cut"                 )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.Trigger,             "trigger"                )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodTaus,            "tau objects"            )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodDiTau,           "ditau pair"             )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.Nocut_GT,            "no cut, GM"             )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.Trigger_GT,          "trigger, GM"            )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodTaus_GT,         "tau objects, GM"        )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.GoodDiTau_GT,        "ditau pair, GM"         )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.TotalWeighted,       "no cut, weighted"       )
        self.out.cutflow.GetXaxis().SetBinLabel(1+self.TotalWeighted_no0PU, "no cut, weighted, PU>0" )
        self.out.cutflow.GetXaxis().SetLabelSize(0.041)
        
    def beginJob(self):
        pass
        
    def endJob(self):
        #check = self.out.outputfile.mkdir('check')
        #self.pileup.SetDirectory(check)
        if not self.isData:
          self.btagTool.setDirectory(self.out.outputfile,'btag')
          self.btagTool_deep.setDirectory(self.out.outputfile,'btag')
        self.out.outputfile.Write()
        self.out.outputfile.Close()
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):        
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #electrons = Collection(event, "Electron")
        #muons = Collection(event, "Muon")
        
        ##print '-'*80
        ngentauhads = 0
        ngentaus = 0
        #if not self.isData:            
        #    for igp in range(event.nGenPart):
        #        if abs(event.GenPart_pdgId[igp])==15 and event.GenPart_status[igp]==2:
        #            genflag = event.GenPart_statusFlags[igp]
        #            binary = format(genflag,'b').zfill(15)
        #            # 0 : isPrompt
        #            # 1 : isDecayedLeptonHadron
        #            # 2 : isTauDecayProduct
        #            # 3 : isPromptTauDecayProduct
        #            # 4 : isDirectTauDecayProduct
        #            # 5 : isDirectPromptTauDecayProduct
        #            # 6 : isDirectHadronDecayProduct
        #            # 7 : isHardProcess
        #            # 8 : fromHardProcess
        #            # 9 : isHardProcessTauDecayProduct
        #            # 10 : isDirectHardProcessTauDecayProduct
        #            # 11 : fromHardProcessBeforeFSR
        #            # 12 : isFirstCopy
        #            # 13 : isLastCopy
        #            # 14 : isLastCopyBeforeFSR
        #            
        #            if int(binary[14])==0: continue
        #            if int(binary[6])==0: continue
        #            #print 'Tau found with status = 2 (pt, eta) = ', event.GenPart_pt[igp], event.GenPart_eta[igp], event.GenPart_statusFlags[igp]
        #            
        #            ngentaus += 1
        #            _pdg_ = -1
        #            _idx_ = event.GenPart_genPartIdxMother[igp]
        #            #_status_ = -1
        #            flag_resonance = False
        #            
        #            while abs(_pdg_) not in [9000002, 9000006]:
        #                if _idx_==-1: break
        #                _pdg_ = event.GenPart_pdgId[_idx_]
        #                # _status_ = event.GenPart_status[_idx_]
        #                _idx_ = event.GenPart_genPartIdxMother[_idx_]
        #                if abs(_pdg_) > 30 and abs(_pdg_) not in [9000002, 9000006]: 
        #                    flag_resonance = True
        #                #print '\t (pdg, mother id) = ', _pdg_, _status_, _idx_
        #            if flag_resonance: continue
        #            _dr_ = 100.
        #            for igvt in range(event.nGenVisTau):
        #                dr = deltaR(event.GenPart_eta[igp], event.GenPart_phi[igp], event.GenVisTau_eta[igvt], event.GenVisTau_phi[igvt])
        #                #print dr, _dr_, event.GenPart_eta[igp], event.GenPart_phi[igp], event.GenVisTau_eta[igvt], event.GenVisTau_phi[igvt]
        #                if _dr_ > dr:
        #                    _dr_ = dr
        #            #print 'match !',_pdg_, event.nGenVisTau,  _dr_
        #            if _dr_ < 0.1:
        #                ngentauhads += 1
        #    
        #    #for igvt in range(event.nGenVisTau):
        #    #    print 'status = ', event.GenVisTau_status[igvt], 'mother ID = ', event.GenVisTau_genPartIdxMother[igvt], 'pt = ', event.GenVisTau_pt[igvt], ', eta = ', event.GenVisTau_eta[igvt]
        #    #    ngentauhads += 1
        #    
        #    if ngentaus != 2:
        #       print 'WOW!!! ngentaus = %d != 2'%(ngentaus)
        
        
        #####################################
        self.out.cutflow.Fill(self.Nocut)
        #if ngentauhads == 2:
        #   self.out.cutflow.Fill(self.Nocut_GT)
        if self.isData:
          self.out.cutflow.Fill(self.TotalWeighted, 1.)
          if event.PV_npvs>0:
            self.out.cutflow.Fill(self.TotalWeighted_no0PU, 1.)
          else:
            return False
        else:
          self.out.cutflow.Fill(self.TotalWeighted, event.genWeight)
          self.out.pileup.Fill(event.Pileup_nTrueInt)
          if event.Pileup_nTrueInt>0:
            self.out.cutflow.Fill(self.TotalWeighted_no0PU, event.genWeight)
          else:
            return False
        #####################################        
        
        
        #if event.HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg or : 
        #print 'trig = ', event.HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg, event.HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg, event.HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg
        #if event.HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg > 0.5 or event.HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg > 0.5 or event.HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg > 0.5:
        if not self.trigger(event):
            return False
        
        #####################################
        self.out.cutflow.Fill(self.Trigger)
        #if ngentauhads == 2:
        #    self.out.cutflow.Fill(self.Trigger_GT)
        #####################################
        
        
        idx_goodtaus = [ ]
        for itau in range(event.nTau):
            #print 'pt=', event.Tau_pt[itau], abs(event.Tau_eta[itau])
            if event.Tau_pt[itau] < self.tauCutPt: continue
            if abs(event.Tau_eta[itau]) > 2.1: continue
            if abs(event.Tau_dz[itau]) > 0.2: continue
            if event.Tau_decayMode[itau] not in [0,1,10]: continue
            if abs(event.Tau_charge[itau])!=1: continue
            #print itau, 'decay mode = ', event.Tau_decayMode[itau] 
            if not self.vlooseIso(event,itau): continue
            idx_goodtaus.append(itau)
        
        if len(idx_goodtaus)<2:
            return False
        
        #####################################
        self.out.cutflow.Fill(self.GoodTaus)
        #if ngentauhads == 2:
        #    self.out.cutflow.Fill(self.GoodTaus_GT)
        #####################################
        
        
        taus = Collection(event, 'Tau')
        ditaus = [ ]
        for idx1 in idx_goodtaus:
          for idx2 in idx_goodtaus:
              if idx1 >= idx2: continue
              dR = taus[idx1].p4().DeltaR(taus[idx2].p4())
              if dR < 0.5: continue
              ditau = DiTauPair(idx1, event.Tau_pt[idx1], event.Tau_rawMVAoldDM[idx1],
                                idx2, event.Tau_pt[idx2], event.Tau_rawMVAoldDM[idx2])
              ditaus.append(ditau)
        
        if len(ditaus)==0:
            return False
        
        ditau = bestDiLepton(ditaus)
        tau1  = taus[ditau.id1].p4()
        tau2  = taus[ditau.id2].p4()
        #print 'chosen tau1 (idx, pt) = ', ditau.id1, ditau.tau1_pt, 'check', tau1.p4().Pt()
        #print 'chosen tau2 (idx, pt) = ', ditau.id2, ditau.tau2_pt, 'check', tau2.p4().Pt()
        
        #####################################
        self.out.cutflow.Fill(self.GoodDiTau)
        #if ngentauhads == 2:
        #    self.out.cutflow.Fill(self.GoodDiTau_GT)
        #####################################
        
        
        jetIds  = [ ]
        bjetIds = [ ]
        jets    = Collection(event, 'Jet')
        #jets = filter(self.jetSel,jets):
        nfjets  = 0
        ncjets  = 0
        nbtag   = 0
        for ijet in range(event.nJet):
            if event.Jet_pt[ijet] < 30: continue
            if abs(event.Jet_eta[ijet]) > 4.7: continue
            if tau1.DeltaR(jets[ijet].p4()) < 0.5: continue
            if tau2.DeltaR(jets[ijet].p4()) < 0.5: continue
            jetIds.append(ijet)
            
            if abs(event.Jet_eta[ijet]) > 2.4:
              nfjets += 1
            else:
              ncjets += 1
            
            if event.Jet_btagDeepB[ijet] > self.deepcsv_wp.medium:
              nbtag += 1
              bjetIds.append(ijet)
        
        if not self.isData and self.vlooseIso(event,ditau.id1) and self.vlooseIso(event,ditau.id2):
          self.btagTool.fillEfficiencies(event,jetIds)
          self.btagTool_deep.fillEfficiencies(event,jetIds)
        
        #eventSum = ROOT.TLorentzVector()
        #
        #for lep in muons :
        #    eventSum += lep.p4()
        #for lep in electrons :
        #    eventSum += lep.p4()
        #for j in filter(self.jetSel,jets):
        #    eventSum += j.p4()
        
        
        # TAU 1
        self.out.pt_1[0]                       = event.Tau_pt[ditau.id1]
        self.out.eta_1[0]                      = event.Tau_eta[ditau.id1]
        self.out.phi_1[0]                      = event.Tau_phi[ditau.id1]
        self.out.mass_1[0]                     = event.Tau_mass[ditau.id1]
        self.out.dxy_1[0]                      = event.Tau_dxy[ditau.id1]
        self.out.dz_1[0]                       = event.Tau_dz[ditau.id1]         
        self.out.leadTkPtOverTauPt_1[0]        = event.Tau_leadTkPtOverTauPt[ditau.id1]
        self.out.chargedIso_1[0]               = event.Tau_chargedIso[ditau.id1]
        self.out.neutralIso_1[0]               = event.Tau_neutralIso[ditau.id1]
        self.out.photonsOutsideSignalCone_1[0] = event.Tau_photonsOutsideSignalCone[ditau.id1]
        self.out.puCorr_1[0]                   = event.Tau_puCorr[ditau.id1]
        self.out.rawAntiEle_1[0]               = event.Tau_rawAntiEle[ditau.id1]
        self.out.rawIso_1[0]                   = event.Tau_rawIso[ditau.id1]
        self.out.rawMVAnewDM2017v2_1[0]        = event.Tau_rawMVAnewDM2017v2[ditau.id1]
        self.out.rawMVAoldDM_1[0]              = event.Tau_rawMVAoldDM[ditau.id1]
        self.out.rawMVAoldDM2017v1_1[0]        = event.Tau_rawMVAoldDM2017v1[ditau.id1]
        self.out.rawMVAoldDM2017v2_1[0]        = event.Tau_rawMVAoldDM2017v2[ditau.id1]
        self.out.q_1[0]                        = event.Tau_charge[ditau.id1]
        self.out.decayMode_1[0]                = event.Tau_decayMode[ditau.id1]
        ###self.out.rawAntiEleCat_1[0]            = event.Tau_rawAntiEleCat[ditau.id1]
        self.out.idAntiEle_1[0]                = ord(event.Tau_idAntiEle[ditau.id1])
        self.out.idAntiMu_1[0]                 = ord(event.Tau_idAntiMu[ditau.id1])
        self.out.idDecayMode_1[0]              = event.Tau_idDecayMode[ditau.id1]
        self.out.idDecayModeNewDMs_1[0]        = event.Tau_idDecayModeNewDMs[ditau.id1]
        self.out.idMVAnewDM2017v2_1[0]         = ord(event.Tau_idMVAnewDM2017v2[ditau.id1])
        self.out.idMVAoldDM_1[0]               = ord(event.Tau_idMVAoldDM[ditau.id1])
        self.out.idMVAoldDM2017v1_1[0]         = ord(event.Tau_idMVAoldDM2017v1[ditau.id1])
        self.out.idMVAoldDM2017v2_1[0]         = ord(event.Tau_idMVAoldDM2017v2[ditau.id1])
        
        
        # GENERATOR 1
        if not self.isData:
          self.out.genPartFlav_1[0]            = ord(event.Tau_genPartFlav[ditau.id1])
          genvistau = Collection(event, 'GenVisTau')
          _drmax_ = 1000
          gendm   = -1
          genpt   = -1
          geneta  = -1
          genphi  = -1
          for igvt in range(event.nGenVisTau):
            _dr_ = genvistau[igvt].p4().DeltaR(tau1)
            if _dr_ < 0.5 and _dr_ < _drmax_:
              _drmax_ = _dr_
              gendm   = event.GenVisTau_status[igvt]
              genpt   = event.GenVisTau_pt[igvt]
              geneta  = event.GenVisTau_eta[igvt]
              genphi  = event.GenVisTau_phi[igvt]
          self.out.gendecayMode_1[0]           = gendm
          self.out.genvistaupt_1[0]            = genpt
          self.out.genvistaueta_1[0]           = geneta
          self.out.genvistauphi_1[0]           = genphi
        
        
        # TAU 2
        self.out.pt_2[0]                       = event.Tau_pt[ditau.id2]
        self.out.eta_2[0]                      = event.Tau_eta[ditau.id2]
        self.out.phi_2[0]                      = event.Tau_phi[ditau.id2]
        self.out.mass_2[0]                     = event.Tau_mass[ditau.id2]
        self.out.dxy_2[0]                      = event.Tau_dxy[ditau.id2]
        self.out.dz_2[0]                       = event.Tau_dz[ditau.id2]         
        self.out.leadTkPtOverTauPt_2[0]        = event.Tau_leadTkPtOverTauPt[ditau.id2]
        self.out.chargedIso_2[0]               = event.Tau_chargedIso[ditau.id2]
        self.out.neutralIso_2[0]               = event.Tau_neutralIso[ditau.id2]
        self.out.photonsOutsideSignalCone_2[0] = event.Tau_photonsOutsideSignalCone[ditau.id2]
        self.out.puCorr_2[0]                   = event.Tau_puCorr[ditau.id2]
        self.out.rawAntiEle_2[0]               = event.Tau_rawAntiEle[ditau.id2]
        self.out.rawIso_2[0]                   = event.Tau_rawIso[ditau.id2]
        self.out.q_2[0]                        = event.Tau_charge[ditau.id2]
        self.out.decayMode_2[0]                = event.Tau_decayMode[ditau.id2]
        ###self.out.rawAntiEleCat_2[0]            = event.Tau_rawAntiEleCat[ditau.id2]
        self.out.idAntiEle_2[0]                = ord(event.Tau_idAntiEle[ditau.id2])
        self.out.idAntiMu_2[0]                 = ord(event.Tau_idAntiMu[ditau.id2])
        self.out.idDecayMode_2[0]              = event.Tau_idDecayMode[ditau.id2]
        self.out.idDecayModeNewDMs_2[0]        = event.Tau_idDecayModeNewDMs[ditau.id2]
        self.out.rawMVAoldDM_2[0]              = event.Tau_rawMVAoldDM[ditau.id2]
        self.out.rawMVAoldDM2017v1_2[0]        = event.Tau_rawMVAoldDM2017v1[ditau.id2]
        self.out.rawMVAoldDM2017v2_2[0]        = event.Tau_rawMVAoldDM2017v2[ditau.id2]
        self.out.rawMVAnewDM2017v2_2[0]        = event.Tau_rawMVAnewDM2017v2[ditau.id2]
        self.out.idMVAoldDM_2[0]               = ord(event.Tau_idMVAoldDM[ditau.id2])
        self.out.idMVAoldDM2017v1_2[0]         = ord(event.Tau_idMVAoldDM2017v1[ditau.id2])
        self.out.idMVAoldDM2017v2_2[0]         = ord(event.Tau_idMVAoldDM2017v2[ditau.id2])
        self.out.idMVAnewDM2017v2_2[0]         = ord(event.Tau_idMVAnewDM2017v2[ditau.id2])
        
        
        # GENERATOR 2
        #print type(event.Tau_genPartFlav[ditau.id2])
        if not self.isData:
          self.out.genPartFlav_2[0]            = ord(event.Tau_genPartFlav[ditau.id2])
          genvistau = Collection(event, 'GenVisTau')
          dRmax  = 1000
          gendm  = -1
          genpt  = -1
          geneta = -1
          genphi = -1
          for igvt in range(event.nGenVisTau):
            dR = genvistau[igvt].p4().DeltaR(tau2)
            if dR < 0.5 and dR < dRmax:
              dRmax  = dR
              gendm  = event.GenVisTau_status[igvt]
              genpt  = event.GenVisTau_pt[igvt]
              geneta = event.GenVisTau_eta[igvt]
              genphi = event.GenVisTau_phi[igvt]
          
          self.out.gendecayMode_2[0]           = gendm
          self.out.genvistaupt_2[0]            = genpt
          self.out.genvistaueta_2[0]           = geneta
          self.out.genvistauphi_2[0]           = genphi
        
        
        # EVENT
        self.out.isData[0]                     = self.isData
        self.out.run[0]                        = event.run
        self.out.luminosityBlock[0]            = event.luminosityBlock
        #print 'event =', event.event & 0xffffffffffffffff, 'original = ', event.event
        self.out.event[0]                      = event.event & 0xffffffffffffffff
        self.out.MET_pt[0]                     = event.MET_pt
        self.out.MET_phi[0]                    = event.MET_phi
        self.out.PuppiMET_pt[0]                = event.PuppiMET_pt
        self.out.PuppiMET_phi[0]               = event.PuppiMET_phi
        ###self.out.MET_significance[0]           = event.MET_significance
        ###self.out.MET_covXX[0]                  = event.MET_covXX
        ###self.out.MET_covXY[0]                  = event.MET_covXY
        ###self.out.MET_covYY[0]                  = event.MET_covYY
        ###self.out.fixedGridRhoFastjetAll[0]     = event.fixedGridRhoFastjetAll
        self.out.npvs[0]                       = event.PV_npvs
        self.out.npvsGood[0]                   = event.PV_npvsGood
        
        if not self.isData:
          self.out.ngentauhads[0]              = ngentauhads
          self.out.ngentaus[0]                 = ngentaus
          self.out.GenMET_pt[0]                = event.GenMET_pt
          self.out.GenMET_phi[0]               = event.GenMET_phi
          self.out.nPU[0]                      = event.Pileup_nPU
          self.out.nTrueInt[0]                 = event.Pileup_nTrueInt
          try:
            self.out.LHE_Njets[0]              = event.LHE_Njets
          except RuntimeError:
            self.out.LHE_Njets[0]              = -1
        #print 'check (LO)', event.LHE_NpLO, type(event.LHE_NpLO)
        #print 'check (NLO)', event.LHE_NpNLO, type(event.LHE_NpNLO)        
        #self.out.LHE_NpLO[0]                   = event.LHE_NpLO
        #self.out.LHE_NpNLO[0]                  = event.LHE_NpNLO
        #print self.out.LHE_Njets[0], event.LHE_Njets 
        #print self.out.event[0], event.event, (event.event & 0xffffffffffffffff)
        #print self.out.LHE_Njets[0], event.LHE_Njets, int(event.LHE_Njets)
        #print event.LHE_NpNLO
        #print self.out.Pileup_nPU, event.Pileup_nPU
        
        
        # JETS
        self.out.njets[0]                      = len(jetIds)
        self.out.njets50[0]                    = len([j for j in jetIds if event.Jet_pt[j]>50])
        self.out.nfjets[0]                     = nfjets
        self.out.ncjets[0]                     = ncjets
        self.out.nbtag[0]                      = nbtag
        
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
        
        if len(bjetIds)>0:
          self.out.bpt_1[0]                    = event.Jet_pt[bjetIds[0]]
          self.out.beta_1[0]                   = event.Jet_eta[bjetIds[0]]
        else:
          self.out.bpt_1[0]                    = -9.
          self.out.beta_1[0]                   = -9.
        
        if len(bjetIds)>1:
          self.out.bpt_2[0]                    = event.Jet_pt[bjetIds[1]]
          self.out.beta_2[0]                   = event.Jet_eta[bjetIds[1]]
        else:
          self.out.bpt_2[0]                    = -9.
          self.out.beta_2[0]                   = -9.
        
        self.out.pfmt_1[0]                     = math.sqrt( 2 * self.out.pt_1[0] * self.out.MET_pt[0] * ( 1 - math.cos(deltaPhi(self.out.phi_1[0], self.out.MET_phi[0])) ) );
        self.out.pfmt_2[0]                     = math.sqrt( 2 * self.out.pt_2[0] * self.out.MET_pt[0] * ( 1 - math.cos(deltaPhi(self.out.phi_2[0], self.out.MET_phi[0])) ) );
        
        self.out.m_vis[0]                      = (tau1 + tau2).M()
        self.out.pt_ll[0]                      = (tau1 + tau2).Pt()
        self.out.dR_ll[0]                      = tau1.DeltaR(tau2)
        self.out.dphi_ll[0]                    = deltaPhi(self.out.phi_1[0], self.out.phi_2[0])
        
        
        # PZETA
        leg1     = ROOT.TVector3(tau1.Px(), tau1.Py(), 0.)
        leg2     = ROOT.TVector3(tau2.Px(), tau2.Py(), 0.)
        met_tlv  = ROOT.TLorentzVector()
        met_tlv.SetPxPyPzE(self.out.MET_pt[0]*math.cos(self.out.MET_phi[0]), 
                           self.out.MET_pt[0]*math.sin(self.out.MET_phi[0]),
                           0,
                           self.out.MET_pt[0])
        metleg   = met_tlv.Vect()
        zetaAxis = ROOT.TVector3(leg1.Unit() + leg2.Unit()).Unit()
        pzetaVis = leg1*zetaAxis + leg2*zetaAxis
        pzetaMET = metleg*zetaAxis
        self.out.pzetamiss[0]  = pzetaMET
        self.out.pzetavis[0]   = pzetaVis
        self.out.pzeta_disc[0] = pzetaMET - 0.5*pzetaVis
        
        
        # VETOS
        self.out.extramuon_veto[0], self.out.extraelec_veto[0], self.out.dilepton_veto[0]  = extraLeptonVetos(event, [-1], [-1], self.name)
        
        
        # WEIGHTS
        if not self.isData:
          if self.doZpt:
            zboson = getZPTMass(event)
            self.out.m_genboson[0]    = zboson.M()
            self.out.pt_genboson[0]   = zboson.Pt()
            self.out.zptweight[0]     = self.recoilTool.getZptWeight(zboson.Pt(),zboson.M())
          diTauLeg1SF   = self.tauSFs.getTriggerSF(   self.out.pt_1, self.out.eta_1, self.out.phi_1 )
          diTauLeg2SF   = self.tauSFs.getTriggerSF(   self.out.pt_2, self.out.eta_2, self.out.phi_2 )
          diTauLeg1SFVT = self.tauSFsVT.getTriggerSF( self.out.pt_1, self.out.eta_1, self.out.phi_1 )
          diTauLeg2SFVT = self.tauSFsVT.getTriggerSF( self.out.pt_2, self.out.eta_2, self.out.phi_2 )
          self.out.genWeight[0]       = event.genWeight
          self.out.trigweight[0]      = diTauLeg1SF*diTauLeg2SF
          self.out.trigweightVT[0]    = diTauLeg1SFVT*diTauLeg2SFVT
          self.out.puweight[0]        = self.puTool.getWeight(event.Pileup_nTrueInt)
          self.out.idisoweight_1[0]   = self.ltfSFs.getSF(self.out.genPartFlav_1[0],self.out.eta_1[0])
          self.out.idisoweight_2[0]   = self.ltfSFs.getSF(self.out.genPartFlav_2[0],self.out.eta_2[0])
          self.out.btagweight[0]      = self.btagTool.getWeight(event,jetIds)
          self.out.btagweight_deep[0] = self.btagTool_deep.getWeight(event,jetIds)
          self.out.weight[0]          = self.out.genWeight[0]*self.out.puweight[0]*self.out.trigweight[0]*self.out.idisoweight_1[0]*self.out.idisoweight_2[0]
        
        
        self.out.tree.Fill() 
        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

#TauTauModule = lambda : TauTauProducer(jetSelection= lambda j : j.pt > 30) 
