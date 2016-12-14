//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec  9 10:09:17 2016 by ROOT version 6.06/08
// from TTree ProcessedTree/ProcessedTree
// found on file: Ntuples-Data-2016RunG-PromptReco-80Xpart4.root
//////////////////////////////////////////////////////////

#ifndef HeaderSMPJ_h
#define HeaderSMPJ_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "Math/GenVector/PxPyPzE4D.h"

class HeaderSMPJ {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxfilterIdList = 1;
   const Int_t kMaxCaloMet__et = 1;
   const Int_t kMaxCaloMet__CaloMetPt = 1;
   const Int_t kMaxCaloMet__sumEt = 1;
   const Int_t kMaxCaloMet__phi = 1;
   const Int_t kMaxPFMet__et = 1;
   const Int_t kMaxPFMet__CaloMetPt = 1;
   const Int_t kMaxPFMet__sumEt = 1;
   const Int_t kMaxPFMet__phi = 1;
   const Int_t kMaxMvaMet__et = 1;
   const Int_t kMaxMvaMet__CaloMetPt = 1;
   const Int_t kMaxMvaMet__sumEt = 1;
   const Int_t kMaxMvaMet__phi = 1;
   const Int_t kMaxTriggerDecision = 1;
   const Int_t kMaxtriggerList = 1;
   const Int_t kMaxL1Prescale = 1;
   const Int_t kMaxHLTPrescale = 1;
   const Int_t kMaxHLTObj = 1;
   const Int_t kMaxL1Obj = 1;
   const Int_t kMaxGenJets_ = 1;
   const Int_t kMaxCaloJets_ = 1;
   const Int_t kMaxCaloJets__genR = 1;
   const Int_t kMaxCaloJets__cor = 1;
   const Int_t kMaxCaloJets__jecLabels = 1;
   const Int_t kMaxCaloJets__unc = 1;
   const Int_t kMaxCaloJets__uncSrc = 1;
   const Int_t kMaxCaloJets__area = 1;
   const Int_t kMaxCaloJets__looseID = 1;
   const Int_t kMaxCaloJets__tightID = 1;
   const Int_t kMaxCaloJets__emf = 1;
   const Int_t kMaxCaloJets__fHPD = 1;
   const Int_t kMaxCaloJets__fRBX = 1;
   const Int_t kMaxCaloJets__n90hits = 1;
   const Int_t kMaxCaloJets__nTrkCalo = 1;
   const Int_t kMaxCaloJets__nTrkVtx = 1;
   const Int_t kMaxPFJets_ = 35;
   const Int_t kMaxPFJets__genR = 35;
   const Int_t kMaxPFJets__cor = 35;
   const Int_t kMaxPFJets__jecLabels = 35;
   const Int_t kMaxPFJets__unc = 35;
   const Int_t kMaxPFJets__uncSrc = 35;
   const Int_t kMaxPFJets__area = 35;
   const Int_t kMaxPFJets__looseID = 35;
   const Int_t kMaxPFJets__tightID = 35;
   const Int_t kMaxPFJets__TCHE = 35;
   const Int_t kMaxPFJets__TCHP = 35;
   const Int_t kMaxPFJets__TCHEpf = 35;
   const Int_t kMaxPFJets__TCHPpf = 35;
   const Int_t kMaxPFJets__SoftMuonTagByIP = 35;
   const Int_t kMaxPFJets__SoftElectronTagByIP = 35;
   const Int_t kMaxPFJets__SoftMuonTag = 35;
   const Int_t kMaxPFJets__SoftElectronTag = 35;
   const Int_t kMaxPFJets__SimpleSecVertexHE = 35;
   const Int_t kMaxPFJets__SimpleSecVertexHP = 35;
   const Int_t kMaxPFJets__SimpleSecVertexHEpf = 35;
   const Int_t kMaxPFJets__SimpleSecVertexHPpf = 35;
   const Int_t kMaxPFJets__CSV = 35;
   const Int_t kMaxPFJets__CSVpf = 35;
   const Int_t kMaxPFJets__CinclSVpf = 35;
   const Int_t kMaxPFJets__CMVApf = 35;
   const Int_t kMaxPFJets__CSVSoftLeptonpf = 35;
   const Int_t kMaxPFJets__CSVpfPositive = 35;
   const Int_t kMaxPFJets__CSVpfNegative = 35;
   const Int_t kMaxPFJets__QGtagger = 35;
   const Int_t kMaxPFJets__partonFlavour = 35;
   const Int_t kMaxPFJets__hadronFlavour = 35;
   const Int_t kMaxPFJets__recommend1 = 35;
   const Int_t kMaxPFJets__recommend2 = 35;
   const Int_t kMaxPFJets__recommend3 = 35;
   const Int_t kMaxPFJets__chf = 35;
   const Int_t kMaxPFJets__nhf = 35;
   const Int_t kMaxPFJets__nemf = 35;
   const Int_t kMaxPFJets__cemf = 35;
   const Int_t kMaxPFJets__muf = 35;
   const Int_t kMaxPFJets__hf_hf = 35;
   const Int_t kMaxPFJets__hf_phf = 35;
   const Int_t kMaxPFJets__hf_hm = 35;
   const Int_t kMaxPFJets__hf_phm = 35;
   const Int_t kMaxPFJets__chm = 35;
   const Int_t kMaxPFJets__nhm = 35;
   const Int_t kMaxPFJets__phm = 35;
   const Int_t kMaxPFJets__elm = 35;
   const Int_t kMaxPFJets__mum = 35;
   const Int_t kMaxPFJets__ncand = 35;
   const Int_t kMaxPFJets__beta = 35;
   const Int_t kMaxPFJets__betaStar = 35;
   const Int_t kMaxPFJets__mpuTrk = 35;
   const Int_t kMaxPFJets__mlvTrk = 35;
   const Int_t kMaxPFJets__mjtTrk = 35;
   const Int_t kMaxPFJets__hof = 35;
   const Int_t kMaxPFJets__pujid = 35;
   const Int_t kMaxPFJets__calojetpt = 35;
   const Int_t kMaxPFJets__calojetef = 35;
   const Int_t kMaxPFJetsCHS_ = 61;
   const Int_t kMaxPFJetsCHS__genR = 61;
   const Int_t kMaxPFJetsCHS__cor = 61;
   const Int_t kMaxPFJetsCHS__jecLabels = 61;
   const Int_t kMaxPFJetsCHS__unc = 61;
   const Int_t kMaxPFJetsCHS__uncSrc = 61;
   const Int_t kMaxPFJetsCHS__area = 61;
   const Int_t kMaxPFJetsCHS__looseID = 61;
   const Int_t kMaxPFJetsCHS__tightID = 61;
   const Int_t kMaxPFJetsCHS__TCHE = 61;
   const Int_t kMaxPFJetsCHS__TCHP = 61;
   const Int_t kMaxPFJetsCHS__TCHEpf = 61;
   const Int_t kMaxPFJetsCHS__TCHPpf = 61;
   const Int_t kMaxPFJetsCHS__SoftMuonTagByIP = 61;
   const Int_t kMaxPFJetsCHS__SoftElectronTagByIP = 61;
   const Int_t kMaxPFJetsCHS__SoftMuonTag = 61;
   const Int_t kMaxPFJetsCHS__SoftElectronTag = 61;
   const Int_t kMaxPFJetsCHS__SimpleSecVertexHE = 61;
   const Int_t kMaxPFJetsCHS__SimpleSecVertexHP = 61;
   const Int_t kMaxPFJetsCHS__SimpleSecVertexHEpf = 61;
   const Int_t kMaxPFJetsCHS__SimpleSecVertexHPpf = 61;
   const Int_t kMaxPFJetsCHS__CSV = 61;
   const Int_t kMaxPFJetsCHS__CSVpf = 61;
   const Int_t kMaxPFJetsCHS__CinclSVpf = 61;
   const Int_t kMaxPFJetsCHS__CMVApf = 61;
   const Int_t kMaxPFJetsCHS__CSVSoftLeptonpf = 61;
   const Int_t kMaxPFJetsCHS__CSVpfPositive = 61;
   const Int_t kMaxPFJetsCHS__CSVpfNegative = 61;
   const Int_t kMaxPFJetsCHS__QGtagger = 61;
   const Int_t kMaxPFJetsCHS__partonFlavour = 61;
   const Int_t kMaxPFJetsCHS__hadronFlavour = 61;
   const Int_t kMaxPFJetsCHS__recommend1 = 61;
   const Int_t kMaxPFJetsCHS__recommend2 = 61;
   const Int_t kMaxPFJetsCHS__recommend3 = 61;
   const Int_t kMaxPFJetsCHS__chf = 61;
   const Int_t kMaxPFJetsCHS__nhf = 61;
   const Int_t kMaxPFJetsCHS__nemf = 61;
   const Int_t kMaxPFJetsCHS__cemf = 61;
   const Int_t kMaxPFJetsCHS__muf = 61;
   const Int_t kMaxPFJetsCHS__hf_hf = 61;
   const Int_t kMaxPFJetsCHS__hf_phf = 61;
   const Int_t kMaxPFJetsCHS__hf_hm = 61;
   const Int_t kMaxPFJetsCHS__hf_phm = 61;
   const Int_t kMaxPFJetsCHS__chm = 61;
   const Int_t kMaxPFJetsCHS__nhm = 61;
   const Int_t kMaxPFJetsCHS__phm = 61;
   const Int_t kMaxPFJetsCHS__elm = 61;
   const Int_t kMaxPFJetsCHS__mum = 61;
   const Int_t kMaxPFJetsCHS__ncand = 61;
   const Int_t kMaxPFJetsCHS__beta = 61;
   const Int_t kMaxPFJetsCHS__betaStar = 61;
   const Int_t kMaxPFJetsCHS__mpuTrk = 61;
   const Int_t kMaxPFJetsCHS__mlvTrk = 61;
   const Int_t kMaxPFJetsCHS__mjtTrk = 61;
   const Int_t kMaxPFJetsCHS__hof = 61;
   const Int_t kMaxPFJetsCHS__pujid = 61;
   const Int_t kMaxPFJetsCHS__calojetpt = 61;
   const Int_t kMaxPFJetsCHS__calojetef = 61;
   const Int_t kMaxgenFlavour = 1;
   const Int_t kMaxgenFlavourHadron = 1;
   const Int_t kMaxmMuon_ = 1;
   const Int_t kMaxmMuon__genR = 1;
   const Int_t kMaxmMuon__cor = 1;
   const Int_t kMaxmMuon__jecLabels = 1;
   const Int_t kMaxmMuon__unc = 1;
   const Int_t kMaxmMuon__uncSrc = 1;
   const Int_t kMaxmMuon__area = 1;
   const Int_t kMaxmMuon__DxyVertex = 1;
   const Int_t kMaxmMuon__DzVertex = 1;
   const Int_t kMaxmMuon__PDGID = 1;
   const Int_t kMaxmMuon__PfIso = 1;
   const Int_t kMaxmMuon__PuChargedHadronIso = 1;
   const Int_t kMaxmMuon__ChargedHadronIso = 1;
   const Int_t kMaxmMuon__NeutralHadronIso = 1;
   const Int_t kMaxmMuon__PhotonIso = 1;
   const Int_t kMaxmElectron_ = 1;
   const Int_t kMaxmElectron__genR = 1;
   const Int_t kMaxmElectron__cor = 1;
   const Int_t kMaxmElectron__jecLabels = 1;
   const Int_t kMaxmElectron__unc = 1;
   const Int_t kMaxmElectron__uncSrc = 1;
   const Int_t kMaxmElectron__area = 1;
   const Int_t kMaxmElectron__DxyVertex = 1;
   const Int_t kMaxmElectron__DzVertex = 1;
   const Int_t kMaxmElectron__PDGID = 1;
   const Int_t kMaxmElectron__PfIso = 1;
   const Int_t kMaxmElectron__PuChargedHadronIso = 1;
   const Int_t kMaxmElectron__ChargedHadronIso = 1;
   const Int_t kMaxmElectron__NeutralHadronIso = 1;
   const Int_t kMaxmElectron__PhotonIso = 1;

   // Declaration of leaf types
 //QCDEvent        *events;
   vector<vector<int> > filterIdList_;
   Bool_t          EvtHdr__mIsPVgood;
   Bool_t          EvtHdr__mHCALNoise;
   Bool_t          EvtHdr__mHCALNoiseNoMinZ;
   Int_t           EvtHdr__mRun;
   Int_t           EvtHdr__mEvent;
   Int_t           EvtHdr__mLumi;
   Int_t           EvtHdr__mBunch;
   Int_t           EvtHdr__mNVtx;
   Int_t           EvtHdr__mNVtxGood;
   Int_t           EvtHdr__mOOTPUEarly;
   Int_t           EvtHdr__mOOTPULate;
   Int_t           EvtHdr__mINTPU;
   Int_t           EvtHdr__mNBX;
   Float_t         EvtHdr__mPVndof;
   Float_t         EvtHdr__mTrPu;
   Float_t         EvtHdr__mPVx;
   Float_t         EvtHdr__mPVy;
   Float_t         EvtHdr__mPVz;
   Float_t         EvtHdr__mBSx;
   Float_t         EvtHdr__mBSy;
   Float_t         EvtHdr__mBSz;
   Float_t         EvtHdr__mPthat;
   Float_t         EvtHdr__mWeight;
   Float_t         EvtHdr__mCaloRho;
   Float_t         EvtHdr__mPFRho;
   Float_t         CaloMet__et_;
   Float_t         CaloMet__CaloMetPt_;
   Float_t         CaloMet__sumEt_;
   Float_t         CaloMet__phi_;
   Float_t         PFMet__et_;
   Float_t         PFMet__CaloMetPt_;
   Float_t         PFMet__sumEt_;
   Float_t         PFMet__phi_;
   Float_t         MvaMet__et_;
   Float_t         MvaMet__CaloMetPt_;
   Float_t         MvaMet__sumEt_;
   Float_t         MvaMet__phi_;
   vector<int>     TriggerDecision_;
   vector<string>  triggerList_;
   vector<int>     L1Prescale_;
   vector<int>     HLTPrescale_;
 //vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > HLTObj_;
 //vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > L1Obj_;
   Int_t           GenJets__;
   Double_t        GenJets__fCoordinates_fX[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fY[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fZ[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fT[kMaxGenJets_];   //[GenJets__]
   Int_t           CaloJets__;
   Double_t        CaloJets__P4__fCoordinates_fX[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__P4__fCoordinates_fY[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__P4__fCoordinates_fZ[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__P4__fCoordinates_fT[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fX[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fY[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fZ[kMaxCaloJets_];   //[CaloJets__]
   Double_t        CaloJets__genP4__fCoordinates_fT[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__genR_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__cor_[kMaxCaloJets_];   //[CaloJets__]
   vector<double>  CaloJets__jecLabels_[kMaxCaloJets_];
   Float_t         CaloJets__unc_[kMaxCaloJets_];   //[CaloJets__]
   vector<float>   CaloJets__uncSrc_[kMaxCaloJets_];
   Float_t         CaloJets__area_[kMaxCaloJets_];   //[CaloJets__]
   Bool_t          CaloJets__looseID_[kMaxCaloJets_];   //[CaloJets__]
   Bool_t          CaloJets__tightID_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__emf_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__fHPD_[kMaxCaloJets_];   //[CaloJets__]
   Float_t         CaloJets__fRBX_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           CaloJets__n90hits_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           CaloJets__nTrkCalo_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           CaloJets__nTrkVtx_[kMaxCaloJets_];   //[CaloJets__]
   Int_t           PFJets__;
   Double_t        PFJets__P4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__P4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__P4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__P4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fX[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fY[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fZ[kMaxPFJets_];   //[PFJets__]
   Double_t        PFJets__genP4__fCoordinates_fT[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__genR_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__cor_[kMaxPFJets_];   //[PFJets__]
   vector<double>  PFJets__jecLabels_[kMaxPFJets_];
   Float_t         PFJets__unc_[kMaxPFJets_];   //[PFJets__]
   vector<float>   PFJets__uncSrc_[kMaxPFJets_];
   Float_t         PFJets__area_[kMaxPFJets_];   //[PFJets__]
   Bool_t          PFJets__looseID_[kMaxPFJets_];   //[PFJets__]
   Bool_t          PFJets__tightID_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__TCHE_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__TCHP_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__TCHEpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__TCHPpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SoftMuonTagByIP_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SoftElectronTagByIP_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SoftMuonTag_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SoftElectronTag_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SimpleSecVertexHE_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SimpleSecVertexHP_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SimpleSecVertexHEpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__SimpleSecVertexHPpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CSV_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CSVpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CinclSVpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CMVApf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CSVSoftLeptonpf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CSVpfPositive_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__CSVpfNegative_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__QGtagger_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__partonFlavour_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hadronFlavour_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__recommend1_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__recommend2_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__recommend3_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__chf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__nhf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__nemf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__cemf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__muf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hf_hf_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hf_phf_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__hf_hm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__hf_phm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__chm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__nhm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__phm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__elm_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mum_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__ncand_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__beta_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__betaStar_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mpuTrk_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mlvTrk_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJets__mjtTrk_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__hof_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__pujid_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__calojetpt_[kMaxPFJets_];   //[PFJets__]
   Float_t         PFJets__calojetef_[kMaxPFJets_];   //[PFJets__]
   Int_t           PFJetsCHS__;
   Double_t        PFJetsCHS__P4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__genR_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__cor_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   vector<double>  PFJetsCHS__jecLabels_[kMaxPFJetsCHS_];
   Float_t         PFJetsCHS__unc_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   vector<float>   PFJetsCHS__uncSrc_[kMaxPFJetsCHS_];
   Float_t         PFJetsCHS__area_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Bool_t          PFJetsCHS__looseID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Bool_t          PFJetsCHS__tightID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__TCHE_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__TCHP_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__TCHEpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__TCHPpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SoftMuonTagByIP_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SoftElectronTagByIP_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SoftMuonTag_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SoftElectronTag_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SimpleSecVertexHE_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SimpleSecVertexHP_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SimpleSecVertexHEpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__SimpleSecVertexHPpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CSV_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CSVpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CinclSVpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CMVApf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CSVSoftLeptonpf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CSVpfPositive_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__CSVpfNegative_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__QGtagger_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__partonFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hadronFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__recommend1_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__recommend2_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__recommend3_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__chf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__nhf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__nemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__cemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__muf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hf_hf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hf_phf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hf_hm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hf_phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__chm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__nhm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__elm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mum_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__ncand_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__beta_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__betaStar_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mpuTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mlvTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mjtTrk_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hof_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__pujid_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__calojetpt_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__calojetef_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   vector<float>   genFlavour_;
   vector<float>   genFlavourHadron_;
   Int_t           mMuon__;
   Double_t        mMuon__P4__fCoordinates_fX[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__P4__fCoordinates_fY[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__P4__fCoordinates_fZ[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__P4__fCoordinates_fT[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__genP4__fCoordinates_fX[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__genP4__fCoordinates_fY[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__genP4__fCoordinates_fZ[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__genP4__fCoordinates_fT[kMaxmMuon_];   //[mMuon__]
   Float_t         mMuon__genR_[kMaxmMuon_];   //[mMuon__]
   Float_t         mMuon__cor_[kMaxmMuon_];   //[mMuon__]
   vector<double>  mMuon__jecLabels_[kMaxmMuon_];
   Float_t         mMuon__unc_[kMaxmMuon_];   //[mMuon__]
   vector<float>   mMuon__uncSrc_[kMaxmMuon_];
   Float_t         mMuon__area_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__DxyVertex_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__DzVertex_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__PDGID_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__PfIso_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__PuChargedHadronIso_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__ChargedHadronIso_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__NeutralHadronIso_[kMaxmMuon_];   //[mMuon__]
   Double_t        mMuon__PhotonIso_[kMaxmMuon_];   //[mMuon__]
   Int_t           mElectron__;
   Double_t        mElectron__P4__fCoordinates_fX[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__P4__fCoordinates_fY[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__P4__fCoordinates_fZ[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__P4__fCoordinates_fT[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__genP4__fCoordinates_fX[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__genP4__fCoordinates_fY[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__genP4__fCoordinates_fZ[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__genP4__fCoordinates_fT[kMaxmElectron_];   //[mElectron__]
   Float_t         mElectron__genR_[kMaxmElectron_];   //[mElectron__]
   Float_t         mElectron__cor_[kMaxmElectron_];   //[mElectron__]
   vector<double>  mElectron__jecLabels_[kMaxmElectron_];
   Float_t         mElectron__unc_[kMaxmElectron_];   //[mElectron__]
   vector<float>   mElectron__uncSrc_[kMaxmElectron_];
   Float_t         mElectron__area_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__DxyVertex_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__DzVertex_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__PDGID_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__PfIso_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__PuChargedHadronIso_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__ChargedHadronIso_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__NeutralHadronIso_[kMaxmElectron_];   //[mElectron__]
   Double_t        mElectron__PhotonIso_[kMaxmElectron_];   //[mElectron__]

   // List of branches
   TBranch        *b_events_filterIdList_;   //!
   TBranch        *b_events_EvtHdr__mIsPVgood;   //!
   TBranch        *b_events_EvtHdr__mHCALNoise;   //!
   TBranch        *b_events_EvtHdr__mHCALNoiseNoMinZ;   //!
   TBranch        *b_events_EvtHdr__mRun;   //!
   TBranch        *b_events_EvtHdr__mEvent;   //!
   TBranch        *b_events_EvtHdr__mLumi;   //!
   TBranch        *b_events_EvtHdr__mBunch;   //!
   TBranch        *b_events_EvtHdr__mNVtx;   //!
   TBranch        *b_events_EvtHdr__mNVtxGood;   //!
   TBranch        *b_events_EvtHdr__mOOTPUEarly;   //!
   TBranch        *b_events_EvtHdr__mOOTPULate;   //!
   TBranch        *b_events_EvtHdr__mINTPU;   //!
   TBranch        *b_events_EvtHdr__mNBX;   //!
   TBranch        *b_events_EvtHdr__mPVndof;   //!
   TBranch        *b_events_EvtHdr__mTrPu;   //!
   TBranch        *b_events_EvtHdr__mPVx;   //!
   TBranch        *b_events_EvtHdr__mPVy;   //!
   TBranch        *b_events_EvtHdr__mPVz;   //!
   TBranch        *b_events_EvtHdr__mBSx;   //!
   TBranch        *b_events_EvtHdr__mBSy;   //!
   TBranch        *b_events_EvtHdr__mBSz;   //!
   TBranch        *b_events_EvtHdr__mPthat;   //!
   TBranch        *b_events_EvtHdr__mWeight;   //!
   TBranch        *b_events_EvtHdr__mCaloRho;   //!
   TBranch        *b_events_EvtHdr__mPFRho;   //!
   TBranch        *b_events_CaloMet__et_;   //!
   TBranch        *b_events_CaloMet__CaloMetPt_;   //!
   TBranch        *b_events_CaloMet__sumEt_;   //!
   TBranch        *b_events_CaloMet__phi_;   //!
   TBranch        *b_events_PFMet__et_;   //!
   TBranch        *b_events_PFMet__CaloMetPt_;   //!
   TBranch        *b_events_PFMet__sumEt_;   //!
   TBranch        *b_events_PFMet__phi_;   //!
   TBranch        *b_events_MvaMet__et_;   //!
   TBranch        *b_events_MvaMet__CaloMetPt_;   //!
   TBranch        *b_events_MvaMet__sumEt_;   //!
   TBranch        *b_events_MvaMet__phi_;   //!
   TBranch        *b_events_TriggerDecision_;   //!
   TBranch        *b_events_triggerList_;   //!
   TBranch        *b_events_L1Prescale_;   //!
   TBranch        *b_events_HLTPrescale_;   //!
   TBranch        *b_events_GenJets__;   //!
   TBranch        *b_GenJets__fCoordinates_fX;   //!
   TBranch        *b_GenJets__fCoordinates_fY;   //!
   TBranch        *b_GenJets__fCoordinates_fZ;   //!
   TBranch        *b_GenJets__fCoordinates_fT;   //!
   TBranch        *b_events_CaloJets__;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fX;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fY;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fZ;   //!
   TBranch        *b_CaloJets__P4__fCoordinates_fT;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fX;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fY;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fZ;   //!
   TBranch        *b_CaloJets__genP4__fCoordinates_fT;   //!
   TBranch        *b_CaloJets__genR_;   //!
   TBranch        *b_CaloJets__cor_;   //!
   TBranch        *b_CaloJets__jecLabels_;   //!
   TBranch        *b_CaloJets__unc_;   //!
   TBranch        *b_CaloJets__uncSrc_;   //!
   TBranch        *b_CaloJets__area_;   //!
   TBranch        *b_CaloJets__looseID_;   //!
   TBranch        *b_CaloJets__tightID_;   //!
   TBranch        *b_CaloJets__emf_;   //!
   TBranch        *b_CaloJets__fHPD_;   //!
   TBranch        *b_CaloJets__fRBX_;   //!
   TBranch        *b_CaloJets__n90hits_;   //!
   TBranch        *b_CaloJets__nTrkCalo_;   //!
   TBranch        *b_CaloJets__nTrkVtx_;   //!
   TBranch        *b_events_PFJets__;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fX;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fY;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fZ;   //!
   TBranch        *b_PFJets__P4__fCoordinates_fT;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fX;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fY;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fZ;   //!
   TBranch        *b_PFJets__genP4__fCoordinates_fT;   //!
   TBranch        *b_PFJets__genR_;   //!
   TBranch        *b_PFJets__cor_;   //!
   TBranch        *b_PFJets__jecLabels_;   //!
   TBranch        *b_PFJets__unc_;   //!
   TBranch        *b_PFJets__uncSrc_;   //!
   TBranch        *b_PFJets__area_;   //!
   TBranch        *b_PFJets__looseID_;   //!
   TBranch        *b_PFJets__tightID_;   //!
   TBranch        *b_PFJets__TCHE_;   //!
   TBranch        *b_PFJets__TCHP_;   //!
   TBranch        *b_PFJets__TCHEpf_;   //!
   TBranch        *b_PFJets__TCHPpf_;   //!
   TBranch        *b_PFJets__SoftMuonTagByIP_;   //!
   TBranch        *b_PFJets__SoftElectronTagByIP_;   //!
   TBranch        *b_PFJets__SoftMuonTag_;   //!
   TBranch        *b_PFJets__SoftElectronTag_;   //!
   TBranch        *b_PFJets__SimpleSecVertexHE_;   //!
   TBranch        *b_PFJets__SimpleSecVertexHP_;   //!
   TBranch        *b_PFJets__SimpleSecVertexHEpf_;   //!
   TBranch        *b_PFJets__SimpleSecVertexHPpf_;   //!
   TBranch        *b_PFJets__CSV_;   //!
   TBranch        *b_PFJets__CSVpf_;   //!
   TBranch        *b_PFJets__CinclSVpf_;   //!
   TBranch        *b_PFJets__CMVApf_;   //!
   TBranch        *b_PFJets__CSVSoftLeptonpf_;   //!
   TBranch        *b_PFJets__CSVpfPositive_;   //!
   TBranch        *b_PFJets__CSVpfNegative_;   //!
   TBranch        *b_PFJets__QGtagger_;   //!
   TBranch        *b_PFJets__partonFlavour_;   //!
   TBranch        *b_PFJets__hadronFlavour_;   //!
   TBranch        *b_PFJets__recommend1_;   //!
   TBranch        *b_PFJets__recommend2_;   //!
   TBranch        *b_PFJets__recommend3_;   //!
   TBranch        *b_PFJets__chf_;   //!
   TBranch        *b_PFJets__nhf_;   //!
   TBranch        *b_PFJets__nemf_;   //!
   TBranch        *b_PFJets__cemf_;   //!
   TBranch        *b_PFJets__muf_;   //!
   TBranch        *b_PFJets__hf_hf_;   //!
   TBranch        *b_PFJets__hf_phf_;   //!
   TBranch        *b_PFJets__hf_hm_;   //!
   TBranch        *b_PFJets__hf_phm_;   //!
   TBranch        *b_PFJets__chm_;   //!
   TBranch        *b_PFJets__nhm_;   //!
   TBranch        *b_PFJets__phm_;   //!
   TBranch        *b_PFJets__elm_;   //!
   TBranch        *b_PFJets__mum_;   //!
   TBranch        *b_PFJets__ncand_;   //!
   TBranch        *b_PFJets__beta_;   //!
   TBranch        *b_PFJets__betaStar_;   //!
   TBranch        *b_PFJets__mpuTrk_;   //!
   TBranch        *b_PFJets__mlvTrk_;   //!
   TBranch        *b_PFJets__mjtTrk_;   //!
   TBranch        *b_PFJets__hof_;   //!
   TBranch        *b_PFJets__pujid_;   //!
   TBranch        *b_PFJets__calojetpt_;   //!
   TBranch        *b_PFJets__calojetef_;   //!
   TBranch        *b_events_PFJetsCHS__;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fX;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fY;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fZ;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fT;   //!
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fX;   //!
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fY;   //!
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fZ;   //!
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fT;   //!
   TBranch        *b_PFJetsCHS__genR_;   //!
   TBranch        *b_PFJetsCHS__cor_;   //!
   TBranch        *b_PFJetsCHS__jecLabels_;   //!
   TBranch        *b_PFJetsCHS__unc_;   //!
   TBranch        *b_PFJetsCHS__uncSrc_;   //!
   TBranch        *b_PFJetsCHS__area_;   //!
   TBranch        *b_PFJetsCHS__looseID_;   //!
   TBranch        *b_PFJetsCHS__tightID_;   //!
   TBranch        *b_PFJetsCHS__TCHE_;   //!
   TBranch        *b_PFJetsCHS__TCHP_;   //!
   TBranch        *b_PFJetsCHS__TCHEpf_;   //!
   TBranch        *b_PFJetsCHS__TCHPpf_;   //!
   TBranch        *b_PFJetsCHS__SoftMuonTagByIP_;   //!
   TBranch        *b_PFJetsCHS__SoftElectronTagByIP_;   //!
   TBranch        *b_PFJetsCHS__SoftMuonTag_;   //!
   TBranch        *b_PFJetsCHS__SoftElectronTag_;   //!
   TBranch        *b_PFJetsCHS__SimpleSecVertexHE_;   //!
   TBranch        *b_PFJetsCHS__SimpleSecVertexHP_;   //!
   TBranch        *b_PFJetsCHS__SimpleSecVertexHEpf_;   //!
   TBranch        *b_PFJetsCHS__SimpleSecVertexHPpf_;   //!
   TBranch        *b_PFJetsCHS__CSV_;   //!
   TBranch        *b_PFJetsCHS__CSVpf_;   //!
   TBranch        *b_PFJetsCHS__CinclSVpf_;   //!
   TBranch        *b_PFJetsCHS__CMVApf_;   //!
   TBranch        *b_PFJetsCHS__CSVSoftLeptonpf_;   //!
   TBranch        *b_PFJetsCHS__CSVpfPositive_;   //!
   TBranch        *b_PFJetsCHS__CSVpfNegative_;   //!
   TBranch        *b_PFJetsCHS__QGtagger_;   //!
   TBranch        *b_PFJetsCHS__partonFlavour_;   //!
   TBranch        *b_PFJetsCHS__hadronFlavour_;   //!
   TBranch        *b_PFJetsCHS__recommend1_;   //!
   TBranch        *b_PFJetsCHS__recommend2_;   //!
   TBranch        *b_PFJetsCHS__recommend3_;   //!
   TBranch        *b_PFJetsCHS__chf_;   //!
   TBranch        *b_PFJetsCHS__nhf_;   //!
   TBranch        *b_PFJetsCHS__nemf_;   //!
   TBranch        *b_PFJetsCHS__cemf_;   //!
   TBranch        *b_PFJetsCHS__muf_;   //!
   TBranch        *b_PFJetsCHS__hf_hf_;   //!
   TBranch        *b_PFJetsCHS__hf_phf_;   //!
   TBranch        *b_PFJetsCHS__hf_hm_;   //!
   TBranch        *b_PFJetsCHS__hf_phm_;   //!
   TBranch        *b_PFJetsCHS__chm_;   //!
   TBranch        *b_PFJetsCHS__nhm_;   //!
   TBranch        *b_PFJetsCHS__phm_;   //!
   TBranch        *b_PFJetsCHS__elm_;   //!
   TBranch        *b_PFJetsCHS__mum_;   //!
   TBranch        *b_PFJetsCHS__ncand_;   //!
   TBranch        *b_PFJetsCHS__beta_;   //!
   TBranch        *b_PFJetsCHS__betaStar_;   //!
   TBranch        *b_PFJetsCHS__mpuTrk_;   //!
   TBranch        *b_PFJetsCHS__mlvTrk_;   //!
   TBranch        *b_PFJetsCHS__mjtTrk_;   //!
   TBranch        *b_PFJetsCHS__hof_;   //!
   TBranch        *b_PFJetsCHS__pujid_;   //!
   TBranch        *b_PFJetsCHS__calojetpt_;   //!
   TBranch        *b_PFJetsCHS__calojetef_;   //!
   TBranch        *b_events_genFlavour_;   //!
   TBranch        *b_events_genFlavourHadron_;   //!
   TBranch        *b_events_mMuon__;   //!
   TBranch        *b_mMuon__P4__fCoordinates_fX;   //!
   TBranch        *b_mMuon__P4__fCoordinates_fY;   //!
   TBranch        *b_mMuon__P4__fCoordinates_fZ;   //!
   TBranch        *b_mMuon__P4__fCoordinates_fT;   //!
   TBranch        *b_mMuon__genP4__fCoordinates_fX;   //!
   TBranch        *b_mMuon__genP4__fCoordinates_fY;   //!
   TBranch        *b_mMuon__genP4__fCoordinates_fZ;   //!
   TBranch        *b_mMuon__genP4__fCoordinates_fT;   //!
   TBranch        *b_mMuon__genR_;   //!
   TBranch        *b_mMuon__cor_;   //!
   TBranch        *b_mMuon__jecLabels_;   //!
   TBranch        *b_mMuon__unc_;   //!
   TBranch        *b_mMuon__uncSrc_;   //!
   TBranch        *b_mMuon__area_;   //!
   TBranch        *b_mMuon__DxyVertex_;   //!
   TBranch        *b_mMuon__DzVertex_;   //!
   TBranch        *b_mMuon__PDGID_;   //!
   TBranch        *b_mMuon__PfIso_;   //!
   TBranch        *b_mMuon__PuChargedHadronIso_;   //!
   TBranch        *b_mMuon__ChargedHadronIso_;   //!
   TBranch        *b_mMuon__NeutralHadronIso_;   //!
   TBranch        *b_mMuon__PhotonIso_;   //!
   TBranch        *b_events_mElectron__;   //!
   TBranch        *b_mElectron__P4__fCoordinates_fX;   //!
   TBranch        *b_mElectron__P4__fCoordinates_fY;   //!
   TBranch        *b_mElectron__P4__fCoordinates_fZ;   //!
   TBranch        *b_mElectron__P4__fCoordinates_fT;   //!
   TBranch        *b_mElectron__genP4__fCoordinates_fX;   //!
   TBranch        *b_mElectron__genP4__fCoordinates_fY;   //!
   TBranch        *b_mElectron__genP4__fCoordinates_fZ;   //!
   TBranch        *b_mElectron__genP4__fCoordinates_fT;   //!
   TBranch        *b_mElectron__genR_;   //!
   TBranch        *b_mElectron__cor_;   //!
   TBranch        *b_mElectron__jecLabels_;   //!
   TBranch        *b_mElectron__unc_;   //!
   TBranch        *b_mElectron__uncSrc_;   //!
   TBranch        *b_mElectron__area_;   //!
   TBranch        *b_mElectron__DxyVertex_;   //!
   TBranch        *b_mElectron__DzVertex_;   //!
   TBranch        *b_mElectron__PDGID_;   //!
   TBranch        *b_mElectron__PfIso_;   //!
   TBranch        *b_mElectron__PuChargedHadronIso_;   //!
   TBranch        *b_mElectron__ChargedHadronIso_;   //!
   TBranch        *b_mElectron__NeutralHadronIso_;   //!
   TBranch        *b_mElectron__PhotonIso_;   //!

   HeaderSMPJ(TTree *tree=0);
   virtual ~HeaderSMPJ();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HeaderSMPJ_cxx
HeaderSMPJ::HeaderSMPJ(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ntuples-Data-2016RunG-PromptReco-80Xpart4.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Ntuples-Data-2016RunG-PromptReco-80Xpart4.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Ntuples-Data-2016RunG-PromptReco-80Xpart4.root:/ak4");
      dir->GetObject("ProcessedTree",tree);

   }
   Init(tree);
}

HeaderSMPJ::~HeaderSMPJ()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HeaderSMPJ::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HeaderSMPJ::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HeaderSMPJ::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("filterIdList_", &filterIdList_, &b_events_filterIdList_);
   fChain->SetBranchAddress("EvtHdr_.mIsPVgood", &EvtHdr__mIsPVgood, &b_events_EvtHdr__mIsPVgood);
   fChain->SetBranchAddress("EvtHdr_.mHCALNoise", &EvtHdr__mHCALNoise, &b_events_EvtHdr__mHCALNoise);
   fChain->SetBranchAddress("EvtHdr_.mHCALNoiseNoMinZ", &EvtHdr__mHCALNoiseNoMinZ, &b_events_EvtHdr__mHCALNoiseNoMinZ);
   fChain->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun, &b_events_EvtHdr__mRun);
   fChain->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent, &b_events_EvtHdr__mEvent);
   fChain->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi, &b_events_EvtHdr__mLumi);
   fChain->SetBranchAddress("EvtHdr_.mBunch", &EvtHdr__mBunch, &b_events_EvtHdr__mBunch);
   fChain->SetBranchAddress("EvtHdr_.mNVtx", &EvtHdr__mNVtx, &b_events_EvtHdr__mNVtx);
   fChain->SetBranchAddress("EvtHdr_.mNVtxGood", &EvtHdr__mNVtxGood, &b_events_EvtHdr__mNVtxGood);
   fChain->SetBranchAddress("EvtHdr_.mOOTPUEarly", &EvtHdr__mOOTPUEarly, &b_events_EvtHdr__mOOTPUEarly);
   fChain->SetBranchAddress("EvtHdr_.mOOTPULate", &EvtHdr__mOOTPULate, &b_events_EvtHdr__mOOTPULate);
   fChain->SetBranchAddress("EvtHdr_.mINTPU", &EvtHdr__mINTPU, &b_events_EvtHdr__mINTPU);
   fChain->SetBranchAddress("EvtHdr_.mNBX", &EvtHdr__mNBX, &b_events_EvtHdr__mNBX);
   fChain->SetBranchAddress("EvtHdr_.mPVndof", &EvtHdr__mPVndof, &b_events_EvtHdr__mPVndof);
   fChain->SetBranchAddress("EvtHdr_.mTrPu", &EvtHdr__mTrPu, &b_events_EvtHdr__mTrPu);
   fChain->SetBranchAddress("EvtHdr_.mPVx", &EvtHdr__mPVx, &b_events_EvtHdr__mPVx);
   fChain->SetBranchAddress("EvtHdr_.mPVy", &EvtHdr__mPVy, &b_events_EvtHdr__mPVy);
   fChain->SetBranchAddress("EvtHdr_.mPVz", &EvtHdr__mPVz, &b_events_EvtHdr__mPVz);
   fChain->SetBranchAddress("EvtHdr_.mBSx", &EvtHdr__mBSx, &b_events_EvtHdr__mBSx);
   fChain->SetBranchAddress("EvtHdr_.mBSy", &EvtHdr__mBSy, &b_events_EvtHdr__mBSy);
   fChain->SetBranchAddress("EvtHdr_.mBSz", &EvtHdr__mBSz, &b_events_EvtHdr__mBSz);
   fChain->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat, &b_events_EvtHdr__mPthat);
   fChain->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight, &b_events_EvtHdr__mWeight);
   fChain->SetBranchAddress("EvtHdr_.mCaloRho", &EvtHdr__mCaloRho, &b_events_EvtHdr__mCaloRho);
   fChain->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho, &b_events_EvtHdr__mPFRho);
   fChain->SetBranchAddress("CaloMet_.et_", &CaloMet__et_, &b_events_CaloMet__et_);
   fChain->SetBranchAddress("CaloMet_.CaloMetPt_", &CaloMet__CaloMetPt_, &b_events_CaloMet__CaloMetPt_);
   fChain->SetBranchAddress("CaloMet_.sumEt_", &CaloMet__sumEt_, &b_events_CaloMet__sumEt_);
   fChain->SetBranchAddress("CaloMet_.phi_", &CaloMet__phi_, &b_events_CaloMet__phi_);
   fChain->SetBranchAddress("PFMet_.et_", &PFMet__et_, &b_events_PFMet__et_);
   fChain->SetBranchAddress("PFMet_.CaloMetPt_", &PFMet__CaloMetPt_, &b_events_PFMet__CaloMetPt_);
   fChain->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_, &b_events_PFMet__sumEt_);
   fChain->SetBranchAddress("PFMet_.phi_", &PFMet__phi_, &b_events_PFMet__phi_);
   fChain->SetBranchAddress("MvaMet_.et_", &MvaMet__et_, &b_events_MvaMet__et_);
   fChain->SetBranchAddress("MvaMet_.CaloMetPt_", &MvaMet__CaloMetPt_, &b_events_MvaMet__CaloMetPt_);
   fChain->SetBranchAddress("MvaMet_.sumEt_", &MvaMet__sumEt_, &b_events_MvaMet__sumEt_);
   fChain->SetBranchAddress("MvaMet_.phi_", &MvaMet__phi_, &b_events_MvaMet__phi_);
   fChain->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
   fChain->SetBranchAddress("triggerList_", &triggerList_, &b_events_triggerList_);
   fChain->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);
   fChain->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
   fChain->SetBranchAddress("GenJets_", &GenJets__, &b_events_GenJets__);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fX", &GenJets__fCoordinates_fX, &b_GenJets__fCoordinates_fX);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fY", &GenJets__fCoordinates_fY, &b_GenJets__fCoordinates_fY);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fZ", &GenJets__fCoordinates_fZ, &b_GenJets__fCoordinates_fZ);
   fChain->SetBranchAddress("GenJets_.fCoordinates.fT", &GenJets__fCoordinates_fT, &b_GenJets__fCoordinates_fT);
   fChain->SetBranchAddress("CaloJets_", &CaloJets__, &b_events_CaloJets__);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fX", &CaloJets__P4__fCoordinates_fX, &b_CaloJets__P4__fCoordinates_fX);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fY", &CaloJets__P4__fCoordinates_fY, &b_CaloJets__P4__fCoordinates_fY);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fZ", &CaloJets__P4__fCoordinates_fZ, &b_CaloJets__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("CaloJets_.P4_.fCoordinates.fT", &CaloJets__P4__fCoordinates_fT, &b_CaloJets__P4__fCoordinates_fT);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fX", &CaloJets__genP4__fCoordinates_fX, &b_CaloJets__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fY", &CaloJets__genP4__fCoordinates_fY, &b_CaloJets__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fZ", &CaloJets__genP4__fCoordinates_fZ, &b_CaloJets__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("CaloJets_.genP4_.fCoordinates.fT", &CaloJets__genP4__fCoordinates_fT, &b_CaloJets__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("CaloJets_.genR_", &CaloJets__genR_, &b_CaloJets__genR_);
   fChain->SetBranchAddress("CaloJets_.cor_", &CaloJets__cor_, &b_CaloJets__cor_);
   fChain->SetBranchAddress("CaloJets_.jecLabels_", &CaloJets__jecLabels_, &b_CaloJets__jecLabels_);
   fChain->SetBranchAddress("CaloJets_.unc_", &CaloJets__unc_, &b_CaloJets__unc_);
   fChain->SetBranchAddress("CaloJets_.uncSrc_", &CaloJets__uncSrc_, &b_CaloJets__uncSrc_);
   fChain->SetBranchAddress("CaloJets_.area_", &CaloJets__area_, &b_CaloJets__area_);
   fChain->SetBranchAddress("CaloJets_.looseID_", &CaloJets__looseID_, &b_CaloJets__looseID_);
   fChain->SetBranchAddress("CaloJets_.tightID_", &CaloJets__tightID_, &b_CaloJets__tightID_);
   fChain->SetBranchAddress("CaloJets_.emf_", &CaloJets__emf_, &b_CaloJets__emf_);
   fChain->SetBranchAddress("CaloJets_.fHPD_", &CaloJets__fHPD_, &b_CaloJets__fHPD_);
   fChain->SetBranchAddress("CaloJets_.fRBX_", &CaloJets__fRBX_, &b_CaloJets__fRBX_);
   fChain->SetBranchAddress("CaloJets_.n90hits_", &CaloJets__n90hits_, &b_CaloJets__n90hits_);
   fChain->SetBranchAddress("CaloJets_.nTrkCalo_", &CaloJets__nTrkCalo_, &b_CaloJets__nTrkCalo_);
   fChain->SetBranchAddress("CaloJets_.nTrkVtx_", &CaloJets__nTrkVtx_, &b_CaloJets__nTrkVtx_);
   fChain->SetBranchAddress("PFJets_", &PFJets__, &b_events_PFJets__);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fX", PFJets__P4__fCoordinates_fX, &b_PFJets__P4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fY", PFJets__P4__fCoordinates_fY, &b_PFJets__P4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fZ", PFJets__P4__fCoordinates_fZ, &b_PFJets__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJets_.P4_.fCoordinates.fT", PFJets__P4__fCoordinates_fT, &b_PFJets__P4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fX", PFJets__genP4__fCoordinates_fX, &b_PFJets__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fY", PFJets__genP4__fCoordinates_fY, &b_PFJets__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fZ", PFJets__genP4__fCoordinates_fZ, &b_PFJets__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJets_.genP4_.fCoordinates.fT", PFJets__genP4__fCoordinates_fT, &b_PFJets__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJets_.genR_", PFJets__genR_, &b_PFJets__genR_);
   fChain->SetBranchAddress("PFJets_.cor_", PFJets__cor_, &b_PFJets__cor_);
   fChain->SetBranchAddress("PFJets_.jecLabels_", PFJets__jecLabels_, &b_PFJets__jecLabels_);
   fChain->SetBranchAddress("PFJets_.unc_", PFJets__unc_, &b_PFJets__unc_);
   fChain->SetBranchAddress("PFJets_.uncSrc_", PFJets__uncSrc_, &b_PFJets__uncSrc_);
   fChain->SetBranchAddress("PFJets_.area_", PFJets__area_, &b_PFJets__area_);
   fChain->SetBranchAddress("PFJets_.looseID_", PFJets__looseID_, &b_PFJets__looseID_);
   fChain->SetBranchAddress("PFJets_.tightID_", PFJets__tightID_, &b_PFJets__tightID_);
   fChain->SetBranchAddress("PFJets_.TCHE_", PFJets__TCHE_, &b_PFJets__TCHE_);
   fChain->SetBranchAddress("PFJets_.TCHP_", PFJets__TCHP_, &b_PFJets__TCHP_);
   fChain->SetBranchAddress("PFJets_.TCHEpf_", PFJets__TCHEpf_, &b_PFJets__TCHEpf_);
   fChain->SetBranchAddress("PFJets_.TCHPpf_", PFJets__TCHPpf_, &b_PFJets__TCHPpf_);
   fChain->SetBranchAddress("PFJets_.SoftMuonTagByIP_", PFJets__SoftMuonTagByIP_, &b_PFJets__SoftMuonTagByIP_);
   fChain->SetBranchAddress("PFJets_.SoftElectronTagByIP_", PFJets__SoftElectronTagByIP_, &b_PFJets__SoftElectronTagByIP_);
   fChain->SetBranchAddress("PFJets_.SoftMuonTag_", PFJets__SoftMuonTag_, &b_PFJets__SoftMuonTag_);
   fChain->SetBranchAddress("PFJets_.SoftElectronTag_", PFJets__SoftElectronTag_, &b_PFJets__SoftElectronTag_);
   fChain->SetBranchAddress("PFJets_.SimpleSecVertexHE_", PFJets__SimpleSecVertexHE_, &b_PFJets__SimpleSecVertexHE_);
   fChain->SetBranchAddress("PFJets_.SimpleSecVertexHP_", PFJets__SimpleSecVertexHP_, &b_PFJets__SimpleSecVertexHP_);
   fChain->SetBranchAddress("PFJets_.SimpleSecVertexHEpf_", PFJets__SimpleSecVertexHEpf_, &b_PFJets__SimpleSecVertexHEpf_);
   fChain->SetBranchAddress("PFJets_.SimpleSecVertexHPpf_", PFJets__SimpleSecVertexHPpf_, &b_PFJets__SimpleSecVertexHPpf_);
   fChain->SetBranchAddress("PFJets_.CSV_", PFJets__CSV_, &b_PFJets__CSV_);
   fChain->SetBranchAddress("PFJets_.CSVpf_", PFJets__CSVpf_, &b_PFJets__CSVpf_);
   fChain->SetBranchAddress("PFJets_.CinclSVpf_", PFJets__CinclSVpf_, &b_PFJets__CinclSVpf_);
   fChain->SetBranchAddress("PFJets_.CMVApf_", PFJets__CMVApf_, &b_PFJets__CMVApf_);
   fChain->SetBranchAddress("PFJets_.CSVSoftLeptonpf_", PFJets__CSVSoftLeptonpf_, &b_PFJets__CSVSoftLeptonpf_);
   fChain->SetBranchAddress("PFJets_.CSVpfPositive_", PFJets__CSVpfPositive_, &b_PFJets__CSVpfPositive_);
   fChain->SetBranchAddress("PFJets_.CSVpfNegative_", PFJets__CSVpfNegative_, &b_PFJets__CSVpfNegative_);
   fChain->SetBranchAddress("PFJets_.QGtagger_", PFJets__QGtagger_, &b_PFJets__QGtagger_);
   fChain->SetBranchAddress("PFJets_.partonFlavour_", PFJets__partonFlavour_, &b_PFJets__partonFlavour_);
   fChain->SetBranchAddress("PFJets_.hadronFlavour_", PFJets__hadronFlavour_, &b_PFJets__hadronFlavour_);
   fChain->SetBranchAddress("PFJets_.recommend1_", PFJets__recommend1_, &b_PFJets__recommend1_);
   fChain->SetBranchAddress("PFJets_.recommend2_", PFJets__recommend2_, &b_PFJets__recommend2_);
   fChain->SetBranchAddress("PFJets_.recommend3_", PFJets__recommend3_, &b_PFJets__recommend3_);
   fChain->SetBranchAddress("PFJets_.chf_", PFJets__chf_, &b_PFJets__chf_);
   fChain->SetBranchAddress("PFJets_.nhf_", PFJets__nhf_, &b_PFJets__nhf_);
   fChain->SetBranchAddress("PFJets_.nemf_", PFJets__nemf_, &b_PFJets__nemf_);
   fChain->SetBranchAddress("PFJets_.cemf_", PFJets__cemf_, &b_PFJets__cemf_);
   fChain->SetBranchAddress("PFJets_.muf_", PFJets__muf_, &b_PFJets__muf_);
   fChain->SetBranchAddress("PFJets_.hf_hf_", PFJets__hf_hf_, &b_PFJets__hf_hf_);
   fChain->SetBranchAddress("PFJets_.hf_phf_", PFJets__hf_phf_, &b_PFJets__hf_phf_);
   fChain->SetBranchAddress("PFJets_.hf_hm_", PFJets__hf_hm_, &b_PFJets__hf_hm_);
   fChain->SetBranchAddress("PFJets_.hf_phm_", PFJets__hf_phm_, &b_PFJets__hf_phm_);
   fChain->SetBranchAddress("PFJets_.chm_", PFJets__chm_, &b_PFJets__chm_);
   fChain->SetBranchAddress("PFJets_.nhm_", PFJets__nhm_, &b_PFJets__nhm_);
   fChain->SetBranchAddress("PFJets_.phm_", PFJets__phm_, &b_PFJets__phm_);
   fChain->SetBranchAddress("PFJets_.elm_", PFJets__elm_, &b_PFJets__elm_);
   fChain->SetBranchAddress("PFJets_.mum_", PFJets__mum_, &b_PFJets__mum_);
   fChain->SetBranchAddress("PFJets_.ncand_", PFJets__ncand_, &b_PFJets__ncand_);
   fChain->SetBranchAddress("PFJets_.beta_", PFJets__beta_, &b_PFJets__beta_);
   fChain->SetBranchAddress("PFJets_.betaStar_", PFJets__betaStar_, &b_PFJets__betaStar_);
   fChain->SetBranchAddress("PFJets_.mpuTrk_", PFJets__mpuTrk_, &b_PFJets__mpuTrk_);
   fChain->SetBranchAddress("PFJets_.mlvTrk_", PFJets__mlvTrk_, &b_PFJets__mlvTrk_);
   fChain->SetBranchAddress("PFJets_.mjtTrk_", PFJets__mjtTrk_, &b_PFJets__mjtTrk_);
   fChain->SetBranchAddress("PFJets_.hof_", PFJets__hof_, &b_PFJets__hof_);
   fChain->SetBranchAddress("PFJets_.pujid_", PFJets__pujid_, &b_PFJets__pujid_);
   fChain->SetBranchAddress("PFJets_.calojetpt_", PFJets__calojetpt_, &b_PFJets__calojetpt_);
   fChain->SetBranchAddress("PFJets_.calojetef_", PFJets__calojetef_, &b_PFJets__calojetef_);
   fChain->SetBranchAddress("PFJetsCHS_", &PFJetsCHS__, &b_events_PFJetsCHS__);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fX", PFJetsCHS__P4__fCoordinates_fX, &b_PFJetsCHS__P4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fY", PFJetsCHS__P4__fCoordinates_fY, &b_PFJetsCHS__P4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fZ", PFJetsCHS__P4__fCoordinates_fZ, &b_PFJetsCHS__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fT", PFJetsCHS__P4__fCoordinates_fT, &b_PFJetsCHS__P4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fX", PFJetsCHS__genP4__fCoordinates_fX, &b_PFJetsCHS__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fY", PFJetsCHS__genP4__fCoordinates_fY, &b_PFJetsCHS__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fZ", PFJetsCHS__genP4__fCoordinates_fZ, &b_PFJetsCHS__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("PFJetsCHS_.genP4_.fCoordinates.fT", PFJetsCHS__genP4__fCoordinates_fT, &b_PFJetsCHS__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("PFJetsCHS_.genR_", PFJetsCHS__genR_, &b_PFJetsCHS__genR_);
   fChain->SetBranchAddress("PFJetsCHS_.cor_", PFJetsCHS__cor_, &b_PFJetsCHS__cor_);
   fChain->SetBranchAddress("PFJetsCHS_.jecLabels_", PFJetsCHS__jecLabels_, &b_PFJetsCHS__jecLabels_);
   fChain->SetBranchAddress("PFJetsCHS_.unc_", PFJetsCHS__unc_, &b_PFJetsCHS__unc_);
   fChain->SetBranchAddress("PFJetsCHS_.uncSrc_", PFJetsCHS__uncSrc_, &b_PFJetsCHS__uncSrc_);
   fChain->SetBranchAddress("PFJetsCHS_.area_", PFJetsCHS__area_, &b_PFJetsCHS__area_);
   fChain->SetBranchAddress("PFJetsCHS_.looseID_", PFJetsCHS__looseID_, &b_PFJetsCHS__looseID_);
   fChain->SetBranchAddress("PFJetsCHS_.tightID_", PFJetsCHS__tightID_, &b_PFJetsCHS__tightID_);
   fChain->SetBranchAddress("PFJetsCHS_.TCHE_", PFJetsCHS__TCHE_, &b_PFJetsCHS__TCHE_);
   fChain->SetBranchAddress("PFJetsCHS_.TCHP_", PFJetsCHS__TCHP_, &b_PFJetsCHS__TCHP_);
   fChain->SetBranchAddress("PFJetsCHS_.TCHEpf_", PFJetsCHS__TCHEpf_, &b_PFJetsCHS__TCHEpf_);
   fChain->SetBranchAddress("PFJetsCHS_.TCHPpf_", PFJetsCHS__TCHPpf_, &b_PFJetsCHS__TCHPpf_);
   fChain->SetBranchAddress("PFJetsCHS_.SoftMuonTagByIP_", PFJetsCHS__SoftMuonTagByIP_, &b_PFJetsCHS__SoftMuonTagByIP_);
   fChain->SetBranchAddress("PFJetsCHS_.SoftElectronTagByIP_", PFJetsCHS__SoftElectronTagByIP_, &b_PFJetsCHS__SoftElectronTagByIP_);
   fChain->SetBranchAddress("PFJetsCHS_.SoftMuonTag_", PFJetsCHS__SoftMuonTag_, &b_PFJetsCHS__SoftMuonTag_);
   fChain->SetBranchAddress("PFJetsCHS_.SoftElectronTag_", PFJetsCHS__SoftElectronTag_, &b_PFJetsCHS__SoftElectronTag_);
   fChain->SetBranchAddress("PFJetsCHS_.SimpleSecVertexHE_", PFJetsCHS__SimpleSecVertexHE_, &b_PFJetsCHS__SimpleSecVertexHE_);
   fChain->SetBranchAddress("PFJetsCHS_.SimpleSecVertexHP_", PFJetsCHS__SimpleSecVertexHP_, &b_PFJetsCHS__SimpleSecVertexHP_);
   fChain->SetBranchAddress("PFJetsCHS_.SimpleSecVertexHEpf_", PFJetsCHS__SimpleSecVertexHEpf_, &b_PFJetsCHS__SimpleSecVertexHEpf_);
   fChain->SetBranchAddress("PFJetsCHS_.SimpleSecVertexHPpf_", PFJetsCHS__SimpleSecVertexHPpf_, &b_PFJetsCHS__SimpleSecVertexHPpf_);
   fChain->SetBranchAddress("PFJetsCHS_.CSV_", PFJetsCHS__CSV_, &b_PFJetsCHS__CSV_);
   fChain->SetBranchAddress("PFJetsCHS_.CSVpf_", PFJetsCHS__CSVpf_, &b_PFJetsCHS__CSVpf_);
   fChain->SetBranchAddress("PFJetsCHS_.CinclSVpf_", PFJetsCHS__CinclSVpf_, &b_PFJetsCHS__CinclSVpf_);
   fChain->SetBranchAddress("PFJetsCHS_.CMVApf_", PFJetsCHS__CMVApf_, &b_PFJetsCHS__CMVApf_);
   fChain->SetBranchAddress("PFJetsCHS_.CSVSoftLeptonpf_", PFJetsCHS__CSVSoftLeptonpf_, &b_PFJetsCHS__CSVSoftLeptonpf_);
   fChain->SetBranchAddress("PFJetsCHS_.CSVpfPositive_", PFJetsCHS__CSVpfPositive_, &b_PFJetsCHS__CSVpfPositive_);
   fChain->SetBranchAddress("PFJetsCHS_.CSVpfNegative_", PFJetsCHS__CSVpfNegative_, &b_PFJetsCHS__CSVpfNegative_);
   fChain->SetBranchAddress("PFJetsCHS_.QGtagger_", PFJetsCHS__QGtagger_, &b_PFJetsCHS__QGtagger_);
   fChain->SetBranchAddress("PFJetsCHS_.partonFlavour_", PFJetsCHS__partonFlavour_, &b_PFJetsCHS__partonFlavour_);
   fChain->SetBranchAddress("PFJetsCHS_.hadronFlavour_", PFJetsCHS__hadronFlavour_, &b_PFJetsCHS__hadronFlavour_);
   fChain->SetBranchAddress("PFJetsCHS_.recommend1_", PFJetsCHS__recommend1_, &b_PFJetsCHS__recommend1_);
   fChain->SetBranchAddress("PFJetsCHS_.recommend2_", PFJetsCHS__recommend2_, &b_PFJetsCHS__recommend2_);
   fChain->SetBranchAddress("PFJetsCHS_.recommend3_", PFJetsCHS__recommend3_, &b_PFJetsCHS__recommend3_);
   fChain->SetBranchAddress("PFJetsCHS_.chf_", PFJetsCHS__chf_, &b_PFJetsCHS__chf_);
   fChain->SetBranchAddress("PFJetsCHS_.nhf_", PFJetsCHS__nhf_, &b_PFJetsCHS__nhf_);
   fChain->SetBranchAddress("PFJetsCHS_.nemf_", PFJetsCHS__nemf_, &b_PFJetsCHS__nemf_);
   fChain->SetBranchAddress("PFJetsCHS_.cemf_", PFJetsCHS__cemf_, &b_PFJetsCHS__cemf_);
   fChain->SetBranchAddress("PFJetsCHS_.muf_", PFJetsCHS__muf_, &b_PFJetsCHS__muf_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_hf_", PFJetsCHS__hf_hf_, &b_PFJetsCHS__hf_hf_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_phf_", PFJetsCHS__hf_phf_, &b_PFJetsCHS__hf_phf_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_hm_", PFJetsCHS__hf_hm_, &b_PFJetsCHS__hf_hm_);
   fChain->SetBranchAddress("PFJetsCHS_.hf_phm_", PFJetsCHS__hf_phm_, &b_PFJetsCHS__hf_phm_);
   fChain->SetBranchAddress("PFJetsCHS_.chm_", PFJetsCHS__chm_, &b_PFJetsCHS__chm_);
   fChain->SetBranchAddress("PFJetsCHS_.nhm_", PFJetsCHS__nhm_, &b_PFJetsCHS__nhm_);
   fChain->SetBranchAddress("PFJetsCHS_.phm_", PFJetsCHS__phm_, &b_PFJetsCHS__phm_);
   fChain->SetBranchAddress("PFJetsCHS_.elm_", PFJetsCHS__elm_, &b_PFJetsCHS__elm_);
   fChain->SetBranchAddress("PFJetsCHS_.mum_", PFJetsCHS__mum_, &b_PFJetsCHS__mum_);
   fChain->SetBranchAddress("PFJetsCHS_.ncand_", PFJetsCHS__ncand_, &b_PFJetsCHS__ncand_);
   fChain->SetBranchAddress("PFJetsCHS_.beta_", PFJetsCHS__beta_, &b_PFJetsCHS__beta_);
   fChain->SetBranchAddress("PFJetsCHS_.betaStar_", PFJetsCHS__betaStar_, &b_PFJetsCHS__betaStar_);
   fChain->SetBranchAddress("PFJetsCHS_.mpuTrk_", PFJetsCHS__mpuTrk_, &b_PFJetsCHS__mpuTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.mlvTrk_", PFJetsCHS__mlvTrk_, &b_PFJetsCHS__mlvTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.mjtTrk_", PFJetsCHS__mjtTrk_, &b_PFJetsCHS__mjtTrk_);
   fChain->SetBranchAddress("PFJetsCHS_.hof_", PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
   fChain->SetBranchAddress("PFJetsCHS_.pujid_", PFJetsCHS__pujid_, &b_PFJetsCHS__pujid_);
   fChain->SetBranchAddress("PFJetsCHS_.calojetpt_", PFJetsCHS__calojetpt_, &b_PFJetsCHS__calojetpt_);
   fChain->SetBranchAddress("PFJetsCHS_.calojetef_", PFJetsCHS__calojetef_, &b_PFJetsCHS__calojetef_);
   fChain->SetBranchAddress("genFlavour_", &genFlavour_, &b_events_genFlavour_);
   fChain->SetBranchAddress("genFlavourHadron_", &genFlavourHadron_, &b_events_genFlavourHadron_);
   fChain->SetBranchAddress("mMuon_", &mMuon__, &b_events_mMuon__);
   fChain->SetBranchAddress("mMuon_.P4_.fCoordinates.fX", &mMuon__P4__fCoordinates_fX, &b_mMuon__P4__fCoordinates_fX);
   fChain->SetBranchAddress("mMuon_.P4_.fCoordinates.fY", &mMuon__P4__fCoordinates_fY, &b_mMuon__P4__fCoordinates_fY);
   fChain->SetBranchAddress("mMuon_.P4_.fCoordinates.fZ", &mMuon__P4__fCoordinates_fZ, &b_mMuon__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("mMuon_.P4_.fCoordinates.fT", &mMuon__P4__fCoordinates_fT, &b_mMuon__P4__fCoordinates_fT);
   fChain->SetBranchAddress("mMuon_.genP4_.fCoordinates.fX", &mMuon__genP4__fCoordinates_fX, &b_mMuon__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("mMuon_.genP4_.fCoordinates.fY", &mMuon__genP4__fCoordinates_fY, &b_mMuon__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("mMuon_.genP4_.fCoordinates.fZ", &mMuon__genP4__fCoordinates_fZ, &b_mMuon__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("mMuon_.genP4_.fCoordinates.fT", &mMuon__genP4__fCoordinates_fT, &b_mMuon__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("mMuon_.genR_", &mMuon__genR_, &b_mMuon__genR_);
   fChain->SetBranchAddress("mMuon_.cor_", &mMuon__cor_, &b_mMuon__cor_);
   fChain->SetBranchAddress("mMuon_.jecLabels_", &mMuon__jecLabels_, &b_mMuon__jecLabels_);
   fChain->SetBranchAddress("mMuon_.unc_", &mMuon__unc_, &b_mMuon__unc_);
   fChain->SetBranchAddress("mMuon_.uncSrc_", &mMuon__uncSrc_, &b_mMuon__uncSrc_);
   fChain->SetBranchAddress("mMuon_.area_", &mMuon__area_, &b_mMuon__area_);
   fChain->SetBranchAddress("mMuon_.DxyVertex_", &mMuon__DxyVertex_, &b_mMuon__DxyVertex_);
   fChain->SetBranchAddress("mMuon_.DzVertex_", &mMuon__DzVertex_, &b_mMuon__DzVertex_);
   fChain->SetBranchAddress("mMuon_.PDGID_", &mMuon__PDGID_, &b_mMuon__PDGID_);
   fChain->SetBranchAddress("mMuon_.PfIso_", &mMuon__PfIso_, &b_mMuon__PfIso_);
   fChain->SetBranchAddress("mMuon_.PuChargedHadronIso_", &mMuon__PuChargedHadronIso_, &b_mMuon__PuChargedHadronIso_);
   fChain->SetBranchAddress("mMuon_.ChargedHadronIso_", &mMuon__ChargedHadronIso_, &b_mMuon__ChargedHadronIso_);
   fChain->SetBranchAddress("mMuon_.NeutralHadronIso_", &mMuon__NeutralHadronIso_, &b_mMuon__NeutralHadronIso_);
   fChain->SetBranchAddress("mMuon_.PhotonIso_", &mMuon__PhotonIso_, &b_mMuon__PhotonIso_);
   fChain->SetBranchAddress("mElectron_", &mElectron__, &b_events_mElectron__);
   fChain->SetBranchAddress("mElectron_.P4_.fCoordinates.fX", &mElectron__P4__fCoordinates_fX, &b_mElectron__P4__fCoordinates_fX);
   fChain->SetBranchAddress("mElectron_.P4_.fCoordinates.fY", &mElectron__P4__fCoordinates_fY, &b_mElectron__P4__fCoordinates_fY);
   fChain->SetBranchAddress("mElectron_.P4_.fCoordinates.fZ", &mElectron__P4__fCoordinates_fZ, &b_mElectron__P4__fCoordinates_fZ);
   fChain->SetBranchAddress("mElectron_.P4_.fCoordinates.fT", &mElectron__P4__fCoordinates_fT, &b_mElectron__P4__fCoordinates_fT);
   fChain->SetBranchAddress("mElectron_.genP4_.fCoordinates.fX", &mElectron__genP4__fCoordinates_fX, &b_mElectron__genP4__fCoordinates_fX);
   fChain->SetBranchAddress("mElectron_.genP4_.fCoordinates.fY", &mElectron__genP4__fCoordinates_fY, &b_mElectron__genP4__fCoordinates_fY);
   fChain->SetBranchAddress("mElectron_.genP4_.fCoordinates.fZ", &mElectron__genP4__fCoordinates_fZ, &b_mElectron__genP4__fCoordinates_fZ);
   fChain->SetBranchAddress("mElectron_.genP4_.fCoordinates.fT", &mElectron__genP4__fCoordinates_fT, &b_mElectron__genP4__fCoordinates_fT);
   fChain->SetBranchAddress("mElectron_.genR_", &mElectron__genR_, &b_mElectron__genR_);
   fChain->SetBranchAddress("mElectron_.cor_", &mElectron__cor_, &b_mElectron__cor_);
   fChain->SetBranchAddress("mElectron_.jecLabels_", &mElectron__jecLabels_, &b_mElectron__jecLabels_);
   fChain->SetBranchAddress("mElectron_.unc_", &mElectron__unc_, &b_mElectron__unc_);
   fChain->SetBranchAddress("mElectron_.uncSrc_", &mElectron__uncSrc_, &b_mElectron__uncSrc_);
   fChain->SetBranchAddress("mElectron_.area_", &mElectron__area_, &b_mElectron__area_);
   fChain->SetBranchAddress("mElectron_.DxyVertex_", &mElectron__DxyVertex_, &b_mElectron__DxyVertex_);
   fChain->SetBranchAddress("mElectron_.DzVertex_", &mElectron__DzVertex_, &b_mElectron__DzVertex_);
   fChain->SetBranchAddress("mElectron_.PDGID_", &mElectron__PDGID_, &b_mElectron__PDGID_);
   fChain->SetBranchAddress("mElectron_.PfIso_", &mElectron__PfIso_, &b_mElectron__PfIso_);
   fChain->SetBranchAddress("mElectron_.PuChargedHadronIso_", &mElectron__PuChargedHadronIso_, &b_mElectron__PuChargedHadronIso_);
   fChain->SetBranchAddress("mElectron_.ChargedHadronIso_", &mElectron__ChargedHadronIso_, &b_mElectron__ChargedHadronIso_);
   fChain->SetBranchAddress("mElectron_.NeutralHadronIso_", &mElectron__NeutralHadronIso_, &b_mElectron__NeutralHadronIso_);
   fChain->SetBranchAddress("mElectron_.PhotonIso_", &mElectron__PhotonIso_, &b_mElectron__PhotonIso_);
   Notify();
}

Bool_t HeaderSMPJ::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HeaderSMPJ::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HeaderSMPJ::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HeaderSMPJ_cxx
