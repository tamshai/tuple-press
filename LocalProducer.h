//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 23 18:09:40 2016 by ROOT version 6.04/14
// from TTree ProcessedTree/ProcessedTree
// found on file: ProcessedTree_data.root
//////////////////////////////////////////////////////////

#ifndef LocalProducer_h
#define LocalProducer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <iostream>
#include <TLorentzVector.h>
#include <THashList.h>

// JEC headers
#include "interface/FactorizedJetCorrector.h"
#include "interface/JetCorrectorParameters.h"

// Header file for the classes stored in the TTree if any.
#include "Math/GenVector/PxPyPzE4D.h"

class LocalProducer {
public :
   TTree          *fChain_ak4;   //!pointer to the analyzed TTree or TChain
   TTree          *fChain_ak7;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxGenJets_ = 64;
   static const Int_t kMaxPFJetsCHS_ = 64;


   // Declaration of leaf types
   TH1F*           TriggerNames;

   Int_t           EvtHdr__mRun;
   Int_t           EvtHdr__mEvent;
   Int_t           EvtHdr__mLumi;

   Float_t         EvtHdr__mPthat;
   Float_t         EvtHdr__mWeight;
   Float_t         EvtHdr__mPFRho;

   Float_t         PFMet__et_;
   Float_t         PFMet__sumEt_;

   vector<int>     TriggerDecision_;
   vector<int>     L1Prescale_;
   vector<int>     HLTPrescale_;

   Int_t           GenJets__;
   Double_t        GenJets__fCoordinates_fX[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fY[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fZ[kMaxGenJets_];   //[GenJets__]
   Double_t        GenJets__fCoordinates_fT[kMaxGenJets_];   //[GenJets__]

   Int_t           PFJetsCHS__;
   Double_t        PFJetsCHS__P4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   
   /*
   Double_t        PFJetsCHS__genP4__fCoordinates_fX[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fY[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fZ[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__genP4__fCoordinates_fT[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   */
   
   Float_t         PFJetsCHS__cor_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__area_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Bool_t          PFJetsCHS__tightID_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

   Float_t         PFJetsCHS__chf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__nhf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__nemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__cemf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__muf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hf_hf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hf_phf_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hf_hm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__hf_phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]*/
   Int_t           PFJetsCHS__chm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__nhm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__phm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__elm_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Int_t           PFJetsCHS__mum_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

   Float_t         PFJetsCHS__QGtagger_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__beta_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__betaStar_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hof_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   
   Float_t         PFJetsCHS__CSV_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

   Float_t         PFJetsCHS__partonFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__hadronFlavour_[kMaxPFJetsCHS_];   //[PFJetsCHS__]

   Int_t           EvtHdr__mRun_ak7;
   Int_t           EvtHdr__mEvent_ak7;
   Int_t           EvtHdr__mLumi_ak7;

   Int_t           PFJetsCHS_ak7__;
   Double_t        PFJetsCHS__P4__fCoordinates_fX_ak7[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fY_ak7[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fZ_ak7[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Double_t        PFJetsCHS__P4__fCoordinates_fT_ak7[kMaxPFJetsCHS_];   //[PFJetsCHS__]
      
   Float_t         PFJetsCHS__cor_ak7_[kMaxPFJetsCHS_];   //[PFJetsCHS__]
   Float_t         PFJetsCHS__area_ak7_[kMaxPFJetsCHS_];   //[PFJetsCHS__]


   // List of branches
   TBranch        *b_events_EvtHdr__mRun;   //!
   TBranch        *b_events_EvtHdr__mEvent;   //!
   TBranch        *b_events_EvtHdr__mLumi;   //!

   TBranch        *b_events_EvtHdr__mPthat;   //!
   TBranch        *b_events_EvtHdr__mWeight;   //!

   TBranch        *b_events_EvtHdr__mPFRho;   //!
   TBranch        *b_events_PFMet__et_;   //!
   TBranch        *b_events_PFMet__sumEt_;   //!

   TBranch        *b_events_TriggerDecision_;   //!
   TBranch        *b_events_L1Prescale_;   //!
   TBranch        *b_events_HLTPrescale_;   //!

   TBranch        *b_events_GenJets__;   //!
   TBranch        *b_GenJets__fCoordinates_fX;   //!
   TBranch        *b_GenJets__fCoordinates_fY;   //!
   TBranch        *b_GenJets__fCoordinates_fZ;   //!
   TBranch        *b_GenJets__fCoordinates_fT;   //!

   TBranch        *b_events_PFJetsCHS__;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fX;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fY;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fZ;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fT;   //!
   /*
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fX;   //!
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fY;   //!
   TBranch        *b_PFJetsCHS__genP4__fCoordinates_fZ;   //!
   TBranch        *b_PFJets"__genP4__fCoordinates_fT;   //!
   */
   TBranch        *b_PFJetsCHS__cor_;   //!
   
   TBranch        *b_PFJetsCHS__area_;   //!
   TBranch        *b_PFJetsCHS__tightID_;   //!


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

   TBranch        *b_PFJetsCHS__QGtagger_;   //!
   TBranch        *b_PFJetsCHS__beta_;   //!
   TBranch        *b_PFJetsCHS__betaStar_;   //!
   TBranch        *b_PFJetsCHS__hof_;   //!

   TBranch        *b_PFJetsCHS__CSV_;   //!

   TBranch        *b_PFJetsCHS__partonFlavour_;   //!
   TBranch        *b_PFJetsCHS__hadronFlavour_;   //!

   TBranch        *b_events_EvtHdr__mRun_ak7;   //!
   TBranch        *b_events_EvtHdr__mEvent_ak7;   //!
   TBranch        *b_events_EvtHdr__mLumi_ak7;   //!

   TBranch        *b_events_PFJetsCHS_ak7__;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fX_ak7;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fY_ak7;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fZ_ak7;   //!
   TBranch        *b_PFJetsCHS__P4__fCoordinates_fT_ak7;   //!

   TBranch        *b_PFJetsCHS__cor_ak7_;  
   TBranch        *b_PFJetsCHS__area_ak7_;   //!


   LocalProducer(TTree *tree_ak4 = 0, TTree *tree_ak7 = 0, TH1F* trgNames = 0, Bool_t mc = false);
   virtual ~LocalProducer();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree_ak4 = 0, TTree *tree_ak7 = 0);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   std::vector<std::string> _pfTriggers = {};

private:

   FactorizedJetCorrector *_JEC;
   Bool_t isMC ;
};

#endif

#ifdef LocalProducer_cxx
LocalProducer::LocalProducer(TTree *tree_ak4, TTree *tree_ak7, TH1F* trgNames, Bool_t mc) : fChain_ak4(0), TriggerNames(trgNames), isMC(mc)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree_ak4, tree_ak7);
   Loop();
}

LocalProducer::~LocalProducer()
{

}

Int_t LocalProducer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain_ak4) return 0;
   return fChain_ak4->GetEntry(entry);
}

Long64_t LocalProducer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain_ak4) 
      return -5;
   Long64_t centry = fChain_ak4->LoadTree(entry);

   if (centry < 0) 
      return centry;

   if (fChain_ak4->GetTreeNumber() != fCurrent) {
      fCurrent = fChain_ak4->GetTreeNumber();
      
   }
   return centry;
}

void LocalProducer::Init(TTree *tree_ak4, TTree *tree_ak7)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree_ak4) return;
   fChain_ak4 = tree_ak4;
   fCurrent = -1;
   fChain_ak4->SetMakeClass(1);


   // Event identification
   fChain_ak4->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun, &b_events_EvtHdr__mRun);
   fChain_ak4->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent, &b_events_EvtHdr__mEvent);
   fChain_ak4->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi, &b_events_EvtHdr__mLumi);
   fChain_ak4->SetBranchAddress("EvtHdr_.mPthat", &EvtHdr__mPthat, &b_events_EvtHdr__mPthat);
   fChain_ak4->SetBranchAddress("EvtHdr_.mWeight", &EvtHdr__mWeight, &b_events_EvtHdr__mWeight);

   fChain_ak4->SetBranchAddress("EvtHdr_.mPFRho", &EvtHdr__mPFRho, &b_events_EvtHdr__mPFRho);
   fChain_ak4->SetBranchAddress("PFMet_.et_", &PFMet__et_, &b_events_PFMet__et_);
   fChain_ak4->SetBranchAddress("PFMet_.sumEt_", &PFMet__sumEt_, &b_events_PFMet__sumEt_);

   // Triggers and prescale values
   fChain_ak4->SetBranchAddress("TriggerDecision_", &TriggerDecision_, &b_events_TriggerDecision_);
   fChain_ak4->SetBranchAddress("L1Prescale_", &L1Prescale_, &b_events_L1Prescale_);   
   fChain_ak4->SetBranchAddress("HLTPrescale_", &HLTPrescale_, &b_events_HLTPrescale_);
   
   // Monte Carlo Generated jets
   fChain_ak4->SetBranchAddress("GenJets_", &GenJets__, &b_events_GenJets__);
   fChain_ak4->SetBranchAddress("GenJets_.fCoordinates.fX", GenJets__fCoordinates_fX, &b_GenJets__fCoordinates_fX);
   fChain_ak4->SetBranchAddress("GenJets_.fCoordinates.fY", GenJets__fCoordinates_fY, &b_GenJets__fCoordinates_fY);
   fChain_ak4->SetBranchAddress("GenJets_.fCoordinates.fZ", GenJets__fCoordinates_fZ, &b_GenJets__fCoordinates_fZ);
   fChain_ak4->SetBranchAddress("GenJets_.fCoordinates.fT", GenJets__fCoordinates_fT, &b_GenJets__fCoordinates_fT);
   
   // PF jets
   fChain_ak4->SetBranchAddress("PFJetsCHS_", &PFJetsCHS__, &b_events_PFJetsCHS__);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fX", PFJetsCHS__P4__fCoordinates_fX, &b_PFJetsCHS__P4__fCoordinates_fX);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fY", PFJetsCHS__P4__fCoordinates_fY, &b_PFJetsCHS__P4__fCoordinates_fY);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fZ", PFJetsCHS__P4__fCoordinates_fZ, &b_PFJetsCHS__P4__fCoordinates_fZ);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fT", PFJetsCHS__P4__fCoordinates_fT, &b_PFJetsCHS__P4__fCoordinates_fT);

   // JEC correction
   fChain_ak4->SetBranchAddress("PFJetsCHS_.cor_", PFJetsCHS__cor_, &b_PFJetsCHS__cor_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.area_", PFJetsCHS__area_, &b_PFJetsCHS__area_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.tightID_", PFJetsCHS__tightID_, &b_PFJetsCHS__tightID_);

   // Composition quantities
   fChain_ak4->SetBranchAddress("PFJetsCHS_.chf_", PFJetsCHS__chf_, &b_PFJetsCHS__chf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.nhf_", PFJetsCHS__nhf_, &b_PFJetsCHS__nhf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.nemf_", PFJetsCHS__nemf_, &b_PFJetsCHS__nemf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.cemf_", PFJetsCHS__cemf_, &b_PFJetsCHS__cemf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.muf_", PFJetsCHS__muf_, &b_PFJetsCHS__muf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.hf_hf_", PFJetsCHS__hf_hf_, &b_PFJetsCHS__hf_hf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.hf_phf_", PFJetsCHS__hf_phf_, &b_PFJetsCHS__hf_phf_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.hf_hm_", PFJetsCHS__hf_hm_, &b_PFJetsCHS__hf_hm_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.hf_phm_", PFJetsCHS__hf_phm_, &b_PFJetsCHS__hf_phm_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.chm_", PFJetsCHS__chm_, &b_PFJetsCHS__chm_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.nhm_", PFJetsCHS__nhm_, &b_PFJetsCHS__nhm_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.phm_", PFJetsCHS__phm_, &b_PFJetsCHS__phm_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.elm_", PFJetsCHS__elm_, &b_PFJetsCHS__elm_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.mum_", PFJetsCHS__mum_, &b_PFJetsCHS__mum_);


   fChain_ak4->SetBranchAddress("PFJetsCHS_.QGtagger_", PFJetsCHS__QGtagger_, &b_PFJetsCHS__QGtagger_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.beta_", PFJetsCHS__beta_, &b_PFJetsCHS__beta_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.betaStar_", PFJetsCHS__betaStar_, &b_PFJetsCHS__betaStar_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.hof_", PFJetsCHS__hof_, &b_PFJetsCHS__hof_);
      
   fChain_ak4->SetBranchAddress("PFJetsCHS_.CSV_", PFJetsCHS__CSV_, &b_PFJetsCHS__CSV_);
   
   // Generator-level flavour
   fChain_ak4->SetBranchAddress("PFJetsCHS_.partonFlavour_", PFJetsCHS__partonFlavour_, &b_PFJetsCHS__partonFlavour_);
   fChain_ak4->SetBranchAddress("PFJetsCHS_.hadronFlavour_", PFJetsCHS__hadronFlavour_, &b_PFJetsCHS__hadronFlavour_);


   // AK7 variables
   fChain_ak7 = tree_ak7;
   if (fChain_ak7) fChain_ak7->SetMakeClass(1);

   fChain_ak7->SetBranchAddress("EvtHdr_.mRun", &EvtHdr__mRun_ak7, &b_events_EvtHdr__mRun_ak7);
   fChain_ak7->SetBranchAddress("EvtHdr_.mEvent", &EvtHdr__mEvent_ak7, &b_events_EvtHdr__mEvent_ak7);
   fChain_ak7->SetBranchAddress("EvtHdr_.mLumi", &EvtHdr__mLumi_ak7, &b_events_EvtHdr__mLumi_ak7);

   fChain_ak7->SetBranchAddress("PFJetsCHS_", &PFJetsCHS_ak7__, &b_events_PFJetsCHS_ak7__);
   fChain_ak7->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fX", PFJetsCHS__P4__fCoordinates_fX_ak7, &b_PFJetsCHS__P4__fCoordinates_fX_ak7);
   fChain_ak7->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fY", PFJetsCHS__P4__fCoordinates_fY_ak7, &b_PFJetsCHS__P4__fCoordinates_fY_ak7);
   fChain_ak7->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fZ", PFJetsCHS__P4__fCoordinates_fZ_ak7, &b_PFJetsCHS__P4__fCoordinates_fZ_ak7);
   fChain_ak7->SetBranchAddress("PFJetsCHS_.P4_.fCoordinates.fT", PFJetsCHS__P4__fCoordinates_fT_ak7, &b_PFJetsCHS__P4__fCoordinates_fT_ak7);

   fChain_ak7->SetBranchAddress("PFJetsCHS_.cor_", PFJetsCHS__cor_ak7_, &b_PFJetsCHS__cor_ak7_);
   fChain_ak7->SetBranchAddress("PFJetsCHS_.area_", PFJetsCHS__area_ak7_, &b_PFJetsCHS__area_ak7_);
   

   Notify();
}

Bool_t LocalProducer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LocalProducer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain_ak4) return;
   fChain_ak4->Show(entry);
}
Int_t LocalProducer::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LocalProducer_cxx
