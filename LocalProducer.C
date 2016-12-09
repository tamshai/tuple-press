#define LocalProducer_cxx
#include "LocalProducer.h"
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <regex>


// Compactify (Use UChar_t since Char_t == signed char)
inline UChar_t mapFloatToUChar(Float_t frac) {
// Convert floats in [0,1] to uchar in [0,255] to save space
// https://stackoverflow.com/questions/599976/convert-quantize-float-range-to-integer-range
  return (UChar_t)(min(255, (int) (frac * 255.)));
}

inline UChar_t mapIntToUChar(Int_t mult) {
// Convert 4 byte integer to 1 byte
// Return value capped at 255
  return (UChar_t)(min(255, mult));
}

inline Char_t mapIntToChar(Int_t num) {
// Convert 4 byte integer to 1 byte
// Return value capped at 255
  return (Char_t)num;
}

// De-compactify
inline Float_t mapUCharToFloat(UChar_t cfrac) {
// Convert char back to float
  return (Float_t)(cfrac/255.);
}

inline Int_t mapUCharToInt(UChar_t cmult) {
// Convert 1 byte to 4 byte integer
// Not absolutely necessary, only for consistency
  return (Int_t)cmult;
}



void LocalProducer::Loop()
{
//   In a ROOT session, you can do:
//      root> .L LocalProducer.C
//      root> LocalProducer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


  assert(fChain_ak4);

  // Save input directory
  TDirectory *currDir = gDirectory;

  std::string outName = ( isMC ? "./FlatTuple_MC.root" : "./FlatTuple_DATA.root");
  TFile *fout = new TFile(outName.c_str(), "RECREATE");

  std::cout << "Output will be saved in file: " << outName << std::endl;

  // Create new tuple/tree into output file
  TTree *tree = new TTree("FlatTree", "FlatTree");

  // Change back to input directory
  currDir->cd();


  // Constants
  const Int_t kMaxNjet = 64;
  const Int_t kMaxNtrg = 64;

  // Number of leading jets for which composition variables are saved
  const Int_t kMaxNComp = 3;


  // PF AK5 jets
  UChar_t njet;
  Float_t pt[kMaxNjet];
  Float_t eta[kMaxNjet];
  Float_t phi[kMaxNjet];
  Float_t E[kMaxNjet];
  Bool_t tightID[kMaxNjet];
  Float_t area[kMaxNjet];
  Float_t jec[kMaxNjet];
  Int_t igen[kMaxNjet];

  // PF AK7 jets
  Int_t njet_ak7;
  Float_t pt_ak7[kMaxNjet];
  Float_t eta_ak7[kMaxNjet];
  Float_t phi_ak7[kMaxNjet];
  Float_t E_ak7[kMaxNjet];
  Float_t area_ak7[kMaxNjet];
  Float_t jec_ak7[kMaxNjet];
  Int_t ak7_to_ak4[kMaxNjet];

  // Jet composition
  UChar_t ncomp;
  UChar_t chf[kMaxNComp];
  UChar_t nhf[kMaxNComp];
  UChar_t phf[kMaxNComp];
  UChar_t elf[kMaxNComp];
  UChar_t muf[kMaxNComp];
  UChar_t hf_hf[kMaxNComp];
  UChar_t hf_phf[kMaxNComp];
  UChar_t hf_hm[kMaxNComp];
  UChar_t hf_phm[kMaxNComp];
  UChar_t chm[kMaxNComp];
  UChar_t nhm[kMaxNComp];
  UChar_t phm[kMaxNComp];
  UChar_t elm[kMaxNComp];
  UChar_t mum[kMaxNComp];   
  UChar_t beta[kMaxNComp];   
  UChar_t bstar[kMaxNComp];
  UChar_t hof[kMaxNComp];  
  UChar_t qgl[kMaxNComp];

  UChar_t csv[kMaxNComp];
  Char_t pfl[kMaxNComp];
  Char_t hfl[kMaxNComp];


  // Generated jets
  UChar_t ngen;
  Float_t gen_pt[kMaxNjet];
  Float_t gen_eta[kMaxNjet];
  Float_t gen_phi[kMaxNjet];
  Float_t gen_E[kMaxNjet];


  // Event identification
  UInt_t run;
  UInt_t lumi;
  ULong64_t event;

  // Triggers
  UInt_t ntrg;
  //std::vector<std::string> triggernames;
  std::vector<std::string> trgfired;
  //char triggernames[64][kMaxNtrg];
  std::vector<int> prescales;

  // MET, SuMET, rho, eventfilter
  Float_t met;
  Float_t sumet;
  Float_t rho;

  // MC variables
  Float_t pthat;
  Float_t mcweight;


  TBranch *b_njet = tree->Branch("njet", &njet, "njet/b");
  TBranch *b_pt = tree->Branch("pt", pt, "pt[njet]/F");
  TBranch *b_eta = tree->Branch("eta", eta, "eta[njet]/F");
  TBranch *b_phi = tree->Branch("phi", phi, "phi[njet]/F");
  TBranch *b_E = tree->Branch("E", E, "E[njet]/F");   
  TBranch *b_tightID = tree->Branch("tightID", tightID, "tightID[njet]/O");
  TBranch *b_area = tree->Branch("area", area, "area[njet]/F");
  TBranch *b_jec = tree->Branch("jec", jec, "jec[njet]/F");
  
  TBranch *b_njet_ak7 = tree->Branch("njet_ak7", &njet_ak7, "njet_ak7/i");
  TBranch *b_pt_ak7 = tree->Branch("pt_ak7", pt_ak7, "pt_ak7[njet_ak7]/F");
  TBranch *b_eta_ak7 = tree->Branch("eta_ak7", eta_ak7, "eta_ak7[njet_ak7]/F");
  TBranch *b_phi_ak7 = tree->Branch("phi_ak7", phi_ak7, "phi_ak7[njet_ak7]/F");
  TBranch *b_E_ak7 = tree->Branch("E_ak7", E_ak7, "E_ak7[njet_ak7]/F");   
  TBranch *b_area_ak7 = tree->Branch("area_ak7", area_ak7, "area_ak7[njet_ak7]/F");
  TBranch *b_jec_ak7 = tree->Branch("jec_ak7", jec_ak7, "jec_ak7[njet_ak7]/F");
  TBranch *b_ak7_to_ak4 = tree->Branch("ak7_to_ak4", ak7_to_ak4, "ak7_to_ak4[njet_ak7]/I");


  TBranch *b_ncomp = tree->Branch("ncomp", &ncomp, "ncomp/b"); // Between 1 and kMaxNComp
  TBranch *b_chf = tree->Branch("chf", chf, "chf[ncomp]/b");   
  TBranch *b_nhf = tree->Branch("nhf", nhf, "nhf[ncomp]/b");   
  TBranch *b_phf = tree->Branch("phf", phf, "phf[ncomp]/b");   
  TBranch *b_elf = tree->Branch("elf", elf, "elf[ncomp]/b");   
  TBranch *b_muf = tree->Branch("muf", muf, "muf[ncomp]/b");   

  TBranch *b_hf_hf  = tree->Branch("hf_hf", hf_hf, "hf_hf[ncomp]/b");   
  TBranch *b_hf_phf = tree->Branch("hf_phf", hf_phf, "hf_phf[ncomp]/b");   
  TBranch *b_hf_hm  = tree->Branch("hf_hm", hf_hm, "hf_hm[ncomp]/b");    
  TBranch *b_hf_phm = tree->Branch("hf_phm", hf_phm, "hf_phm[ncomp]/b");
  
  TBranch *b_chm = tree->Branch("chm", chm, "chm[ncomp]/b");   
  TBranch *b_nhm = tree->Branch("nhm", nhm, "nhm[ncomp]/b");   
  TBranch *b_phm = tree->Branch("phm", phm, "phm[ncomp]/b");   
  TBranch *b_elm = tree->Branch("elm", elm, "elm[ncomp]/b");   
  TBranch *b_mum = tree->Branch("mum", mum, "mum[ncomp]/b");
   
  TBranch *b_hof      = tree->Branch("hof", hof, "hof[ncomp]/b");   
  TBranch *b_beta     = tree->Branch("beta", beta, "beta[ncomp]/b");   
  TBranch *b_bstar    = tree->Branch("bstar", bstar, "bstar[ncomp]/b");
  TBranch *b_qgl      = tree->Branch("qgl", qgl, "qgl[ncomp]/b");

  TBranch *b_csv    = tree->Branch("csv", csv, "csv[ncomp]/b");
  TBranch *b_pfl    = tree->Branch("pfl", pfl, "pfl[ncomp]/B");
  TBranch *b_hfl    = tree->Branch("hfl", hfl, "hfl[ncomp]/B");


  
  if (isMC) {
    TBranch *b_ngen     = tree->Branch("ngen", &ngen, "ngen/b");
    TBranch *b_gen_pt   = tree->Branch("gen_pt", gen_pt, "gen_pt[ngen]/F");
    TBranch *b_gen_eta  = tree->Branch("gen_eta", gen_eta, "gen_eta[ngen]/F");
    TBranch *b_gen_phi  = tree->Branch("gen_phi", gen_phi, "gen_phi[ngen]/F");
    TBranch *b_gen_E    = tree->Branch("gen_E", gen_E, "gen_E[ngen]/F");
  }

  TBranch *b_run      = tree->Branch("run", &run, "run/i");
  TBranch *b_lumi     = tree->Branch("lumi", &lumi, "lumi/i");
  TBranch *b_event    = tree->Branch("event", &event, "event/l");

  TBranch *b_ntrg         = tree->Branch("ntrg", &ntrg, "ntrg/i");
  TBranch *b_trgfired     = tree->Branch("trgfired", &trgfired);
  TBranch *b_prescales    = tree->Branch("prescales", &prescales);

  TBranch *b_met      = tree->Branch("met", &met, "met/F");
  TBranch *b_sumet    = tree->Branch("sumet", &sumet, "sumet/F");
  TBranch *b_rho      = tree->Branch("rho", &rho, "rho/F");
  TBranch *b_pthat    = tree->Branch("pthat", &pthat, "pthat/F");
  TBranch *b_mcweight = tree->Branch("mcweight", &mcweight, "mcweight/F");
  



  assert(fChain_ak4 && "AK4 tree invalid!" );
  assert(fChain_ak7 && "AK7 tree invalid!" );
  
  fChain_ak4->SetBranchStatus("*",0);  // disable all branches
  fChain_ak7->SetBranchStatus("*",0);


  // To begin, read only these 3 variables (AK7) Does it work?
  fChain_ak7->SetBranchStatus("EvtHdr_.mRun",1); // run
  fChain_ak7->SetBranchStatus("EvtHdr_.mLumi",1); // lumi
  fChain_ak7->SetBranchStatus("EvtHdr_.mEvent",1); // event



  // Enable the remaining variables
  fChain_ak4->SetBranchStatus("PFJetsCHS_",1); // njet
  fChain_ak4->SetBranchStatus("PFJetsCHS_.P4_.fCoordinates.f*",1); // Four-momentum
  
  fChain_ak4->SetBranchStatus("PFJetsCHS_.tightID_",1); // tightID
  fChain_ak4->SetBranchStatus("PFJetsCHS_.area_",1); // area
  fChain_ak4->SetBranchStatus("PFJetsCHS_.cor_",1); // jec
  
  // Composition values
  fChain_ak4->SetBranchStatus("PFJetsCHS_.chf_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.nhf_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.nemf_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.cemf_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.muf_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.hf_*",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.chm_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.nhm_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.phm_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.elm_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.mum_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.QGtagger_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.beta_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.betaStar_",1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.hof_",1);

  fChain_ak4->SetBranchStatus("PFJetsCHS_.CSV_", 1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.partonFlavour_", 1);
  fChain_ak4->SetBranchStatus("PFJetsCHS_.hadronFlavour_", 1);
  
  if (isMC) { // MC
    fChain_ak4->SetBranchStatus("GenJets_",1); // ngen
    fChain_ak4->SetBranchStatus("GenJets_.fCoordinates.f*",1);
  }
  
  fChain_ak4->SetBranchStatus("TriggerDecision_",1);
  fChain_ak4->SetBranchStatus("L1Prescale_",1);
  fChain_ak4->SetBranchStatus("HLTPrescale_",1);

  fChain_ak4->SetBranchStatus("EvtHdr_.mRun",1); // run
  fChain_ak4->SetBranchStatus("EvtHdr_.mLumi",1); // lumi
  fChain_ak4->SetBranchStatus("EvtHdr_.mEvent",1); // event
  fChain_ak4->SetBranchStatus("PFMet_.et_",1); // met
  fChain_ak4->SetBranchStatus("PFMet_.sumEt_",1); // sumet
  fChain_ak4->SetBranchStatus("EvtHdr_.mPFRho",1); // rho
  fChain_ak4->SetBranchStatus("EvtHdr_.mPthat",1); // pthat
  fChain_ak4->SetBranchStatus("EvtHdr_.mWeight",1); // mcweight


  fChain_ak7->SetBranchStatus("PFJetsCHS_", 1);
  fChain_ak7->SetBranchStatus("PFJetsCHS_.P4_.fCoordinates.f*", 1);
  fChain_ak7->SetBranchStatus("PFJetsCHS_.cor_", 1);
  fChain_ak7->SetBranchStatus("PFJetsCHS_.area_", 1);



  //std::cout << "Mapping AK7 events to AK4 events for data" << std::endl;
  map<Int_t, map<Int_t, map<UInt_t, Long64_t> > > ak7entry;

  /*
  Long64_t nentries7 = fChain_ak7->GetEntries();
  for (Long64_t jentry7 = 0; jentry7 < nentries7; jentry7++) {
    fChain_ak7->GetEntry(jentry7);

    // Check for duplicates
    assert(ak7entry[EvtHdr__mRun_ak7][EvtHdr__mLumi_ak7][EvtHdr__mEvent_ak7] == 0);

    // Save entry number
    ak7entry[EvtHdr__mRun_ak7][EvtHdr__mLumi_ak7][EvtHdr__mEvent_ak7] = jentry7;
  } 
  std::cout << "Found mapping for " << ak7entry.size() << " runs" << std::endl;
  */

  // Helper variables
  TLorentzVector p4, p4_ak4, p4_ak7, p4gen;
  
  // Total number of events
  Long64_t nentries = fChain_ak4->GetEntries(); 
  std::cout << "Total entries: " << nentries << std::endl;

  // DEBUG:
  // Change number of events here
  nentries = 1000000;

  // Process triggers
  // Shorten the name of PF triggers
  if (!isMC) {    

    // Trigger names
    assert(TriggerNames);

    auto trgAxis = TriggerNames->GetXaxis();

    for (int i = trgAxis->GetFirst(); i != trgAxis->GetLast(); ++i) {

      std::regex pfstr("HLT_PFJet([0-9]+)_v[0-9]+");
      std::string shortTrig;

      std::string trgName = trgAxis->GetBinLabel(i);
      if (trgName.size() == 0) 
        continue;
      else if (std::regex_match(trgName, pfstr)) 
        shortTrig = std::regex_replace(trgName, pfstr, "jt$1", std::regex_constants::format_no_copy);
      else
        shortTrig = "";

      _pfTriggers.push_back(shortTrig);
    }
  }


  {
    // Full redoing of JEC    
    const char *s;
    const char *p = "/afs/cern.ch/user/m/mhaapale/work/public/TuplePress/CondFormats/JetMETObjects/data/";
    std::cout << "Path: " << p << std::endl;
    std::string gt("Spring16_25nsV8p2_");
    gt = gt + ( isMC ? "MC_" : "DATA_" );
    const char *t = gt.c_str();
    std::cout << t << std::endl;
    const char *a = "AK4";

    _JEC = 0;

    cout << "Loading "<<p<<"PFchs JEC" << endl;
    s = Form("%s%sL1FastJet_%sPFchs.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l1 = new JetCorrectorParameters(s);
    s = Form("%s%sL2Relative_%sPFchs.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2 = new JetCorrectorParameters(s);
    s = Form("%s%sL3Absolute_%sPFchs.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l3 = new JetCorrectorParameters(s);
    s = Form("%s%sL2L3Residual_%sPFchs.txt",p,t,a); cout<<s<<endl<<flush;
    JetCorrectorParameters *par_l2l3res = new JetCorrectorParameters(s);


    vector<JetCorrectorParameters> vpar;
    vpar.push_back(*par_l1);
    vpar.push_back(*par_l2);
    vpar.push_back(*par_l3);
    if (true) vpar.push_back(*par_l2l3res); // Only for data
    _JEC = new FactorizedJetCorrector(vpar);

  } // JEC redone


  // Iterating over the events
  for (Long64_t jentry = 0; jentry != nentries; ++jentry) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain_ak4->GetEntry(jentry);

    // Event identification
    run = EvtHdr__mRun;
    event = EvtHdr__mEvent;
    lumi = EvtHdr__mLumi;

    // MET, SuMET, rho
    met = PFMet__et_;
    sumet = PFMet__sumEt_;
    rho = EvtHdr__mPFRho;




    // Sort jets wrt. pT ??

    /*

    */

    // Number of leading jets for which composition variables are saved
    ncomp = min(kMaxNComp, PFJetsCHS__);

    // Jet index in the output arrays (after pT cut)
    Int_t i_out = 0;

    for (Int_t i = 0; i != PFJetsCHS__; ++i) {

      p4.SetPxPyPzE(  PFJetsCHS__P4__fCoordinates_fX[i], PFJetsCHS__P4__fCoordinates_fY[i],
              PFJetsCHS__P4__fCoordinates_fZ[i], PFJetsCHS__P4__fCoordinates_fT[i]);

      Bool_t cur_id = PFJetsCHS__tightID_[i];
      Float_t cur_area = PFJetsCHS__area_[i];
      Float_t cur_jec = PFJetsCHS__cor_[i]; 



      // Update JEC from text file
      if (true) {
        _JEC->setRho(rho);
        _JEC->setJetA(cur_area);
        _JEC->setJetPt(p4.Pt());
        _JEC->setJetE(p4.E());
        _JEC->setJetEta(p4.Eta());

        Float_t new_jec =  _JEC->getCorrection();

        // Redo JEC
        p4 *= new_jec/cur_jec;
        cur_jec = new_jec;

      }
    

      // pT selection after JEC
      const Float_t minPt = 15;
      if (p4.Pt() > minPt) { // Break first jet with pT<15 GeV ? 

        assert(i_out < kMaxNjet && "Number of jets to be saved exceeds 'kMaxNjet' !");

        pt[i_out] = p4.Pt();
        eta[i_out] = p4.Eta();
        phi[i_out] = p4.Phi();
        E[i_out] = p4.E();

        tightID[i_out] = cur_id;
        area[i_out] = cur_area;
        jec[i_out] = cur_jec; 

        // Jet composition, only for ncomp leading jets
        if (i_out < ncomp) {

          assert( i_out < kMaxNComp && "Number of composition variables to be saved exceeds 'kMaxNComp'" );

          chf[i_out]  = mapFloatToUChar(PFJetsCHS__chf_[i]);
          //std::cout << PFJetsCHS__chf_[i] << ' ' << chf[i_out] << ' ' << mapUCharToFloat(chf[i_out]) << std::endl;
          
          nhf[i_out]  = mapFloatToUChar(PFJetsCHS__nhf_[i]);
          phf[i_out]  = mapFloatToUChar(PFJetsCHS__nemf_[i]);
          elf[i_out]  = mapFloatToUChar(PFJetsCHS__cemf_[i]);
          muf[i_out]  = mapFloatToUChar(PFJetsCHS__muf_[i]);
          hf_hf[i_out]    = mapFloatToUChar(PFJetsCHS__hf_hf_[i]);
          hf_phf[i_out]   = mapFloatToUChar(PFJetsCHS__hf_phf_[i]);

          hf_hm[i_out]    = mapIntToUChar(PFJetsCHS__hf_hm_[i]);
          hf_phm[i_out]   = mapIntToUChar(PFJetsCHS__hf_phm_[i]);
          chm[i_out]      = mapIntToUChar(PFJetsCHS__chm_[i]);
          //std::cout << PFJetsCHS__chm_[i] << ' ' << chm[i_out] << ' ' << mapUCharToInt(chm[i_out]) << std::endl;
    
          nhm[i_out]  = mapIntToUChar(PFJetsCHS__nhm_[i]);
          phm[i_out]  = mapIntToUChar(PFJetsCHS__phm_[i]);
          elm[i_out]  = mapIntToUChar(PFJetsCHS__elm_[i]);
          mum[i_out]  = mapIntToUChar(PFJetsCHS__mum_[i]);

          qgl[i_out]   = mapFloatToUChar(PFJetsCHS__QGtagger_[i]);
          beta[i_out]  = mapFloatToUChar(PFJetsCHS__beta_[i]);
          bstar[i_out] = mapFloatToUChar(PFJetsCHS__betaStar_[i]);
          hof[i_out]   = mapFloatToUChar(PFJetsCHS__hof_[i]);

          csv[i_out] = mapFloatToUChar(PFJetsCHS__CSV_[i]);
 
          hfl[i_out] = mapIntToChar(PFJetsCHS__hadronFlavour_[i]);
          pfl[i_out] = mapIntToChar(PFJetsCHS__partonFlavour_[i]);
    
        }

        // Increment index
        ++i_out;
      }
    }

    // Ugly hack to save the number of jets saved
    njet = mapIntToUChar(i_out); 


    /*
    // Read corresponding AK7 event using the entry mapping
    fChain_ak7->GetEntry(ak7entry[EvtHdr__mRun][EvtHdr__mLumi][EvtHdr__mEvent]);

    // Sanity check
    if (EvtHdr__mRun != EvtHdr__mRun_ak7 || 
      EvtHdr__mEvent != EvtHdr__mEvent_ak7 ||
      EvtHdr__mLumi != EvtHdr__mLumi_ak7) {
      
      assert(false && "Mismatch between AK4 and AK7 events!!");
    }

    
    // Keep only four leading jets
    njet_ak7 = min(PFJetsCHS_ak7__, 3);
    for (Int_t i = 0; i != njet_ak7; ++i) {

      p4_ak7.SetPxPyPzE(  PFJetsCHS__P4__fCoordinates_fX_ak7[i], PFJetsCHS__P4__fCoordinates_fY_ak7[i],
                PFJetsCHS__P4__fCoordinates_fZ_ak7[i], PFJetsCHS__P4__fCoordinates_fT_ak7[i]);

      // 4-momentum (corrected!)
      // Do the correction again??
      pt_ak7[i]   = p4_ak7.Pt();   
      eta_ak7[i]  = p4_ak7.Eta();
      phi_ak7[i]  = p4_ak7.Phi();
      E_ak7[i]    = p4_ak7.E();

      // Area and jet energy correction
      area_ak7[i] = PFJetsCHS__area_ak7_[i];
      jec_ak7[i]  = PFJetsCHS__cor_ak7_[i]; 

      // Matching AK7 jet to AK4
      ak7_to_ak4[i] = -1;

      // Search AK4 jet with minimum distance to this PFjet   
      float rMin(999);
      for (Int_t ak4_i = 0; ak4_i != njet; ++ak4_i) {
        
        // Initialize AK4 jet
        p4_ak4.SetPtEtaPhiE(  pt[ak4_i],  eta[ak4_i],
                    phi[ak4_i], E[ak4_i]);  
        
        // Distance between jets
        double deltaR = p4_ak7.DeltaR(p4_ak4);

        if (deltaR < rMin) {
          rMin = deltaR;
          ak7_to_ak4[i] = ak4_i;
        }
      }
    }
    */    

    // MC truth jets
    if (isMC) {

      ngen = mapIntToUChar(GenJets__);
      pthat = EvtHdr__mPthat;
      mcweight = EvtHdr__mWeight;

      for (UChar_t i = 0; i != ngen; ++i) {
        p4gen.SetPxPyPzE(   GenJets__fCoordinates_fX[i], GenJets__fCoordinates_fY[i],
                  GenJets__fCoordinates_fZ[i], GenJets__fCoordinates_fT[i]);
        
        gen_pt[i] = p4gen.Pt(); 
        gen_eta[i] = p4gen.Eta();
        gen_phi[i] = p4gen.Phi();
        gen_E[i] = p4gen.E();

      }
    }
    

    // Save trigger bits for data
    // Partially ripped from fillHistos.C
    if (!isMC) {

      trgfired.clear();
      prescales.clear();
      for (unsigned int itrg = 0; itrg != _pfTriggers.size(); ++itrg ) {
        
        std::string strg = _pfTriggers[itrg];
        bool pass = ( TriggerDecision_[itrg] == 1 ) && ( strg != "" ); // -1, 0 or 1
        
        if (pass) {
          int prsc  = HLTPrescale_[itrg]*L1Prescale_[itrg];
          if (prsc > 0) {
            trgfired.push_back(strg);
            prescales.push_back(prsc);

            //std::cout << prsc << std::endl;
          }
          else {
            std::cout << "Error for trigger " << strg << " prescales: "
                      << "L1  =" << L1Prescale_[itrg]
                      << "HLT =" << HLTPrescale_[itrg] << std::endl;
          }
        }
      }
    }


    // Save event into file
    tree->Fill();

  }

  fout->cd();
  tree->Print();
  tree->Write();
  fout->Close();
}
