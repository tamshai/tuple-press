
#ifndef JetP4_h
#define JetP4_h

class JetP4
{

public:
	JetP4(Int_t njet) : fNjet(njet) {}
	~JetP4() {}

  setPt(Double_t p, Int_t index);
  setEta(Double_t e, Int_t index);
  setPhi()

/*
	// Setters

  // PF AK5 jets
  UChar_t njet;
  Float_t pt[kMaxNjet];
  Float_t eta[kMaxNjet];
  Float_t phi[kMaxNjet];
  Float_t E[kMaxNjet];
  Bool_t tightID[kMaxNjet];
  Float_t area[kMaxNjet];
  Float_t jec[kMaxNjet];

*/


	// PF AK5 jets, return TLorentzVector??
	//UChar_t njet;
	//Float_t pt[kMaxNjet];
	//Float_t eta[kMaxNjet];
	//Float_t phi[kMaxNjet];
	//Float_t E[kMaxNjet];
	//Bool_t tightID[kMaxNjet];
	//Float_t area[kMaxNjet];
	//Float_t jec[kMaxNjet];
	//Int_t igen[kMaxNjet];



private:


  // PF AK5 jets
  enum { kMaxNjet = 32 };
  UChar_It fNjet;
  UShort_t fPt[kMaxNjet];
  UChar_t fEta[kMaxNjet];
  UChar_t fPhi[kMaxNjet];
  Float_t fE[kMaxNjet];
  /*
  Bool_t tightID[kMaxNjet];
  Float_t area[kMaxNjet];
  Float_t jec[kMaxNjet];
  Int_t igen[kMaxNjet];
  */
};




#endif