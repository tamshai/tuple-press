
// Run in ROOT: root mk_LocalProducer.C+

{
	
	// Make this global variable?? 
	// (Otherwise has to be defined again in the actual ROOT script)
	const TString JECpath( "/afs/cern.ch/user/m/mhaapale/work/public/TuplePress/CondFormats/JetMETObjects" );
	

	// Include path to local JEC libraries 
	gInterpreter->AddIncludePath(JECpath);
	
	TString  curPath( gROOT->GetMacroPath() );
	gROOT->SetMacroPath( curPath + ":" + JECpath);


	// Compile macros
    gROOT->ProcessLine(".L src/Utilities.cc+");
    gROOT->ProcessLine(".L src/JetCorrectorParameters.cc+");
    gROOT->ProcessLine(".L src/SimpleJetCorrector.cc+");
    gROOT->ProcessLine(".L src/FactorizedJetCorrector.cc+");

	gROOT->ProcessLine(".L LocalProducer.C+g");



	// MC or DATA
	Bool_t isMC = false;

	// Name of the tree
	const std::string treeName = "ProcessedTree";

	// AK4 and AK7 trees
	TChain *chain_ak4 = new TChain(("ak4/" + treeName).c_str());
	TChain *chain_ak7 = new TChain(("ak7/" + treeName).c_str());

	// Here are the input files!
	std::vector<std::string> fileNames;

	if (!isMC) {

		std::vector<std::string> names =
			{
				"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart1.root",
              	"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart2.root",
               	"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart3.root",
               	"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016RunG-PromptReco-80Xpart4.root"
            };

	    fileNames = names;

	}
	else if (isMC) { 
		std::vector<std::string> names =
			{ 	
				//"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/MC/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp_80X.root" ---> Corrupted, do not use!
				"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/MC/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp_80X_v2.root"
            	//"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/MC/Ntuples-MC-Pythia8-Flat15to7000-25ns-CUETM1-13TeV.root"
            };
	    fileNames = names;
	}


	for (auto name : fileNames) {
		chain_ak4->Add(name.c_str());
		chain_ak7->Add(name.c_str());
	}

	// Retrieve trigger names
	TFile *f = TFile::Open(fileNames[0].c_str());
	TDirectory *dir = (TDirectory*)f->Get("ak4");
	TH1F* trgNames;
	dir->GetObject("TriggerNames", trgNames);    

	// Process the trees
	LocalProducer(chain_ak4, chain_ak7, trgNames, isMC);

}
