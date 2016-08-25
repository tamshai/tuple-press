{

	// Compile to get faster execution

	// Run in ROOT: root mk_LocalProducer.C+

	gROOT->ProcessLine(".L LocalOpenDataTreeProducer.C+g");

	// Input tuples
	const std::string treeName = "ProcessedTree";

	// AK4 and AK7 trees
	TChain *chain_ak4 = new TChain(("ak4/" + treeName).c_str());
	TChain *chain_ak7 = new TChain(("ak7/" + treeName).c_str());


	// Here are the input files! (https://twiki.cern.ch/twiki/bin/viewauth/CMS/InclusiveJetsLegacy)
	std::vector<std::string> fileNames  = 
		{ 	//"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016-RunC-v2-part1.root",
        	//"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016-RunC-v2-part2.root",
        	//"root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016-RunC-v2-part3.root",
    		"tuples0.root"
		};

	for (std::string name : fileNames) {
		chain_ak4->Add(name.c_str());
		chain_ak7->Add(name.c_str());
	}

	// Retrieve list of trigger names
	TFile *f = TFile::Open(fileNames[0].c_str());
	TDirectory *dir = (TDirectory*)f->Get("ak4");
	TH1F* trgNames;
	dir->GetObject("TriggerNames", trgNames);    

	// Process the trees
	LocalOpenDataTreeProducer(chain_ak4, chain_ak7, trgNames);

}
