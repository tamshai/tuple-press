{

	// Compile to get faster execution

	// Run in ROOT: root mk_LocalProducer.C+

	gROOT->ProcessLine(".L LocalOpenDataTreeProducer.C+g");

	// Input tuples
	const std::string treeName = "ProcessedTree";
	const int numFiles = 3;

	// Here are the input files!
	const char *fileNames[numFiles] = { "root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016-RunC-v2-part1.root",
	                                   "root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016-RunC-v2-part2.root",
	                                   "root://eoscms.cern.ch//store/group/phys_smp/Multijet/13TeV/Data/2016/Ntuples-Data-2016-RunC-v2-part3.root",
	                                };

	// AK4 and AK7 trees
	TChain *chain_ak4 = new TChain(("ak4/" + treeName).c_str());
	TChain *chain_ak7 = new TChain(("ak7/" + treeName).c_str());

	for (int i = 0; i < numFiles; ++i) {
		chain_ak4->Add(fileNames[i]);
		chain_ak7->Add(fileNames[i]);
	}

	// Retrieve trigger names
	TFile *f = TFile::Open(fileNames[0]);
	TDirectory *dir = (TDirectory*)f->Get("ak4");
	TH1F* trgNames;
	dir->GetObject("TriggerNames", trgNames);    

	// Process the trees
	LocalOpenDataTreeProducer(chain_ak4, chain_ak7, trgNames);

}
