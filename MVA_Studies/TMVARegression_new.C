// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVARegression                                                     *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l TMVARegression.C\(\"LD,MLP\"\)                                      *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVAReg.root" can be analysed with the use of dedicated       *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 **********************************************************************************/

#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVARegGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#endif

using namespace TMVA;

void TMVARegression_new( TString myMethodList = "" ){

  //#! Parameters
  myMethodList = "BDT";
  TString outfileName( "TMVAReg.root" );  
  TString fname = "/fdata/hepx/store/user/lpernie/Hhh_For_Regression/TT_1k.root";
  TString tree_name="DiHiggsWWBBAna/evtree";  
  TCut mycut = "b1jet_leptonPtRel>0 && b1jet_leptonPtRel<5 && met_pt>20 && muon1_pt>20 && TMath::Abs(muon1_eta)<2.4 && muon2_pt>10 && TMath::Abs(muon2_eta)<2.4 && mass_l1l2>12 && b1jet_pt>20 && TMath::Abs(b1jet_eta)<2.4 && b2jet_pt>20 && TMath::Abs(b2jet_eta)<2.4 && mass_trans>10 && mass_l1l2<76 && mass_b1b2>30 && dR_l1l2b1b2>0.2 && dR_l1l2b1b2<4 && dR_l1l2<3.8 && dR_b1b2<3.8";
  // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
  // if you use your private .rootrc, or run from a different directory, please copy the 
  // corresponding lines from .rootrc

  // methods to be processed can be given as an argument; use format:
  //
  // mylinux~> root -l TMVARegression.C\(\"myMethod1,myMethod2,myMethod3\"\)
  //

  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 0;
  Use["PDEFoam"]         = 1; 
  Use["KNN"]             = 1;
  // --- Linear Discriminant Analysis
  Use["LD"]              = 1;
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 1;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  // --- Neural Network
  Use["MLP"]             = 1; 
  // --- Support Vector Machine 
  Use["SVM"]             = 0;
  // --- Boosted Decision Trees
  Use["BDT"]             = 0;
  Use["BDTG"]            = 1;
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVARegression" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
	std::string regMethod(mlist[i]);

	if (Use.find(regMethod) == Use.end()) {
	  std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	  for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	  std::cout << std::endl;
	  return;
	}
	Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------
  // --- Here the preparation phase begins
  // Create a new root output file
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory will
  // then run the performance analysis for you.
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/ 
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in 
  // front of the "Silent" argument in the option string
  TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile, 
	"!V:!Silent:Color:DrawProgressBar" );

  // If you wish to modify default settings 
  // (please check "src/Config.h" to see all available global options)
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

  //#!Define the input variables that shall be used for the MVA training
  //  note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  //  [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  factory->AddVariable( "numOfVertices",      "numOfVertices", "units", 'I' );
  factory->AddVariable( "b1jet_pt",           "b1jet_pt", "units", 'F' );
  factory->AddVariable( "b1jet_eta",          "b1jet_eta", "units", 'F' );
  factory->AddVariable( "b1jet_mt",           "b1jet_mt", "units", 'F' );
  factory->AddVariable( "b1jet_leadTrackPt",  "b1jet_leadTrackPt", "units", 'F' );
  factory->AddVariable( "b1jet_leptonDeltaR", "b1jet_leptonDeltaR", "units", 'F' );
  //factory->AddVariable( "b1jet_leptonPtRel",  "b1jet_leptonPtRel", "units", 'F' );
  factory->AddVariable( "b1jet_leptonPt",     "b1jet_leptonPt", "units", 'F' );
  factory->AddVariable( "b1jet_vtxPt",        "b1jet_vtxPt", "units", 'F' );
  factory->AddVariable( "b1jet_vtxMass",      "b1jet_vtxMass", "units", 'F' );
  factory->AddVariable( "b1jet_vtxNtracks",   "b1jet_vtxNtracks", "units", 'F' );
  factory->AddVariable( "b1jet_neHEF",        "b1jet_neHEF", "units", 'F' );
  factory->AddVariable( "b1jet_neEmEF",       "b1jet_neEmEF", "units", 'F' );
  factory->AddVariable( "b1jet_vtx3DSig",     "b1jet_vtx3DSig", "units", 'F' );
  factory->AddVariable( "b1jet_vtx3DVal",     "b1jet_vtx3DVal", "units", 'F' );

  // You can add so-called "Spectator variables", which are not used in the MVA training, 
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
  // input variables, the response values of all trained MVAs, and the spectator variables
  factory->AddSpectator( "spec1:=mass_b1b2",  "mass_b1b2", "units", 'F' );

  // Add the variable carrying the regression target
  factory->AddTarget( "b1genjet_pt/b1jet_pt" ); 

  // It is also possible to declare additional targets for multi-dimensional regression, ie:
  // -- factory->AddTarget( "fvalue2" );
  // BUT: this is currently ONLY implemented for MLP

  // Read training and test data (see TMVAClassification for reading ASCII files)
  // load the signal and background event samples from ROOT trees
  TFile *input(0);
  if (!gSystem->AccessPathName( fname ))
    input = TFile::Open( fname ); // check if file in local directory exists
  if (!input) {
    std::cout << "ERROR: could not open data file" << std::endl;
    exit(1);
  }
  std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;

  // --- Register the regression tree
  TTree *regTree = (TTree*)input->Get(tree_name.Data());

  // global event weights per tree (see below for setting event-wise weights)
  Double_t regWeight  = 1.0;   

  // You can add an arbitrary number of regression trees
  factory->AddRegressionTree( regTree, regWeight );

  // This would set individual event weights (the variables defined in the 
  // expression need to exist in the original TTree)
  factory->SetWeightExpression( "reweighting", "MuonPtreweighting" );

  // Apply additional cuts on the signal and background samples (can be different)
  // tell the factory to use all remaining events in the trees after training for testing:
  factory->PrepareTrainingAndTestTree( mycut, 
	"nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );

  // If no numbers of events are given, half of the events in the tree are used 
  // for training, and the other half for testing:
  //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  

  // ---- Book MVA methods
  //
  // please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

  // PDE - RS method
  if (Use["PDERS"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERS", 
	  "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=40:NEventsMax=60:VarTransform=None" );
  // And the options strings for the MinMax and RMS methods, respectively:
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   

  if (Use["PDEFoam"])
    factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam", 
	  "!H:!V:MultiTargetRegression=F:TargetSelection=Mpv:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Compress=T:Kernel=None:Nmin=10:VarTransform=None" );

  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"])
    factory->BookMethod( TMVA::Types::kKNN, "KNN", 
	  "nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

  // Linear discriminant
  if (Use["LD"])
    factory->BookMethod( TMVA::Types::kLD, "LD", 
	  "!H:!V:VarTransform=None" );

  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"]) 
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
	  "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=MC:SampleSize=100000:Sigma=0.1:VarTransform=D" );

  if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options) .. the formula of this example is good for parabolas
    factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
	  "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:PopSize=100:Cycles=3:Steps=30:Trim=True:SaveBestGen=1:VarTransform=Norm" );

  if (Use["FDA_MT"]) 
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
	  "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

  if (Use["FDA_GAMT"]) 
    factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
	  "!H:!V:Formula=(0)+(1)*x0+(2)*x1:ParRanges=(-100,100);(-100,100);(-100,100):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

  // Neural network (MLP)
  if (Use["MLP"])
    factory->BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

  // Support Vector Machine
  if (Use["SVM"])
    factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

  // Boosted Decision Trees
  if (Use["BDT"])
    factory->BookMethod( TMVA::Types::kBDT, "BDT",
	  "!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );

  if (Use["BDTG"])
    factory->BookMethod( TMVA::Types::kBDT, "BDTG",
	  "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );
  // --------------------------------------------------------------------------------------------------

  // ---- Now you can tell the factory to train, test, and evaluate the MVAs

  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();    

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVARegression is done!" << std::endl;      

  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVARegGui( outfileName );
}
