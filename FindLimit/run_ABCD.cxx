// COntinuous Likelihood-ABCD
// author: Stefano Giagu <stefano.giagu@cern.ch>
// Version: 4.0 (cleaned up as demo) June 27th, 2016
// Modified by G. Watts for the CalRatio analysis

/*
 Two variables x,y create 4 regions: A,B,C,D with A signal dominated region and
 B,C,D control regions (sidebands).

 ^ y (sumpt)
 |
 |--------------------------+
 |    C         |     B     |
 |              |           |
 |              |           |
 |              |           |
 |              |           |
 |--------------+-----------|
 |    D         |     A     |
 |              |           |
 |              |  Signal   |
 |              |           |
 +-----------------------------> x (deltaPhi)
 ^

 ABCD Ansatz A:B = D:C -->  A = B * D/C

 */

#include <iostream>
#include <fstream>
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooRealSumPdf.h"
#include "RooDataSet.h"
#include "TTree.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "RooRandom.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/NumberCountingUtils.h"
#include "TROOT.h"
#include "TSystem.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumEventsTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;


// StandardHypoTestInvDemo: modified (ouput & some settings) standard rooStat macro for performing an inverted  hypothesis test for computing an interval from L.Moneta / K.Cramer

bool plotHypoTestResult = true;          // plot test statistic result at each point
bool writeResult = true;                 // write HypoTestInverterResult in a file 
TString resultFileName;                  // file with results (by default is built automatically using the workspace input file name)
bool optimize = true;                    // optmize evaluation of test statistic 
bool useVectorStore = true;              // convert data to use new roofit data store 
bool generateBinned = false;             // generate binned data sets 
bool noSystematics = false;              // force all systematics to be off (i.e. set all nuisance parameters as constat
                                         // to their nominal values)
double nToysRatio = 2;                   // ratio Ntoys S+b/ntoysB
double maxPOI = -1;                      // max value used of POI (in case of auto scan) 
bool useProof = false;                    // use Proof Light when using toys (for freq or hybrid)
int nworkers = 21;                        // number of worker for Proof
bool rebuild = false;                    // re-do extra toys for computing expected limits and rebuild test stat
                                         // distributions (N.B this requires much more CPU (factor is equivalent to nToyToRebuild)
int nToyToRebuild = 100;                 // number of toys used to rebuild 
int initialFit = -1;                     // do a first  fit to the model (-1 : default, 0 skip fit, 1 do always fit) 
int randomSeed = -1;                     // random seed (if = -1: use default value, if = 0 always random )
                                         // NOTE: Proof uses automatically a random seed

int nAsimovBins = 0;                     // number of bins in observables used for Asimov data sets (0 is the default and it is given by workspace, typically is 100)

bool reuseAltToys = false;                // reuse same toys for alternate hypothesis (if set one gets more stable bands)

std::string massValue = "";              // extra string to tag output file of result 
std::string  minimizerType = "";                  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
int   printLevel = 0;                    // print level for debugging PL test statistics and calculators  

// internal class to run the inverter and more

namespace RooStats { 

   class HypoTestInvTool{

   public:
      HypoTestInvTool();
      ~HypoTestInvTool(){};

      HypoTestInverterResult * 
      RunInverter(Int_t enne, TString esse, RooWorkspace * w, 
                  const char * modelSBName, const char * modelBName, 
                  const char * dataName,
                  int type,  int testStatType, 
                  bool useCLs, 
                  int npoints, double poimin, double poimax, int ntoys, 
                  bool useNumberCounting = false, 
                  const char * nuisPriorName = 0);



      Double_t 
      AnalyzeResult( HypoTestInverterResult * r,
                     int calculatorType,
                     int testStatType, 
                     bool useCLs,  
                     int npoints,
                     const char * fileNameBase = 0 );

      void SetParameter(const char * name, const char * value);
      void SetParameter(const char * name, bool value);
      void SetParameter(const char * name, int value);
      void SetParameter(const char * name, double value);

   private:

      bool mPlotHypoTestResult;
      bool mWriteResult;
      bool mOptimize;
      bool mUseVectorStore;
      bool mGenerateBinned;
      bool mUseProof;
      bool mRebuild;
      bool mReuseAltToys; 
      int     mNWorkers;
      int     mNToyToRebuild;
      int     mPrintLevel;
      int     mInitialFit; 
      int     mRandomSeed; 
      double  mNToysRatio;
      double  mMaxPoi;
      int mAsimovBins;
      std::string mMassValue;
      std::string mMinimizerType;                  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
      TString     mResultFileName; 
   };

} // end namespace RooStats

RooStats::HypoTestInvTool::HypoTestInvTool() : mPlotHypoTestResult(true),
                                               mWriteResult(false),
                                               mOptimize(true),
                                               mUseVectorStore(true),
                                               mGenerateBinned(false),
                                               mUseProof(false),
                                               mRebuild(false),
                                               mReuseAltToys(false),
                                               mNWorkers(4),
                                               mNToyToRebuild(100),
                                               mPrintLevel(0),
                                               mInitialFit(-1),
                                               mRandomSeed(-1),
                                               mNToysRatio(2),
                                               mMaxPoi(-1),
                                               mAsimovBins(0),
                                               mMassValue(""),
                                               mMinimizerType(""),
                                               mResultFileName() {
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, bool value){
   //
   // set boolean parameters
   //

   std::string s_name(name);

   if (s_name.find("PlotHypoTestResult") != std::string::npos) mPlotHypoTestResult = value;
   if (s_name.find("WriteResult") != std::string::npos) mWriteResult = value;
   if (s_name.find("Optimize") != std::string::npos) mOptimize = value;
   if (s_name.find("UseVectorStore") != std::string::npos) mUseVectorStore = value;
   if (s_name.find("GenerateBinned") != std::string::npos) mGenerateBinned = value;
   if (s_name.find("UseProof") != std::string::npos) mUseProof = value;
   if (s_name.find("Rebuild") != std::string::npos) mRebuild = value;
   if (s_name.find("ReuseAltToys") != std::string::npos) mReuseAltToys = value;

   return;
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, int value){
   //
   // set integer parameters
   //

   std::string s_name(name);

   if (s_name.find("NWorkers") != std::string::npos) mNWorkers = value;
   if (s_name.find("NToyToRebuild") != std::string::npos) mNToyToRebuild = value;
   if (s_name.find("PrintLevel") != std::string::npos) mPrintLevel = value;
   if (s_name.find("InitialFit") != std::string::npos) mInitialFit = value;
   if (s_name.find("RandomSeed") != std::string::npos) mRandomSeed = value;
   if (s_name.find("AsimovBins") != std::string::npos) mAsimovBins = value;

   return;
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, double value){
   //
   // set double precision parameters
   //

   std::string s_name(name);

   if (s_name.find("NToysRatio") != std::string::npos) mNToysRatio = value;
   if (s_name.find("MaxPOI") != std::string::npos) mMaxPoi = value;

   return;
}



void
RooStats::HypoTestInvTool::SetParameter(const char * name, const char * value){
   //
   // set string parameters
   //

   std::string s_name(name);

   if (s_name.find("MassValue") != std::string::npos) mMassValue.assign(value);
   if (s_name.find("MinimizerType") != std::string::npos) mMinimizerType.assign(value);
   if (s_name.find("ResultFileName") != std::string::npos) mResultFileName = value;

   return;
}



Double_t 
StandardHypoTestInvDemo(int enne, TString esse, const char * infile = 0,
                        const char * wsName = "combined",
                        const char * modelSBName = "ModelConfig",
                        const char * modelBName = "",
                        const char * dataName = "obsData",                 
                        int calculatorType = 0,
                        int testStatType = 0, 
                        bool useCLs = true ,  
                        int npoints = 6,   
                        double poimin = 0,  
                        double poimax = 5, 
                        int ntoys=1000,
                        bool useNumberCounting = false,
                        const char * nuisPriorName = 0){

			RooWorkspace::autoImportClassCode(kTRUE);
/*

  Other Parameter to pass in tutorial
  apart from standard for filename, ws, modelconfig and data

  type = 0 Freq calculator 
  type = 1 Hybrid calculator
  type = 2 Asymptotic calculator  
  type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)

  testStatType = 0 LEP
  = 1 Tevatron 
  = 2 Profile Likelihood
  = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
  = 4 Profiel Likelihood signed ( pll = -pll if mu < mu_hat) 
  = 5 Max Likelihood Estimate as test statistic
  = 6 Number of observed event as test statistic

  useCLs          scan for CLs (otherwise for CLs+b)    

  npoints:        number of points to scan , for autoscan set npoints = -1 

  poimin,poimax:  min/max value to scan in case of fixed scans 
  (if min >  max, try to find automatically)                           

  ntoys:         number of toys to use 

  useNumberCounting:  set to true when using number counting events 

  nuisPriorName:   name of prior for the nnuisance. This is often expressed as constraint term in the global model
  It is needed only when using the HybridCalculator (type=1)
  If not given by default the prior pdf from ModelConfig is used. 

  extra options are available as global paramwters of the macro. They major ones are: 

  plotHypoTestResult   plot result of tests at each point (TS distributions) (defauly is true)
  useProof             use Proof   (default is true) 
  writeResult          write result of scan (default is true)
  rebuild              rebuild scan for expected limits (require extra toys) (default is false)
  generateBinned       generate binned data sets for toys (default is false) - be careful not to activate with 
  a too large (>=3) number of observables 
  nToyRatio            ratio of S+B/B toys (default is 2)
    

*/

  
  
   TString fileName(infile);
   if (fileName.IsNull()) { 
      fileName = "results/example_combined_GaussExample_model.root";
      std::cout << "Use standard file generated with HistFactory : " << fileName << std::endl;
   }
  
   // open file and check if input file exists
   TFile * file = TFile::Open(fileName); 
  
   // if input file was specified but not found, quit
   if(!file && !TString(infile).IsNull()){
      cout <<"file " << fileName << " not found" << endl;
      return -999;
   } 
  
   // if default file not found, try to create it
   if(!file ){
      // Normally this would be run on the command line
      cout <<"will run standard hist2workspace example"<<endl;
      gROOT->ProcessLine(".! prepareHistFactory .");
      gROOT->ProcessLine(".! hist2workspace config/example.xml");
      cout <<"\n\n---------------------"<<endl;
      cout <<"Done creating example input"<<endl;
      cout <<"---------------------\n\n"<<endl;
    
      // now try to access the file again
      file = TFile::Open(fileName);
    
   }
  
   if(!file){
      // if it is still not there, then we can't continue
      cout << "Not able to run hist2workspace to create example input" <<endl;
      return -999;
   }
  


   HypoTestInvTool calc;

   // set parameters
   calc.SetParameter("PlotHypoTestResult", plotHypoTestResult);
   calc.SetParameter("WriteResult", writeResult);
   calc.SetParameter("Optimize", optimize);
   calc.SetParameter("UseVectorStore", useVectorStore);
   calc.SetParameter("GenerateBinned", generateBinned);
   calc.SetParameter("NToysRatio", nToysRatio);
   calc.SetParameter("MaxPOI", maxPOI);
   calc.SetParameter("UseProof", useProof);
   calc.SetParameter("NWorkers", nworkers);
   calc.SetParameter("Rebuild", rebuild);
   calc.SetParameter("ReuseAltToys", reuseAltToys);
   calc.SetParameter("NToyToRebuild", nToyToRebuild);
   calc.SetParameter("MassValue", massValue.c_str());
   calc.SetParameter("MinimizerType", minimizerType.c_str());
   calc.SetParameter("PrintLevel", printLevel);
   calc.SetParameter("InitialFit",initialFit);
   calc.SetParameter("ResultFileName",resultFileName);
   calc.SetParameter("RandomSeed",randomSeed);
   calc.SetParameter("AsimovBins",nAsimovBins);


   RooWorkspace * w = dynamic_cast<RooWorkspace*>( file->Get(wsName) );
   HypoTestInverterResult * r = 0;  
   std::cout << w << "\t" << fileName << std::endl;
   if (w != NULL) {
      r = calc.RunInverter(enne, esse, w, modelSBName, modelBName,
                           dataName, calculatorType, testStatType, useCLs,
                           npoints, poimin, poimax,  
                           ntoys, useNumberCounting, nuisPriorName );    
      if (!r) { 
         std::cerr << "Error running the HypoTestInverter - Exit " << std::endl;
         return -999;          
      }
   }
   else { 
      // case workspace is not present look for the inverter result
      std::cout << "Reading an HypoTestInverterResult with name " << wsName << " from file " << fileName << std::endl;
      r = dynamic_cast<HypoTestInverterResult*>( file->Get(wsName) ); //
      if (!r) { 
         std::cerr << "File " << fileName << " does not contain a workspace or an HypoTestInverterResult - Exit " 
                   << std::endl;
         file->ls();
         return -999; 
      }
   }		
  
   Double_t expLimitSG = calc.AnalyzeResult( r, calculatorType, testStatType, useCLs, npoints, infile );
  
   return expLimitSG;
}



Double_t 
RooStats::HypoTestInvTool::AnalyzeResult( HypoTestInverterResult * r,
                                          int calculatorType,
                                          int testStatType, 
                                          bool useCLs,  
                                          int npoints,
                                          const char * fileNameBase ){

   // analyze result produced by the inverter, optionally save it in a file 
   
  
   double lowerLimit = 0;
   double llError = 0;
#if defined ROOT_SVN_VERSION &&  ROOT_SVN_VERSION >= 44126
   if (r->IsTwoSided()) {
      lowerLimit = r->LowerLimit();
      llError = r->LowerLimitEstimatedError();
   }
#else
   lowerLimit = r->LowerLimit();
   llError = r->LowerLimitEstimatedError();
#endif

   double upperLimit = r->UpperLimit();
   double ulError = r->UpperLimitEstimatedError();

   //std::cout << "DEBUG : [ " << lowerLimit << " , " << upperLimit << "  ] " << std::endl;
      
   if (lowerLimit < upperLimit*(1.- 1.E-4) && lowerLimit != 0) 
      std::cout << "The computed lower limit is: " << lowerLimit << " +/- " << llError << std::endl;
   std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
  
   // compute expected limit
   std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
   std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
   std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
   std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
   std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
   std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;

   //formatted ouput for cut&paste into code
   //std::cout << std::endl;
   //std::cout << "Double_t muobs   = " << upperLimit << ";" << std::endl;
   //std::cout << "Double_t muexp   = " << r->GetExpectedUpperLimit(0) << ";" << std::endl;
   //std::cout << "Double_t muexpm2 = " << r->GetExpectedUpperLimit(-2) << ";" << std::endl;
   //std::cout << "Double_t muexpm1 = " << r->GetExpectedUpperLimit(-1) << ";" << std::endl;
   //std::cout << "Double_t muexpp1 = " << r->GetExpectedUpperLimit(1) << ";" << std::endl;
   //std::cout << "Double_t muexpp2 = " << r->GetExpectedUpperLimit(2) << ";" << std::endl;
   //std::cout << std::endl;

   Double_t SGexpected = r->GetExpectedUpperLimit(0);
  
   // write result in a file 
   if (r != NULL && mWriteResult) {
    
      // write to a file the results
      const char *  calcType = (calculatorType == 0) ? "Freq" : (calculatorType == 1) ? "Hybr" : "Asym";
      const char *  limitType = (useCLs) ? "CLs" : "Cls+b";
      const char * scanType = (npoints < 0) ? "auto" : "grid";
      if (mResultFileName.IsNull()) {
         mResultFileName = TString::Format("%s_%s_%s_ts%d_",calcType,limitType,scanType,testStatType);      
         //strip the / from the filename
         if (mMassValue.size()>0) {
            mResultFileName += mMassValue.c_str();
            mResultFileName += "_";
         }
    
         TString name = fileNameBase; 
         name.Replace(0, name.Last('/')+1, "");
         mResultFileName += name;
      }

      TFile * fileOut = new TFile(mResultFileName,"RECREATE");
      r->Write();
      fileOut->Close();                                                                     
   }   
  
  
   // plot the result ( p values vs scan points) 
   std::string typeName = "";
   if (calculatorType == 0 )
      typeName = "Frequentist";
   if (calculatorType == 1 )
      typeName = "Hybrid";   
   else if (calculatorType == 2 || calculatorType == 3) { 
      typeName = "Asymptotic";
      mPlotHypoTestResult = false; 
   }
  
   const char * resultName = r->GetName();
   TString plotTitle = TString::Format("%s CL Scan for workspace %s",typeName.c_str(),resultName);
   HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);

   // plot in a new canvas with style
   TString c1Name = TString::Format("%s_Scan",typeName.c_str());
   TCanvas * c1 = new TCanvas(c1Name); 
   c1->SetLogy(false);

   //plot->Draw("CLb 2CL");  // plot all and Clb
   //plot->Draw("EXP");  // plot all and Clb
   plot->Draw("");  // plot all and Clb
   c1->SaveAs("brasilianFlag.pdf");
   c1->SaveAs("brasilianFlag.C");

   // if (useCLs) 
   //    plot->Draw("CLb 2CL");  // plot all and Clb
   // else 
   //    plot->Draw("");  // plot all and Clb
  
   const int nEntries = r->ArraySize();
  
   // plot test statistics distributions for the two hypothesis 
   if (mPlotHypoTestResult) { 
      TCanvas * c2 = new TCanvas();
      if (nEntries > 1) { 
         int ny = TMath::CeilNint(TMath::Sqrt(nEntries));
         int nx = TMath::CeilNint(double(nEntries)/ny);
         c2->Divide( nx,ny);
      }
      for (int i=0; i<nEntries; i++) {
         if (nEntries > 1) c2->cd(i+1);
         SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
         pl->SetLogYaxis(true);
         pl->Draw();
      }
   }

   return SGexpected;
}



// internal routine to run the inverter
HypoTestInverterResult *
RooStats::HypoTestInvTool::RunInverter(Int_t enne, TString esse, RooWorkspace * w,
                                       const char * modelSBName, const char * modelBName, 
                                       const char * dataName, int type,  int testStatType, 
                                       bool useCLs, int npoints, double poimin, double poimax, 
                                       int ntoys,
                                       bool useNumberCounting,
                                       const char * nuisPriorName ){

   std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;
  
   w->Print();
  
  
   RooAbsData * data = w->data(dataName); 
   if (!data) { 
      Error("StandardHypoTestDemo","Not existing data %s",dataName);
      return 0;
   }
   else 
      std::cout << "Using data set " << dataName << std::endl;
  
   if (mUseVectorStore) { 
      RooAbsData::setDefaultStorageType(RooAbsData::Vector);
      data->convertToVectorStore() ;
   }
  
  
   // get models from WS
   // get the modelConfig out of the file
   ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
   ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);
  
   if (!sbModel) {
      Error("StandardHypoTestDemo","Not existing ModelConfig %s",modelSBName);
      return 0;
   }
   // check the model 
   if (!sbModel->GetPdf()) { 
      Error("StandardHypoTestDemo","Model %s has no pdf ",modelSBName);
      return 0;
   }
   if (!sbModel->GetParametersOfInterest()) {
      Error("StandardHypoTestDemo","Model %s has no poi ",modelSBName);
      return 0;
   }
   if (!sbModel->GetObservables()) {
      Error("StandardHypoTestInvDemo","Model %s has no observables ",modelSBName);
      return 0;
   }
   if (!sbModel->GetSnapshot() ) { 
      Info("StandardHypoTestInvDemo","Model %s has no snapshot  - make one using model poi",modelSBName);
      sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
   }
  
   // case of no systematics
   // remove nuisance parameters from model
   if (noSystematics) { 
      const RooArgSet * nuisPar = sbModel->GetNuisanceParameters();
      if (nuisPar && nuisPar->getSize() > 0) { 
         std::cout << "StandardHypoTestInvDemo" << "  -  Switch off all systematics by setting them constant to their initial values" << std::endl;
         RooStats::SetAllConstant(*nuisPar);
      }
      if (bModel) { 
         const RooArgSet * bnuisPar = bModel->GetNuisanceParameters();
         if (bnuisPar) 
            RooStats::SetAllConstant(*bnuisPar);
      }
   }
  
   if (!bModel || bModel == sbModel) {
      Info("StandardHypoTestInvDemo","The background model %s does not exist",modelBName);
      Info("StandardHypoTestInvDemo","Copy it from ModelConfig %s and set POI to zero",modelSBName);
      bModel = (ModelConfig*) sbModel->Clone();
      bModel->SetName(TString(modelSBName)+TString("_with_poi_0"));      
      RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
      if (!var) return 0;
      double oldval = var->getVal();
      var->setVal(0);
      bModel->SetSnapshot( RooArgSet(*var)  );
      var->setVal(oldval);
   }
   else { 
      if (!bModel->GetSnapshot() ) { 
         Info("StandardHypoTestInvDemo","Model %s has no snapshot  - make one using model poi and 0 values ",modelBName);
         RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
         if (var) { 
            double oldval = var->getVal();
            var->setVal(0);
            bModel->SetSnapshot( RooArgSet(*var)  );
            var->setVal(oldval);
         }
         else { 
            Error("StandardHypoTestInvDemo","Model %s has no valid poi",modelBName);
            return 0;
         }         
      }
   }

   // check model  has global observables when there are nuisance pdf
   // for the hybrid case the globobs are not needed
   if (type != 1 ) { 
      bool hasNuisParam = (sbModel->GetNuisanceParameters() && sbModel->GetNuisanceParameters()->getSize() > 0);
      bool hasGlobalObs = (sbModel->GetGlobalObservables() && sbModel->GetGlobalObservables()->getSize() > 0);
      if (hasNuisParam && !hasGlobalObs ) {  
         // try to see if model has nuisance parameters first 
         RooAbsPdf * constrPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisanceConstraintPdf_sbmodel");
         if (constrPdf) { 
            Warning("StandardHypoTestInvDemo","Model %s has nuisance parameters but no global observables associated",sbModel->GetName());
            Warning("StandardHypoTestInvDemo","\tThe effect of the nuisance parameters will not be treated correctly ");
         }
      }
   }


  
   // run first a data fit 
  
   const RooArgSet * poiSet = sbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar*)poiSet->first();
  
   std::cout << "StandardHypoTestInvDemo : POI initial value:   " << poi->GetName() << " = " << poi->getVal()   << std::endl;  
  
   // fit the data first (need to use constraint )
   TStopwatch tw; 

   bool doFit = initialFit;
   if (testStatType == 0 && initialFit == -1) doFit = false;  // case of LEP test statistic
   if (type == 3  && initialFit == -1) doFit = false;         // case of Asymptoticcalculator with nominal Asimov
   double poihat = 0;

   if (minimizerType.size()==0) minimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
   else 
      ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minimizerType.c_str());
    
   Info("StandardHypoTestInvDemo","Using %s as minimizer for computing the test statistic",
        ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str() );
   
   if (doFit)  { 

      // do the fit : By doing a fit the POI snapshot (for S+B)  is set to the fit value
      // and the nuisance parameters nominal values will be set to the fit value. 
      // This is relevant when using LEP test statistics

      Info( "StandardHypoTestInvDemo"," Doing a first fit to the observed data ");
      RooArgSet constrainParams;
      if (sbModel->GetNuisanceParameters() ) constrainParams.add(*sbModel->GetNuisanceParameters());
      RooStats::RemoveConstantParameters(&constrainParams);
      tw.Start(); 
      RooFitResult * fitres = sbModel->GetPdf()->fitTo(*data,InitialHesse(false), Hesse(false),
                                                       Minimizer(minimizerType.c_str(),"Migrad"), Strategy(0), PrintLevel(mPrintLevel+1), Constrain(constrainParams), Save(true) );
      if (fitres->status() != 0) { 
         Warning("StandardHypoTestInvDemo","Fit to the model failed - try with strategy 1 and perform first an Hesse computation");
         fitres = sbModel->GetPdf()->fitTo(*data,InitialHesse(true), Hesse(false),Minimizer(minimizerType.c_str(),"Migrad"), Strategy(1), PrintLevel(mPrintLevel+1), Constrain(constrainParams), Save(true) );
      }
      if (fitres->status() != 0) 
         Warning("StandardHypoTestInvDemo"," Fit still failed - continue anyway.....");
  
  
      poihat  = poi->getVal();
      std::cout << "StandardHypoTestInvDemo - Best Fit value : " << poi->GetName() << " = "  
                << poihat << " +/- " << poi->getError() << std::endl;
      std::cout << "Time for fitting : "; tw.Print(); 
  
      //save best fit value in the poi snapshot 
      sbModel->SetSnapshot(*sbModel->GetParametersOfInterest());
      std::cout << "StandardHypoTestInvo: snapshot of S+B Model " << sbModel->GetName() 
                << " is set to the best fit value" << std::endl;


      TFile pippo("FitResults.root","RECREATE");
      fitres->Write();
      pippo.Close();
   }

   // print a message in case of LEP test statistics because it affects result by doing or not doing a fit 
   if (testStatType == 0) {
      if (!doFit) 
         Info("StandardHypoTestInvDemo","Using LEP test statistic - an initial fit is not done and the TS will use the nuisances at the model value");
      else 
         Info("StandardHypoTestInvDemo","Using LEP test statistic - an initial fit has been done and the TS will use the nuisances at the best fit value");
   }


   // build test statistics and hypotest calculators for running the inverter 
  
   SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(),*bModel->GetPdf());

   // null parameters must includes snapshot of poi plus the nuisance values 
   RooArgSet nullParams(*sbModel->GetSnapshot());
   if (sbModel->GetNuisanceParameters()) nullParams.add(*sbModel->GetNuisanceParameters());
   if (sbModel->GetSnapshot()) slrts.SetNullParameters(nullParams);
   RooArgSet altParams(*bModel->GetSnapshot());
   if (bModel->GetNuisanceParameters()) altParams.add(*bModel->GetNuisanceParameters());
   if (bModel->GetSnapshot()) slrts.SetAltParameters(altParams);
  
   // ratio of profile likelihood - need to pass snapshot for the alt
   RatioOfProfiledLikelihoodsTestStat 
      ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   if (testStatType == 11) ropl.SetSubtractMLE(true);
   ropl.SetPrintLevel(mPrintLevel);
   ropl.SetMinimizer(minimizerType.c_str());
  
   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   if (testStatType == 3) profll.SetOneSided(true);
   if (testStatType == 4) profll.SetSigned(true);
   profll.SetMinimizer(minimizerType.c_str());
   profll.SetPrintLevel(mPrintLevel);

   profll.SetReuseNLL(mOptimize);
   slrts.SetReuseNLL(mOptimize);
   ropl.SetReuseNLL(mOptimize);

   if (mOptimize) { 
      profll.SetStrategy(0);
      ropl.SetStrategy(0);
      ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   }
  
   if (mMaxPoi > 0) poi->setMax(mMaxPoi);  // increase limit
  
   MaxLikelihoodEstimateTestStat maxll(*sbModel->GetPdf(),*poi); 
   NumEventsTestStat nevtts;

   AsymptoticCalculator::SetPrintLevel(mPrintLevel);
  
   // create the HypoTest calculator class 
   HypoTestCalculatorGeneric *  hc = 0;
   if (type == 0) hc = new FrequentistCalculator(*data, *bModel, *sbModel);
   else if (type == 1) hc = new HybridCalculator(*data, *bModel, *sbModel);
   // else if (type == 2 ) hc = new AsymptoticCalculator(*data, *bModel, *sbModel, false, mAsimovBins);
   // else if (type == 3 ) hc = new AsymptoticCalculator(*data, *bModel, *sbModel, true, mAsimovBins);  // for using Asimov data generated with nominal values 
   else if (type == 2 ) hc = new AsymptoticCalculator(*data, *bModel, *sbModel, false );
   else if (type == 3 ) hc = new AsymptoticCalculator(*data, *bModel, *sbModel, true );  // for using Asimov data generated with nominal values 
   else {
      Error("StandardHypoTestInvDemo","Invalid - calculator type = %d supported values are only :\n\t\t\t 0 (Frequentist) , 1 (Hybrid) , 2 (Asymptotic) ",type);
      return 0;
   }
  
   // set the test statistic 
   TestStatistic * testStat = 0;
   if (testStatType == 0) testStat = &slrts;
   if (testStatType == 1 || testStatType == 11) testStat = &ropl;
   if (testStatType == 2 || testStatType == 3 || testStatType == 4) testStat = &profll;
   if (testStatType == 5) testStat = &maxll;
   if (testStatType == 6) testStat = &nevtts;

   if (testStat == 0) { 
      Error("StandardHypoTestInvDemo","Invalid - test statistic type = %d supported values are only :\n\t\t\t 0 (SLR) , 1 (Tevatron) , 2 (PLR), 3 (PLR1), 4(MLE)",testStatType);
      return 0;
   }
  
  
   ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
   if (toymcs && (type == 0 || type == 1) ) { 
      // look if pdf is number counting or extended
      if (sbModel->GetPdf()->canBeExtended() ) { 
         if (useNumberCounting)   Warning("StandardHypoTestInvDemo","Pdf is extended: but number counting flag is set: ignore it ");
      }
      else { 
         // for not extended pdf
         if (!useNumberCounting  )  { 
            int nEvents = data->numEntries();
            Info("StandardHypoTestInvDemo","Pdf is not extended: number of events to generate taken  from observed data set is %d",nEvents);
            toymcs->SetNEventsPerToy(nEvents);
         }
         else {
            Info("StandardHypoTestInvDemo","using a number counting pdf");
            toymcs->SetNEventsPerToy(1);
         }
      }

      toymcs->SetTestStatistic(testStat);
    
      if (data->isWeighted() && !mGenerateBinned) { 
         Info("StandardHypoTestInvDemo","Data set is weighted, nentries = %d and sum of weights = %8.1f but toy generation is unbinned - it would be faster to set mGenerateBinned to true\n",data->numEntries(), data->sumEntries());
      }
      toymcs->SetGenerateBinned(mGenerateBinned);
  
      toymcs->SetUseMultiGen(mOptimize);
    
      if (mGenerateBinned &&  sbModel->GetObservables()->getSize() > 2) { 
         Warning("StandardHypoTestInvDemo","generate binned is activated but the number of ovservable is %d. Too much memory could be needed for allocating all the bins",sbModel->GetObservables()->getSize() );
      }

      // set the random seed if needed
      if (mRandomSeed >= 0) RooRandom::randomGenerator()->SetSeed(mRandomSeed); 
    
   }
  
   // specify if need to re-use same toys
   if (reuseAltToys) {
      hc->UseSameAltToys();
   }
  
   if (type == 1) { 
      HybridCalculator *hhc = dynamic_cast<HybridCalculator*> (hc);
      assert(hhc);
    
      hhc->SetToys(ntoys, (int) (ntoys/mNToysRatio)); // can use less ntoys for b hypothesis 
    
      // remove global observables from ModelConfig (this is probably not needed anymore in 5.32)
      bModel->SetGlobalObservables(RooArgSet() );
      sbModel->SetGlobalObservables(RooArgSet() );
    
    
      // check for nuisance prior pdf in case of nuisance parameters 
      if (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters() ) {

         // fix for using multigen (does not work in this case)
         toymcs->SetUseMultiGen(false);
         ToyMCSampler::SetAlwaysUseMultiGen(false);

         RooAbsPdf * nuisPdf = 0; 
         if (nuisPriorName) nuisPdf = w->pdf(nuisPriorName);
         // use prior defined first in bModel (then in SbModel)
         if (!nuisPdf)  { 
            Info("StandardHypoTestInvDemo","No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the model");
            if (bModel->GetPdf() && bModel->GetObservables() ) 
               nuisPdf = RooStats::MakeNuisancePdf(*bModel,"nuisancePdf_bmodel");
            else 
               nuisPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisancePdf_sbmodel");
         }   
         if (!nuisPdf ) {
            if (bModel->GetPriorPdf())  { 
               nuisPdf = bModel->GetPriorPdf();
               Info("StandardHypoTestInvDemo","No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf->GetName());            
            }
            else { 
               Error("StandardHypoTestInvDemo","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived");
               return 0;
            }
         }
         assert(nuisPdf);
         Info("StandardHypoTestInvDemo","Using as nuisance Pdf ... " );
         nuisPdf->Print();
      
         const RooArgSet * nuisParams = (bModel->GetNuisanceParameters() ) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
         RooArgSet * np = nuisPdf->getObservables(*nuisParams);
         if (np->getSize() == 0) { 
            Warning("StandardHypoTestInvDemo","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
         }
         delete np;
      
         hhc->ForcePriorNuisanceAlt(*nuisPdf);
         hhc->ForcePriorNuisanceNull(*nuisPdf);
      
      
      }
   } 
   else if (type == 2 || type == 3) { 
      if (testStatType == 3) ((AsymptoticCalculator*) hc)->SetOneSided(true);  
      if (testStatType != 2 && testStatType != 3)  
         Warning("StandardHypoTestInvDemo","Only the PL test statistic can be used with AsymptoticCalculator - use by default a two-sided PL");
   }
   else if (type == 0 || type == 1) 
      ((FrequentistCalculator*) hc)->SetToys(ntoys, (int) (ntoys/mNToysRatio)); 

  
   // Get the result
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  
  
  
   HypoTestInverter calc(*hc);
   calc.SetConfidenceLevel(0.95);
   //calc.SetConfidenceLevel(0.90);
  
  
   calc.UseCLs(useCLs);
   calc.SetVerbose(true);
  
   // can speed up using proof-lite
   if (mUseProof && mNWorkers > 1) { 
    //ProofConfig pc(*w, mNWorkers, "ippolitv@atlas-ui-03.roma1.infn.it:21001", kFALSE);
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << "ENNE " << enne << " ESSE " << esse << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
      ProofConfig pc(*w, enne, esse, kFALSE);
    //ProofConfig pc(*w, mNWorkers, "", kFALSE); // FIXME
      toymcs->SetProofConfig(&pc);    // enable proof
   }
  
  
   if (npoints > 0) {
      if (poimin > poimax) { 
         // if no min/max given scan between MLE and +4 sigma 
         poimin = int(poihat);
         poimax = int(poihat +  4 * poi->getError());
      }
      std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
      calc.SetFixedScan(npoints,poimin,poimax);
   }
   else { 
      //poi->setMax(10*int( (poihat+ 10 *poi->getError() )/10 ) );
      std::cout << "Doing an  automatic scan  in interval : " << poi->getMin() << " , " << poi->getMax() << std::endl;
   }
  
   tw.Start();
   HypoTestInverterResult * r = calc.GetInterval();
   std::cout << "Time to perform limit scan \n";
   tw.Print();
  
   if (mRebuild) {
      calc.SetCloseProof(1);
      tw.Start();
      SamplingDistribution * limDist = calc.GetUpperLimitDistribution(true,mNToyToRebuild);
      std::cout << "Time to rebuild distributions " << std::endl;
      tw.Print();
    
      if (limDist) { 
         std::cout << "expected up limit " << limDist->InverseCDF(0.5) << " +/- " 
                   << limDist->InverseCDF(0.16) << "  " 
                   << limDist->InverseCDF(0.84) << "\n"; 
      
         //update r to a new updated result object containing the rebuilt expected p-values distributions
         // (it will not recompute the expected limit)
         if (r) delete r;  // need to delete previous object since GetInterval will return a cloned copy
         r = calc.GetInterval();
      
      }
      else 
         std::cout << "ERROR : failed to re-build distributions " << std::endl; 
   }
  
   return r;
}



void ReadResult(const char * fileName, const char * resultName="", bool useCLs=true) { 
   // read a previous stored result from a file given the result name

   Int_t enne = -1;
   TString esse = "";
   StandardHypoTestInvDemo(enne, esse, fileName, resultName,"","","",0,0,useCLs);
}


#ifdef USE_AS_MAIN
int main() {
    StandardHypoTestInvDemo();
}
#endif


/* Simultanous ABCD code  (S.giagu) */
/* 
 * n[4] = {n_A, n_B, n_C, n_D} <-- number of observed events in regions A, B, C, D 
 * s[4] = {s_A, s_B, s_C, s_D} <-- number of signal events in regions A, B, C, D 
 * b[4] = {b_A, b_B, b_C, b_D} <-- number of BG events (estimated from MC) in regions A, B, C, D 
 * c[4] = {c_A, c_B, c_C, c_D} <-- number of other BG events like cosmics etc.. (indipendently estimated from data) in regions A, B, C, D 
 * useB = kFALSE <-- don't use BG events estimated from MC;  kTRUE <-- use them 
 * useC = kFALSE <-- don't use other BG events (like cosmics etc..) indipendently estimated from data;  kTRUE <-- use them
 * blindA: kTRUE <-- keep signal region blind (i.e. test done assuming n_A = ABCD_exp_A),  kFALSE <-- use obs events in signal region
 *
 */
Double_t simultaneousABCD(Double_t n[4], Double_t s[4], Double_t b[4], Double_t c[4],
                          TString out_filename = "ABCD_ws.root", Bool_t useB = kFALSE, Bool_t useC = kFALSE, Bool_t blindA = kTRUE)
{

    // set RooFit random seed to a fix value for reproducible results
    RooRandom::randomGenerator()->SetSeed(4357);

    // init
    RooWorkspace::autoImportClassCode(kTRUE); // set default behaviour of RooWorkspace when importing new classes

    //Inputs
    // signal
    Double_t ns_A = s[0];
    if (ns_A <= 0) {
      std::cout << "ERROR: 0 signal events in signal region (A)!!! --> Check inputs!  s[0] = " << ns_A << std::endl;
      return -1.0;
    }
    Double_t ns_B = s[1]/ns_A;
    Double_t ns_C = s[2]/ns_A;
    Double_t ns_D = s[3]/ns_A;
    // data
    Double_t nd_A = n[0];
    Double_t nd_B = n[1];
    Double_t nd_C = n[2];
    Double_t nd_D = n[3];
    // MC based BG
    Double_t nb_A = b[0];
    Double_t nb_B = b[1];
    Double_t nb_C = b[2];
    Double_t nb_D = b[3];
    // Independently DATA based BG
    Double_t nc_A = c[0];
    Double_t nc_B = c[1];
    Double_t nc_C = c[2];
    Double_t nc_D = c[3];

    // Some initial printout ...
    Double_t nd_A_expected = 0.0;
    if (nd_C > 0) nd_A_expected = nd_B * nd_D / nd_C;
    std::cout << "Obs events in signal region (A) estimated from control regions using PLAIN ABCD: " << nd_A_expected << std::endl;
    std::cout << "              " << std::endl;

    if (blindA) { //don't use observed data in signal region (for expected yields) but expectation from PLAIN ABCD
      nd_A = nd_A_expected;
    }

    std::cout << "Input yields with signal region " << (blindA?"blinded: ":"unblinded: ") << std::endl;
    std::cout << "Observed A/B/C/D: " << nd_A << " / " << nd_B << " / " << nd_C << " / " << nd_D << std::endl;
    std::cout << "Signal   A/B/C/D: " << ns_A << " / " << ns_B << " / " << ns_C << " / " << ns_D << std::endl;
    if (useB)  
       std::cout << "MC BG A/B/C/D:    " << nb_A << " / " << nb_B << " / " << nb_C << " / " << nb_D << std::endl;
    if (useC)  
      std::cout << "Data BG A/B/C/D:  " << nc_A << " / " << nc_B << " / " << nc_C << " / " << nc_D << std::endl;
    std::cout << "              " << std::endl;

    // guess some initial values (to speedup fit convergence) for the parameters ...
    Double_t sr_B = ns_B / ns_A;
    Double_t sr_C = ns_C / ns_A;
    Double_t sr_D = ns_D / ns_A;
    Double_t mu_guess = 3.0 / ns_A; // starting guess for mu = N_sig_UL / N_sig_exp assume no signal observed (~3 events UL at 95% CL)

    std::cout << "Signal ratios B/A C/A and D/A: " << sr_B << " / " << sr_C << " / " << sr_D << std::endl;
    std::cout << "Guess mu: " << mu_guess << std::endl;
    std::cout << "              " << std::endl;

    Double_t nq_A_guess = nd_A_expected; // starting guess from ABCD events in signal region (from ABCD ansatz + other BG)
    if (useB) nq_A_guess -= nb_A;
    if (useC) nq_A_guess -= nc_A;
    if (nq_A_guess < 0) nq_A_guess = 0.0;
    Double_t nq_B_guess = nd_B;
    Double_t nq_D_guess = nd_D;
    std::cout << "Guess Multijet BG A/B/D: " << nq_A_guess << " / " << nq_B_guess << " / " << nq_D_guess << std::endl;
    std::cout << "              " << std::endl;


    if (nq_A_guess <= 0) {
        nq_A_guess = 3.0;
        std::cout << "WARNING: zero events for nq_A_guess, used 3 events" << std::endl;
    }
    Double_t ta_A = nq_A_guess; // tau parameters guesses (see Likelihood ABCD note in stat forum for definition)

    Double_t ta_B = 3.0;
    Double_t ta_D = 3.0;
    if (nq_B_guess > 0)
        ta_B = nq_B_guess / nq_A_guess;
    else
        std::cout << "WARNING: zero events for nq_B_guess, used 3 events" << std::endl;
    if (nq_D_guess > 0)
        ta_D = nq_D_guess / nq_A_guess;
    else
        std::cout << "WARNING: zero events for nq_D_guess, used 3 events" << std::endl;

    std::cout << "tau Multijet A/B/D: " << ta_A << " / " << ta_B << " / " << ta_D << std::endl;
    std::cout << "              " << std::endl;


    // make model
    RooWorkspace* wspace = new RooWorkspace("wspace", "ABCD workspace");
    wspace->addClassDeclImportDir("."); // add code import paths
    wspace->addClassImplImportDir(".");

    // observed events
    wspace->factory("NA[0,1000]");
    wspace->factory("NB[0,1000]");
    wspace->factory("NC[0,1000]");
    wspace->factory("ND[0,1000]");

    // POI
    wspace->factory(TString::Format("mu[%f,0,1000]", mu_guess));  // mu = NsA/Ns0 (Ns0 = expected events)
    //note: SM means mu=0 (used for the BG only hypotesis for the expected limit

    // pdf parameters
    wspace->factory(TString::Format("lumi[%f]",1.0));    // Luminosity (scale factor wrt the luinosity on dat)

    wspace->factory(TString::Format("Ns0[%f,0,1000]",  ns_A)); //Expected signal in region A
    wspace->factory(TString::Format("effB[%f,0,1000]", sr_B)); //Sig. eff. in region B wrt region A
    wspace->factory(TString::Format("effC[%f,0,1000]", sr_C)); //Sig. eff. in region C wrt region A
    wspace->factory(TString::Format("effD[%f,0,1000]", sr_D)); //Sig. eff. in region D wrt region A

    wspace->factory(TString::Format("Nq[%f,0,1000]",   ta_A)); //number of Multijet events in region A
    wspace->factory(TString::Format("tauB[%f,0,1000]", ta_B)); //Multijet BG eff. in region B wrt region A
    wspace->factory(TString::Format("tauD[%f,0,1000]", ta_D)); //Multijet BG eff. in region D wrt region A

    if (useB) {
       wspace->factory(TString::Format("NbA[%f,0,1000]", nb_A)); //number of MC BG events in region A
       wspace->factory(TString::Format("NbB[%f,0,1000]", nb_B)); // in region B
       wspace->factory(TString::Format("NbC[%f,0,1000]", nb_C)); // C
       wspace->factory(TString::Format("NbD[%f,0,1000]", nb_D)); // ... and D
    }

    if (useC) {
       wspace->factory(TString::Format("NcA[%f,0,1000]", nc_A)); //number of Other data-driven BG events in region A
       wspace->factory(TString::Format("NcB[%f,0,1000]", nc_B)); // in region B
       wspace->factory(TString::Format("NcC[%f,0,1000]", nc_C)); // C
       wspace->factory(TString::Format("NcD[%f,0,1000]", nc_D)); // ... and D
    }

    //Systematic uncertanties' nuisance parameters 
    //lumi
    wspace->factory("alpha_lumi[1, 0, 10]");
    wspace->factory("nom_lumi[1, 0, 10]");
    wspace->factory("nom_sigma_lumi[0.021]");  // <--  2.1% final run2 2015
    wspace->factory("Gaussian::constraint_lumi(nom_lumi, alpha_lumi, nom_sigma_lumi)");

    wspace->factory("alpha_S[1, 0, 2]");  //systematic nuisance on signal (efficiencies etc.) and on MC bg
    wspace->factory("nom_S[1, 0, 10]");
    // trigger eff: muNS 6%  3mu6: 5.8%  caloratio 11% --> weighted with efficiencies in signal MC: 7%
    // reco eff: 20%
    // pileup on SumPt: 5%
    // pt resolution: 10%
    // JES: 1%
    // total: 24%
    wspace->factory("nom_sigma_S[0.24]"); // 24% totale displaced LJ analysis 2016
    wspace->factory("Gaussian::constraint_S(nom_S, alpha_S, nom_sigma_S)");

    wspace->factory("alpha_Q[1, 0, 2]");  //systematic nuisance on Multijet   
    wspace->factory("nom_Q[1, 0, 5]");
    wspace->factory("nom_sigma_Q[0.30]");   //30%% on QCD from ABCD variations and closure tests
    wspace->factory("Gaussian::constraint_Q(nom_Q, alpha_Q, nom_sigma_Q)");

    if (useC) {
       wspace->factory("alpha_C[1, 0, 10]");  //systematic nuisance on other data-driven BG
       wspace->factory("nom_C[1, 0, 10]");
       wspace->factory("nom_sigma_C[0.32]");   //20% on cosmic BG (just as an example)
       wspace->factory("Gaussian::constraint_C(nom_C, alpha_C, nom_sigma_C)");
    }

    if (useB) {
       wspace->factory("alpha_B[1, 0, 10]");  //systematic nuisance on MC BG
       wspace->factory("nom_B[1, 0, 10]");
       wspace->factory("nom_sigma_B[0.32]");   //10% on MC BG (just as an example)
       wspace->factory("Gaussian::constraint_B(nom_B, alpha_B, nom_sigma_B)");
    }

    // PDF
    wspace->factory("prod::NsA(mu,Ns0,lumi,alpha_lumi,alpha_S)");      // expected signal events in A: mu*Ns0*L
    wspace->factory("prod::NsB(mu,Ns0,effB,lumi,alpha_lumi,alpha_S)"); // expected signal events in B: mu*Ns0*L*effB
    wspace->factory("prod::NsC(mu,Ns0,effC,lumi,alpha_lumi,alpha_S)");
    wspace->factory("prod::NsD(mu,Ns0,effD,lumi,alpha_lumi,alpha_S)");

    wspace->factory("prod::NbQA(Nq,lumi, alpha_Q)");    // expected Multijet (i.e. ABCD) BG events in A: Nq*SysQ
    wspace->factory("prod::NbQB(Nq,lumi, tauB)");       // expected Multijet BG events in B: Nq*tauB;
    wspace->factory("prod::NbQC(Nq,lumi, tauB,tauD)");  // ABCD ansatz: expected Multijet BG events in C: Nq*tauB*tauD;
    wspace->factory("prod::NbQD(Nq,lumi, tauD)");       // expected Multijet BG events in D: Nq*tauD;

    if (useB) {
       wspace->factory("prod::NbBA(NbA,lumi,alpha_lumi,alpha_B)");   // expected MC BG events in region A: NbA*sysB*sysLumi (MC)
       wspace->factory("prod::NbBB(NbB,lumi,alpha_lumi,alpha_B)");
       wspace->factory("prod::NbBC(NbC,lumi,alpha_lumi,alpha_B)");
       wspace->factory("prod::NbBD(NbD,lumi,alpha_lumi,alpha_B)");
    }

    if (useC) {
       wspace->factory("prod::NbCA(NcA,lumi, alpha_C)");   // expected cosm. BG events in A: NcA*sysCos 
       wspace->factory("prod::NbCB(NcB,lumi, alpha_C)");
       wspace->factory("prod::NbCC(NcC,lumi, alpha_C)");
       wspace->factory("prod::NbCD(NcD,lumi, alpha_C)");
    }

    if (useC && useB) {
       wspace->factory("sum::NexpA(NsA,NbQA,NbBA,NbCA)");
       wspace->factory("sum::NexpB(NsB,NbQB,NbBB,NbCB)");
       wspace->factory("sum::NexpC(NsC,NbQC,NbBC,NbCC)");
       wspace->factory("sum::NexpD(NsD,NbQD,NbBD,NbCD)");
    } else if (useC && !useB) {
       wspace->factory("sum::NexpA(NsA,NbQA,NbCA)");
       wspace->factory("sum::NexpB(NsB,NbQB,NbCB)");
       wspace->factory("sum::NexpC(NsC,NbQC,NbCC)");
       wspace->factory("sum::NexpD(NsD,NbQD,NbCD)");
    } else if (!useC && useB) {
       wspace->factory("sum::NexpA(NsA,NbQA,NbBA)");
       wspace->factory("sum::NexpB(NsB,NbQB,NbBB)");
       wspace->factory("sum::NexpC(NsC,NbQC,NbBC)");
       wspace->factory("sum::NexpD(NsD,NbQD,NbBD)");
    } else {
       wspace->factory("sum::NexpA(NsA,NbQA)");
       wspace->factory("sum::NexpB(NsB,NbQB)");
       wspace->factory("sum::NexpC(NsC,NbQC)");
       wspace->factory("sum::NexpD(NsD,NbQD)");
    }

    wspace->factory("Poisson::obsA(NA,NexpA)");
    wspace->factory("Poisson::obsB(NB,NexpB)");
    wspace->factory("Poisson::obsC(NC,NexpC)");
    wspace->factory("Poisson::obsD(ND,NexpD)");


    if (useC && useB) {
      wspace->factory("PROD::model(obsA,obsB,obsC,obsD,constraint_lumi,constraint_Q,constraint_S, constraint_C, constraint_B)");
    } else if (useC && !useB) {
      wspace->factory("PROD::model(obsA,obsB,obsC,obsD,constraint_lumi,constraint_Q,constraint_S, constraint_C)");
    } else if (!useC && useB) {
      wspace->factory("PROD::model(obsA,obsB,obsC,obsD,constraint_lumi,constraint_Q,constraint_S, constraint_B)");
    } else {
      wspace->factory("PROD::model(obsA,obsB,obsC,obsD,constraint_lumi,constraint_Q,constraint_S)");
    }

    // sets
    TString the_poi  = "mu";
    TString the_nuis = ",Nq,tauB,tauD,alpha_lumi,alpha_Q,alpha_S";
    TString the_glob = ",nom_lumi,nom_Q,nom_S";
    if (useC) {
        the_nuis += ",alpha_C";
        the_glob += ",nom_C";
    }
    if (useB) {
        the_nuis += ",alpha_B";
        the_glob += ",nom_B";
    }
    std::cout << "poi:        " << the_poi << std::endl;
    std::cout << "nuisances:  " << the_nuis << std::endl;
    std::cout << "global var: " << the_glob << std::endl;

    wspace->defineSet("obs", "NA,NB,NC,ND");
    wspace->defineSet("poi",  the_poi);
    wspace->defineSet("nuis", the_nuis);
    wspace->defineSet("glob", the_glob);

    //fix needed parameters

    TString    interesting = "," + the_poi + "," + the_nuis; // note leading comma

    TIterator *itr;
    itr = wspace->pdf("model")->getParameters(*wspace->set("obs"))->createIterator();
    TObject *obj(0);
    RooRealVar *rrv(0);
    RooCategory *rc(0);
    while ((obj = itr->Next())) {
        if (interesting.Contains(TString::Format(",%s", obj->GetName())) == kFALSE) {
            // try RooRealVar
            rrv = dynamic_cast<RooRealVar*>(obj);
            rc = dynamic_cast<RooCategory*>(obj);

            if (rrv) {
                rrv -> setConstant(kTRUE);
                std::cout << "setting " << rrv->GetName() << " const(val = " << rrv->getVal() << ")" << std::endl;
            } else if (rc) {
                rc -> setConstant(kTRUE);
                std::cout << "setting " << rc->GetName() << " const(val = " << rc->getLabel() << ")" << std::endl;
            } else
                std::cout << "unable to set " << obj->GetName() << " const(not RRV)" << std::endl;
        }
    } // loop on pdf parameters

    // inspect workspace
    wspace->Print();

    ////////////////////////////////////////////////////////////
    // Generate toy data
    // generate toy data assuming current value of the parameters
    // add Verbose() to see how it's being generated
    //RooDataSet* data = wspace->pdf("model")->generate(*wspace->set("obs"),1);
    //data->Print("v");
    // or input real data
    // add Verbose() to see how it's being generated
    RooDataSet* data = new RooDataSet("data", "obsData", *wspace->set("obs"));
    wspace->var("NA")->setVal(nd_A);
    wspace->var("NB")->setVal(nd_B);
    wspace->var("NC")->setVal(nd_C);
    wspace->var("ND")->setVal(nd_D);
    data->add(*wspace->set("obs"));
    data->Print("v");

    // import into workspace.
    wspace->import(*data, RooFit::Rename("obsData"));

    /////////////////////////////////////////////////////
    // Now the statistical tests
    // model config
    ModelConfig* mc = new ModelConfig("mc");
    mc->SetWorkspace(*wspace);
    mc->SetPdf(*wspace->pdf("model"));
    mc->SetObservables(*wspace->set("obs"));
    mc->SetParametersOfInterest(*wspace->set("poi"));
    mc->SetNuisanceParameters(*wspace->set("nuis"));
    mc->SetGlobalObservables(RooArgSet());
    wspace->import(*mc);

    wspace->Print();

    std::cout << "Writing on " << out_filename << std::endl;
    wspace->writeToFile(out_filename);
    std::cout << "All OK." << std::endl;

    // CLs test

    Int_t type = 0;     
// type = 0 Freq calculator   (for toys)
// type = 1 Hybrid calculator (don't use this)
// type = 2 Asymptotic calculator  (for asymptotic approximation)
// type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
//
    Int_t testStatType = 3;   
// testStatType = 0 LEP
//              = 1 Tevatron 
//              = 2 Profile Likelihood two sided
//              = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)  <-- for Upper Limits (use this)
//              = 4 Profile Likelihood signed ( pll = -pll if mu < mu_hat) 
//              = 5 Max Likelihood Estimate as test statistic
//              = 6 Number of observed event as test statistic

    Double_t par_poi_min = 0.0;   // mu scanned from par_poi_min to par_poi_max with par_npointscan steps
    Double_t par_poi_max = 2.0;
    Int_t    par_npointscan = 100; // default: 100
    Int_t    par_ntoys = 2500; // number of events in Asimov sample in case of type 2 or 3, number of events in each toys for type 0; defaul: 50000

    Double_t score = StandardHypoTestInvDemo(0,"", out_filename, "wspace", "mc", "mc", "obsData", type, testStatType, true, par_npointscan, par_poi_min, par_poi_max, par_ntoys);

    return score;
}

// driver (blind: kTRUE --> blinded signal region,  kFALSE: unblinded
void run_ABCD(Bool_t blind=kTRUE) {

  // inputs Lepton Jets
   Double_t xs[4] = {94.2786, 20.5344, 0.263938, 1.79478}; //  signal
  //Double_t xs[4] = {94.2786/100, 20.5344/100, 0.263938/100, 1.79478/100}; // scaled signal
  Double_t xn[4] = {46.,21.,6.,10.};  // Lepton Jets obs data
  
  //  CalRatio signal
  // Cristiano's Xsec; regions as defined in our note: A = B*C/D. Here it's defined A = B*D/C so changing C <--> D as def in the note
  //  Double_t xs[4] = {7.4715,1.12939,0.0939429,0.494767}; // CalRatio signal, 1000_150_9m; Cristiano's Xsec = 0.1233 pb
  //Double_t xs[4] = {119.659,17.4703,1.00324,7.26485}; // CalRatio signal, 600_150_9m; Cristiano's Xsec = 2pb: 
  // Double_t xs[4] = {415.528,62.7199,3.68,25.76}; // CalRatio signal, 400_100_9m; Cristiano's Xsec = 9.5 pb
  
  // 1pb Xsec
  //Double_t xs[4] = {7.4715/0.1233,1.12939/0.1233,0.0939429/0.1233,0.494767/0.1233}; // CalRatio signal, 1000_150_9m; 1pb Xsec
  //Double_t xs[4] = {119.659/2,17.4703/2,1.00324/2,7.26485/2}; // CalRatio signal, 600_150_9m; 1pb Xsec
  //Double_t xs[4] = {415.528/9.5,62.7199/9.5,3.68/9.5,25.76/9.5}; // CalRatio signal, 400_100_9m; 1pb Xsec
  

  //Double_t xn[4] = {1.,21.,12.,10.};  // CalRatio obs data: regions as defined in our note: A = B*C/D. Here it's defined A = B*D/C so changing C <--> D as def in the note. Setting A to 1 (blinded)

    Double_t xb[4] = {0.,0.,0.,0.}; // MC BG not used in this example
    Double_t xc[4] = {0.,0.,0.,0.}; // other data-driven BG not used in this example

    Double_t z = simultaneousABCD(xn, xs, xb, xc,"ABCD_ws_test.root", kFALSE, kFALSE, blind);

    std::cout << "Input: obs  A/B/C/D: " << xn[0] << " / " << xn[1] << " / " << xn[2] << " / " << xn[3] << endl;
    std::cout << "Input: sig  A/B/C/D: " << xs[0] << " / " << xs[1] << " / " << xs[2] << " / " << xs[3] << endl;
    std::cout << "Result: score: " << z << std::endl;

}
