#ifndef __HypoTestInvTool__
#define __HypoTestInvTool__

#include "RooStats/HypoTestInverterResult.h"

//
// Tool to run XXXX
//
class HypoTestInvTool {

public:
	HypoTestInvTool();
	~HypoTestInvTool() {};

	RooStats::HypoTestInverterResult *
		RunInverter(Int_t enne, TString esse, RooWorkspace * w,
			const char * modelSBName, const char * modelBName,
			const char * dataName,
			int type, int testStatType,
			bool useCLs,
			int npoints, double poimin, double poimax, int ntoys,
			bool useNumberCounting = false,
			const char * nuisPriorName = 0);



	Double_t
		AnalyzeResult(RooStats::HypoTestInverterResult * r,
			int calculatorType,
			int testStatType,
			bool useCLs,
			int npoints,
			const char * fileNameBase = 0);

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
	bool mNoSystematics;
	std::string mMassValue;
	std::string mMinimizerType;                  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
	TString     mResultFileName;
};

#endif
