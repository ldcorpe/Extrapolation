/////////////////////////////////////////////////////////////////////////
//   This class has been automatically generated 
//   (at Tue Aug 23 14:49:37 2016 by ROOT version 5.34/36)
//   from TTree extrapTree/Used as input for the extrapolation
//   found on file: C:\Users\gordo\Downloads\LLPExtrapolationMCTree-200pi50lt5mW.root
/////////////////////////////////////////////////////////////////////////


#ifndef ntuple_extrapTree_h
#define ntuple_extrapTree_h

// System Headers needed by the proxy
#if defined(__CINT__) && !defined(__MAKECINT__)
   #define ROOT_Rtypes
   #define ROOT_TError
#endif
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TPad.h>
#include <TH1.h>
#include <TSelector.h>
#include <TBranchProxy.h>
#include <TBranchProxyDirector.h>
#include <TBranchProxyTemplate.h>
#include <TFriendProxy.h>
using namespace ROOT;

// forward declarations needed by this particular proxy


// Header needed by this particular proxy


class junk_macro_parsettree_extrapTree_Interface {
   // This class defines the list of methods that are directly used by ntuple_extrapTree,
   // and that can be overloaded in the user's script
public:
   void junk_macro_parsettree_extrapTree_Begin(TTree*) {}
   void junk_macro_parsettree_extrapTree_SlaveBegin(TTree*) {}
   Bool_t junk_macro_parsettree_extrapTree_Notify() { return kTRUE; }
   Bool_t junk_macro_parsettree_extrapTree_Process(Long64_t) { return kTRUE; }
   void junk_macro_parsettree_extrapTree_SlaveTerminate() {}
   void junk_macro_parsettree_extrapTree_Terminate() {}
};


class ntuple_extrapTree : public TSelector, public junk_macro_parsettree_extrapTree_Interface {
public :
   TTree          *fChain;         //!pointer to the analyzed TTree or TChain
   TH1            *htemp;          //!pointer to the histogram
   TBranchProxyDirector fDirector; //!Manages the proxys

   // Optional User methods
   TClass         *fClass;    // Pointer to this class's description

   // Proxy for each of the branches, leaves and friends of the tree
   TIntProxy      eventNumber;
   TBoolProxy     PassedCalRatio;
   TDoubleProxy   llp1_pt;
   TDoubleProxy   llp1_eta;
   TDoubleProxy   llp1_phi;
   TDoubleProxy   llp1_E;
   TDoubleProxy   llp1_Lxy;
   TDoubleProxy   llp2_pt;
   TDoubleProxy   llp2_eta;
   TDoubleProxy   llp2_phi;
   TDoubleProxy   llp2_E;
   TDoubleProxy   llp2_Lxy;
   TDoubleProxy   event_weight;
   TBoolProxy     RegionA;
   TBoolProxy     RegionB;
   TBoolProxy     RegionC;
   TBoolProxy     RegionD;


   ntuple_extrapTree(TTree *tree=0) : 
      fChain(0),
      htemp(0),
      fDirector(tree,-1),
      fClass                (TClass::GetClass("ntuple_extrapTree")),
      eventNumber   (&fDirector,"eventNumber"),
      PassedCalRatio(&fDirector,"PassedCalRatio"),
      llp1_pt       (&fDirector,"llp1_pt"),
      llp1_eta      (&fDirector,"llp1_eta"),
      llp1_phi      (&fDirector,"llp1_phi"),
      llp1_E        (&fDirector,"llp1_E"),
      llp1_Lxy      (&fDirector,"llp1_Lxy"),
      llp2_pt       (&fDirector,"llp2_pt"),
      llp2_eta      (&fDirector,"llp2_eta"),
      llp2_phi      (&fDirector,"llp2_phi"),
      llp2_E        (&fDirector,"llp2_E"),
      llp2_Lxy      (&fDirector,"llp2_Lxy"),
      event_weight  (&fDirector,"event_weight"),
      RegionA       (&fDirector,"RegionA"),
      RegionB       (&fDirector,"RegionB"),
      RegionC       (&fDirector,"RegionC"),
      RegionD       (&fDirector,"RegionD")
      { }
   ~ntuple_extrapTree();
   Int_t   Version() const {return 1;}
   void    Begin(::TTree *tree);
   void    SlaveBegin(::TTree *tree);
   void    Init(::TTree *tree);
   Bool_t  Notify();
   Bool_t  Process(Long64_t entry);
   void    SlaveTerminate();
   void    Terminate();

   ClassDef(ntuple_extrapTree,0);


//inject the user's code
#include "junk_macro_parsettree_extrapTree.C"
};

#endif


#ifdef __MAKECINT__
#pragma link C++ class ntuple_extrapTree;
#endif


inline ntuple_extrapTree::~ntuple_extrapTree() {
   // destructor. Clean up helpers.

}

inline void ntuple_extrapTree::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fDirector.SetTree(fChain);
   if (htemp == 0) {
      htemp = fDirector.CreateHistogram(GetOption());
      htemp->SetTitle("junk_macro_parsettree_extrapTree.C");
      fObject = htemp;
   }
}

Bool_t ntuple_extrapTree::Notify()
{
   // Called when loading a new file.
   // Get branch pointers.
   fDirector.SetTree(fChain);
   junk_macro_parsettree_extrapTree_Notify();
   
   return kTRUE;
}
   

inline void ntuple_extrapTree::Begin(TTree *tree)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   junk_macro_parsettree_extrapTree_Begin(tree);

}

inline void ntuple_extrapTree::SlaveBegin(TTree *tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   Init(tree);

   junk_macro_parsettree_extrapTree_SlaveBegin(tree);

}

inline Bool_t ntuple_extrapTree::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.

   // WARNING when a selector is used with a TChain, you must use
   //  the pointer to the current TTree to call GetEntry(entry).
   //  The entry is always the local entry number in the current tree.
   //  Assuming that fChain is the pointer to the TChain being processed,
   //  use fChain->GetTree()->GetEntry(entry).


   fDirector.SetReadEntry(entry);
   junk_macro_parsettree_extrapTree();
   junk_macro_parsettree_extrapTree_Process(entry);
   return kTRUE;

}

inline void ntuple_extrapTree::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   junk_macro_parsettree_extrapTree_SlaveTerminate();
}

inline void ntuple_extrapTree::Terminate()
{
   // Function called at the end of the event loop.
   htemp = (TH1*)fObject;
   Int_t drawflag = (htemp && htemp->GetEntries()>0);
   
   if (gPad && !drawflag && !fOption.Contains("goff") && !fOption.Contains("same")) {
      gPad->Clear();
   } else {
      if (fOption.Contains("goff")) drawflag = false;
      if (drawflag) htemp->Draw(fOption);
   }
   junk_macro_parsettree_extrapTree_Terminate();
}
