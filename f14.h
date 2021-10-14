//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep 26 01:31:23 2021 by ROOT version 5.34/36
// from TTree f14/This is F14RoutineS tree
// found on file: /data2/zbfu/Dat/UrQMD_Data/ROOT/f14_ecm007.70_imp12.48_3.dat.root
//////////////////////////////////////////////////////////

#ifndef f14_cxx
#define f14_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>
// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxfTracksOut=9999;
class f14
{
public:
   TTree *fChain;  //!pointer to the analyzed TTree or TChain
   Int_t fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   //F14RoutineS     *collision;
   UInt_t fUniqueID;
   UInt_t fBits;
   Int_t fNtrackOut;
   Int_t fEvtHdr_judge;
   Int_t fEvtHdr_eventId;
   Double_t fEvtHdr_imp;
   Int_t fColHdr_collisionId;
   Int_t fColHdr_Ntrack;
   Int_t fColHdr_timestep;
   Int_t fColHdr_Ncoll;
   Int_t fColHdr_Nelastic;
   Int_t fColHdr_Ninelastic;
   Int_t fColHdr_NpauliB;
   Int_t fColHdr_Ndecays;
   Int_t fColHdr_Nhard;
   Int_t fColHdr_Nsoft;
   Int_t fColHdr_Nanother;
   Int_t fTracksOut_;
   UInt_t fTracksOut_fUniqueID[kMaxfTracksOut]; //[fTracksOut_]
   UInt_t fTracksOut_fBits[kMaxfTracksOut];     //[fTracksOut_]
   Double_t fTracksOut_t[kMaxfTracksOut];       //[fTracksOut_]
   Double_t fTracksOut_rx[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_ry[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_rz[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_p0[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_px[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_py[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_pz[kMaxfTracksOut];      //[fTracksOut_]
   Double_t fTracksOut_mass[kMaxfTracksOut];    //[fTracksOut_]
   Int_t fTracksOut_ityp[kMaxfTracksOut];       //[fTracksOut_]
   Int_t fTracksOut_iso3[kMaxfTracksOut];       //[fTracksOut_]
   Int_t fTracksOut_chg[kMaxfTracksOut];        //[fTracksOut_]
   Int_t fTracksOut_Npc[kMaxfTracksOut];        //[fTracksOut_]
   Int_t fTracksOut_Nc[kMaxfTracksOut];         //[fTracksOut_]
   Int_t fTracksOut_pptype[kMaxfTracksOut];     //[fTracksOut_]

   // List of branches
   TBranch *b_collision_fUniqueID;           //!
   TBranch *b_collision_fBits;               //!
   TBranch *b_collision_fNtrackOut;          //!
   TBranch *b_collision_fEvtHdr_judge;       //!
   TBranch *b_collision_fEvtHdr_eventId;     //!
   TBranch *b_collision_fEvtHdr_imp;         //!
   TBranch *b_collision_fColHdr_collisionId; //!
   TBranch *b_collision_fColHdr_Ntrack;      //!
   TBranch *b_collision_fColHdr_timestep;    //!
   TBranch *b_collision_fColHdr_Ncoll;       //!
   TBranch *b_collision_fColHdr_Nelastic;    //!
   TBranch *b_collision_fColHdr_Ninelastic;  //!
   TBranch *b_collision_fColHdr_NpauliB;     //!
   TBranch *b_collision_fColHdr_Ndecays;     //!
   TBranch *b_collision_fColHdr_Nhard;       //!
   TBranch *b_collision_fColHdr_Nsoft;       //!
   TBranch *b_collision_fColHdr_Nanother;    //!
   TBranch *b_collision_fTracksOut_;         //!
   TBranch *b_fTracksOut_fUniqueID;          //!
   TBranch *b_fTracksOut_fBits;              //!
   TBranch *b_fTracksOut_t;                  //!
   TBranch *b_fTracksOut_rx;                 //!
   TBranch *b_fTracksOut_ry;                 //!
   TBranch *b_fTracksOut_rz;                 //!
   TBranch *b_fTracksOut_p0;                 //!
   TBranch *b_fTracksOut_px;                 //!
   TBranch *b_fTracksOut_py;                 //!
   TBranch *b_fTracksOut_pz;                 //!
   TBranch *b_fTracksOut_mass;               //!
   TBranch *b_fTracksOut_ityp;               //!
   TBranch *b_fTracksOut_iso3;               //!
   TBranch *b_fTracksOut_chg;                //!
   TBranch *b_fTracksOut_Npc;                //!
   TBranch *b_fTracksOut_Nc;                 //!
   TBranch *b_fTracksOut_pptype;             //!

   f14(TTree *tree = 0);
   virtual ~f14();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef f14_cxx
f14::f14(TTree *tree) : fChain(0)
{
   // if parameter tree is not specified (or zero), connect the file
   // used to generate this class and read the Tree.
   if (tree == 0)
   {
      TChain * chain = new TChain("f14");
      chain->Add("../dataSet/f14_ecm007.70_imp12.48_3.dat.root/f14");
      // chain->Add("../dataSet/f14_ecm007.70_imp12.48_4.dat.root/f14");
      tree = chain;
   }
   Init(tree);
}

f14::~f14()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t f14::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t f14::LoadTree(Long64_t entry)
{
   // Set the environment to read one entry
   if (!fChain)
      return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0)
      return centry;
   if (fChain->GetTreeNumber() != fCurrent)
   {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void f14::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree)
      return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_collision_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_collision_fBits);
   fChain->SetBranchAddress("fNtrackOut", &fNtrackOut, &b_collision_fNtrackOut);
   fChain->SetBranchAddress("fEvtHdr.judge", &fEvtHdr_judge, &b_collision_fEvtHdr_judge);
   fChain->SetBranchAddress("fEvtHdr.eventId", &fEvtHdr_eventId, &b_collision_fEvtHdr_eventId);
   fChain->SetBranchAddress("fEvtHdr.imp", &fEvtHdr_imp, &b_collision_fEvtHdr_imp);
   fChain->SetBranchAddress("fColHdr.collisionId", &fColHdr_collisionId, &b_collision_fColHdr_collisionId);
   fChain->SetBranchAddress("fColHdr.Ntrack", &fColHdr_Ntrack, &b_collision_fColHdr_Ntrack);
   fChain->SetBranchAddress("fColHdr.timestep", &fColHdr_timestep, &b_collision_fColHdr_timestep);
   fChain->SetBranchAddress("fColHdr.Ncoll", &fColHdr_Ncoll, &b_collision_fColHdr_Ncoll);
   fChain->SetBranchAddress("fColHdr.Nelastic", &fColHdr_Nelastic, &b_collision_fColHdr_Nelastic);
   fChain->SetBranchAddress("fColHdr.Ninelastic", &fColHdr_Ninelastic, &b_collision_fColHdr_Ninelastic);
   fChain->SetBranchAddress("fColHdr.NpauliB", &fColHdr_NpauliB, &b_collision_fColHdr_NpauliB);
   fChain->SetBranchAddress("fColHdr.Ndecays", &fColHdr_Ndecays, &b_collision_fColHdr_Ndecays);
   fChain->SetBranchAddress("fColHdr.Nhard", &fColHdr_Nhard, &b_collision_fColHdr_Nhard);
   fChain->SetBranchAddress("fColHdr.Nsoft", &fColHdr_Nsoft, &b_collision_fColHdr_Nsoft);
   fChain->SetBranchAddress("fColHdr.Nanother", &fColHdr_Nanother, &b_collision_fColHdr_Nanother);
   fChain->SetBranchAddress("fTracksOut", &fTracksOut_, &b_collision_fTracksOut_);
   fChain->SetBranchAddress("fTracksOut.fUniqueID", fTracksOut_fUniqueID, &b_fTracksOut_fUniqueID);
   fChain->SetBranchAddress("fTracksOut.fBits", fTracksOut_fBits, &b_fTracksOut_fBits);
   fChain->SetBranchAddress("fTracksOut.t", fTracksOut_t, &b_fTracksOut_t);
   fChain->SetBranchAddress("fTracksOut.rx", fTracksOut_rx, &b_fTracksOut_rx);
   fChain->SetBranchAddress("fTracksOut.ry", fTracksOut_ry, &b_fTracksOut_ry);
   fChain->SetBranchAddress("fTracksOut.rz", fTracksOut_rz, &b_fTracksOut_rz);
   fChain->SetBranchAddress("fTracksOut.p0", fTracksOut_p0, &b_fTracksOut_p0);
   fChain->SetBranchAddress("fTracksOut.px", fTracksOut_px, &b_fTracksOut_px);
   fChain->SetBranchAddress("fTracksOut.py", fTracksOut_py, &b_fTracksOut_py);
   fChain->SetBranchAddress("fTracksOut.pz", fTracksOut_pz, &b_fTracksOut_pz);
   fChain->SetBranchAddress("fTracksOut.mass", fTracksOut_mass, &b_fTracksOut_mass);
   fChain->SetBranchAddress("fTracksOut.ityp", fTracksOut_ityp, &b_fTracksOut_ityp);
   fChain->SetBranchAddress("fTracksOut.iso3", fTracksOut_iso3, &b_fTracksOut_iso3);
   fChain->SetBranchAddress("fTracksOut.chg", fTracksOut_chg, &b_fTracksOut_chg);
   fChain->SetBranchAddress("fTracksOut.Npc", fTracksOut_Npc, &b_fTracksOut_Npc);
   fChain->SetBranchAddress("fTracksOut.Nc", fTracksOut_Nc, &b_fTracksOut_Nc);
   fChain->SetBranchAddress("fTracksOut.pptype", fTracksOut_pptype, &b_fTracksOut_pptype);
   Notify();
}

Bool_t f14::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void f14::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t f14::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}

void f14::Loop()
{
   if (fChain == 0)
      return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0)
         break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
#endif // #ifdef f14_cxx
