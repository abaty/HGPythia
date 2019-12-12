#include "Pythia8/Pythia.h"
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFolder.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>

using namespace Pythia8;
Int_t EvalMult(Pythia * pythia,TClonesArray &arr, Double_t etamin1=2.5, Double_t etamax1=5., Double_t etamin2=100, Double_t etamax2=100);

// 1 mb = 0.1 fm^2
// 1 b  = 100 fm^2
// 1 fm^2 = 0.01 b
// 1 fm^2 = 10 mb

void glauberPythia(const int RNGSeedOffset = 0,
                   const std::string outputName = "gtree.root",
                   const Int_t nev        = 100,
		   const Double_t energy  = 5020,
                   const Bool_t dotree    = 1,
                   const Double_t sigHin  = -1,		   
		   const Double_t sigS    = 57.,
                   const Int_t A          = 208,
                   const Int_t B          = 208,
		   const Double_t bmin    = 0., 
		   const Double_t bmax    = 20., 
		   const Bool_t shadow    = kFALSE, 
		   const Bool_t elastic   = kFALSE) 
{
  // libraries needed ...
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");

  //initialize pythia
  Pythia pythia;
  pythia.readString("Random:SetSeed = " + std::to_string( 35789 + RNGSeedOffset));
  pythia.readString("Beams:eCM = " + std::to_string( (double) energy) );
  pythia.readString("SoftQCD:inelastic = on");
  pythia.init();

  TStopwatch swatch;
  swatch.Start();

  Double_t cetamin=2.5;
  Double_t cetamax=5;

  // hard cross-section to determine MPI in mb
  Double_t sigH = sigHin;
  if (sigHin<0) { // all for ptmin=2
    if (shadow) {
      if (energy<201)
	sigH = 7.71968269;
      else if (energy<2800)
	sigH = 39.5451393;
      else if (energy<5100)
	sigH = 54.2068253;
      else if (energy<8100)
	sigH = 69.7713470;
    } else {
      if (energy<20) {
	sigH = 0.161440969;
	cetamin=1.5;
	cetamax=2.5;
      } else if (energy<40){
	sigH = 1.07414019;
	cetamin=1.5;
	cetamax=2.5;
      } else if (energy<64) {
	sigH = 2.32993174;
	cetamin=2;
	cetamax=3;
      } else if (energy<201) {
	sigH = 11.6641903;
	cetamin=3; //phenix
	cetamax=4;
      } else if (energy<2800)
	sigH = 85.2298813;
      else if (energy<5100)
	sigH = 124.296341;
      else if (energy<5500)
	sigH = 130.82;
      else if (energy<6400)
	sigH = 144.17;
      else if (energy<8100)
	sigH = 166.184998;
    }
  }

  //some distributions for heavy ions
  TF1* rwsF[2];
  if (A==208)
    rwsF[0] =  new TF1("wsPb", "7.208e-4*4.*TMath::Pi()*x^2/(1+exp((x-6.62)/0.546))", 0., 20.);
  else if (A==197)
    rwsF[0] =  new TF1("wsAu", "8.596e-04*4.*TMath::Pi()*x^2/(1+exp((x-6.38)/0.535))", 0., 20.);
  else if (A==129) {
    TF2 *f = new TF2("wsXe2a","x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",0,15,0.0,TMath::Pi());
    rwsF[0]=f; 
    f->SetNpx(120);
    f->SetNpy(120);
    f->SetParameters(5.36,0.59,0.18,0);
  }
  else if (A==40)
    rwsF[0] =  new TF1("wsAr", "1.*TMath::Pi()*x^2/(1+exp((x-3.53)/0.542))", 0., 20.);
  else if (A==16) {
    rwsF[0] =  new TF1("wsO", "x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0, 20.);
    rwsF[0]->SetParameters(2.608,0.513,-0.051);
  } else if (A==6)
    rwsF[0] =  new TF1("wsC", "7.208e-4*4.*TMath::Pi()*x^2*(1.-0.149*(x/2.46)**2)/(1+exp((x-2.46)/0.522))", 0., 20.);
  else if (A==3)
    rwsF[0] =  new TF1("wsHe", "7.208e-4*4.*TMath::Pi()*x^2*(1.+0.517*(x/0.964)**2)/(1+exp((x-0.964)/0.322))", 0., 20.);
  else 
    rwsF[0] = new TF1("prot","x*x*exp(-x/0.234)",0,5);

  if (B==208)
    rwsF[1] =  new TF1("wsPb", "7.208e-4*4.*TMath::Pi()*x^2/(1+exp((x-6.62)/0.546))", 0., 20.);
  else if (B==197)
    rwsF[1] =  new TF1("wsAu", "8.596e-04*4.*TMath::Pi()*x^2/(1+exp((x-6.38)/0.535))", 0., 20.);
  else if (B==129) {
    TF2 *f = new TF2("wsXe2a","x*x*TMath::Sin(y)/(1+exp((x-[0]*(1+[2]*0.315*(3*pow(cos(y),2)-1.0)+[3]*0.105*(35*pow(cos(y),4)-30*pow(cos(y),2)+3)))/[1]))",0,15,0.0,TMath::Pi());
    rwsF[1]=f; 
    f->SetNpx(120);
    f->SetNpy(120);
    f->SetParameters(5.36,0.59,0.18,0);
  }
  else if (B==40)
    rwsF[1] =  new TF1("wsAr", "1.*TMath::Pi()*x^2/(1+exp((x-3.53)/0.542))", 0., 20.);
  else if (B==16) {
    rwsF[1] =  new TF1("wsO", "x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0, 20.);
    rwsF[1]->SetParameters(2.608,0.513,-0.051);
  } else if (B==6)
    rwsF[1] =  new TF1("wsC", "7.208e-4*4.*TMath::Pi()*x^2*(1.-0.149*(x/2.46)**2)/(1+exp((x-2.46)/0.522))", 0., 20.);
  else if (B==3)
    rwsF[1] =  new TF1("wsHe", "7.208e-4*4.*TMath::Pi()*x^2*(1.+0.517*(x/0.964)**2)/(1+exp((x-0.964)/0.322))", 0., 20.);
  else if (B==1)
    rwsF[1] = new TF1("prot","x*x*exp(-x/0.234)",0,5);

  // matter distribution in proton
  TF1* eik = new TF1("eik", "[0]^2/[1] * ([0]*x)^3*TMath::BesselK(3, [0] * x)", 0., 10.);
  eik->SetParameter(0, 3.9);
  eik->SetParameter(1, 96);
 
  // Ranges for histos
  Float_t npartMax, ncolMax;
  Float_t npartMin, ncolMin;
  Int_t   ncolB, npartB;
  
  if (B > 4) {
    npartMax = A+B;
    ncolMax  = 2500.;
    npartMin =    0.;
    ncolMin  =    0.;
    ncolB    = 2500.;
    npartB   =  100.;
  } else if (B == 4) {
    npartMax = 100;
    ncolMax  = 100;
    npartMin =    0.;
    ncolMin  =    0.;
    ncolB    =  100.;
    npartB   =  100.;
  } else {
    npartMax =  30.;
    ncolMax  =  40.;
    npartMin = -0.5;
    ncolMin  = -0.5;
    ncolB    =  40.;
    npartB   =  30.;
  }

  Double_t bins[20000]={0};
  Int_t counter=0;
  Double_t be=0;
  while(be<250) {
    bins[counter]=be;
    be++;
    counter++;
  }
  while(be<1000) {
    bins[counter]=be;
    be+=5;
    counter++;
  }
  while(be<10000) {
    bins[counter]=be;
    be+=10;
    counter++;
  }
  while(be<20000) {
    bins[counter]=be;
    be+=100;
    counter++;
  }
  while(be<50000) {
    bins[counter]=be;
    be+=250;
    counter++;
  }
  bins[counter]=100000;  

  Int_t ncollbins=2500;
  Int_t nhardbins=4500;
  Int_t ntrackbins=2000;
  Int_t ntrackmax=100000;
  if (0) { //pp
    ncollbins=1;
    nhardbins=10;
    ntrackbins=200;
    ntrackmax=2000;
  }

  //define some histograms
  TH2F *hNpartVsB = new TH2F("hNpartVsB",";Npart;B (fm)",420,0,420,200,0,20);
  hNpartVsB->Sumw2();
  TH2F *hNcollVsB = new TH2F("hNcollVsB",";Ncoll;B (fm)",ncollbins,0,ncollbins,200,0,20);
  hNcollVsB->Sumw2();
  TH2F *hNhardVsB = new TH2F("hNhardVsB",";Nhard;B (fm)",nhardbins,0,nhardbins,200,0,20);
  hNhardVsB->Sumw2();
  TH2F *hNhardVsNcoll = new TH2F("hNhardVsNcoll",";Nhard;Ncoll",nhardbins,0,nhardbins,ncollbins,0,ncollbins);
  hNhardVsNcoll->Sumw2();
  TH2F *hNcollVsMult = new TH2F("hNcollVsMult",";Ncoll;Mult",ncollbins,0,ncollbins,counter,bins);
  hNcollVsMult->Sumw2();
  TH2F *hNhardVsMult = new TH2F("hNhardVsMult",";Nhard;Mult",nhardbins,0,nhardbins,counter,bins);
  hNhardVsMult->Sumw2();
  TH2F *hNtracksVsMult = new TH2F("hNtracksVsMult",";Ntracks;Mult",ntrackbins,0,ntrackmax,counter,bins);
  hNtracksVsMult->Sumw2();
  TH2F *hPtVsMult = new TH2F("hPtVsMult",";p_{T} (GeV);Mult",250,0,100,counter,bins);
  hPtVsMult->Sumw2();
  TH2F *hPtVsB = new TH2F("hPtVsB",";p_{T} (GeV);B (fm)",250,0,100,200,0,20);
  hPtVsB->Sumw2();
  TProfile *hNhcratVsB = new TProfile("hNhcratVsB",";B (fm);Nhard/Ncoll",200,0,20);
  TProfile *hNhcratVsN = new TProfile("hNhcratVsN",";Ncoll;Nhard/Ncoll",ncollbins,0,ncollbins);
  TH2F *hPtVsMPIpp = new TH2F("hPtVsMPIpp",";p_{T} (GeV);#MPI",250,0,100,30,0,30);
  hPtVsMPIpp->Sumw2();
  TH2F *hMultVsMPIpp = new TH2F("hMultVsMPIpp",";Mult;#MPI",200,0,200,30,0,30);
  hMultVsMPIpp->Sumw2();

  //variables needed for Glauber
  Double_t x,y;
  Double_t xx[2][208];
  Double_t yy[2][208];
  Int_t   wounded[2][208];
  Int_t   woundedT[2][208];
  Int_t   woundedJ[2][208];
  Int_t   nMPI[100];
  const Double_t dmax = 1.43;
  Int_t nPsiT = 0;
  Int_t nPsi2C = 0;

  // Excentricity
  TH1F* epsSTH = new TH1F("epsSTH", "", 100, -1., 1.);
  TH1F* epsNPH = new TH1F("epsNPH", "", 100, -1., 1.);

  // Ncol and Npart
  TH1F* npaH  = new TH1F("npaH", "", npartB, npartMin, npartMax);
  TH1F* ncolH = new TH1F("ncolH", "", ncolB, ncolMin,  ncolMax);
  TH2F* ncoH2 = new TH2F("ncoH2", "", ncolB, ncolMin,  ncolMax, npartB, npartMin, npartMax);
  TH1F* bH    = new TH1F("bH", "",    100, 0.,  20.);

  // NN Impact parameter
  TProfile* bnnH  = new TProfile("bnnH",  "bNN vs ncol ",  ncolB,  ncolMin,  ncolMax);
  TProfile* bnnNH = new TProfile("bnnNH", "bNN vs npart",  npartB, npartMin, npartMax);
  TProfile* bnnBH = new TProfile("bnnBH", "", 100, 0., 20.);

  // Minijets
  TProfile* jetH   = new TProfile("jetH",  "#jets per collision vs ncol", ncolB, ncolMin, ncolMax);
  TH1F*     jetTH  = new TH1F("jetTH",  "Total number of jets (normalized)", 60, ncolMin, 59.5);
  TProfile* jetBH  = new TProfile("jetBH", "#jets per collisions vs b", 100, 0, 20.);

  // Comparison hit and miss to eikonal
  TProfile* ratioH1 = new TProfile("ratioH1",  "", npartB, npartMin, npartMax);
  TProfile* ratioH2 = new TProfile("ratioH2",  "", 200, 0., 20.);
  TProfile* ncolMC  = new TProfile("ncolMC",  "", 100, 0., 20.);
  TProfile* ncolEI  = new TProfile("ncolEI",  "", 100, 0., 20.);

  // 
  TH1F* jpsiH = new TH1F("jpsiH",  " ", 100, 0., 100.);
  TH1F* wH    = new TH1F("wH",     "woundedness ", 100, -0.5, 99.5);

  // open a file for an event tree 
  TFile hfile(outputName.c_str(),"RECREATE","glauber tree");
  hfile.SetCompressionLevel(5);
  typedef struct 
  {
    Int_t mpi[20]; //no longer stored in tree
    Int_t ncoll;
    Int_t npart;
    Float_t b;
    Float_t bnn;
    Int_t mult;
    Int_t mul1;
    Int_t muv0;
    Int_t nch;
    Int_t nch1;
    Int_t nch2;
    Int_t nch3;
    Int_t nch4;
    Int_t nch5;
    Int_t nch8;
    Int_t nmpi;
  } EVENT;
    
  static EVENT gevent;
  TTree *tree = 0;
  if (dotree) {
    tree = new TTree("T","Glauber event tree");
    if (0) /*individual mpi not really needed */
      tree->Branch("mpi",    gevent.mpi,   "mpi[20]/I");
    tree->Branch("ncoll", &gevent.ncoll, "ncoll/I");
    tree->Branch("npart", &gevent.npart, "npart/I");
    tree->Branch("b",     &gevent.b,     "b/F");
    tree->Branch("bnn",   &gevent.bnn,   "bnn/F");
    tree->Branch("mult",  &gevent.mult,  "mult/I"); // forward mult symmetric
    tree->Branch("mul1",  &gevent.mul1,  "mul1/I"); // mult 1<eta<2 
    tree->Branch("muv0",  &gevent.muv0,  "muv0/I"); // forward mult v0
    tree->Branch("nch",   &gevent.nch,   "nch/I");  // midrapity mult
    tree->Branch("nch1",  &gevent.nch1,  "nch1/I"); // midrap pt > 1 
    tree->Branch("nch2",  &gevent.nch2,  "nch2/I"); // midrap pt > 2 
    tree->Branch("nch3",  &gevent.nch3,  "nch3/I"); // midrap pt > 3 
    tree->Branch("nch4",  &gevent.nch4,  "nch4/I"); // midrap pt > 4 
    tree->Branch("nch5",  &gevent.nch5,  "nch5/I"); // midrap pt > 5
    tree->Branch("nch8",  &gevent.nch8,  "nch8/I"); // midrap pt > 8
    tree->Branch("nmpi",  &gevent.nmpi,  "nmpi/I");
  }
  ncolMC->Sumw2();
  ncolEI->Sumw2();

  TClonesArray arr("TLorentzVector",99999);
  TClonesArray arrpp("TLorentzVector",9999);

  // Event loop
  Int_t    nAcc1 = 0;
  Int_t    nAcc2 = 0;
  Float_t  dsum0 = 0;
  Int_t    ncol0 = 0;

  Int_t nevrun=0;
  for (Int_t iev = 0; iev < nev; nevrun++) {
    for (Int_t i = 0; i < 20; i++) gevent.mpi[i] = 0;
    // Position the nucleons
    Float_t b = TMath::Sqrt(bmin * bmin + gRandom->Rndm()*(bmax*bmax - bmin * bmin));
    // initialise counters
    Int_t nPsi2 = 0;
    Int_t njetT = 0;
    
    for (Int_t i = 0; i < A+B; i++) { 
      Int_t j      = (i >= A)? 1 : 0;
      Int_t k      = (i >= A)? i - A : i;
      Double_t dx  = (i >= A)? -b/2. : b/2.;
      if (j == 0 || (j == 1 && B != 1)) {
	TF2 *f2 = dynamic_cast<TF2*>(rwsF[j]);
	if (f2) {
          Double_t r;
          Double_t theta;
          f2->GetRandom2(r,theta);
          Double_t phi = 2*TMath::Pi()*gRandom->Rndm();
          x = r * TMath::Sin(phi) * TMath::Sin(theta);
          y = r * TMath::Cos(phi) * TMath::Sin(theta);
	} else {
	  Double_t r      = rwsF[j]->GetRandom();
	  Double_t phi    = 2. * TMath::Pi() * gRandom->Rndm();
	  Double_t costh  = 2. * gRandom->Rndm() - 1.;
	  Double_t costh2 = costh * costh;
	  Double_t sinth  = 0.;
	  if (costh2 < 1.) {
	    sinth = TMath::Sqrt(1. - costh2);
	  }
	  x = r * sinth * TMath::Cos(phi);
	  y = r * sinth * TMath::Sin(phi);      
	}
      }
      xx[j][k] = x + dx;
      yy[j][k] = y;
    }
      
    for (Int_t i = 0; i < 208; i++) {
      wounded[0][i]  = 0;
      wounded[1][i]  = 0;
      woundedT[0][i] = 0;
      woundedT[1][i] = 0;
      woundedJ[0][i] = 0;
      woundedJ[1][i] = 0;
    }

    Int_t ncol  = 0;
    Int_t ncolT = 0;
    Int_t nPsi  = 0;
    // flag the wounded nucleons
    Float_t dsum = 0;
    for (Int_t i = 0; i < A; i++) {
      for (Int_t j = 0; j < B; j++) {
	Double_t r2 = (xx[0][i] - xx[1][j]) *  (xx[0][i] - xx[1][j]) +
                      (yy[0][i] - yy[1][j]) *  (yy[0][i] - yy[1][j]);
	if ((r2 < dmax * dmax)) {
	  woundedT[0][i] += 1;
	  woundedT[1][j] += 1;
	  ncolT++;
	  if (gRandom->Rndm() < 3.7e-3) {
	    nPsiT++;
	    nPsi++;
	    if (woundedJ[0][i] > 0 || woundedJ[1][j] > 0)
	      nPsi2++;
	    woundedJ[0][i] += 1;
	    woundedJ[1][j] += 1;
	  }
	}
	      
	Double_t rrb   = TMath::Min(1., b * b / 35.2 / 1.44);
	Double_t aphx  = 0.1 * 4./3. * 4.92 * TMath::Sqrt(1. - rrb);
	Double_t sigHS = sigH - aphx * 103.65;
	if (!shadow) sigHS = sigH;
	// Interaction probability
	Double_t d = TMath::Sqrt(r2);
	if (d > 5) continue;
	Double_t b02  = 0.5 * sigS * 0.1 / TMath::Pi();
	r2 /= b02;
	Double_t gstot0 = 1.;
	if (elastic) gstot0 = 2.*(1.-TMath::Exp(-(sigS+sigHS)/sigS*eik->Eval(0.001)));
	r2 /= gstot0;
	Double_t chi    = eik->Eval(sqrt(r2));
	Double_t gs     = (1. - TMath::Exp(-2.* (sigS+sigHS)/sigS*chi));
	Double_t gstot  = 2.*(1.-TMath::Sqrt(1-gs));
	Double_t rantot = gRandom->Rndm() * gstot0;
	if (rantot  > gstot && elastic)           continue;
	if (rantot  > gs)                         continue;
	wounded[0][i] = 1;
	wounded[1][j] = 1;
	dsum+=d;
	dsum0+=d;
	ncol++;
	ncol0++;
	// minijets
	Double_t tt = 2. * chi * sigHS/sigS;
	Double_t ts = 2. * chi;
	
	Int_t njet = 0;
	//cout << rantot << " " << tt << " " << ts << endl;
	if (rantot < (TMath::Exp(-tt) * (1.-TMath::Exp(-ts)))) {
	  gevent.mpi[0]++;
	  continue;
	}

	Double_t xr = - TMath::Log(TMath::Exp(-tt) + gRandom->Rndm()*(1.-TMath::Exp(-tt)));
	while(1) {
	  njet++;
	  xr-=TMath::Log(gRandom->Rndm());
	  //cout << "xr " << njet << " " <<  xr << " " << tt << endl;
	  if (xr > tt) break;
	}
	gevent.mpi[njet]++;
	njetT+=njet;
      } // target
    } // projectlile
  
    //Glauber is finished here
    
    if (ncol<1) {
      --nevrun;
      continue;
    }
       
 
    //
    // Generate Pythia events
    //      
    Int_t mpiT[20];
      for (Int_t j = 0; j < 20; j++) {
        mpiT[j] = gevent.mpi[j];
      }
    Int_t multP   = 0;
    Int_t multP1  = 0;
    Int_t multPv0 = 0;
    arr.Clear();
    while (true) {
      //count how many more pp events are needed and store as nc
      Int_t nc = 0;
      for (Int_t j = 0; j < 20; j++) {
	    nc += mpiT[j];
      }

      //if nc is 0 we are done....
      if (nc == 0) break;
      arrpp.Clear();
     
      //otherwise get a PYTHIA event...
      while(true){
        if (pythia.next()) break;
      }
      Int_t mpic = pythia.info.nMPI();
     
      //get multiplicities for this event 
      Int_t mult0 = EvalMult(&pythia, arrpp,cetamin,cetamax);
      Int_t mulv0 = EvalMult(&pythia, arrpp,-3.7,-1.7,2.8,5.1);
      Int_t mult1 = EvalMult(&pythia, arrpp,1,2);
      if (1) { //fill pp
	    Int_t npp = arrpp.GetEntries();
	    hMultVsMPIpp->Fill(npp,mpic);
	    for (Int_t k=0; k<npp; ++k) {
	      TVector3 *track = static_cast<TVector3*>(arrpp.At(k));
	      Double_t pt = track->Pt();
	      hPtVsMPIpp->Fill(pt,mpic);
	    }
      }
      if (mpic > 20) continue;
      //add this pp event mulitplicites to the total counters for a HI event, and subtract off the event from the number of events needed at that given number of MPIs
      if (mpiT[mpic] > 0) {
	    mpiT[mpic]--;
	    multP   += mult0;
	    multPv0 += mulv0;
	    multP1  += mult1;
        arr.AbsorbObjects(&arrpp);
      } else {
	    continue;
      }
    }

    printf("%5d %5d %13.3f %6d\n", iev, ncol, b, multP);

    //event generation is done, below is some bookkeeping histograms

    Double_t xd1[208];
    Double_t yd1[208];  
    Double_t xd2[208];
    Double_t yd2[208];
    
    Int_t iw1    = 0;
    Int_t iw2    = 0;
    Int_t npartT = 0;
    Double_t mx2   = 0.;
    Double_t my2   = 0.;
    Double_t mx    = 0.;
    Double_t my    = 0.;

    Double_t mxy = 0.;
    // Calculate excentricity
    // (standard and participant)
    for (Int_t i = 0; i < A+B; i++) { 
      Int_t j = (i >= A)? 1 : 0;
      Int_t k = (i >= A)? i - A : i;
      if (woundedT[j][k]) npartT++;
      if (ncolT > 0) wH->Fill(woundedT[j][k]);
      if (!wounded[j][k]) continue;
      if (j == 0) {
	xd1[iw1] = xx[j][k];
	yd1[iw1] = yy[j][k];	
	iw1++;
      } else {
	xd2[iw2] = xx[j][k];
	yd2[iw2] = yy[j][k];	
	iw2++;
      }
      mx  += xx[j][k];
      my  += yy[j][k];
      
      mx2 += (xx[j][k] * xx[j][k]);
      my2 += (yy[j][k] * yy[j][k]);
      mxy += (xx[j][k] * yy[j][k]);
    }

    if (ncolT > 0)  nAcc1++;
    if (npartT < 416 && npartT > 0) nPsi2C+=nPsi2;
    if ((iw1+iw2) == 0) 
      continue;
    Double_t iw = iw1+iw2;
    mx2 -= (mx*mx/iw);
    my2 -= (my*my/iw);
    mxy -= (mx*my/iw);
    
    Double_t epsST = (my2-mx2) / (mx2+my2);
    Double_t epsNP = TMath::Sqrt((my2-mx2) * (my2-mx2) + 4. * mxy * mxy) / (mx2+my2);
    epsNPH->Fill(epsNP);
    epsSTH->Fill(epsST);

    if (npartT < 416 && npartT > 0) jpsiH->Fill(Float_t(nPsi));
    if (npartT > 0) {
      npaH->Fill(Float_t(npartT));
      ncoH2->Fill(Float_t(ncolT), Float_t(npartT));
    }

    if (ncol > 0) {
      nAcc2++;
      ncolH->Fill(Float_t(ncol));
      bH->Fill(b);
      bnnH->Fill(Float_t(ncol), dsum/ncol);
      bnnNH->Fill(iw, dsum/ncol);
      bnnBH->Fill(b,  dsum/ncol);
      jetH->Fill(ncol, Float_t(njetT)/ncol * 71./125.);
      jetTH->Fill(Float_t(njetT), 1.);
      jetBH->Fill(b, Float_t(njetT)/ncol);
      
      if (npartT > 0) ratioH1->Fill(npartT, iw/npartT);
      ratioH2->Fill(b, Float_t(ncolT)/Float_t(ncol));
      if (ncol   > 0) ncolEI->Fill(b, Float_t(ncol));
      if (ncolT  > 0) ncolMC->Fill(b, Float_t(ncolT));

      Int_t vsym=multP;
      Int_t Nch1=0;
      Int_t Nch2=0;
      Int_t Nch3=0;
      Int_t Nch4=0;
      Int_t Nch5=0;
      Int_t Nch8=0;
      Int_t Nch=arr.GetEntries();
      for (Int_t i=0; i<Nch; ++i) {
	TVector3 *track = static_cast<TVector3*>(arr.At(i));
	Double_t pt = track->Pt();
	hPtVsMult->Fill(pt,vsym);
	hPtVsB->Fill(pt,b);
	if (pt>1)
	  Nch1++;
	if (pt>2)
	  Nch2++;
	if (pt>3)
	  Nch3++;
	if (pt>4)
	  Nch4++;
	if (pt>5)
	  Nch5++;
	if (pt>8)
	  Nch8++;
      }

      hNpartVsB->Fill(npartT,b);
      hNcollVsB->Fill(ncol,b);
      hNhardVsB->Fill(njetT,b);
      hNhardVsNcoll->Fill(njetT,ncol);
      hNhcratVsB->Fill(b,njetT/ncol);
      hNhcratVsN->Fill(ncol,njetT/ncol);
      hNcollVsMult->Fill(ncol,vsym);
      hNhardVsMult->Fill(njetT,vsym);
      hNtracksVsMult->Fill(Nch,vsym);

      gevent.mult  = vsym;
      gevent.mul1  = multP1;
      gevent.muv0  = multPv0;
      gevent.ncoll = ncol;
      gevent.b     = b;
      gevent.bnn   = dsum/ncol;
      gevent.npart = npartT;
      gevent.nch   = Nch;
      gevent.nch1  = Nch1;
      gevent.nch2  = Nch2;
      gevent.nch3  = Nch3;
      gevent.nch4  = Nch4;
      gevent.nch5  = Nch5;
      gevent.nch8  = Nch8;
      gevent.nmpi  = njetT;
      if (tree)
	tree->Fill();
      iev++;
    }
  } // event loop
  
  if (tree) {
    tree->Print();
  }
      
  hNpartVsB->Write();
  hNcollVsB->Write();
  hNhardVsB->Write();
  hNhardVsNcoll->Write();
  hNcollVsMult->Write();
  hNhardVsMult->Write();
  hNtracksVsMult->Write();
  hPtVsMult->Write();
  hPtVsB->Write();
  hNhcratVsB->Write();
  hNhcratVsN->Write();
  hfile.Write();
  hPtVsMPIpp->Write();
  hMultVsMPIpp->Write();
  hfile.Close();

  swatch.Stop();
  swatch.Print();

  printf("Total x-section %13.3f (barn) \n", 0.01 * TMath::Pi() * bmax * bmax * Float_t(nAcc1)/Float_t(nevrun));
  printf("Total x-section %13.3f (barn) \n", 0.01 * TMath::Pi() * bmax * bmax * Float_t(nAcc2)/Float_t(nevrun));
  printf("Mean impact parameter %13.3f \n", dsum0/ncol0);
 }
    


Int_t EvalMult(Pythia * pythia, TClonesArray &arr, Double_t etamin1, Double_t etamax1, Double_t etamin2, Double_t etamax2)
{
    Int_t nparticle = pythia->event.size();
    Int_t mult = 0;
    Int_t arri = arr.GetEntries();
    for (Int_t part=0; part<nparticle; part++) {
        if( !(pythia->event[part].isFinal())) continue;
        if( !(pythia->event[part].isCharged())) continue;
        Float_t eta = pythia->event[part].eta();
        Float_t etaabs = TMath::Abs(eta);
        if (etamin2>=99) {
            if (etaabs > etamin1 && etaabs < etamax1) {
                mult++;
            }
        } else {
            if ((eta>etamin1 && eta<etamax1) ||
                (eta>etamin2 && eta<etamax2)) {
                mult++;
            }
        }
        if (etaabs<1) {
            new(arr[arri]) TVector3( pythia->event[part].px(), pythia->event[part].py(), pythia->event[part].pz());
            arri++;
        }
    } // particle loop
    return mult;
}

int main(int argc, const char* argv[]){

  if(argc != 4)
  {
    std::cout << "Usage: HGPythia <nEvents> <jobNumber> <outputName>" << std::endl;
    return 1;
  }  

  int nEvents = std::atoi(argv[1]); 
  int RNGoffset = std::atoi(argv[2]); 
  std::string outputName = argv[3];

  //glauberPythia(RNGoffset, outputName, nEvents);//PbPb (default)
  glauberPythia( RNGoffset,outputName, nEvents, 5020, 1 , -1 , 57, 1, 1, 0 , 5);  //pp settings

  return 1;
}
