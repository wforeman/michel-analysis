/////////////////////////////////////////////////////////////////////
// 
//  AnaMichel.cxx
//
//  This is a ROOT analysis macro which takes in files produced as 
//  output from the analysis/filter module "MichelAna_module" in 
//  lariatsoft v06_61_00_03.
//  
//  The input ROOT files should be kept in a directory called
//  "files" located in the same area in which the macro is run.
//  
//  Author: W. Foreman (wforeman_at_uchicago.edu)
//  Last updated Mar 2019
//
//////////////////////////////////////////////////////////////////// 

#include "anamichel.h"

// Default mode:
//bool  fDoLikelihoodFit  = true;
//bool  fDoBareElectrons  = true;
//int   fMaxMCEvts        = -200000;
//int   fMaxMCEvts_IsoEl  = -100000;

// Fast mode
bool  fDoLikelihoodFit  = false;
bool  fDoBareElectrons  = true;
int   fMaxMCEvts        = 400000;
int   fMaxMCEvts_IsoEl  = 100000;

bool  fSavePlots        = true;
int   fRunMode          = 2;

bool fHybridRes_UseChisquareWeighting = true;

float fChargeSignFilter = 0; // allow only electrons (-1), positrons (1)


float fChargeRes;
float fPhelRes;

TH1D* hElectronLifetime;
TCanvas* cdefault;

bool fUseTrueEnergy = false;
bool fRequireContainment = false;

bool fSmearTruePE = false;
bool fCorrectQuenching = true;
std::default_random_engine generator;

float fFidMarginX;
float fFidMarginY;
float fFidMarginZ;


bool fDebugFlag = false;

// ----------------------------------------
// Functions for describing parameterization
// of Q resolution as double-gaussian for Michel electron
// sample:
//TF1* funcP_resQ;
//TF1* fFuncP_resQ;
///TF1* funcP_resL;
//TF1* fFuncP_resL;
TF1* fFuncP_Q_Gaus;
TF1* funcP_Q_Gaus;
TF1* fFuncP_L_Gaus;
TF1* funcP_L_Gaus;
TF1* fFuncP_Q;
TF1* funcP_Q;
TF1* fFuncP_L;
TF1* funcP_L;
std::vector< TF1* > funcP_Q_fp;
std::vector< TF1* > funcP_L_fp;
std::vector< TF1* > funcP_Q_Gaus_fp;
std::vector< TF1* > funcP_L_Gaus_fp;
//TF1* funcP_L;
//TF1* fFuncP_L;
//std::vector< TF1* > funcP_resQ_fp;
//std::vector< TF1* > funcP_resL_fp;



// Bethe formula for energy deposition in matter from theoretical calculation
// for muons: 
TF1 *f_dEdx_mu;
TF1 *f_MPV_mu;

// ##########################################################################

TRandom2* fRand;
float refQPerCm = 0.;
float refLPerCm = 0.;
float fEField;
TF1* f_ModBoxRecomb;
TF1* f_BirksRecomb;
TF1* f_ModBoxRecomb_Efield;
TF1* f_BirksRecomb_Efield;
TF1* f_BirksRecombArgoNeuT_Efield;
TF1* f_muStopPower; // MEAN dEdx as function of residual range
TF1* f_dADCdx_from_dEdx;
TF1* f_dADCdx_from_rr;
TF1* f_dQdx_from_dEdx;
float fCalAreaConstants[2][2];
int fCalPlane=1;
int fCalType=1;
float fCalRR1 = 3.;
float fCalRR2 = 25.;
bool fCalFixB = true;
bool fCalFixConst = false;
bool fUseModBoxRecomb = true;

float fFastRatioNom = 0.469;
float fChargeSmearFactor = 0.0;
bool fUsePmt[2];

float fChargeResFromMC = 0.9;
float fPEResFromMC[2];

float fQCorrFactor = 1.0;
float fLCorrFactor = 1.0;

float fMaxPE[2];

float fWgt;

//#############################################################################
void AnaMichel(){

  if (!(gInterpreter->IsLoaded("vector")))
    gInterpreter->ProcessLine("#include <vector>");

//  gSystem->Load("libdict.so");

  fRand = new TRandom2(10);
  fSmearMode = 0;

  gROOT->Reset();
  

  fFastRatio          = fFastRatioNom; 
  fLArDensity               = 1.370; // This is what LArSoft uses (ave Run II)
  fEField                   = 0.484;
  
  fApplyCalibCorr           = false;
  
  fdTcut                    = 1800;

  fMinMuClusterSize         = 8;
  fMinElClusterSize         = 4;
  fMaxElClusterSize         = 9999; 
  fMinBraggSlope            = -9999.; //20.; //20 // ADC
  fMinMuLinearity           = 0.7; // .6
  fMinMuClusterHitsEndFit   = 5;
  fMinDecayAngle2D          = 15.;
  fMaxDecayAngle2D          = 165.;  
  fMinElShowerSize          = 4; // 5
  fMaxElShowerSize          = 500; // 60
  fMinNumPts3D              = 3; // 3 in Run 2?
  fMinFracHits3D            = 0.15;
  fFidMarginX               = 4.;
  fFidMarginY               = 4.;
  fFidMarginZ               = 4.;
  
  // phasing out these cuts...
  fMinProjDist              = -20.;
  fMinShowerFrac            = 0.0;
  fMaxCovAtBnd              = -0.8;
  fMinFracMuHitsLinear      = 0.0; // 0.5

 
  fPE_Require3DShower=false;
  fSetMinimizer = "simplex";
  pmttag[0] ="HMM PMT";
  pmttag[1] ="ETL PMT";
  diamtag[0]="#oslash76mm";
  diamtag[1]="#oslash51mm";
 
  fSetStrategy  = 2; 
  fFixScale = 0;
  fFixP[1]  = 0;
  fFixK[1]  = 0;
  fFixP[0]  = 0;
  fFixK[0]  = 0;
  fFixFastRatio = 1;
  fFixSmearPrompt = 1;
  fFixSmearTotal = 0;
  fMinimizeBothPMTs = 0;

  fDecayFit_MinDecayTime    = 350.;
  fDecayFit_MaxDecayTime    = 7150.; 
  fDecayFit_T1              = fDecayFit_MinDecayTime;
  fDecayFit_T2              = fDecayFit_MaxDecayTime;
  fDecayFit_BG              = 0.1;
  fDecayFit_TauPlus         = 2197.;
  fDecayFit_TauMinus        = 615.; 
  fDecayFit_FixFreeLifetime = true;
  fDecayFit_FixBG           = false;
  fDecayFit_ReqShower1      = false;
  fDecayFit_ReqShower0      = false;


  // Initialize histograms
  Init(fRunMode);
//  MakeWeightMap();
  //TimePlots();

}

void Init() { Init(2); }

void Init(int mode) {
 
  fRunMode = mode; 
   
  fCalAreaConstants[1][1] = 0.05525;
  fCalAreaConstants[1][0] = 0.027;
  
  fMaxPE[0]                 = -900;
  fMaxPE[1]                 = -350;
  
  if( fRunMode == 1 ) {
    runtag                    = "Run I";
    fRunStart                 = 6245;
    fRunEnd                   = 6345;
    fCalAreaConstants[0][0]   = 0.0273;
    fCalAreaConstants[0][1]   = 0.0910;
    fTimeStart                = 1434672090; // unix time for Run 6245
    fTimeEnd                  = 1435614552+60; // unix time for Run 6344
    fMaxRun                   = -1;
    fGateDelay                = 350.;
    filenameData              = "files/MichelAna_run1.root";
    filenameMC                = "files/MichelAna_mc1.root";
    filenameMC_ns             = "files/MichelAna_mc1_nosmear.root";
//    filenameData              = "files/MichelAna_run1_default.root";
//    filenameMC                = "files/MichelAna_mc2_default.root";

    
    fChargeResFromMC          = 1.0;
    fPEResFromMC[0]           = 1.0;
    fPEResFromMC[1]           = 1.0;


    fUsePmt[0]                = false;
    fUsePmt[1]                = true;

    // Last optimized 10/1/2018
    fScaleFactor[1]           =  1.016;
    fTrigEff_P[1]             =  25.74;
    fTrigEff_K[1]             =  9.098; //8.8615; 
    fSmearFactor[1]           =  0.152; //1.842;
    fSmearFactorPr[1]         =  0.0;
    fAmpCut[1]                = 80.;
    fScaleFactor[0]           = 0.;
    fSmearFactor[0]           = 0.;
    fTrigEff_P[0]             = 0.;
    fTrigEff_K[0]             = 0.;

    fMinMuPromptPE[0]         = 0.;
    fMinMuPromptPE[1]         = 0.; //30.;

  }

  if( fRunMode == 2 ) {
    runtag                    = "Run IIB";
    fUsePmt[0]                = true;
    fUsePmt[1]                = true;
    fRunStart                 = 9438;
    fRunEnd                   = 9981;
    fCalAreaConstants[0][0]      = 0.0224;
    fCalAreaConstants[0][1]      = 0.0480;
    fTimeStart                = 1465344332; // unix time for Run 9447
    fTimeEnd                  = 1469462754+60; // unix time for Run 9980
    fMinRun                   = -9448;
    fMaxRun                   = -9563; //9900;
    fGateDelay                = 350.;
    filenameData              = "files/MichelAna_run2b.root";
    filenameMC                = "files/MichelAna_mc2b.root";
    filenameMC_ns             = "files/MichelAna_mc2b_nosmear.root";
    //filenameData              = "files/MichelAna_run2b_default.root";
    //filenameMC                = "files/MichelAna_mc2_default.root";
    
    fChargeResFromMC          = 1.0;
    fPEResFromMC[0]           = 1.0;
    fPEResFromMC[1]           = 1.0;

    // Last optimized 10/1/2018
    fScaleFactor[0]           = 1.021;
    fTrigEff_P[0]             = 81.31; //82.48;
    fTrigEff_K[0]             = 8.271; //5.779;
    fSmearFactor[0]           = 0.002; //2.626; //2.08562;
    fSmearFactorPr[0]         = 0.;
    
    fScaleFactor[1]           = 0.957; //1.0112; 
    fTrigEff_P[1]             = 22.72; //15.8668; //20.5889;
    fTrigEff_K[1]             = 4.506; //4.95885; 
    fSmearFactor[1]           = 0.238; //2.548;
    fSmearFactorPr[1]         =  0.;
    
    fAmpCut[0]                = 130.;
    fAmpCut[1]                = 40.;
    
    fMinMuPromptPE[0]         = 0.; //100.;//80.;
    fMinMuPromptPE[1]         = 0.;//30.;
    
  }
    
  fDayEnd                   = float( fTimeEnd-fTimeStart)/86400.;
 
  f_muStopPower = new TF1("muStopPwr","[0]*x^[1]+[2]",0.01,100.);
    f_muStopPower->SetParameter(0, 5.44);
    f_muStopPower->SetParameter(1,-0.550);
    f_muStopPower->SetParameter(2, 1.446);
  
  f_ModBoxRecomb  = new TF1("_ModBoxRecomb",_ModBoxRecomb,0.5,100.,3);
  f_ModBoxRecomb  ->FixParameter(0,fBoxModelA);
  f_ModBoxRecomb  ->FixParameter(1,fBoxModelB);
  f_ModBoxRecomb  ->FixParameter(2,fEField);
  f_ModBoxRecomb   ->SetNpx(1000);
  
  f_ModBoxRecomb_Efield  = new TF1("_ModBoxRecomb_Efield",_ModBoxRecomb_Efield,0.5,100.,3);
  f_ModBoxRecomb_Efield  ->FixParameter(0,fBoxModelA);
  f_ModBoxRecomb_Efield  ->FixParameter(1,fBoxModelB);
  f_ModBoxRecomb_Efield  ->FixParameter(2,2.3); // dE/dx
  f_ModBoxRecomb_Efield   ->SetParName(0,"#alpha");
  f_ModBoxRecomb_Efield   ->SetParName(1,"#beta [(kV/cm)(g/cm^{2})/MeV]");
  f_ModBoxRecomb_Efield   ->SetParName(2,"dE/dx [MeV/cm]");
  f_ModBoxRecomb_Efield   ->SetLineColor(kBlue);

  f_BirksRecomb_Efield  = new TF1("_BirksRecomb_Efield",_BirksRecomb_Efield,0.5,100.,3);
  f_BirksRecomb_Efield  ->FixParameter(0,fBirksModelA);
  f_BirksRecomb_Efield  ->FixParameter(1,fBirksModelk);
  f_BirksRecomb_Efield  ->FixParameter(2,2.3); // dE/dx
  f_BirksRecomb_Efield   ->SetParName(0,"A_{B}");
  f_BirksRecomb_Efield   ->SetParName(1,"k_{B} [(kV/cm)(g/cm^{2})/MeV]");
  f_BirksRecomb_Efield   ->SetParName(2,"dE/dx [MeV/cm]");
  f_BirksRecomb_Efield   ->SetLineColor(kRed);
  
  f_BirksRecombArgoNeuT_Efield  = new TF1("_BirksRecomb_Efield",_BirksRecomb_Efield,0.5,100.,3);
  f_BirksRecombArgoNeuT_Efield  ->FixParameter(0,0.806);
  f_BirksRecombArgoNeuT_Efield  ->FixParameter(1,0.052);
  f_BirksRecombArgoNeuT_Efield  ->FixParameter(2,2.3); // dE/dx
  f_BirksRecombArgoNeuT_Efield   ->SetParName(0,"A_{B}");
  f_BirksRecombArgoNeuT_Efield   ->SetParName(1,"k_{B} [(kV/cm)(g/cm^{2})/MeV]");
  f_BirksRecombArgoNeuT_Efield   ->SetParName(2,"dE/dx [MeV/cm]");
  f_BirksRecombArgoNeuT_Efield   ->SetLineColor(kOrange+7);
  
  f_BirksRecomb  = new TF1("_BirksRecomb",_BirksRecomb,0.5,100.,3);
  f_BirksRecomb  ->FixParameter(0,fBirksModelA);
  f_BirksRecomb  ->FixParameter(1,fBirksModelk);
  f_BirksRecomb  ->FixParameter(2,fEField);
  f_BirksRecomb   ->SetNpx(1000);
 
  f_dADCdx_from_rr    = new TF1("dADCdx_from_rr",dADCdx_from_rr,1.,50.,1);
  
  f_dADCdx_from_dEdx  = new TF1("dADCdx_from_dEdx",dADCdx_from_dEdx,1.,50.,4);
  f_dADCdx_from_dEdx  ->SetParameter(1,fBoxModelA);
  f_dADCdx_from_dEdx  ->SetParameter(2,fBoxModelB);
  f_dADCdx_from_dEdx  ->SetParameter(3,fEField);
  f_dADCdx_from_dEdx  ->SetParName(0,"C_{cal} (ADC/e^{-})");
  f_dADCdx_from_dEdx  ->SetParName(1,"Box Model #alpha");
  f_dADCdx_from_dEdx  ->SetParName(2,"Box Model #beta");
  f_dADCdx_from_dEdx  ->SetParName(3,"Electric field");
  if( !fUseModBoxRecomb ) {
    f_dADCdx_from_dEdx  ->SetParameter(1,fBirksModelA);
    f_dADCdx_from_dEdx  ->SetParameter(2,fBirksModelk);
  }
  
  f_dQdx_from_dEdx  = new TF1("dQdx_from_dEdx",dQdx_from_dEdx,1.,50.,3);
  f_dQdx_from_dEdx  ->SetParameter(0,fBoxModelA);
  f_dQdx_from_dEdx  ->SetParameter(1,fBoxModelB);
  f_dQdx_from_dEdx  ->SetParameter(2,fEField);
  if( !fUseModBoxRecomb ) {
    f_dQdx_from_dEdx  ->SetParameter(0,fBirksModelA);
    f_dQdx_from_dEdx  ->SetParameter(1,fBirksModelk);
  }

  // Q+L energy liklihood (given Q and L measured)
  f_LogLEnergyQL  = new TF1("LogLEnergyQL",_logLEnergyQL,5.,150.,13);
  f_LogLEnergyQL  ->SetNpx(50);
  f_LogLEnergyQL  ->SetParName(0,"Q_ion");
  f_LogLEnergyQL  ->SetParName(1,"Q_ph");
  f_LogLEnergyQL  ->SetParName(2,"L");
  f_LogLEnergyQL  ->SetParName(3,"excRatio");
  f_LogLEnergyQL  ->SetParName(4,"recombFrac");
  f_LogLEnergyQL  ->SetParName(5,"W");
  f_LogLEnergyQL  ->SetParName(6,"aveDriftTime");
  f_LogLEnergyQL  ->SetParName(7,"elecLifetime");
  f_LogLEnergyQL  ->SetParName(8,"chargeRes");
  f_LogLEnergyQL  ->SetParName(9,"fracVis");
  f_LogLEnergyQL  ->SetParName(10,"peRes");
  f_LogLEnergyQL  ->SetParameter("excRatio",fExcRatio);
  f_LogLEnergyQL  ->SetParameter("recombFrac",fRecomb);
  f_LogLEnergyQL  ->SetParameter("W",fWion);
  f_LogLEnergyQL  ->SetParName(11,"useResFuncQ");
  f_LogLEnergyQL  ->SetParName(12,"useResFuncL");
          f_LogLEnergyQL->SetParameter("useResFuncQ",1);
          f_LogLEnergyQL->SetParameter("useResFuncL",1);



  // the Sternheimer parameterization for the density effect is used
  // values are obtained from Mokhov (2001) in atomic data
  double  x0=0.2;
  double  x1=3.0;
  double  Cbar=5.2146;
  double  a=0.19559;
  double  k=3.0;
  
  // most probable energy loss and bethe-bloch  constants
  double  K=0.307075;
  int     z=1;            // number of electrons
  int     Z=18;           // atomic number
  double  A=39.948;       // atomic mass
  double  m_e=0.511;      // electron mass [MeV]
  double  M_mu=105.6;     // muon mass [MeV]
  double  M_pi=139.6;     // pion mass [MeV]
  double  M_ka=493.7;     // kaon mass [MeV]
  double  M_proton=938.3;
  double  I=188e-6;       // mean exc. energy for LAr [MeV]
  double  width=0.67;  // 0.67    // wire pitch, detector width [cm]
  double  widthErr = 0.0002; // 0.0002
  double  rho=fLArDensity;// density of LAr [g/cm^3]
  double  j=0.2;          // Landau's constant
 
  f_dEdx_mu = new TF1("dEdx_mu","([0]*([5]*[5]/x/x+1)*(0.5*log([1]*pow(x/[5],4)/(2*sqrt(x*x/([5]*[5])+1)*[6]+1))-(1/([5]*[5]/x/x+1))-((x<[5]*exp(0.2*log(10)))*0+([5]*exp(0.2*log(10))<=x && x<[5]*exp(3.0*log(10)))*(2*log(x/[5])-[2]+[3]*pow([4]-log(x/[5])/log(10),3.))+(x>=[5]*exp(3.0*log(10)))*(2*log(x/[5])-[2]))/2)) * [7]",10,10000);
  f_dEdx_mu->SetLineColor(kGreen+2);
  f_dEdx_mu->SetParameter(0,K*Z/A);
  f_dEdx_mu->SetParameter(1,4*m_e*m_e/I/I);
  f_dEdx_mu->SetParameter(2,Cbar);
  f_dEdx_mu->SetParameter(3,a);
  f_dEdx_mu->SetParameter(4,x1);
  f_dEdx_mu->SetParameter(5,M_mu);
  f_dEdx_mu->SetParameter(6,m_e/M_mu);
  f_dEdx_mu->SetParameter(7,rho);
  f_dEdx_mu->SetNpx(100000);
 
  f_MPV_mu = new TF1("MPV_mu","([0]*([6]*[6]/x/x+1)*(log([1]*(x*x/[6]/[6]))+log([0]/[7]*([6]*[6]/x/x+1))+0.2-(1/([6]*[6]/x/x+1))-((x<[6]*exp(0.2*log(10)))*0+([6]*exp(0.2*log(10))<=x && x<[6]*exp(3.0*log(10)))*(2*log(x/[6])-[3]+[4]*pow([5]-log(x/[6])/log(10),3.))+(x>=[6]*exp(3.0*log(10)))*(2*log(x/[6])-[3])))/[8]) * [9]",10,10000);
  f_MPV_mu->SetLineColor(kGreen+1);
  f_MPV_mu->SetLineStyle(7);
  f_MPV_mu->SetParameter(0,(K/2)*(Z/A)*width*rho);
  f_MPV_mu->SetParameter(1,2*m_e/I);
  f_MPV_mu->SetParameter(3,Cbar);
  f_MPV_mu->SetParameter(4,a);
  f_MPV_mu->SetParameter(5,x1);
  f_MPV_mu->SetParameter(6,M_mu);
  f_MPV_mu->SetParameter(7,I);
  f_MPV_mu->SetParameter(8,width*rho);
  f_MPV_mu->SetParameter(9,rho);
  f_MPV_mu->SetNpx(100000);

  

  /*
  // Parameterizing Q resolution distribution for Michel MC sample (MC2)
  float Emin = 5.;
  float Emax = 55.;
  funcP_resQ = new TF1("funcP_resQ","(1./ ((1.+[1])*sqrt(2*3.14159))) * ( exp(-0.5*(x/[0])^2)/[0] + [1]*exp(-0.5*((x-[2])/[3])^2)/[3] )");
  fFuncP_resQ = new TF1("fFuncP_resQ",_funcP_resQ,0.,Emax,1);
  fFuncP_resQ ->SetParName(0,"energy");
  funcP_resQ_fp.resize(4);
  funcP_resQ_fp[0]   = new TF1("funcP_resQ_peakSig","[0]/x + [1]",0.,Emax);
    funcP_resQ_fp[0] ->SetParameter(0,1.90); // +/- 0.05
    funcP_resQ_fp[0] ->SetParameter(1,0.009558); // +/- 0.0015
  funcP_resQ_fp[1]   = new TF1("funcP_resQ_BGr","[0]+[1]*x+[2]*x^2+[3]*x^3",0,Emax);
    funcP_resQ_fp[1] ->SetParameter(0,-0.7937); // +/- 0.21
    funcP_resQ_fp[1] ->SetParameter(1,0.1448); // +/- 0.025
    funcP_resQ_fp[1] ->SetParameter(2,-0.00424); // +/- 0.0009
    funcP_resQ_fp[1] ->SetParameter(3,3.833e-05); // +/- 1e-05
  funcP_resQ_fp[2]   = new TF1("funcP_resQ_BGmean","[0]/x^3 + [1]",0.,Emax);
    funcP_resQ_fp[2]  ->SetParameter(0,609.4); // +/- 42
    funcP_resQ_fp[2]  ->SetParameter(1,-0.0957); // +/- 0.0036
  funcP_resQ_fp[3]   = new TF1("funcP_resQ_BGsig","[0]/x^3 + [1]/x + [2] + [3]*x",0.,Emax);
    funcP_resQ_fp[3]  ->SetParameter(0,963.5); // +/- 170
    funcP_resQ_fp[3]  ->SetParameter(1,-9.154); // +/- 2.2
    funcP_resQ_fp[3]  ->SetParameter(2,0.903); // +/- 0.11
    funcP_resQ_fp[3]  ->SetParameter(3,-0.009523); // +/- 0.0017
 */

  // Parameterizing Q distribution for Michel MC sample (MC2)
  float Emin = 5.;
  float Emax = 55.;
  funcP_Q = new TF1("funcP_Q","(1./ ((1.+[2])*sqrt(2*3.14159))) * ( exp(-0.5*((x-[0])/[1])^2)/[1] + [2]*exp(-0.5*((x-[3])/[4])^2)/[4] )");
  funcP_Q_fp.resize(5);
  funcP_Q_fp[0]   = new TF1("funcP_Q_muPeak","[0] + [1]*x",Emin,Emax);
    funcP_Q_fp[0] ->SetParameter(0,1630.); // +/- 4286
    funcP_Q_fp[0] ->SetParameter(1,2.77e4); // +/- 135
  funcP_Q_fp[1]   = new TF1("funcP_Q_rsigPeak","[0]/x^2 + [1]/x + [2]",Emin,Emax);
    funcP_Q_fp[1] ->SetParameter(0,2.8); // +/- 3.7
    funcP_Q_fp[1] ->SetParameter(1,1.77); // +/- 0.367
    funcP_Q_fp[1] ->SetParameter(2,0.0197); // +/- 0.0074
  funcP_Q_fp[2]   = new TF1("funcP_Q_rBG","[0]+[1]*x+[2]*x^2+[3]*x^3",Emin,Emax);
    funcP_Q_fp[2] ->SetParameter(0,-0.5937); // +/- 0.08
    funcP_Q_fp[2] ->SetParameter(1,0.0986); 
    funcP_Q_fp[2] ->SetParameter(2,-0.002144);
    funcP_Q_fp[2] ->SetParameter(3,1.211e-5); 
  funcP_Q_fp[3]   = new TF1("funcP_Q_rmuBG","[0]/x^[1] + [2]",Emin,Emax);
    funcP_Q_fp[3]  ->SetParameter(0,87.44); 
    funcP_Q_fp[3]  ->SetParameter(1,1.956);
    funcP_Q_fp[3]  ->SetParameter(2,0.8012);
  funcP_Q_fp[4]   = new TF1("funcP_Q_rsigBG","[0] + [1]*x + [2]*x^2",Emin,Emax);
    funcP_Q_fp[4]  ->SetParameter(0,2.931); 
    funcP_Q_fp[4]  ->SetParameter(1,0.00300); 
    funcP_Q_fp[4]  ->SetParameter(2,0.00156); 
  fFuncP_Q = new TF1("fFuncP_Q",_funcP_Q,-500e3,2500e3,2);
  fFuncP_Q ->SetParName(0,"energy");
  fFuncP_Q ->SetParName(1,"scale");
  fFuncP_Q ->SetParameter("scale",1);
  
  // Parameterizing L distribution for Michel MC sample (MC2)
  //float Emax = 100.;
  funcP_L = new TF1("funcP_L","exp(-0.5*((x-[0])/[1])^2)/([1]*sqrt(2*3.14159))");
  funcP_L_fp.resize(2);
  funcP_L_fp[0]   = new TF1("funcP_L_muPeak","[0] + [1]*x",Emin,Emax);
    funcP_L_fp[0] ->SetParameter(0,2800.); // +/- 2572
    funcP_L_fp[0] ->SetParameter(1,2.34e4); // +/- 83
  funcP_L_fp[1]   = new TF1("funcP_Q_rsigPeak","[0]/x^[1] + [2]",Emin,Emax);
    funcP_L_fp[1] ->SetParameter(0,0.958); // +/- 3.7
    funcP_L_fp[1] ->SetParameter(1,1.012); // +/- 0.367
    funcP_L_fp[1] ->SetParameter(2,0.0593); // +/- 0.367
  fFuncP_L = new TF1("fFuncP_L",_funcP_L,-500e3,3000e3,2);
  fFuncP_L ->SetParName(0,"energy");
  fFuncP_L ->SetParName(1,"scale");
  fFuncP_L ->SetParameter("scale",1);

  // Generic Gaussian for isolated el showers
  funcP_Q_Gaus = new TF1("funcP_Q_Gaus","exp(-0.5*((x-[0])/[1])^2)/([1]*sqrt(2*3.14159))");
  funcP_Q_Gaus_fp.resize(2);
  funcP_Q_Gaus_fp[0]   = new TF1("funcP_Q_Gaus_muPeak","[0] + [1]*x",Emin,Emax);
    funcP_Q_Gaus_fp[0] ->SetParameter(0,-9136.); // +/- 2572
    funcP_Q_Gaus_fp[0] ->SetParameter(1,2.772e4); // +/- 83
  funcP_Q_Gaus_fp[1]    = new TF1("funcP_Q_Gaus_rsigPeak","[0]/x^[1] + [2]",Emin,Emax);
    funcP_Q_Gaus_fp[1]  ->SetParameter(0,0.3488); // +/- 3.7
    funcP_Q_Gaus_fp[1]  ->SetParameter(1,0.9447); // +/- 0.367
    funcP_Q_Gaus_fp[1]  ->SetParameter(2,0.01283); // +/- 0.367
    funcP_Q_Gaus_fp[1]  ->SetParLimits(2,0.,1.);
  fFuncP_Q_Gaus = new TF1("fFuncP_Q_Gaus",_funcP_Q_Gaus,-500e3,3000e3,2);
  fFuncP_Q_Gaus ->SetParName(0,"energy");
  fFuncP_Q_Gaus ->SetParName(1,"scale");
  fFuncP_Q_Gaus ->SetParameter("scale",1);
  
  funcP_L_Gaus = new TF1("funcP_L_Gaus","exp(-0.5*((x-[0])/[1])^2)/([1]*sqrt(2*3.14159))");
  funcP_L_Gaus_fp.resize(2);
  funcP_L_Gaus_fp[0]   = new TF1("funcP_L_Gaus_muPeak","[0] + [1]*x",Emin,Emax);
    funcP_L_Gaus_fp[0] ->SetParameter(0,-8735); // +/- 2572
    funcP_L_Gaus_fp[0] ->SetParameter(1,2.556e4); // +/- 83
  funcP_L_Gaus_fp[1]   = new TF1("funcP_L_Gaus_rsigPeak","[0]/x^[1] + [2]",Emin,Emax);
    funcP_L_Gaus_fp[1] ->SetParameter(0,0.3682); // +/- 3.7
    funcP_L_Gaus_fp[1] ->SetParameter(1,0.5246); // +/- 0.367
    funcP_L_Gaus_fp[1] ->SetParameter(2,-0.008254); // +/- 0.367
    funcP_L_Gaus_fp[1]  ->SetParLimits(2,0.,1.);
  fFuncP_L_Gaus = new TF1("fFuncP_L_Gaus",_funcP_L_Gaus,-500e3,3000e3,2);
  fFuncP_L_Gaus ->SetParName(0,"energy");
  fFuncP_L_Gaus ->SetParName(1,"scale");
  fFuncP_L_Gaus ->SetParameter("scale",1);
  fFuncP_L_Gaus = new TF1("fFuncP_L_Gaus",_funcP_L_Gaus,-500e3,3000e3,2);
  fFuncP_L_Gaus ->SetParName(0,"energy");
  fFuncP_L_Gaus ->SetParName(1,"scale");
  fFuncP_L_Gaus ->SetParameter("scale",1);
  /*
  funcP_resL = new TF1("funcP_resL","exp(-0.5*(x/[0])^2) / ([0]*sqrt(2*3.14159))");
  fFuncP_resL = new TF1("fFuncP_resL",_funcP_resL,0.,Emax,1);
  fFuncP_resL ->SetParName(0,"energy");
  funcP_resL_fp.resize(1);
  funcP_resL_fp[0]   = new TF1("funcP_L_sig","sqrt([0]^2/x + [1]^2)",1.,Emax);
    funcP_resL_fp[0] ->SetParameter(0,36.89/100.); // +/- 0.82
    funcP_resL_fp[0] ->SetParameter(0,5.479/100.); // +/- 0.17
  */
  
  
  Init(filenameData, filenameMC, filenameMC_ns);

}


//################################################################################
void Init(std::string filenameData, std::string filenameMC, std::string filenameMC_nosmear){

  std::cout<<"Initializing analysis for "<<runtag.c_str()<<"\n";
  outFile = new TFile("plots.root","RECREATE");
  std::cout<<filenameData.c_str()<<"\n";
  fFile    = new TFile(filenameData.c_str(), "read");
  fTree[0]   = (TTree*)fFile->Get("michelana/anatree");
  setBranches(fTree[0]);
  fFileMC  = new TFile(filenameMC.c_str(), "read");
  fTree[1]   = (TTree*)fFileMC->Get("michelana/anatree");
  setBranches(fTree[1]);
  std::cout<<filenameMC_nosmear.c_str()<<"\n";
  fFileMC_ns = new TFile(filenameMC_nosmear.c_str(),"read");
  fTreeMC_ns   = (TTree*)fFileMC_ns->Get("michelana/anatree");
  setBranches(fTreeMC_ns);
  
  // =================================================================================
  // HISTOGRAMS
  // =================================================================================
  
  
  // Misc. diagnostic plots--------------------------
  int nbins = fRunEnd - fRunStart;
  hRunVsEnergy = new TH2D("RunVsEnergy",Form("%s;Run number;Q-based energy [MeV]",runtag.c_str()),nbins,fRunStart,fRunEnd,40,0,80); //100,0,2e6);
  hRunVsEnergy->SetOption("colz");
  hRunVsLight = new TH2D("RunVsLight",Form("%s;Run number; Reconstructed shower photons [#gamma]",runtag.c_str()),nbins,fRunStart,fRunEnd,40,0,2e6); //100,0,2e6);
  hRunVsLight->SetOption("colz");
  hRunVsCharge = new TH2D("RunVsCharge",Form("%s;Run number; Reconstructed shower charge [e-]",runtag.c_str()),nbins,fRunStart,fRunEnd,40,0,4e6); //100,0,2e6);
  hRunVsCharge->SetOption("colz");
  hRunVsElectronLifetime = new TH2D("RunVsElectronLifetime",Form("%s;Run number; Electron lifetime from database [#mus]",runtag.c_str()),nbins,fRunStart,fRunEnd,80,0,1600); //100,0,2e6);
  hRunVsElectronLifetime->SetOption("colz");

  hElShowerCentroid_X[0]  = new TH1D("ElShowerCentroid_X",";X [cm]",24,0.,47.5);
  hElShowerCentroid_Y[0]  = new TH1D("ElShowerCentroid_Y",";Y [cm]",20,-20.,20);
  hElShowerCentroid_Z[0]  = new TH1D("ElShowerCentroid_Z",";Z [cm]",18,0.,90.);
  hElShowerCentroid_X[1]  = (TH1D*)hElShowerCentroid_X[0]->Clone("ElShowerCentroid_X_mc"); 
  hElShowerCentroid_Y[1]  = (TH1D*)hElShowerCentroid_Y[0]->Clone("ElShowerCentroid_Y_mc"); 
  hElShowerCentroid_Z[1]  = (TH1D*)hElShowerCentroid_Z[0]->Clone("ElShowerCentroid_Z_mc"); 
  hElShowerVis[0]         = new TH1D("ElShowerVis",";Fractional visibility (#times 10^{-3})",75,0.,1.5);
  hElShowerVis[1]         = (TH1D*)hElShowerVis[0]->Clone("ElShowerVis_mc");
  
  hTrueMuEndpt            = (TH3D*)fFileMC->Get("michelana/true/TrueMuEndpt");
  hMCSpatialWgt           = (TH3D*)hTrueMuEndpt->Clone("MCSpatialWgt");
  hTrueMuEndpt_ZX         = (TH2D*)fFileMC->Get("michelana/true/TrueMuEndpt_ZX");
  hMCSpatialWgt_ZX        = (TH2D*)hTrueMuEndpt_ZX->Clone("MCSpatialWgt_ZX");
  
  hElShowerCentroid[0]          = new TH3D("ElShowerCentroid","",3,0.,47.5,3,-20.,20.,5,0.,90.);
  hElShowerCentroid[1]          = (TH3D*)hElShowerCentroid[0]->Clone("ElShowerCentroid_mc");
  //hMCSpatialWgt                 = (TH3D*)hElShowerCentroid[0]->Clone("MCSpatialWgt");
  
  // 12-hour bins
  nbins = (fTimeEnd-fTimeStart)/(24*60*60);
  float daytotal = float(fTimeEnd-fTimeStart)/86400.;
  hTimeVsLight = new TH2D("TimeVsLight",";Time (days);Photons #it{L} [#gamma]",nbins,0.,daytotal, 40,0.,2e6); //100,0,2e6);
  hTimeVsLight->SetOption("colz");
  hTimeVsCharge = new TH2D("TimeVsCharge",";Time (days);Charge #it{Q} [e-]",nbins,0.,daytotal,40,0,4e6); //100,0,2e6);
  hTimeVsCharge->SetOption("colz");
  hTimeVsElectronLifetime = new TH2D("TimeVsElectronLifetime",";Time (days); Electron lifetime from database [#mus]",nbins,0.,daytotal,80,0,1800); //100,0,2e6);
  hTimeVsElectronLifetime->SetOption("colz");
  hChargePerCm        = new TH1D("ChargePerCm",";Charge per unit length [e-/cm]",40,0,2e5);
  hLightPerCm        = new TH1D("LightPerCm",";Photons per unit length [#gamma/cm]",40,0,1e5);
  hTimeVsChargePerCm  = new TH2D("TimeVsChargePerCm",";Time (days);Charge per unit length [e-/cm]",nbins,0.,daytotal,40,0,2e5);
  hTimeVsChargePerCm->SetOption("colz");
  hTimeVsLightPerCm  = new TH2D("TimeVsLightPerCm",";Time (days);Photons per unit length [#gamma/cm]",nbins,0.,daytotal,40,0,1e5);
  hTimeVsLightPerCm->SetOption("colz");

  
  hDecayTimeVs = new TH2D("DecayTimeVs","DecayTimeVs",150,0,7500,80,0,80);
  hQvsL[0] = new TH2D("QvsL","Q vs. L",100,0.,4e6,100,0,2e6);
  hQvsL[1] = (TH2D*)hQvsL[0]->Clone("QvsL_MC");
  
  hElectronLifetime = new TH1D("ElectronLifetime",";Electron lifetime [#mus]",100,0.,2000.);
 
  // Light plots ------------------------------------ 
  float x1_S[2]             = {   0.,     0.    };
  float x2_S[2]             = {   1500.,  500.  };
  int nbins_S[2]            = {   50,     50    };
  float x1_S100[2]          = {   0.,    0.    };
  float x2_S100[2]          = {   400.,  140.  };
  int nbins_S100[2]         = {   40,    35    };
  hPE_prompt[0][0]          = new TH1D("0_PE_prompt","HMM Prompt Light;Prompt light detected by HMM PMT [pe];Events / 10pe",nbins_S100[0],x1_S100[0],x2_S100[0]);
  hPE_prompt[0][1]          = new TH1D("1_PE_prompt","ETL Prompt Light;Prompt light detected by ETL PMT [pe];Events / 4pe",nbins_S100[1],x1_S100[1],x2_S100[1]);
  hPE_prompt[1][0]          = (TH1D*)hPE_prompt[0][0]->Clone("0_PE_prompt_mc");
  hPE_prompt[1][1]          = (TH1D*)hPE_prompt[0][1]->Clone("1_PE_prompt_mc");
  hPE_total[0][0]           = new TH1D("0_PE_total","HMM Total Light;Total light detected by HMM PMT [pe];Events / 30pe",nbins_S[0],x1_S[0],x2_S[0]);
  hPE_total[0][1]           = new TH1D("1_PE_total","ETL Total Light;Total light detected by ETL PMT [pe];Events / 10pe",nbins_S[1],x1_S[1],x2_S[1]);
  hPE_total[1][0]           = (TH1D*)hPE_total[0][0]->Clone("0_PE_total_mc");
  hPE_total[1][1]           = (TH1D*)hPE_total[0][1]->Clone("1_PE_total_mc");
  hPE_total_qc[0][0]           = new TH1D("0_PE_total_qc","HMM Total Light;Total light detected by HMM PMT [pe];Events / 30pe",nbins_S[0],x1_S[0],x2_S[0]);
  hPE_total_qc[0][1]           = new TH1D("1_PE_total_qc","ETL Total Light;Total light detected by ETL PMT [pe];Events / 10pe",nbins_S[1],x1_S[1],x2_S[1]);
  hPE_total_qc[1][0]           = (TH1D*)hPE_total[0][0]->Clone("0_PE_total_qc_mc");
  hPE_total_qc[1][1]           = (TH1D*)hPE_total[0][1]->Clone("1_PE_total_qc_mc");
  hPE[0]                    = new TH1D("PE_Data","Total photoelectrons;Total combined light [pe];Events / 10pe",150,0,1500);
  hPE[1]                    = (TH1D*)hPE[0]->Clone("PE_MC");
  hPE_prompt_pretrigMC[0]   = (TH1D*)hPE_prompt[1][0]->Clone("0_PE_prompt_pretrigmc");
  hPE_prompt_pretrigMC[1]   = (TH1D*)hPE_prompt[1][1]->Clone("1_PE_prompt_pretrigmc");
  hPE_prompt_posttrigMC[0]  = (TH1D*)hPE_prompt[1][0]->Clone("0_PE_prompt_posttrigmc");
  hPE_prompt_posttrigMC[1]  = (TH1D*)hPE_prompt[1][1]->Clone("1_PE_prompt_posttrigmc");
  hPE_prompt_nosmearMC[0]   = (TH1D*)hPE_prompt[1][0]->Clone("0_PE_prompt_nosmear");
  hPE_prompt_nosmearMC[1]   = (TH1D*)hPE_prompt[1][1]->Clone("1_PE_prompt_nosmear");
  hPE_total_nosmearMC[0]    = (TH1D*)hPE_total[1][0]->Clone("0_PE_total_nosmear");
  hPE_total_nosmearMC[1]    = (TH1D*)hPE_total[1][1]->Clone("1_PE_total_nosmear");
  hLightYield[0]            = new TH1D("LightYield_Data",";Light yield [pe/MeV];Events",300,0.,30.);
  hLightYield[1]            = (TH1D*)hLightYield[0]->Clone("LightYield_MC");
 
  hTrue_LY                  = new TH1D("True_LY","True LY;Light yield per event [pe/MeV]",500,0.,50.);
  
  hLY_ZX[0]                 = new TH2D("LY_ZX_Data","Measured light yield;Z [cm];X [cm]", 9,0.,90.,5,0.,47.5);
  hLY_ZX[0]                 ->SetOption("colz");
  hLY_ZX[1]                 = (TH2D*)hLY_ZX[0]->Clone("LY_ZX_MC");
  hE_ZX[0]                  = new TH2D("E_ZX_Data","Deposited energy;Z [cm];X [cm]", 9,0.,90.,5,0.,47.5);
  hE_ZX[1]                  = (TH2D*)hE_ZX[0]->Clone("E_ZX_MC");
  hLY_entries_ZX[0]         = (TH2D*)hLY_ZX[0]->Clone("hLY_entries_ZX_Data");
  hLY_entries_ZX[1]         = (TH2D*)hLY_ZX[0]->Clone("hLY_entries_ZX_MC");
  hLY_entries_ZY[0]         = (TH2D*)hLY_ZX[0]->Clone("hLY_entries_ZY_Data");
  hLY_entries_ZY[1]         = (TH2D*)hLY_ZX[0]->Clone("hLY_entries_ZY_MC");

  hLY_ZY[0]                 = new TH2D("LY_ZY_Data","Measured light yield;Z [cm];Y [cm]", 9,0.,90.,5,-20.,20.);
  hLY_ZY[0]->SetOption("colz");
  hLY_ZY[1]                 = (TH2D*)hLY_ZY[0]->Clone("LY_ZY_MC");
  hE_ZY[0]                  = new TH2D("E_ZY_Data","Deposited energy;Z [cm];Y [cm]", 9,0.,90.,5,-20.,20.);
  hE_ZY[1]                  = (TH2D*)hE_ZY[0]->Clone("E_ZY_MC");

  hPE_total_compare         = new TH2D("PE_total_compare",";Total light in HMM [pe];Total light in ETL [pe]",
    nbins_S[0]*2., x1_S[0], x2_S[0]*2.,
    nbins_S[1]*2., x1_S[1], x2_S[1]*2.);

  // Prompt light fraction  ------------------------
  hPromptFrac[0]        = new TH1D("PromptFrac_Data","Prompt fraction",100,0.,1.);
  hPromptFrac[1]        = (TH1D*)hPromptFrac[0]->Clone("PromptFrac_MC");
  hPromptFracRaw[0]        = new TH1D("PromptFracRaw_Data","Prompt fraction (prior to mu late-light corr.)",100,0.,1.);
  hPromptFracRaw[1]        = (TH1D*)hPromptFrac[0]->Clone("PromptFracRaw_MC");

  // Shower anglular metric ------------------------
  hShowerAngleMean[0]   = new TH1D("ShowerAngleMean_Data","Mean angle of shower products relative to electron dir.;#theta [rad]",100,0.,60.);
  hShowerAngleRMS[0]   = new TH1D("ShowerAngleRMS_Data","Angular RMS of shower products relative to electron dir.;RMS [rad]",100,0.,60.);
  hShowerAngle[0]   = new TH2D("ShowerAngle_Data","#theta [rad];RMS [rad]",100,0.,60.,100,0.,60.);
  hShowerAngle[0]->SetOption("colz"); 
  hShowerAngleMean[1] = (TH1D*)hShowerAngleMean[0]->Clone("ShowerAngleMean_MC"); 
  hShowerAngleRMS[1]  = (TH1D*)hShowerAngleRMS[0]->Clone("ShowerAngleRMS_MC");  
  hShowerAngle[1]   =   (TH2D*)hShowerAngle[0]->Clone("ShowerAngle_MC");

  for(int i=0; i<2; i++){ 
    for(int j=0; j<2; j++){
      hPE_total_vs_prompt[i][j] = new TH2D(Form("PE_total_vs_prompt_%d_pmt%d",i,j),Form("PMT %d (MC=%d);Total light [pe];Prompt light [pe]",j,i),
        nbins_S[j],x1_S[j],x2_S[j], nbins_S100[j],x1_S100[j],x2_S100[j]);
      hPE_total_vs_prompt[i][j] ->SetOption("colz");
    }
  }

  // Energy plots ----------------------------------
  float E_x2 = 100;
  int E_bins = 50;
  int Eres_bins = 60;
  float Eres_x1   = 2.5; //2.5;
  float Eres_x2   = 62.5; //62.5;
  float res_bins = 300; // 150
  float res_max = 3.0;

  hTrue_NContainShwr    = new TH1D("True_NContainShwr","Number of fully contained showers ;True electron energy [MeV];Number of events",12, Eres_x1, Eres_x2);
  hQ[0]                 = new TH1D("Q_Data","Charge of Michel shower;Reconstructed charge Q [e-];Events",30,0,3e6);
  hQ[1]                 = (TH1D*)hQ[0]->Clone("Q_MC");
  hL[0]                 = new TH1D("L_Data","Photons of Michel shower;Reconstructed light L [#gamma];Events",30,0,3e6);
  hL[1]                 = (TH1D*)hL[0]->Clone("L_MC");
  hTrue_Q               = new TH1D("TrueQ_Data","True deposited charge of Michel shower;True deposited charge Q [e-];Events",30,0,3e6);
  hTrue_L               = new TH1D("TrueL_Data","True photons produced in Michel shower;True produced light L [#gamma];Events",30,0,3e6);
  hQRes                 = new TH1D("QRes","Michel electron shower charge resolution;Q_{reco} / Q_{true}^{dep}",120,-1.2,1.2);
  hLRes                 = new TH1D("LRes","Michel electron shower light resolution;L_{reco} / L_{true}^{dep}",120,-1.2,1.2);
  
  hEnergyTrk[0]         = new TH1D("EnergyTrk_Data","Charge-based Michel Ionization Energy;Reconstructed ionization energy [MeV];Events / 2 MeV",E_bins,0,E_x2);
  hEnergyTrk[1]         = (TH1D*)hEnergyTrk[0]->Clone("EnergyTrk_MC");
  //hEnergyTrkErr         = (TH1D*)hEnergyTrk[0]->Clone("EnergyTrk_MC_Err");
  
  hEnergyQ[0]           = new TH1D("EnergyQ_Data","E_{Q} = ( Q / R ) #times W_{ion} (assume uniform R);Reconstructed shower energy [MeV];Events / 2 MeV",E_bins, 0, E_x2);
  hEnergyQ[1]           = (TH1D*)hEnergyQ[0]     ->Clone("EnergyQ_MC");
  hEnergyL[0]           = new TH1D("EnergyL_Data","E_{L} = L(#alpha + 1 - 1/R) #times W_{ion} (assume uniform R);Reconstructed shower energy [MeV];Events / 2 MeV",E_bins, 0, E_x2);
  hEnergyL[1]           = (TH1D*)hEnergyL[0]     ->Clone("EnergyL_MC");
  hEnergyQ_2R[0]           = new TH1D("EnergyQ_2R_Data","E_{Q} = ( Q_{e} / R_{e} + Q_{i} / R_{i} ) #times W_{ion} (non-uniform R);Reconstructed shower energy [MeV];Events / 2 MeV",E_bins, 0, E_x2);
  hEnergyQ_2R[1]           = (TH1D*)hEnergyQ_2R[0]     ->Clone("EnergyQ_2R_MC");
  hEnergyQL[0]          = new TH1D("EnergyQL_Data","E_{Q+L} = ( Q + L ) #times W_{ph};Reconstructed shower energy [MeV];Events / 2 MeV",E_bins, 0, E_x2);
  hEnergyQL[1]          = (TH1D*)hEnergyQL[0]    ->Clone("EnergyQL_MC");
  hEnergyQL_LogL[0]     = new TH1D("EnergyQL_LogL_Data","E_{Q+L}^{likelihood} (log-likelihood, uniform R);Reconstructed shower energy [MeV];Events / 2 MeV",E_bins, 0, E_x2);
  hEnergyQL_LogL[1]     = (TH1D*)hEnergyQL_LogL[0]    ->Clone("EnergyQL_LogL_MC");
  hEnergyQL_LogL_Fail[0]     = new TH1D("EnergyQL_LogL_Fail_Data","E_{Q+L} for events failing Log-L;Reconstructed shower energy [MeV];Events / 2 MeV",E_bins, 0, E_x2);
  hEnergyQL_LogL_Fail[1]     = (TH1D*)hEnergyQL_LogL[0]    ->Clone("EnergyQL_LogL_Fail_MC");
 
  hEnergy_QvsL[0]       = new TH2D("Energy_QvsL_Data",";Charge-based energy E_{Q} [MeV];Light-based energy E_{L} [MeV]",E_bins, 0, E_x2, E_bins, 0, E_x2);
  hEnergy_QvsL[1]       = (TH2D*)hEnergy_QvsL[0]->Clone("Energy_QvsL_MC");
  
  hEvs_Qtrue                = new TH2D("Evs_Qtrue","True Q vs. true Michel energy dep.;True energy deposited [MeV];True charge [e^{-}]",Eres_bins,Eres_x1,Eres_x2,75,0.,3000e3);
  hEvs_Ltrue                = new TH2D("Evs_Ltrue","True L vs. true Michel energy dep.;True energy deposited [MeV];True light [#gamma]",Eres_bins,Eres_x1,Eres_x2,75,0.,3000e3);
  hEvs_Q                = new TH2D("Evs_Q","Q vs. true Michel energy dep.;True energy deposited [MeV];Reco charge [e^{-}]",Eres_bins,Eres_x1,Eres_x2,75,0.,3000e3);
  hEvs_L                = new TH2D("Evs_L","L vs. true Michel energy dep.;True energy deposited [MeV];Reco light [#gamma]",Eres_bins,Eres_x1,Eres_x2,75,0.,3000e3);
 
  hEvsRes_Q             = new TH2D("EvsRes_Q","Q resolution vs. true Michel energy dep.;True energy deposited [MeV];( reco - true ) / true",Eres_bins,Eres_x1,Eres_x2,res_bins,-1.*res_max,res_max);
//  hEvsRes_Q             = new TH2D("EvsRes_Q","Q resolution vs. true Michel energy dep.;True energy deposited [MeV];( reco - true ) / true",Eres_bins*5,Eres_x1,Eres_x2,res_bins,-1.*res_max,res_max);
  hEvsRes_L             = new TH2D("EvsRes_L","L resolution vs. true Michel energy dep.;True energy deposited [MeV];( reco - true ) / true",Eres_bins,Eres_x1,Eres_x2,res_bins,-res_max,res_max);
  hEvsRes_E_Trk           = new TH2D("EvsRes_E_Trk","Electron ion. track;True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
  hEvsRes_E_Q             = new TH2D("EvsRes_E_Q","Q-only shower energy (uniform R);True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
  hEvsRes_E_Q_2R          = new TH2D("EvsRes_E_Q_2R","Q-only shower energy (non-uniform R);True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
  hEvsRes_E_QL            = new TH2D("EvsRes_E_QL","Q+L shower energy;True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
  hEvsRes_E_QL_LogL       = new TH2D("EvsRes_E_QL_LogL","Q+L log-likelihood shower energy (uniform R);True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 

  
  hEnergy_QvsQL[0]      = new TH2D("Energy_QvsQL_Data",";Q-only Energy [MeV];Q+L Energy [MeV]",E_bins,0.,E_x2,E_bins,0.,E_x2);
    hEnergy_QvsQL[0]    ->SetOption("colz");
 
  //h_Qcol_vs_QcolRes   = new TH2D("Qcol_vs_QcolRes","Collected charge resolution;True charge at wires, Q^{col}_{true} [e-]; ( Q^{col}_{reco} - Q^{col}_{true} ) / Q^{col}_{true}",50,0.,5000e3,200,-0.5,0.5);
  //h_Pe_vs_PeRes       = new TH2D("Pe_vs_PeRes","Photoelectron resolution;True photoelectrons from PMT, S_{true} [pe]; ( S_{reco} - S_{true} ) / Q_{true}",60,0.,1200.,200,-0.5,0.5);
  hEnergyVsResPE        = new TH2D("EnergyVsResPE","PE resolution vs. true Michel energy dep.;True energy deposited [MeV];PE resolution",10,0.,50.,120,-1.2,1.2);
  
  //hEnergyDepResQ        = new TH1D("EnergyDepResQ","Charge-based Michel Shower Energy Resolution;(E_{reco} - E_{true}^{dep}) / E_{true}^{dep};Events",res_bins, -1.*res_max, res_max);
  //hEnergyDepResQL       = new TH1D("EnergyDepResQL","Combined Charge+Light Michel Shower Energy Resolution;(E_{reco} - E_{true}^{dep}) / E_{true}^{dep};Events",res_bins, -1.*res_max, res_max);
  hEnergy_QvsQL[1]      = (TH2D*)hEnergy_QvsQL[0] ->Clone("Energy_QvsQL_MC");
    hEnergy_QvsQL[1]    ->SetOption("colz");
  hPERes                  = new TH1D("PERes","PE resolution;(reco-true)/true",120,-1.2,1.2);
  hTrue_EnergyDepTrk      = new TH1D("True_EnergyDepTrk","True energy deposited by Michel electron ionization track;True energy deposited (E_{true}^{dep}) [MeV]",E_bins,0,E_x2);
  hTrue_EnergyDep         = new TH1D("True_EnergyDep","True energy deposited by Michel electron shower;True energy deposited (E_{true}^{dep}) [MeV]",E_bins,0,E_x2);
  hTrue_EnergyDep_PreCuts         = new TH1D("True_EnergyDep_PreCuts","True energy deposited by Michel electron shower;True energy deposited (E_{true}^{dep}) [MeV]",E_bins, 0, E_x2);
  hTrue_Energy         = new TH1D("True_Energy","True initial energy of Michel electron shower;Energy [MeV]",E_bins, 0, E_x2);
  hTrue_EnergyVsContainment = new TH2D("True_EnergyVsContainment",";Initial energy #it{E} [MeV];Fraction of energy deposited in LAr",40,0,80,120,0,1.2);
  hTrue_ZX_Vs_Containment = new TH2D("True_ZX_Vs_Containment",";Z-coordinate of #mu endpoint [cm];X-coordinate of #mu endpoint [cm];Containment fraction",45,0.,90.,24,0.,47.5);
  hTrue_ZX_Vs_Energy      = (TH2D*)hTrue_ZX_Vs_Containment->Clone("True_ZX_Vs_Energy");
  hTrue_Energy_vs_RecoEnergy    = new TH2D("True_Energy_vs_RecoEnergy",";Michel electron energy [MeV];Reconstructed energy [MeV]",120,0.,60.,120,0.,60.);
    hTrue_Energy_vs_RecoEnergy->SetOption("colz");
  hTrue_EnergyDep_vs_RecoEnergy = new TH2D("True_EnergyDep_vs_RecoEnergy",";Energy deposited [MeV];Reconstructed energy [MeV]",120,0.,60.,120,0.,60.);
    hTrue_EnergyDep_vs_RecoEnergy->SetOption("colz");
  hTrue_Dist_vs_Containment = new TH2D("True_Dist_vs_Containment",";Projected distance to TPC boundary [cm];Fractional shower containment",100,0.,100.,101,0,1.01);
    hTrue_Dist_vs_Containment->SetOption("colz");
  hTrue_EnergyVsEnergyDepTrk  = new TH2D("True_EnergyVsEnergyDepTrk","Fully-contained electrons;Initial electron energy [MeV];Energy deposited in ionization track [MeV]",
    120,0,60.,
    120,0,60.);
    hTrue_EnergyVsEnergyDepTrk->SetOption("colz");
  
  // Electron dE/dx plots ---------------------------
  hdEdx[0]                  = new TH1D("dEdx_Data","Michel Electron dE/dx;dE/dx [MeV/cm];Events / 0.1 MeV/cm",50,0,5); 
  hdEdx[1]                  = (TH1D*)hdEdx[0]->Clone("dEdx_MC"); 
 
  // Recombination ----------------------------------
  hR[0]                   = new TH1D("R_Data",";Measured survival fraction (#it{R});Events",30,0.4,1.0);
  hR[1]                   = (TH1D*)hR[0]->Clone("R_MC");
  hTrue_R                 = (TH1D*)hR[0]->Clone("R_MCTrue");
  hTrue_RTrack                 = (TH1D*)hR[0]->Clone("RTrack_MCTrue");
  hTrue_RShower                 = (TH1D*)hR[0]->Clone("RShower_MCTrue");
  //hBeta[0]                = new TH1D("Beta_Data",";Q/L;Events",50,0.,5.);
  //hBeta[1]                = (TH1D*)hBeta[0]->Clone("Beta_MC");
//  hTrue_Beta                 = (TH1D*)hBeta[0]->Clone("Beta_MCTrue");
  
  // Decay time plots -------------------------------
  int numdTbins=(fDecayFit_MaxDecayTime-fDecayFit_MinDecayTime)/100;
  hDecayTime[0]       = new TH1D("DecayTime_Data", ";#DeltaT [ns]; Events / 100 ns", numdTbins,fDecayFit_MinDecayTime,fDecayFit_MaxDecayTime);
  hDecayTime[1]       = (TH1D*)hDecayTime[0]->Clone("DecayTime_MC");

  // =================================================================================
  // Trigger efficiency equation
  // =================================================================================
  trigEffCut[0] = new TF1("trigEffCut0","[3]*(1.-[2]/(1.+(x/[0])^[1] ))",0.,2000.);
  trigEffCut[1] = (TF1*)trigEffCut[0]->Clone("trigEffCut1");
  SetTrigEffParams();
  
  hMCSpatialWgt_ZX->Reset();
  hMCSpatialWgt->Reset();

  // ==================================================================================
  // Event cut histograms
  h_evtCut_MuPE_prompt[0][0]= new TH1D("evtCut_MuPE_prompt_HMM",";Muon prompt light in HMM PMT [pe]",60,0.,1200.);
  h_evtCut_MuPE_prompt[1][0]= (TH1D*)h_evtCut_MuPE_prompt[0][0]->Clone("evtCut_MuPE_prompt_HMM_MC");
  h_evtCut_MuPE_prompt[0][1]= new TH1D("evtCut_MuPE_prompt_ETL",";Muon prompt light in ETL PMT [pe]",50,0.,400.);
  h_evtCut_MuPE_prompt[1][1]= (TH1D*)h_evtCut_MuPE_prompt[0][1]->Clone("evtCut_MuPE_prompt_ETL_MC");
  h_evtCut_MuClusterSize[0] = new TH1D("evtCut_MuClusterSize",";Muon cluster size [hits]",80,0,80);
  h_evtCut_MuClusterSize[1] = (TH1D*)h_evtCut_MuClusterSize[0]->Clone("evtCut_MuClusterSize_MC");
  h_evtCut_ElClusterSize[0] = new TH1D("evtCut_ElClusterSize",";Electron cluster size [hits]",60,0,60);
  h_evtCut_ElClusterSize[1] = (TH1D*)h_evtCut_ElClusterSize[0]->Clone("evtCut_ElClusterSize_MC");
  h_evtCut_BraggSlope[0] = new TH1D("evtCut_BraggSlope",";Muon Bragg slope [ADC/cm]",50,-200,1800);
  h_evtCut_BraggSlope[1] = (TH1D*)h_evtCut_BraggSlope[0]->Clone("evtCut_BraggSlope_MC");
  h_evtCut_MuLinearity[0] = new TH1D("evtCut_MuLinearity",";Ave muon cluster linearity",50,0,1.02);
  h_evtCut_MuLinearity[1] = (TH1D*)h_evtCut_MuLinearity[0]->Clone("evtCut_MuLinearity_MC");
  h_evtCut_MuDir[0] = new TH1D("evtCut_MuDir",";Number muon hits for 2D dir. fit",25,0,25);
  h_evtCut_MuDir[1] = (TH1D*)h_evtCut_MuDir[0]->Clone("evtCut_MuDir_MC");
  h_evtCut_DecayAngle[0] = new TH1D("evtCut_DecayAngle",";2D decay angle [deg]",36,0,180);
  h_evtCut_DecayAngle[1] = (TH1D*)h_evtCut_DecayAngle[0]->Clone("evtCut_DecayAngle_MC");
  h_evtCut_ShowerSize[0] = new TH1D("evtCut_ShowerSize",";Electron shower size [hits]",60,0,60);
  h_evtCut_ShowerSize[1] = (TH1D*)h_evtCut_ShowerSize[0]->Clone("evtCut_ShowerSize_MC");
  h_evtCut_NumPts3D[0] = new TH1D("evtCut_NumPts3D",";Number 3D points formed",40,0,40);
  h_evtCut_NumPts3D[1] = (TH1D*)h_evtCut_NumPts3D[0]->Clone("evtCut_NumPts3D_MC");
  h_evtCut_FracPts3D[0] = new TH1D("evtCut_FracPts3D",";Fraction hits made into 3D points",21,0,1.05);
  h_evtCut_FracPts3D[1] = (TH1D*)h_evtCut_FracPts3D[0]->Clone("evtCut_FracPts3D_MC");
  h_evtCut_ProjDist[0] = new TH1D("evtCut_ProjDist",";Projected distance to wireplanes [cm]",40,0.,200.);
  h_evtCut_ProjDist[1] = (TH1D*)h_evtCut_ProjDist[0]->Clone("evtCut_ProjDist_MC");


  // Preliminary loop over data and MC trees
  //fMaxMCEvts = -200000; //-1.;
  std::cout<<"Looping over data tree...\n";
  Loop(0);
  std::cout<<"Looping over MC tree...\n";
  Loop(1);
  //EnergyPlots();
  
  //if( fApplyCalibCorr ) {
  //  MakeCalibrationPlots();
  //  Loop(0);
  //}
  
  // ====================================================
  // Set random seed (for smearing)
  fRand->SetSeed(88);
  generator.seed(88);

}

//###########################################################################
void Print(){
    std::cout
    <<"===================================================\n"
    <<"Data tree                : "<<fTree[0]->GetEntries()<<" entries\n"
    <<"MC tree                  : "<<fTree[1]->GetEntries()<<" entries\n"
    <<"MC tree (no smear)       : "<<fTreeMC_ns->GetEntries()<<" entries\n\n";
  for(size_t i=0; i<2; i++){
    float chi2w_p = GetChi2Weighted(hPE_prompt[0][i],hPE_prompt[1][i]);
    float chi2w_t = GetChi2Weighted(hPE_total[0][i],hPE_total[1][i]);
    float chi2_p  = GetChi2(hPE_prompt[0][i],hPE_prompt[1][i]);
    float chi2_t  = GetChi2(hPE_total[0][i],hPE_total[1][i]);
    std::cout
    <<"---------------------------------------------------\n"
    <<"CH"<<i<<"\n"
    <<"  - Scale factor         : "<<fScaleFactor[i]<<"\n"
    <<"  - Trig eff P           : "<<fTrigEff_P[i]<<"\n"
    <<"  - Trig eff K           : "<<fTrigEff_K[i]<<"\n"
    <<"  - Smear factor         : "<<fSmearFactor[i]<<"\n"
    <<"  - Prompt chi2_w        : "<<chi2w_p<<"\n"
    <<"  - Total chi2_w         : "<<chi2w_t<<"\n"
    <<"  - Prompt chi2          : "<<chi2_p<<"\n"
    <<"  - Total chi2           : "<<chi2_t<<"\n";
  }
    std::cout
    <<"===================================================\n"
    <<"  P fixed?       : ["<<fFixP[0]<<", "<<fFixP[1]<<"]\n"
    <<"  K fixed?       : ["<<fFixK[0]<<", "<<fFixK[1]<<"]\n"
    <<"  scale fixed?   : "<<fFixScale<<"\n"
    <<"  smear fixed?   : "<<fFixSmearTotal<<"\n";
}


//################################################################################
// LOOP
//  This function loops over both the MC and data TTrees and makes all hsitograms
//  and event counts.  Has option to also apply smearing to light data.

void Loop(int isMC=0){ Loop(fTree[isMC], bool(isMC), false );}
void Loop(TTree* tree, bool isMC, bool doSmearing ) {
  
  fRand->SetSeed(89);
  
   
  // ===================================================
  // First off, is this data or MC?
  //  type = 0 --> data
  //  type = 1 --> MC 
  int type= int(isMC);
 
  
  // =====================================================
  // Reset the histograms

  hR[type]->Reset();
  hDecayTime[type]->Reset();
  hEnergyTrk[type]->Reset();
  hEnergyQ[type]->Reset();
  hEnergyL[type]->Reset();
  hEnergyQ_2R[type]->Reset();
  hEnergyQL[type]->Reset();
  hEnergyQL_LogL[type]->Reset();
  hEnergyQL_LogL_Fail[type]->Reset();
  hEnergy_QvsL[type]->Reset();
  hEnergy_QvsQL[type]->Reset();
  hL[type]      ->Reset();
  hQ[type]      ->Reset();
  hPE[type]     ->Reset();
  hLightYield[type]->Reset();
  hdEdx[type]   ->Reset();
  hPromptFrac[type]->Reset();
  hPromptFracRaw[type]->Reset();
  hLY_ZX[type]->Reset();
  hLY_entries_ZX[type]->Reset();
  hLY_entries_ZY[type]->Reset();
  hE_ZX[type]->Reset();
  hLY_ZY[type]->Reset();
  hE_ZY[type]->Reset();
  hElShowerCentroid_X[type]->Reset();
  hElShowerCentroid_Y[type]->Reset();
  hElShowerCentroid_Z[type]->Reset();
  hElShowerCentroid[type]->Reset();
  hElShowerVis[type]->Reset();
  h_evtCut_MuPE_prompt[type][0]->Reset();
  h_evtCut_MuPE_prompt[type][1]->Reset();
  h_evtCut_MuClusterSize[type]->Reset();
  h_evtCut_ElClusterSize[type]->Reset();
//  h_evtCut_BraggSlope[type]->Reset();
  h_evtCut_MuLinearity[type]->Reset();
  h_evtCut_MuDir[type]->Reset();
  h_evtCut_DecayAngle[type]->Reset();
  h_evtCut_ShowerSize[type]->Reset();
  h_evtCut_NumPts3D[type]->Reset();
  h_evtCut_FracPts3D[type]->Reset();
  //h_evtCut_ProjDist[type]->Reset();
  if( !isMC ) {
    hElectronLifetime->Reset();
    hRunVsElectronLifetime->Reset();
    hRunVsEnergy->Reset();
    hRunVsLight->Reset();
    hTimeVsElectronLifetime->Reset();
    hTimeVsLight->Reset();
    hTimeVsCharge->Reset();
    hTimeVsChargePerCm->Reset();
    hTimeVsLightPerCm->Reset();
    hChargePerCm->Reset();
    hLightPerCm->Reset();
    hPE_total_compare->Reset();
  }
  if( isMC ) {
    hQRes->Reset();
    hLRes->Reset();
//    h_Qcol_vs_QcolRes;
//    h_Pe_vs_PeRes;
    hPERes->Reset();
    //hEnergyDepResQ->Reset();
    //hEnergyDepResQL->Reset();
    //hEnergyResQ->Reset();
    //hEnergyResQL->Reset();

    hEvs_Qtrue->Reset();
    hEvs_Q->Reset();    
    hEvs_L->Reset();    
    hEvs_Ltrue->Reset();    
    hEvsRes_E_Trk->Reset();
    hEvsRes_E_Q->Reset();
    hEvsRes_E_Q_2R->Reset();
    hEvsRes_E_QL->Reset();
    hEvsRes_E_QL_LogL->Reset();
    hEvsRes_Q->Reset();
    hEvsRes_L->Reset();

    //hEnergyVsResQ->Reset();
    hEnergyVsResPE->Reset();
    //hEnergyVsResL->Reset();
    hTrue_Q->Reset();
    hTrue_L->Reset();
    hTrue_EnergyVsEnergyDepTrk->Reset();
    hTrue_EnergyDepTrk->Reset();
    hTrue_EnergyDep->Reset();
    hTrue_R->Reset();
    hTrue_RTrack->Reset();
    hTrue_RShower->Reset();
    hTrue_ZX_Vs_Containment->Reset();
    hTrue_ZX_Vs_Energy->Reset();
    hTrue_Energy_vs_RecoEnergy->Reset();
    hTrue_EnergyDep_vs_RecoEnergy->Reset();
    hTrue_Dist_vs_Containment->Reset();
    hTrue_LY->Reset();
    hTrue_NContainShwr->Reset();
  }
  for(int ch=0; ch<2; ch++){
    hPE_prompt[type][ch]->Reset();
    hPE_total[type][ch]->Reset();
    hPE_total_qc[type][ch]->Reset();
    hPE_prompt_pretrigMC[ch]->Reset();
    hPE_prompt_posttrigMC[ch]->Reset();
    hPE_prompt_nosmearMC[ch]->Reset();
    hPE_total_nosmearMC[ch]->Reset();
  } 
  
  // ==================================================
  // Set trig efficiency parameters 
  SetTrigEffParams();

  // ====================================================
  // begin loop
  int kMax = (int)tree->GetEntries();
  if( isMC && fMaxMCEvts > 0 && kMax > fMaxMCEvts ) kMax = fMaxMCEvts;
  for(int i=0; i<kMax; i++){
  
    if( (i % 10000) == 0 ) std::cout<<"...tree entry "<<i<<"\n";
    
    //............
    //if( i != 3800 ) continue;
    //if( i > 3800 ) break;
   
    tree->GetEntry(i);
    
    if( !isMC && ((fMaxRun > 0 && fRunNumber > fMaxRun) || (fMinRun > 0 && fRunNumber < fMinRun ))) continue;

    if( fChargeSignFilter != 0 && fTrue_ChargeSign != fChargeSignFilter ) continue;

    float fWgt = 1.0;
//    if( isMC ) fWgt *= GetMCSpatialWgt( fTrue_MuTrackEnd_X, fTrue_MuTrackEnd_Y, fTrue_MuTrackEnd_Z );
   
    // convert time into days
    float day = float(fEventTime - fTimeStart)/86400.;
   
    if( !isMC ) {
      hRunVsElectronLifetime->Fill(fRunNumber, fElectronLifetime ); 
      hTimeVsElectronLifetime->Fill(day, fElectronLifetime);
      hElectronLifetime->Fill( fElectronLifetime );
    }
   
    
    fTriggered = true;
   
   
    // --------------------------------------------------------
    // The scaling, smearing, and trigger efficiency parameters
    // will be used to adjust visibility and PEs...
    float vPE_prompt_nosmear[2]={-999.};
    float vPE_total_nosmear[2]={-999.};
    
    // Reset shower visibility
    fElShowerVis = 0.;
    
    // ---------------------------------------------
    // Loop over the PMTs
    int   N = 0;
    for(size_t ch=0; ch<2; ch++){
        
      if( (*fQE_ScaleFactor)[ch] <= 0. ) continue;
      N++;
    
      // Revise the shower visibility for this channel
      float scaleFac = fScaleFactor[ch]/(*fQE_ScaleFactor)[ch];
      fTrue_ElShowerVisCh[ch]*= scaleFac;
      fElShowerVisCh[ch]     *= scaleFac;
      fElShowerVis           += fElShowerVisCh[ch];
   
      // ------------------------------------------
      // If this is MC, then scale and smear the light
      if( isMC ) {
        
        // Scale PEs
        float scaleFacPrompt = scaleFac;
        //if( fFastRatio != fFastRatioNom ) scaleFacPrompt *= fFastRatio / fFastRatioNom;
        fMuPE_prompt[ch] *= scaleFacPrompt;
        fPE_prompt[ch] *= scaleFacPrompt;
        fPE_total[ch] *= scaleFac;
        fPE_total_qc[ch] *= scaleFac;
        fTrue_PE_prompt[ch] *= scaleFacPrompt;
        fTrue_PE_total[ch] *= scaleFac;

        // ---------------------------------------------
        if( fSmearTruePE ) {
          
          // From the true L, "fake" the PE based on visibility and quenching.
          // This is the preferred method when scaling the LY up and down, as it
          // properly incorporates changes in counting statistics.

          // TO DO: add quenching
         
          /* 
          std::binomial_distribution<int> distPrompt( fTrue_ElShowerPhotonsPrompt, fTrue_ElShowerVisCh[ch] );
          std::binomial_distribution<int> distTotal( fTrue_ElShowerPhotons, fTrue_ElShowerVisCh[ch] );
          */
          std::binomial_distribution<int> distPrompt( fTrue_ElShowerPhotonsPromptQuenched, fTrue_ElShowerVisCh[ch] );
          std::binomial_distribution<int> distTotal( fTrue_ElShowerPhotonsQuenched, fTrue_ElShowerVisCh[ch] );
          fTrue_PE_prompt[ch] = distPrompt(generator);
          fTrue_PE_total[ch] = distTotal(generator);
          
          fPE_prompt[ch] = fTrue_PE_prompt[ch];
          fPE_total[ch] = fTrue_PE_total[ch];
          
        }
      
        if( doSmearing ) {
    
          // Save unsmeared PE
          vPE_prompt_nosmear[ch] = fPE_prompt[ch];
          vPE_total_nosmear[ch]  = fPE_total[ch];

          // Do smearing on total light integral
            //if( fSmearFactor[ch] >= 0 )
            //  fPE_total[ch]     += fRand->Gaus(0.,fSmearFactor[ch]*std::sqrt((double)fPE_total[ch]));
          if( fSmearFactor[ch] >= 0 ){
            if( fSmearMode == 0 ) fPE_total[ch]     += fRand->Gaus(0.,fSmearFactor[ch]*fPE_total[ch]);
            else                  fPE_total[ch]     += fRand->Gaus(0.,fSmearFactor[ch]*std::sqrt(fPE_total[ch]));
          }

          // Quenching correction
          fPE_total_qc[ch] = fPE_total[ch];
          if( fCorrectQuenching && fPE_total[ch] > 0 ) fPE_total_qc[ch]  = CorrectForQuenching( fPE_total[ch], fPE_prompt[ch], (*fMuContamCorr_EffTau)[1] );
    
          // Apply trigger cut
          if( fRand->Rndm() > trigEffCut[ch]->Eval(fPE_prompt[ch]) ) fTriggered = false;
    
          // Calc prompt fraction
          fPromptFrac[ch] = fPE_prompt[ch] / fPE_total[ch];
      
        }//<-- end do smearing
      }//<-- end isMC
    }//<-- end loop over PMTs
    if( N == 0 ) fTriggered = false;


    //if( fPE_prompt[1] < 0 || fPE_total[1] < 0 ) continue;

    // --------------------------------------------
    // Smear charge
    if( isMC && fChargeSmearFactor > 0. ){
      float q0 = fElShowerCharge;
      fElShowerCharge += fRand->Gaus(0., fElShowerCharge*fChargeSmearFactor);
      fElShowerEnergy *= (fElShowerCharge/q0);
    }
 
    // ---------------------------------------------
    // Some truth level plots about the MC sample (before trig cut)
    if( isMC ) {
      hTrue_Energy            ->Fill(fTrue_ElEnergy,fWgt);
      hTrue_EnergyDep_PreCuts ->Fill(fTrue_ElShowerEnergyDep,fWgt);
      hTrue_EnergyVsContainment->Fill(fTrue_ElEnergy, fTrue_ElShowerEnergyDep/fTrue_ElEnergy,fWgt);
      hTrue_ZX_Vs_Energy        ->Fill(fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_X, fTrue_ElEnergy);
      hTrue_ZX_Vs_Containment   ->Fill(fTrue_MuTrackEnd_Z, fTrue_MuTrackEnd_X, fTrue_ElShowerEnergyDep);
      if( fElShowerEnergy > 0. ) {
        hTrue_Energy_vs_RecoEnergy    ->Fill( fTrue_ElEnergy, fElShowerEnergy);
        hTrue_EnergyDep_vs_RecoEnergy ->Fill( fTrue_ElShowerEnergyDep, fElShowerEnergy);
      }
      if( fTrue_IsElTrkContained ) hTrue_EnergyVsEnergyDepTrk->Fill(fTrue_ElEnergy, fTrue_ElTrackEnergyDep);
        
      /* 
      // Determine projected distance of shower to nearest wall
      TVector3 muEnd(fTrue_MuTrackEnd_X, fTrue_MuTrackEnd_Y, fTrue_MuTrackEnd_Z);
      TVector3 elShowerCen( fTrue_ElShowerCentroid_X, fTrue_ElShowerCentroid_Y, fTrue_ElShowerCentroid_Z );
      TVector3 elShowerDir = (elShowerCen - muEnd);
      elShowerDir.SetMag(1.);
          // calculate distance to each of the TPC's boundary planes
          // and keep track of which one is closet
          float minDist = 9999.;
          float dist;
          dist = ( 0. - muEnd.X() ) / elShowerDir.X();
          if( dist > 0. && dist < minDist ) minDist = dist;
          dist = ( 47.5 - muEnd.X() ) / elShowerDir.X();
          if( dist > 0. && dist < minDist ) minDist = dist;
          dist = ( -20. - muEnd.Y() ) / elShowerDir.Y();
          if( dist > 0. && dist < minDist ) minDist = dist;
          dist = ( 20. - muEnd.Y() ) / elShowerDir.Y();
          if( dist > 0. && dist < minDist ) minDist = dist;
          dist = ( 0. - muEnd.Z() ) / elShowerDir.Z();
          if( dist > 0. && dist < minDist ) minDist = dist;
          dist = ( 90. - muEnd.Z() ) / elShowerDir.Z();
          if( dist > 0. && dist < minDist ) minDist = dist;
            hTrue_Dist_vs_Containment->Fill( minDist, fTrue_ElShowerEnergyDep / fTrue_ElEnergy );
          */
    }
    
    // ---------------------------------------------
    // Fill plots for pre- and post-trig cut comparisons.
    // If event failed the trigger cut, skip this event
    if( fElShowerEnergy > 0 && fDecayTime > fdTcut ) {for(int ch=0; ch<2; ch++) hPE_prompt_pretrigMC[ch]->Fill(fPE_prompt[ch]);}
    if( !fTriggered ) continue;
    if( fElShowerEnergy > 0 && fDecayTime > fdTcut ) {for(int ch=0; ch<2; ch++) hPE_prompt_posttrigMC[ch]->Fill(fPE_prompt[ch]);}
 

    // -----------------------------------------------
    // Add up total PE, and evaluate prompt PE cut, saturation cut
    bool ampCut = true;
    bool saturationCut = true;
    bool nhitsCut = true;
    bool maxPeCut = true;
    bool muPromptCut = true;
    fTrue_ElShowerPhel = 0;
    fElShowerPhel = 0;
    fElShowerPhel_prompt = 0;
    fElShowerPhel_qc = 0;
    float fMuPE_prompt_comb=0;
    for(int ch=0; ch<2; ch++){
      if( fPE_prompt[ch] <= 0 ) continue;
      if( fUsePmt[ch] ) {
        fElShowerPhel  += fPE_total[ch];
        fElShowerPhel_qc += fPE_total_qc[ch];
        fElShowerPhel_prompt += fPE_prompt[ch];
        fTrue_ElShowerPhel += fTrue_PE_total[ch];
      }
      fMuPE_prompt_comb += fMuPE_prompt[ch];
      if( fMaxPE[ch] > 0. && fPE_total[ch] > fMaxPE[ch] ) maxPeCut=false;
      if( fMuPE_prompt[ch] < fMinMuPromptPE[ch] ) muPromptCut = false;
      if( fElPulseSaturated[ch] ) saturationCut = false;
      if( fAmplitude[ch] < fAmpCut[ch] ) ampCut = false;
//      if( fNumOpHits[ch] != fNumOpHits0[ch] ) nhitsCut = false;
    }

    if( isMC ) {
      if( fTrue_ElShowerEnergyDep > 0 ) hTrue_LY  ->Fill( fTrue_ElShowerPhel / fTrue_ElShowerEnergyDep, fWgt );
      if( fTrue_IsElShwrContained ) hTrue_NContainShwr->Fill( fTrue_ElEnergy );
    }


    // Halt processing and move to next event if containment
    // criteria are not met:
    if( !fIsRealData && (fRequireContainment && !fTrue_IsElShwrContained) ) continue;
   
    // ===================================================
    // Begin event cut evaluation
  
    bool shower2D           = false; 
    bool shower2D_Pl0       = false;
    bool goodShower2D       = false;
    bool goodShower2D_Pl0   = false;
    bool goodShower3D       = false;
    bool goodShower3D_dTcut = false;
    
    //nEvt_TotalEvents++;


    // require "good"  shower on induction
    if( fElShowerSize_Pl0 > 0 ) {
      shower2D_Pl0      = true;
      goodShower2D_Pl0  = true;
    }
    if( fElShowerSize > 0 ) shower2D = true;

 
    // Michel op ID
    if( fMichelOpticalID && fDecayTime > fGateDelay ) { 
      nEvt_OpID_timeRange++; 
    
      // Muon prompt light cut
      if( muPromptCut ) {  
        nEvt_OpID_muPromptCut++;
      
        // cluster boundary found
        if( fMuClusterSize > 0 && (fMaxCovAtBnd < 0 || fCovAtBnd <= fMaxCovAtBnd) ) {
          nEvt_Shwr2D_BndFound++;
            
          h_evtCut_MuClusterSize[type]->Fill(fMuClusterSize);
          h_evtCut_ElClusterSize[type]->Fill(fElClusterSize);
          
          // cluster size cuts
          if( fElClusterSize >= fMinElClusterSize &&
              fElClusterSize <= fMaxElClusterSize &&
              fMuClusterSize >= fMinMuClusterSize ) {
            nEvt_Shwr2D_ClusterSize++;
            
                              
            // Bragg peak slope
            if( fMinBraggSlope < 0 || fBraggSlope > fMinBraggSlope ) {
              nEvt_Shwr2D_Slope++;
            
              h_evtCut_MuLinearity[type]->Fill(fMuAveLinearity);

              // muon linearity
              if( fMuAveLinearity > fMinMuLinearity &&
                  fFracMuHitsLinear > fMinFracMuHitsLinear ) {
                nEvt_Shwr2D_MuLinearity++;
              
                h_evtCut_MuDir[type]->Fill( fMuClusterHitsEndFit );

                // muon terminal direction fit
                if( fMuClusterHitsEndFit >= fMinMuClusterHitsEndFit ) {
                  nEvt_Shwr2D_MuDir++;
                
                  h_evtCut_DecayAngle[type]->Fill( fDecayAngle2D );

                  // decay angle cut (req on *either* plane)
                  if( (fDecayAngle2D > fMinDecayAngle2D && 
                       fDecayAngle2D < fMaxDecayAngle2D)
                      ) {
                    nEvt_Shwr2D_DecayAngle++;
                    
                    if( fElShowerSize >= fElClusterSize ) 
                    h_evtCut_ShowerSize[type]->Fill( fElShowerSize );
                  
                    // shower size cuts
                    if( fElShowerSize >= fMinElShowerSize && fElShowerSize <= fMaxElShowerSize ) {
                      nEvt_Shwr2D_ShowerSize++;
                   
                      // shower completeness frac (skip if the total 
                      if( fElShowerFrac > fMinShowerFrac || fElShowerSize*(1./fElShowerFrac-1.) < 10 ) {
                        nEvt_Shwr2D_ShowerFrac++; 
                        
                        goodShower2D = true;
                      
                      // ---------------------------------------
                      // 3D shower boolean  
                      if( !fRequire3DMuEndpt || ( fMuEnd3D_X > 0 ) ) {
                        if( fMuEnd3D_X > 0 ) nEvt_3DMuEnd++;
                      
                        // shower on pl0 (boolean set earlier) 
                        if( shower2D_Pl0 ) {
                          nEvt_Shwr3D_ShowersOnBothPlanes++;
                              
                            h_evtCut_NumPts3D[type]->Fill( fNumPts3D);
                        
                            // num 3D pts
                            if( fNumPts3D >= fMinNumPts3D && fElShowerVis > 0 ) {
                              nEvt_Shwr3D_Npts++;
                               
                              fFracHits3D = float(fNumPts3D)/fElShowerSize; 
                              h_evtCut_FracPts3D[type]->Fill( fFracHits3D );

                              // fraction of hits used
                              if( fFracHits3D > fMinFracHits3D ) {
                                
                                nEvt_Shwr3D_FracHits3D++;
          
                                hElShowerCentroid_X[type]->Fill( fElShowerCentroid_X,fWgt );
                                hElShowerCentroid_Y[type]->Fill( fElShowerCentroid_Y,fWgt );
                                hElShowerCentroid_Z[type]->Fill( fElShowerCentroid_Z,fWgt );
                                hElShowerCentroid[type]->Fill( 
                                  fElShowerCentroid_X, 
                                  fElShowerCentroid_Y, 
                                  fElShowerCentroid_Z, fWgt);

                                // Shower centroid cut
                                if( fElShowerCentroid_X > fFidMarginX &&
                                    fElShowerCentroid_X < 47.5-fFidMarginX &&
                                    fElShowerCentroid_Y > -20.+fFidMarginY && 
                                    fElShowerCentroid_Y < 20.-fFidMarginY &&
                                    fElShowerCentroid_Z > fFidMarginZ && 
                                    fElShowerCentroid_Z < 90.-fFidMarginZ ) {
                                  nEvt_Shwr3D_Centroid++;
                            
                                  //h_evtCut_ProjDist[type]->Fill( fProjDist3DToWires );

                                
                                  // Electron direction cuts (to eliminate events
                                  // where a large portion of charge might be deposited
                                  // on the other side of the wireplanes, or where
                                  // charge may intersect and interfere with PMTs)
                                  if( fProjDist3DToWires < 0 || fProjDist3DToWires > fMinProjDist ) {
                                 
                                    nEvt_Shwr3D_Direction++; 
                                    goodShower3D=true;

                                    // decay time cut
                                    if( fDecayTime > fdTcut ) {
                                      nEvt_dTcut++;
                                      goodShower3D_dTcut=true;
                                      h_evtCut_MuPE_prompt[type][0]->Fill( fMuPE_prompt[0] );
                                      h_evtCut_MuPE_prompt[type][1]->Fill( fMuPE_prompt[1] );
                                    }//<-- dTcut

                                  }// shower dir. cut
                                }/// centroid fid cut
                              }//<-- frac 3D pts
                            }//<-- num 3D pts
                          }//<-- showers on both planes
                        }//<-- (3D muon endpoint -- optional)
                        // -----------------------------------------
                      }//<-- shower frac 
                    }//<-- shower size
                  }//<-- decay angle cut
                }//<-- mu dir
              }//<-- mu linearity
            }//<-- Bragg peak slope
          }//<-- clstr size
        }//<-- clstr bnd found
      }//<-- opID mu prompt cut
    }//<-- opID time range


    
    // =====================================================
    // Calibration using crossing muons
    /*
    if( fCrsMuLength > 10 
        && fCrsMuCharge > 0 
        && fCrsMuPhotons > 0 ) {
      hChargePerCm->Fill( fCrsMuCharge / fCrsMuLength );
      hLightPerCm->Fill( fCrsMuPhotons / fCrsMuLength );
      hTimeVsChargePerCm->Fill( day, fCrsMuCharge / fCrsMuLength );
      hTimeVsLightPerCm->Fill( day, fCrsMuPhotons / fCrsMuLength );
    }
    */

    // ======================================================
    // Case of bare electron showers
    if( fTrue_ElShowerEnergyDep > 0 && fTrue_MuTrackEnd_X < -9. && fTrue_MuTrackEnd_Y < -29. ){
      fTrue_IsBareElShower = true;  
      //h_vis_vs_energy->Fill(fElShowerVis,fElShowerEnergy);
    }
    
    if( fTrue_IsBareElShower && fElShowerSize >= 1 && fNumPts3D >= 1 ) {
      goodShower2D = true;
      goodShower3D = true;
      goodShower3D_dTcut = true;
      fDecayTime = fdTcut+10;
      fMuEnd3D_X = 99.;
    }


    // =====================================================
    // Light and calorimetry plots

    // ======================================================
    // Light diagnostic plots
    //  - 2D shower
    //  - dT cut
    if( goodShower2D &&  fDecayTime >= fdTcut ) { //&&  fMuEnd3D_X > 0.) {
    //if( goodShower2D &&  fDecayTime >= fdTcut &&  fMuEnd3D_X > 0.) {
      
      // require the 2D direction of the shower be going backward
      // (away from wireplanes)
//      if( fElDir2D_X > 0. ) {

        if( !fPE_Require3DShower || (fPE_Require3DShower && goodShower3D) ) {
          if( !isMC ) hPE_total_compare->Fill( fPE_total[0], fPE_total[1] );
          for(int ch=0; ch<2; ch++){
            hPE_prompt[type][ch] ->Fill( fPE_prompt[ch], fWgt);
            hPE_total[type][ch]  ->Fill( fPE_total[ch], fWgt);
            hPE_total_qc[type][ch]  ->Fill( fPE_total_qc[ch], fWgt);
            hPE_prompt_nosmearMC[ch] ->Fill( vPE_prompt_nosmear[ch], fWgt);
            hPE_total_nosmearMC[ch]  ->Fill( vPE_total_nosmear[ch], fWgt);
            if( fDecayTime > 2000 ) hPE_total_vs_prompt[type][ch]->Fill( fPE_total[ch], fPE_prompt[ch] );
          }
          hPromptFrac[type]->Fill( fElShowerPhel_prompt / fElShowerPhel_qc );
        }
      //}
    
    }



    // =====================================================
    // Calorimetry plots
    
    if( goodShower2D ) {
     
      if( goodShower3D ) {
        
        // Calculate charge-based energy for electron ion. track
        float ionQ = fElTrackCharge;
        float energyQtrk  = ( ionQ / fRecombIon ) * fWion; 
        
        hEnergyTrk[type]->Fill( energyQtrk, fWgt );
        if( isMC ) {
          float trueE_trk = fTrue_ElTrackEnergyDep;
          float resEQ_Trk = (energyQtrk - trueE_trk ) / trueE_trk;
          hEvsRes_E_Trk  ->Fill( trueE_trk, resEQ_Trk);
          hTrue_EnergyDepTrk->Fill(trueE_trk);
        }
          
      
        if( goodShower3D_dTcut ) {
          
          // Define charges and photons 
          float Q_ion       = fElTrackCharge;
          float Q_phot      = fElShowerCharge-fElTrackCharge;
          float Q_shower    = fElShowerCharge;
          float L_shower    = fElShowerPhel_qc / fElShowerVis;
 
          // Calculate charge-based energies 
          float energyQ     = ( Q_shower / fRecomb ) * fWion; // shower energy (uniform recomb)
          float energyQ_2R  = ( Q_ion/fRecombIon 
                                + Q_phot/fRecombPh ) * fWion; // shower energy (non-uniform recomb)

          // Calculate Q+L based energy
          float energyL       = -9.;
          float energyQL      = -9.;
          float energyQL_LogL = -9.;
          if( Q_shower > 0 && L_shower > 0 && fElShowerVis > 0 ) {

            // L-only estimator of energy          
            energyL     = ( L_shower / ( fExcRatio + 1. - fRecomb ) ) * fWion;
            
            // simplest case, no recomb assumption needed
            energyQL = ( Q_shower + L_shower )  * fWph; 
          
            // Log-likelihood method.  First set parameters:
            if( fDoLikelihoodFit ) {
              float E_ll = 1.;
              float E_ul = 150.;
              float dE = 0.01;
              f_LogLEnergyQL->SetParameter("aveDriftTime",fAveDriftTime);
              f_LogLEnergyQL->SetParameter("elecLifetime",fElectronLifetime);
              f_LogLEnergyQL->SetParameter("fracVis",fElShowerVis);
              f_LogLEnergyQL->SetParameter("chargeRes", 0.); // assume 0 for now
              f_LogLEnergyQL->SetParameter("peRes",     0.); // assume 0 for now
              f_LogLEnergyQL->SetParameter("L",L_shower);
              f_LogLEnergyQL->SetParameter("Q_ion",Q_shower);
              energyQL_LogL   = f_LogLEnergyQL->GetMinimumX( E_ll, E_ul, 0.01,100,false);
              // Check for cases where minimization failed (ie, the most likely
              // energy is very close to the Emin/Emax bounds)
              if( energyQL_LogL < E_ll + 3*dE || energyQL_LogL > E_ul - 3*dE ) {
                energyQL_LogL = -9;
                hEnergyQL_LogL_Fail[type]    ->Fill( energyQL );
              }
            } else {
              energyQL_LogL = energyQL;
            }

          }
            
          float trueE_dep   = fTrue_ElShowerEnergyDep;
          float trueE       = trueE_dep; if( fUseTrueEnergy ) trueE = fTrue_ElEnergy;
          float trueL_shower = fTrue_ElShowerPhotons;
          float trueQ_shower = fTrue_ElShowerCharge;
           
           // ----------------------------------------------------
          
          // Fill histograms for visibility, Q, L, PE
          if( fElShowerVis > 0 ) hElShowerVis[type]->Fill( fElShowerVis*1000., fWgt);
          hPE[type]       ->Fill(  fElShowerPhel,    fWgt);
          hL[type]        ->Fill( L_shower, fWgt);
          hQ[type]        ->Fill( Q_shower,     fWgt);
          hQvsL[type]     ->Fill( Q_shower, L_shower, fWgt);
          
          // Fill histograms for energy
          hEnergyQ[type]    ->Fill( energyQ );
          hEnergyL[type]    ->Fill( energyL );
          hEnergyQ_2R[type] ->Fill( energyQ_2R );
          hEnergyQL[type]         ->Fill( energyQL, fWgt);
          hEnergyQL_LogL[type]    ->Fill( energyQL_LogL, fWgt);
          hEnergy_QvsL[type]->Fill( energyQ, energyL);
          hEnergy_QvsQL[type]->Fill( energyQ, energyQL);
          if( fIsRealData ) {
            hTimeVsCharge     ->Fill( day, Q_shower );
            hTimeVsLight     ->Fill( day, L_shower );
            hRunVsEnergy      ->Fill(fRunNumber, energyQ );
            hRunVsLight     ->Fill( fRunNumber, L_shower );
            hRunVsCharge     ->Fill( fRunNumber, Q_shower );
          }
          
          // Recombination histogram
//        hR[type]        ->Fill( (1+fExcRatio)/(1+(fLCorrFactor*L_shower)/(fQCorrFactor*Q_shower)));

          // Light yield histogram and map
          hLightYield[type]->Fill( fElShowerPhel / energyQ, fWgt);
          hLY_ZX[type]      ->Fill( fElShowerCentroid_Z, fElShowerCentroid_X, fElShowerPhel );
          hE_ZX[type]      ->Fill( fElShowerCentroid_Z, fElShowerCentroid_X, energyQ );
          hLY_ZY[type]    ->Fill( fElShowerCentroid_Z, fElShowerCentroid_Y, fElShowerPhel );
          hE_ZY[type]     ->Fill( fElShowerCentroid_Z, fElShowerCentroid_Y, energyQ );
          hLY_entries_ZX[type]->Fill( fElShowerCentroid_Z, fElShowerCentroid_X );
          hLY_entries_ZY[type]->Fill( fElShowerCentroid_Z, fElShowerCentroid_Y );
          
          // Truth information and energy resolution histograms
          if( !fIsRealData ) {
          
            float Qcol_shower     = fElShowerChargeCol;
            float Qcol_total      = fTotalChargeCol[1]; // coll. plane
            float trueQcol_total  = fTrue_TotalChargeCol;
            
            // true deposited energy 
            hTrue_EnergyDep ->Fill(trueE_dep, fWgt);
            
            // collected charge + PE resolution distrbutions 
            //h_Qcol_vs_QcolRes ->Fill( trueQcol_total, (Qcol_total - trueQcol_total)/trueQcol_total );
            //h_Pe_vs_PeRes     ->Fill( fTrue_ElShowerPhel, (fElShowerPhel - fTrue_ElShowerPhel)/fTrue_ElShowerPhel );
            hPERes  ->Fill( (fElShowerPhel-fTrue_ElShowerPhel)/fTrue_ElShowerPhel, fWgt);
          
            // charge and light (Q and L) 
            hQRes ->Fill( (Q_shower-trueQ_shower) / trueQ_shower, fWgt );
            hLRes ->Fill( (L_shower-trueL_shower) / trueL_shower, fWgt);
          
          
            hEvsRes_E_Q         ->Fill( trueE, (energyQ-trueE)/trueE );
            hEvsRes_E_Q_2R      ->Fill( trueE, (energyQ_2R-trueE)/trueE );
            if( energyQL > 0 ) 
              hEvsRes_E_QL        ->Fill( trueE, (energyQL-trueE)/trueE );
            if( energyQL_LogL > 0 ) 
              hEvsRes_E_QL_LogL   ->Fill( trueE, (energyQL_LogL-trueE)/trueE );
          
            hEvs_Qtrue          ->Fill( trueE, trueQ_shower ); 
            if( Q_shower > 0 ) hEvs_Q              ->Fill( trueE, Q_shower );
            if( L_shower > 0 ) hEvs_L              ->Fill( trueE, L_shower );
            hEvs_Ltrue          ->Fill( trueE, trueL_shower );
            hEvsRes_Q           ->Fill( trueE, (Q_shower-trueQ_shower)/trueQ_shower);
            hEvsRes_L           ->Fill( trueE, (L_shower-trueL_shower)/trueL_shower);
            /* 
            if( trueE < 7.5 && (energyQL-trueE)/trueE > 0. && fUseTrueEnergy && fDebugFlag ) {
              std::cout
              <<"Weird Q+L resolution: \n"
              <<"  TrueE  = "<<trueE<<"\n"
              <<"  DepE   = "<<fTrue_ElShowerEnergyDep<<"\n" 
              <<"  Q      = "<<Q_shower<<"\n"
              <<"  L      = "<<L_shower<<"\n"
              <<"  PE     = "<<fElShowerPhel_qc<<"\n"
              <<"  vis    = "<<fElShowerVis<<"\n"
              <<"  EQL    = "<<energyQL<<"\n"
              <<"  res    = "<<(energyQL-trueE)/trueE<<"\n";
            }*/

            //hEvsRes_Q           ->Fill( trueE, (Q_shower-hypQ_shower)/trueQ_shower);
            //hEvsRes_L           ->Fill( trueE, (L_shower-hypL_shower)/trueL_shower);
            //hEnergyVsResPE  ->Fill( trueE, resPE, fWgt);

            //hEnergyDepResQ    ->Fill( resEQ, fWgt);
            //hEnergyDepResQL   ->Fill( resEQL, fWgt);

          }// end if MC

        }//<-- goodshower3D && dT cut
      }// good shower3D 
    }//good shower 2D

    
    
    
    // =====================================================
    // Decay time plots
    if( fMichelOpticalID  ) { 

      if( !fDecayFit_ReqShower1 || goodShower2D ) {

        if( fDecayTime > fDecayFit_MinDecayTime && fDecayTime < fDecayFit_MaxDecayTime ) {
        nEvt_DecayFit_TimeRange++;

          // amplitude cuts
          if( isMC || ampCut ) {
            nEvt_DecayFit_AmpCut++;
            hDecayTime[type]->Fill(fDecayTime, fWgt); 
          }//<-- decayfit: amp cut
      
        }
      }
      
    }
      




  }// end loop over events

  ScaleMC();
  hTrue_ZX_Vs_Containment   ->Divide(hTrue_ZX_Vs_Energy);
 
  hLY_ZX[type]->Divide(hE_ZX[type]);
  hLY_ZY[type]->Divide(hE_ZY[type]);
  
}





// #####################################################################################
float GetMCSpatialWgt(float x, float y, float z){
 
  float w = 1.0;
  int xbin, ybin, zbin;
 
  // Map is structured with Z on the X-axis and X on the Y-axis
  //zbin = hTrueMuEndpt_ZX->GetXaxis()->FindBin(z);
  //xbin = hTrueMuEndpt_ZX->GetYaxis()->FindBin(x);
  //w *= hMCSpatialWgt_ZX->GetBinContent(zbin,xbin);

  xbin = hMCSpatialWgt->GetXaxis()->FindBin(x); 
  ybin = hMCSpatialWgt->GetYaxis()->FindBin(y); 
  zbin = hMCSpatialWgt->GetZaxis()->FindBin(z); 
  w   *= hMCSpatialWgt->GetBinContent(xbin,ybin,zbin);

  if( w > 0. ) {
    return w;
  } else {
    return 1.0;
  }

}

  

//######################################################################################3
// Scale the MC distributions to data
void ScaleMC(){
  float promptScale[2]        = {1., 1.};
  float totalScale[2]         = {1., 1.};
  for(int ch=0; ch<2; ch++){
    //if( hPE_prompt[0][ch]->GetEntries() <= 0 || hPE_total[0][ch] <= 0 ) continue;
    if( hPE_prompt[1][ch]->Integral() > 0 ) promptScale[ch] = hPE_prompt[0][ch]->Integral()/hPE_prompt[1][ch]->Integral();
    if( hPE_total[1][ch]->Integral() > 0 ) totalScale[ch] = hPE_total[0][ch]->Integral()/hPE_total[1][ch]->Integral();
    hPE_prompt_pretrigMC[ch]  ->Scale(promptScale[ch]);
    hPE_prompt_posttrigMC[ch]  ->Scale(promptScale[ch]);
    hPE_prompt[1][ch]     ->Scale(promptScale[ch]);
    hPE_total[1][ch]      ->Scale(totalScale[ch]);
    hPE_total_qc[1][ch]      ->Scale(totalScale[ch]);
    hPE_prompt_nosmearMC[ch]     ->Scale(promptScale[ch]);
    hPE_total_nosmearMC[ch]      ->Scale(totalScale[ch]);
  }
  hPE[1]->Scale( hPE[0]->Integral() / hPE[1]->Integral() );
  if( hDecayTime[0]->GetEntries() > 0 && hDecayTime[1]->GetEntries() > 0 ) {
    hDecayTime[1]->Scale( hDecayTime[0]->GetEntries() / hDecayTime[1]->Integral() );
  }
  if( hEnergyTrk[0]->GetEntries() > 0 && hEnergyTrk[1]->GetEntries() > 0 ) {
    hEnergyTrk[1]->Scale( hEnergyTrk[0]->Integral() / hEnergyTrk[1]->Integral() );
  }
  if( hEnergyQ[0]->GetEntries() > 0 && hEnergyQ[1]->GetEntries() > 0 ) {
    hQ[1]      ->Scale( hQ[0]->Integral() / hQ[1]->Integral() );
    hEnergyQ[1]->Scale( hEnergyQ[0]->Integral() / hEnergyQ[1]->Integral() );
    hEnergyL[1]->Scale( hEnergyL[0]->Integral() / hEnergyL[1]->Integral() );
    //hEnergyDepResQ->Scale( hEnergyQ[0]->GetEntries() / hEnergyQ[1]->Integral() );
    //hEnergyResQ->Scale( hEnergyQ[0]->Integral() / hEnergyQ[1]->Integral() );
    ScaleHistoDataMC( hEnergyQ_2R);
  }
  if( hEnergyQL[0]->GetEntries() > 0 && hEnergyQL[1]->GetEntries() > 0 ) {
//    hL[1]       ->Scale( hL[0]->GetEntries() / hL[1]->Integral() );
    ScaleHistoDataMC(hL);
    ScaleHistoDataMC(hEnergyQL_LogL);
    hR[1]       ->Scale( hR[0]->GetEntries() / hR[1]->Integral() );
    hTrue_R     ->Scale( GetHistMax(hR[1]) / GetHistMax(hTrue_R)); 
    //hTrue_RTrack     ->Scale( GetHistMax(hTrue_RTrack[1]) / GetHistMax(hTrue_RTrack)); 
    //hTrue_RShower     ->Scale( GetHistMax(hTrue_RShower[1]) / GetHistMax(hTrue_RShower)); 
    //hBeta[1]       ->Scale( hBeta[0]->GetEntries() / hBeta[1]->Integral() );
//    hTrue_Beta     ->Scale( GetHistMax(hBeta[1]) / GetHistMax(hTrue_Beta)); 
    hEnergyQL[1]->Scale( hEnergyQL[0]->GetEntries() / hEnergyQL[1]->Integral() );
    //hEnergyDepResQL->Scale( hEnergyQL[0]->GetEntries() / hEnergyDepResQL->Integral() );
    //hEnergyResQL->Scale( hEnergyQL[0]->GetEntries() / hEnergyResQL->Integral() );
    hElShowerCentroid_X[1]->Scale( hElShowerCentroid_X[0]->GetEntries() / hElShowerCentroid_X[0]->Integral() );
    hElShowerCentroid_Y[1]->Scale( hElShowerCentroid_Y[0]->GetEntries() / hElShowerCentroid_Y[0]->Integral() );
    hElShowerCentroid_Z[1]->Scale( hElShowerCentroid_Z[0]->GetEntries() / hElShowerCentroid_Z[0]->Integral() );
    hElShowerVis[1]->Scale( hElShowerVis[0]->GetEntries() / hElShowerVis[1]->Integral() );
  }
  hdEdx[1]->Scale( hdEdx[0]->Integral() / hdEdx[1]->Integral() );
}



//##########################################################################
// This function creates the light-based result plots and those used to
// demonstrate smearing and trigger efficiency.
void LightPlots(){
  
  float textSize = 0.035;
  float mar_l  = 0.15;
  float mar_r  = 0.05;
  float mar_t  = 0.05;
  float mar_b  = 0.15;
  float leg_x1 = 0.50;
  float leg_x2 = 1.-mar_r-0.02;
  float leg_y1 = 0.7;
  float leg_y2 = 1.-mar_t-0.02;
  float text_x1 =  mar_l+0.02;
  float text_x2 =  0.4;
  float text_y1 = 0.75;
  float text_y2 = 1.-mar_t-0.03;
  float axisTitleSize = 0.045;

  /*
  // Setting the Latex Header
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextFont(42);
  t->SetTextSize(textSize);
  t->SetTextAlign(13);
  */
  

  gStyle->SetOptStat(0);

  //=======================================================================
  // 1) Trigger efficiency:
  //=======================================================================
  TCanvas*    cTrigEff = new TCanvas("trigeff","trigeff",1000,500);
  TLegend*    lTrigEff[2];
  TPaveText*  tTrigEff[2];

  cTrigEff -> Divide(2,1);
  for(int ch=0; ch<2; ch++){
    cTrigEff -> cd(ch+1);
  
    float hmax = GetHistMax(hPE_prompt_posttrigMC[ch]);
    hPE_prompt_pretrigMC[ch]->Scale(1/hmax);
    hPE_prompt_posttrigMC[ch]->Scale(1/hmax);  
    

    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t ); 
    hPE_prompt_pretrigMC[ch]->SetTitle("");
    hPE_prompt_pretrigMC[ch]->SetMaximum(hPE_prompt_pretrigMC[ch]->GetMaximum()*1.4);
    
    hPE_prompt_pretrigMC[ch]->SetLineColor(kBlack);
    hPE_prompt_pretrigMC[ch]->GetYaxis()->SetTitleOffset(1.6);
    hPE_prompt_pretrigMC[ch]->GetYaxis()->SetTitleSize(axisTitleSize);
    hPE_prompt_pretrigMC[ch]->GetXaxis()->SetTitleOffset(1.3);
    hPE_prompt_pretrigMC[ch]->GetXaxis()->SetTitleSize(axisTitleSize);
    hPE_prompt_pretrigMC[ch]->GetYaxis()->SetTitle("Events (normlized)");
    hPE_prompt_pretrigMC[ch]->Draw("hist");
    hPE_prompt_pretrigMC[ch]->GetXaxis()->SetRangeUser(0., hPE_prompt_pretrigMC[ch]->GetXaxis()->GetXmax()*0.7);
    
    hPE_prompt_posttrigMC[ch]       ->SetLineColor(kBlack);
    hPE_prompt_posttrigMC[ch]       ->SetFillColor(38);
    hPE_prompt_posttrigMC[ch]       ->Draw("hist same");
    trigEffCut[ch]->SetParameter(3, hPE_prompt_posttrigMC[ch]->GetBinContent(hPE_prompt_posttrigMC[ch]->GetMaximumBin()) ); 
    trigEffCut[ch]          ->SetLineWidth(2);
    trigEffCut[ch]          ->SetLineColor(kRed);
    if( fTrigEff_P[ch] > 0 ) trigEffCut[ch]          ->Draw("same");
    hPE_prompt_pretrigMC[ch]->Draw("sameaxis"); // redraw axis
     
    lTrigEff[ch] = MakeLegend( leg_x1, leg_y2, textSize, 5); 
    lTrigEff[ch] ->AddEntry(hPE_prompt_pretrigMC[ch], "MC (all events)",       "L");
    lTrigEff[ch] ->AddEntry(hPE_prompt_posttrigMC[ch],"MC w/ trig. eff. cut", "LF");
    lTrigEff[ch] ->AddEntry(trigEffCut[ch],           "Trig. eff. function","L"); 
    lTrigEff[ch]  ->AddEntry((TObject*)0, Form("  P = %6.2f pe",trigEffCut[ch]->GetParameter(0)),"");
    lTrigEff[ch]  ->AddEntry((TObject*)0, Form("  K = %5.2f",trigEffCut[ch]->GetParameter(1)),"");
    lTrigEff[ch] ->Draw();
  
    tTrigEff[ch] = MakeTextBox( text_x1, leg_y2, textSize, 3);
    tTrigEff[ch] -> AddText("#bf{LArIAT}"); 
    sprintf(buffer,"%s MC",runtag.c_str()); tTrigEff[ch]->AddText(buffer);
    sprintf(buffer,"%s",pmttag[ch].c_str()); tTrigEff[ch]->AddText(buffer);
    tTrigEff[ch]->Draw();

    
    //AddTextLine(t, text_x1, leg_y2, 1, "#bf{LArIAT Preliminary}");
    //sprintf(buffer,"%s MC",runtag.c_str());
    //AddTextLine(t, text_x1, leg_y2, 2, buffer);
    //sprintf(buffer,"%s",pmttag[ch].c_str());
    //AddTextLine(t, text_x1, leg_y2, 3, buffer);

  }

  //=======================================================================
  // 2) Compare smeared and unsmeared total light distributions
  //=======================================================================
  TCanvas* cNoSmear[2];
  TLegend* lNoSmear[2];
  TPaveText* tNoSmear[2];

  gStyle->SetOptStat(0);

  leg_x1=0.6;

  // Make clones of needed histograms to avoid complicating
  // other plots containing these histograms
  TH1D* hcPE_total_nosmearMC[2];
  TH1D* hcPE_total[2][2];

  for(int ch=0; ch<2; ch++){
    
    tNoSmear[ch] = MakeTextBox( text_x1, leg_y2, textSize, 3);
    tNoSmear[ch] -> AddText("#bf{LArIAT}"); 
    sprintf(buffer,"%s MC",runtag.c_str()); tNoSmear[ch]->AddText(buffer);
    sprintf(buffer,"%s",pmttag[ch].c_str()); tNoSmear[ch]->AddText(buffer);
    
    hcPE_total_nosmearMC[ch] = (TH1D*)hPE_total_nosmearMC[ch]->Clone("tns0");
    hcPE_total[1][ch]           = (TH1D*)hPE_total[1][ch]->Clone("t0");
    cNoSmear[ch] = new TCanvas(Form("NoSmear%d",ch),Form("NoSmear%d",ch),500,500);
    gPad->SetLeftMargin(mar_l);
    gPad->SetRightMargin(mar_r);
    gPad->SetTopMargin(mar_t);
    gPad->SetBottomMargin(mar_b);
    hcPE_total_nosmearMC[ch] ->SetFillStyle(0);
    hcPE_total_nosmearMC[ch] ->SetLineColor(kBlack);
    hcPE_total_nosmearMC[ch] ->SetLineStyle(2);
    hcPE_total_nosmearMC[ch] ->GetYaxis()->SetTitleOffset(1.6);
    hcPE_total_nosmearMC[ch] ->GetYaxis()->SetTitleSize(axisTitleSize);
    hcPE_total_nosmearMC[ch] ->GetXaxis()->SetTitleOffset(1.3);
    hcPE_total_nosmearMC[ch] ->GetXaxis()->SetTitleSize(axisTitleSize);
    hcPE_total_nosmearMC[ch] ->SetTitle("");
    hcPE_total_nosmearMC[ch] ->SetMaximum(hcPE_total_nosmearMC[ch]->GetMaximum()*1.2);
    hcPE_total_nosmearMC[ch] ->Draw("hist");
    hcPE_total[1][ch]        ->SetFillStyle(0);
    hcPE_total[1][ch]        ->SetLineColor(kBlack);
    hcPE_total[1][ch]       ->Draw("same hist");
    hcPE_total_nosmearMC[ch]->Draw("sameaxis"); // redraw axis
    
    lNoSmear[ch] = MakeLegend( leg_x1, leg_y2, textSize, 3); 
    lNoSmear[ch] ->AddEntry(hcPE_total_nosmearMC[ch], "MC before smearing",       "L");
    lNoSmear[ch] ->AddEntry(hcPE_total[1][ch],        "MC after smearing", "L");
//    lNoSmear[ch] ->AddEntry((TObject*)0,  Form("#sigma = %5.2f#sqrt{S}",fSmearFactor[ch]),"");
    lNoSmear[ch] ->AddEntry((TObject*)0,  Form("#sigma = %5.2f%%",fSmearFactor[ch]*100.),"");
    lNoSmear[ch] ->Draw();
    tNoSmear[ch]  ->Draw();

  }




  //=======================================================================
  // 3) Official data-MC comparison plots
  //=======================================================================
  TCanvas* cDataMC[2];
  TLegend* lDataMC[2][2];
  TPaveText* header[2];

  leg_x1=0.6;
  leg_y1=0.75;

  // Make clones of needed histograms to avoid complicating
  // other plots containing these histograms
  TH1D* hc2PE_prompt[2][2];
  TH1D* hc2PE_total[2][2];
  TH1D* hc2PE_prompt_mcerr[2];
  TH1D* hc2PE_total_mcerr[2];
  
  for(int ch=0; ch<2; ch++){
  
    hc2PE_prompt[0][ch] = (TH1D*)hPE_prompt[0][ch]->Clone();
    hc2PE_prompt[1][ch] = (TH1D*)hPE_prompt[1][ch]->Clone();
    hc2PE_prompt_mcerr[ch] = (TH1D*)hPE_prompt[1][ch]->Clone();
    
    hc2PE_total[0][ch] = (TH1D*)hPE_total[0][ch]->Clone();
    hc2PE_total[1][ch] = (TH1D*)hPE_total[1][ch]->Clone();
    hc2PE_total_mcerr[ch] = (TH1D*)hPE_total[1][ch]->Clone();
//    hc2PE_total[0][ch] = (TH1D*)hPE_total_qc[0][ch]->Clone();
//    hc2PE_total[1][ch] = (TH1D*)hPE_total_qc[1][ch]->Clone();
//    hc2PE_total_mcerr[ch] = (TH1D*)hPE_total_qc[1][ch]->Clone();
    /*   
    hc2PE_prompt[0][ch] = hPE_prompt[0][ch];
    hc2PE_prompt[1][ch] = hPE_prompt[1][ch];
    hc2PE_total[0][ch] =  hPE_total_qc[0][ch];
    hc2PE_total[1][ch] =  hPE_total_qc[1][ch];
    hc2PE_prompt_mcerr[ch] = hPE_prompt[1][ch];
    hc2PE_total_mcerr[ch] = hPE_total_qc[1][ch];
    */

    if( ch==0 ) cDataMC[ch] = new TCanvas("DataMC0","DataMC0",1000,500);
    if( ch==1 ) cDataMC[ch] = new TCanvas("DataMC1","DataMC1",1000,500);
    cDataMC[ch] -> Divide(2,1);
    int index=0;
    
    // ---------------------------------------
    // Prompt light
    index=0;
    cDataMC[ch] -> cd(index+1);
    gPad->SetLeftMargin(mar_l);
    gPad->SetRightMargin(mar_r);
    gPad->SetTopMargin(mar_t);
    hc2PE_prompt[1][ch] ->GetYaxis()->SetTitleOffset(1.3);
    hc2PE_prompt[1][ch] ->SetTitle("");
    hc2PE_prompt[1][ch] ->SetMaximum(hc2PE_prompt[0][ch]->GetMaximum()*1.3);
    hc2PE_prompt[1][ch] ->SetLineColor(kBlack);
    hc2PE_prompt[1][ch] ->SetFillColor(kYellow-10); //kGray); //(kRed-10); // kGray
    hc2PE_prompt[1][ch]  ->GetYaxis()->SetTitleSize(0.045);
    hc2PE_prompt[1][ch]  ->GetXaxis()->SetTitleSize(0.045);
    hc2PE_prompt[1][ch] ->DrawCopy("hist");
    
    hc2PE_prompt_mcerr[ch]->SetFillColor(kGray+2);
    hc2PE_prompt_mcerr[ch]->SetFillStyle(3002);
    hc2PE_prompt_mcerr[ch]->DrawCopy("same E2");
    
    hc2PE_prompt[0][ch] ->SetMarkerStyle(20);
    hc2PE_prompt[0][ch] ->SetMarkerSize(0.8);
    hc2PE_prompt[0][ch] ->SetMarkerColor(kRed+2);
    hc2PE_prompt[0][ch] ->SetLineColor(kRed+2);
    hc2PE_prompt[0][ch] ->SetLineWidth(2);
    hc2PE_prompt[0][ch] ->DrawCopy("same E P X0");
    
    lDataMC[ch][index] = MakeLegend(leg_x1, leg_y2, textSize, 2);
    lDataMC[ch][index] ->AddEntry(hc2PE_prompt[1][ch], "MC prediction",       "F");
    lDataMC[ch][index] ->AddEntry(hc2PE_prompt[0][ch], "Data", "LPE");
    lDataMC[ch][index] ->Draw();
   
    /* 
    AddTextLine(t, text_x1, leg_y2, 1, "#bf{LArIAT Preliminary}");
    sprintf(buffer,"%s: %s (%s)",runtag.c_str(),pmttag[ch].c_str(),diamtag[ch].c_str());
    AddTextLine(t, text_x1, leg_y2, 2, buffer);
    sprintf(buffer,"%i Events",(int)hc2PE_prompt[0][ch]->GetEntries() );
    AddTextLine(t, text_x1, leg_y2, 3, buffer);
    */
    TPaveText* hd1 = MakeTextBox(text_x1, leg_y2, textSize, 3);
    hd1->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
    hd1->AddText("Cosmic Michel e^{+/-}");
    hd1->AddText(Form("%s",pmttag[ch].c_str()));
//    hd1->AddText(Form("%i Events",(int)hc2PE_prompt[0][ch]->GetEntries() ) );
//    sprintf(buffer,"%s: %s",runtag.c_str(),pmttag[ch].c_str()); header[ch]->AddText(buffer);
//    sprintf(buffer,"%i Events",(int)hc2PE_prompt[0][ch]->GetEntries() ); header[ch]->AddText(buffer);
    hd1->Draw();

    // Add Chi2 info
    TPaveText* pt1 = MakeTextBox(leg_x1, leg_y2-2.*textSize-0.01, textSize, 4);
    pt1->AddText(Form("Data-MC #chi^{2}_{w} = %5.2f", GetChi2Weighted(hc2PE_prompt[0][ch],hc2PE_prompt[1][ch]) ));
    pt1->AddText("Cuts:");
    pt1->AddText(Form("  #DeltaT > %3.1f #mus",fdTcut/1000.));
    if( fPE_Require3DShower ){ 
      pt1->AddText ("  3D shower reco'd");
    } else {
      pt1->AddText ("  2D shower reco'd");
//      pt1->AddText ("  3D #mu endpt reco'd");
    }
    pt1->Draw();
     
    
    // ---------------------------------------
    // Total light
    index=1;
    cDataMC[ch] -> cd(index+1);
    gPad->SetLeftMargin(mar_l);
    gPad->SetRightMargin(mar_r);
    gPad->SetTopMargin(mar_t);
    hc2PE_total[1][ch] ->GetYaxis()->SetTitleOffset(1.2);
    hc2PE_total[1][ch] ->SetTitle("");
    hc2PE_total[1][ch] ->SetMaximum(hc2PE_total[0][ch]->GetMaximum()*1.3);
    hc2PE_total[1][ch] ->SetLineColor(kBlack);
    //hc2PE_total[1][ch] ->SetFillColor(kGray);  //(kRed-10);
    hc2PE_total[1][ch] ->SetFillColor(kYellow-10); //kGray); //(kRed-10); // kGray
    hc2PE_total[1][ch]  ->GetYaxis()->SetTitleSize(0.045);
    hc2PE_total[1][ch]  ->GetXaxis()->SetTitleSize(0.045);
    hc2PE_total[1][ch] ->DrawCopy("hist");
    
    hc2PE_total_mcerr[ch]->SetFillColor(kGray+2);
    hc2PE_total_mcerr[ch]->SetFillStyle(3002);
    hc2PE_total_mcerr[ch]->DrawCopy("same E2");
    
    hc2PE_total[0][ch] ->SetMarkerStyle(20);
    hc2PE_total[0][ch] ->SetMarkerSize(0.8);
    hc2PE_total[0][ch] ->SetMarkerColor(kRed+2);
    hc2PE_total[0][ch] ->SetLineColor(kRed+2);
    hc2PE_total[0][ch] ->SetLineWidth(2);
    hc2PE_total[0][ch] ->DrawCopy("same E P X0");
    
    lDataMC[ch][index] = MakeLegend(leg_x1, leg_y2, textSize, 2);
    lDataMC[ch][index] ->AddEntry(hc2PE_total[1][ch], "MC prediction",       "F");
    lDataMC[ch][index] ->AddEntry(hc2PE_total[0][ch], "Data", "LPE");
    lDataMC[ch][index] ->Draw();
    
    TPaveText* hd2 = MakeTextBox(text_x1, leg_y2, textSize, 3);
    hd2->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
    hd2->AddText("Cosmic Michel e^{+/-}");
    hd2->AddText(Form("%s",pmttag[ch].c_str()));
//    sprintf(buffer,"%s: %s",runtag.c_str(),pmttag[ch].c_str()); header[ch]->AddText(buffer);
//    sprintf(buffer,"%i Events",(int)hc2PE_prompt[0][ch]->GetEntries() ); header[ch]->AddText(buffer);
    hd2->Draw();
    
    // Add Chi2 info
    TPaveText* pt2 = MakeTextBox(leg_x1, leg_y2-2.*textSize-0.01, textSize, 4);
    pt2->AddText(Form("Data-MC #chi^{2}_{w} = %5.2f", GetChi2Weighted(hc2PE_total[0][ch],hc2PE_total[1][ch]) ));
    pt2->AddText("Cuts:");
    pt2->AddText(Form("  #DeltaT > %3.1f #mus",fdTcut/1000.));
    if( fPE_Require3DShower ){ 
      pt2->AddText ("  3D shower reco'd");
    } else {
      pt2->AddText ("  2D shower reco'd");
//      pt2->AddText ("  3D #mu endpt reco'd");
    }
    pt2->Draw();
    
  }

  
  
  //=======================================================================
  // 4) Combined PE distribution
  //=======================================================================
  // skip for now
  
  
  // To add multiple stat boxes:
  //
  // // Force drawing of current histogram stat box
  // c->Update()
  // 
  // TPaveStats *p2 = (TPaveStats*)h2->FindObject("stats");
  // p2->SetY1NDC(0.05);
  // p2->SetY2NDC(0.45);
  // p2->SetTextColor(kRed);
  // c->Modified();

  TCanvas* cextra = new TCanvas("cextra","cextra",100,100);

}





//##########################################################################
// Plot energy and energy resolution histograms
void EnergyPlots(bool doResolutionSlices = false){

//  gStyle              ->SetTitleFont(62,"t");
  
  // Start program timer
  std::time_t startTime = time(0);
  

  std::cout
  <<"====================================================\n"
  <<"Making energy plots...\n";

  Loop(0);
  
  float textSize = 0.035;
//  float textSize = 0.04;
  float mar_l   = 0.15;
  float mar_r  = 0.05;
  float mar_t  = 0.05;
  float mar_b   = 0.15;
  float leg_x1 = 0.60;
  float leg_x2 = 1.-mar_r-0.02;
  float leg_y1 = 0.7;
  float leg_y2 = 1.-mar_t-0.02;
  float text_x1 =  mar_l+0.02;
  float text_x2 =  0.4;
  float text_y1 = 0.75;
  float text_y2 = 1.-mar_t-0.03;
  float axisTitleSize = 0.045;
  float axisLabelSize = 0.045; 
  
  float statH = 0.12;
  float statW = 0.17;
  
  //int ii=0;
  int dataMax, mcMax;
 
  // Scale the true deposited energy distributions 
  hTrue_EnergyDepTrk  ->Scale( hEnergyTrk[0]->GetEntries()/(float)hTrue_EnergyDepTrk->Integral() );
  hTrue_EnergyDep     ->Scale( hEnergyQ[0]->GetEntries() / (float)hTrue_EnergyDep->Integral() );
  
  
  std::cout
  <<"\n-----------------------------------------------\n"
  <<"Q and L distributions...\n";
  
  
  //=======================================================================
  // 1) Q and L
  //=======================================================================
  TCanvas* cQL  = new TCanvas("QL","QL",1000,500);
  cQL  ->Divide(2,1);  
 
  // Q
  cQL  ->cd(1);
  gPad->SetMargin(mar_l, 0.08, mar_b, mar_t);
  gStyle      ->SetOptStat(0);
  dataMax = hQ[0]->GetBinContent(hQ[0]->GetMaximumBin());
  mcMax = hQ[1]->GetBinContent(hQ[1]->GetMaximumBin());
  hQ[1]->SetTitle("");
  hQ[1] ->SetFillColor(kBlue-10); //(kBlue-10);
  hQ[1] ->SetLineWidth(0);
  hQ[1]->SetStats(0);
  FormatAxes( hQ[1], axisTitleSize, axisLabelSize, 1.0, 1.5); 
  hQ[1] ->SetMaximum( std::max(mcMax,dataMax)*1.2 );
  hQ[1] ->DrawCopy("hist");
  hQ[0] ->SetMarkerStyle(20);
  hQ[0] ->SetMarkerSize(0.7);
  hQ[0] ->SetLineColor(kBlack); //(kBlue+2);
  hQ[0] ->SetMarkerColor(kBlack); //(kBlue+2);
  //hQ[0] ->GetYaxis()->SetTitleOffset(1.5);
  hQ[0] ->DrawCopy("same P E X0");
  hQ[1] ->DrawCopy("sameaxis"); // redraw axis
  
  TPaveText* ptq = MakeTextBox(mar_l + 0.02, leg_y2, textSize, 2);
  ptq->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  ptq->AddText("Cosmic Michel e^{+/-}");
  ptq  ->Draw();

  TLegend* lq = MakeLegend(0.55, leg_y2, textSize, 2);
  lq          ->AddEntry(hQ[1],"MC reco","LF");
  lq          ->AddEntry(hQ[0],"Data","LPE");
  lq          ->Draw();
  
  TPaveText* ptq2 = MakeTextBox(0.55, leg_y2-textSize*3, textSize, 3);
  ptq2 ->AddText(Form("#LT Q #GT_{data} = %5.1f #times10^{4} e^{-}",hQ[0]->GetMean()/1.e4));
  ptq2 ->AddText(Form("#LT Q #GT_{MC} = %5.1f #times10^{4} e^{-}",hQ[1]->GetMean()/1.e4));
  ptq2  ->Draw();
  
  // ......................
  
  //TPaveStats *st_q0 = (TPaveStats*)hQ[0]->FindObject("stats");
  //TPaveStats *st_q1 = (TPaveStats*)hQ[1]->FindObject("stats");
  //st_q0->SetOptStat(1110);
  //st_q1->SetOptStat(1110);

  // ........................

  
  
  
  // L
  cQL  ->cd(2);
  gPad->SetMargin(mar_l, 0.08, mar_b, mar_t);
  gStyle      ->SetOptStat(0);
  dataMax = hL[0]->GetBinContent(hL[0]->GetMaximumBin());
  mcMax = hL[1]->GetBinContent(hL[1]->GetMaximumBin());
  hL[1]->SetTitle("");
  hL[1] ->SetFillColor(kOrange-9); //(kRed-10);
  hL[1] ->SetLineWidth(0);
  hL[1]->SetStats(0);
  FormatAxes( hL[1], axisTitleSize, axisLabelSize, 1.0, 1.5); 
  hL[1] ->SetMaximum( std::max(mcMax,dataMax)*1.2 );
  hL[1] ->DrawCopy("hist");
  hL[0] ->SetMarkerStyle(20);
  hL[0] ->SetMarkerSize(0.7);
  hL[0] ->SetLineColor(kBlack); //(kOrange+2);
  hL[0] ->SetMarkerColor(kBlack); //(kOrange+2);
  //hL[0] ->GetYaxis()->SetTitleOffset(1.4);
  hL[0] ->DrawCopy("same P E X0");
  hL[1] ->DrawCopy("sameaxis"); // redraw axis
  
  TPaveText* ptl = MakeTextBox(mar_l + 0.02, leg_y2, textSize, 2);
  ptl->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  ptl->AddText("Cosmic Michel e^{+/-}");
  //ptl ->AddText(Form("%s Dataset",runtag.c_str()));
  ptl  ->Draw();

  TLegend* ll = MakeLegend(0.55, leg_y2, textSize, 2);
  ll          ->AddEntry(hL[1],"MC reco","LF");
  ll          ->AddEntry(hL[0],"Data","LPE");
  ll          ->Draw();
  
  TPaveText* ptl2 = MakeTextBox(0.55, leg_y2-textSize*3, textSize, 3);
  ptl2 ->AddText(Form("#LT L #GT_{data} = %5.1f #times10^{4} #gamma",hL[0]->GetMean()/1.e4));
  ptl2 ->AddText(Form("#LT L #GT_{MC} = %5.1f #times10^{4} #gamma",hL[1]->GetMean()/1.e4));
  ptl2  ->Draw();



     
  


  //======================================================================
  // 2) Q-based ionization-only energy spectrum
  // =====================================================================
  std::cout
  <<"\n-----------------------------------------------\n"
  <<"Ionization-only energy (electron track)...\n";
  TCanvas* cEnergyTrk = new TCanvas("energyTrk","energyTrk",600,600);
  gPad->SetMargin(mar_l,mar_r,mar_b,mar_t);
  gStyle->SetOptStat(0);
  
  hEnergyTrk[1]->SetTitle("");
  hEnergyTrk[1]->SetMaximum( 1.1*std::max(GetHistMax(hEnergyTrk[1]),GetHistMax(hEnergyTrk[0]) ) );
  hEnergyTrk[1]->SetMaximum( 1.1*std::max(float(hEnergyTrk[1]->GetMaximum()),GetHistMax(hTrue_EnergyDepTrk)) );
  hEnergyTrk[1]->SetStats(0);
  hEnergyTrk[1]->SetFillColor(kWhite);
  hEnergyTrk[1]->SetLineColor(kBlack);
  FormatAxes(hEnergyTrk[1], axisTitleSize, axisLabelSize, 1.0, 1.5 );
  hEnergyTrk[1]->DrawCopy("hist");

  TH1D* hEnergyTrkErr = (TH1D*)hEnergyTrk[1]->Clone("energyTrkMCErr");
  //hEnergyTrkErr = (TH1D*)hEnergyTrk[1]->Clone("energyTrk_MC_Err");
  hEnergyTrkErr->SetFillStyle(3001);
  hEnergyTrkErr->SetFillColor(kGray+1);
  hEnergyTrkErr->SetMarkerColor(kGray+1);
  hEnergyTrkErr->SetLineWidth(0);
  //hEnergyTrkMCErr->SetFillStyle(3001);
  //hEnergyTrkMCErr->SetFillColor(kGray);
  //hEnergyTrkMCErr->SetMarkerColor(kGray);
  hEnergyTrkErr->DrawCopy("same E2");
  hEnergyTrk[1]->DrawCopy("same hist");
 
 /* 
  hEnergyTrkMCErr->SetFillStyle(3002);
  hEnergyTrkMCErr->SetFillColor(kGray+2);
  hEnergyTrkMCErr->SetMarkerColor(kGray+2);
  hEnergyTrkMCErr->SetLineWidth(0);
  //hEnergyTrkMCErr->SetFillStyle(3001);
  //hEnergyTrkMCErr->SetFillColor(kGray);
  //hEnergyTrkMCErr->SetMarkerColor(kGray);
  hEnergyTrkMCErr->DrawCopy("same E2");
  */

  hEnergyTrk[0] ->SetMarkerStyle(20);
  hEnergyTrk[0] ->SetMarkerSize(0.7);
  hEnergyTrk[0] ->SetLineColor(kBlue);
  hEnergyTrk[0] ->SetMarkerColor(kBlue);
  hEnergyTrk[0] ->DrawCopy("same P E X0");
  hEnergyTrk[1] ->DrawCopy("same axis"); // redraw axis

  hTrue_EnergyDepTrk ->SetLineStyle(2);
  hTrue_EnergyDepTrk ->SetLineColor(kBlack);
  hTrue_EnergyDepTrk->DrawCopy("hist same");

  TPaveText* headTrk = MakeTextBox(leg_x1, leg_y2, textSize, 4);
  headTrk->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  //headTrk ->AddText(Form("%i Events",(int)hEnergyTrk[0]->GetEntries()));
  headTrk->AddText("Michel Electrons"); //e^{+/-}");
  headTrk ->AddText("Q-only ion. trk energy")->SetTextColor(kBlue);
  headTrk ->AddText("w/recomb. correction")->SetTextColor(kBlue);
  headTrk->Draw();

  //TPaveText* ptTrk  = MakeTextBox(leg_x1, leg_y2-4.5*textSize - 0.02, textSize, 1);
  //ptTrk   ->AddText(Form("Data-MC_{reco} #chi^{2}_{#nu} = %5.2f", GetChi2(hEnergyTrk[0],hEnergyTrk[1])));
  //ptTrk   ->AddText(Form("Data-MC_{true} #chi^{2}_{#nu} = %5.2f", GetChi2(hEnergyTrk[0],hTrue_EnergyDepTrk)));
//  ptTrk   ->AddText("Cuts:");
//  ptTrk   ->AddText(Form("  #DeltaT > %3.1f #mus",fdTcut/1000.));
//  ptTrk   ->AddText (Form("  3D shower (Npts #geq %1d)",fMinNumPts3D ));
  //ptTrk   ->Draw();
  TLegend* legTrk = MakeLegend(leg_x1+0.08, leg_y2 - 5.5*textSize-0.07, textSize, 3);
  legTrk  ->AddEntry(hTrue_EnergyDepTrk, "MC true E_{dep}^{ion}", "L");   
  legTrk  ->AddEntry(hEnergyTrk[1], "MC reco", "L");
  legTrk  ->AddEntry(hEnergyTrk[0], "Data", "LPE");
  legTrk  ->Draw();
  
  
  
  // --------------------------------------------------
  // Format all energy histograms
  axisTitleSize = 0.045;
  axisLabelSize = 0.040; 

  // true E dep ------------------
    hTrue_EnergyDep ->SetLineStyle(2);
    hTrue_EnergyDep ->SetLineColor(kBlack);

  // E_Trk
    hEnergyTrk[1]->SetFillColor(kWhite);
    hEnergyTrk[1]->SetLineColor(kBlack);
    FormatAxes(hEnergyTrk[1], axisTitleSize, axisLabelSize, 1.0, 1.5);
    // MC error bars
    TH1D* hEnergyTrk_mcerr = (TH1D*)hEnergyTrk[1]->Clone("energyTrk_mcerr");
    hEnergyTrk_mcerr->SetFillStyle(3002);
    hEnergyTrk_mcerr->SetFillColor(kGray+2);
    hEnergyTrk_mcerr->SetMarkerColor(kBlack);
    // Data
    hEnergyTrk[0] ->SetMarkerStyle(20);
    hEnergyTrk[0] ->SetMarkerSize(0.7);
    hEnergyTrk[0] ->SetLineColor(kBlue);
    hEnergyTrk[0] ->SetMarkerColor(kBlue);

 
  // E_Q --------------------------
    // MC
    hEnergyQ[1]->SetFillColor(kGray);
    hEnergyQ[1]->SetLineColor(kBlack);
    FormatAxes(hEnergyQ[1], axisTitleSize, axisLabelSize, 1.0, 1.5);
    // MC error bars
    TH1D* hEnergyQ_mcerr = (TH1D*)hEnergyQ[1]->Clone("energyQ_mcerr");
    hEnergyQ_mcerr->SetFillStyle(3002);
    hEnergyQ_mcerr->SetFillColor(kGray+2);
    hEnergyQ_mcerr->SetMarkerColor(kBlack);
    // Data
    hEnergyQ[0] ->SetMarkerStyle(20);
    hEnergyQ[0] ->SetMarkerSize(0.7);
    hEnergyQ[0] ->SetLineColor(kBlue);
    hEnergyQ[0] ->SetMarkerColor(kBlue);
  
  // E_Q varying recomb ---------------
    CopyHistoFormat(hEnergyQ[1],hEnergyQ_2R[1]);
    CopyHistoFormat(hEnergyQ[0],hEnergyQ_2R[0]);
    hEnergyQ_2R[0] ->SetMarkerStyle(4);
    hEnergyQ_2R[0] ->SetMarkerSize(0.7);
 
  // E_L ------------------------------
    CopyHistoFormat(hEnergyQ[1],hEnergyL[1]);
    // MC error bands
    TH1D* hEnergyL_mcerr = (TH1D*)hEnergyL[1]->Clone("energyL_mcerr");
    CopyHistoFormat(hEnergyQ_mcerr,hEnergyL_mcerr);
    // Data
    CopyHistoFormat(hEnergyQ[0],hEnergyL[0]);
    hEnergyL[0] ->SetMarkerStyle(21);
    hEnergyL[0] ->SetLineColor(kRed+2);
    hEnergyL[0] ->SetMarkerColor(kRed+2);
  //
  
  // E_QL --------------------------
    // MC
    CopyHistoFormat(hEnergyQ[1],hEnergyQL[1]);
    // MC error bands
    TH1D* hEnergyQL_mcerr = (TH1D*)hEnergyQL[1]->Clone("energyQL_mcerr");
    CopyHistoFormat(hEnergyQ_mcerr,hEnergyQL_mcerr);
    // Data
    CopyHistoFormat(hEnergyQ[0],hEnergyQL[0]);
    hEnergyQL[0] ->SetMarkerStyle(21);
    hEnergyQL[0] ->SetLineColor(kMagenta+2);
    hEnergyQL[0] ->SetMarkerColor(kMagenta+2);
  
  // E_QL LogL (uniform recomb ) --------------------
    // MC
    CopyHistoFormat(hEnergyQL[1],hEnergyQL_LogL[1]);
    TH1D* hEnergyQL_LogL_mcerr = (TH1D*)hEnergyQL_LogL[1]->Clone("energyQL_LogL_mcerr");
    CopyHistoFormat(hEnergyQ_mcerr,hEnergyQL_LogL_mcerr);
    //hEnergyQL_LogL[1] ->SetTitleSize(0.08);
    //hEnergyQL_LogL[1] ->GetXaxis()->SetTitleSize(0.08);
    //FormatAxes(hEnergyQL_LogL[1], axisTitleSize, axisLabelSize, 1.0, 1.4);
    //FormatAxes(hEnergyQL_LogL[1], 0.05, 0.05, 1.0, 1.2);
    //hEnergyQL_LogL[1] ->SetFillColor(kGray); 
    //hEnergyQL_LogL[1] ->SetLineColor(kGray+2);
  
    // Data
    CopyHistoFormat(hEnergyQL[0],hEnergyQL_LogL[0]);
    hEnergyQL_LogL[0] ->SetMarkerStyle(23);
    hEnergyQL_LogL[0] ->SetMarkerSize(0.8);
    hEnergyQL_LogL[0] ->SetMarkerColor(kGreen+2);
    hEnergyQL_LogL[0]->SetLineColor(kGreen+2);
  
  
  
  
  

  //======================================================================
  // 3) Q vs. Q+L energy spectra
  //====================================================================
  std::cout
  <<"\n-----------------------------------------------\n"
  <<"Q vs. Q+L shower energy...\n";
  TCanvas* cEnergy = new TCanvas("energy","energy",1000,500);
  cEnergy->Divide(2,1);
  
  // ------------------------------------------
  // Q-only energy
  cEnergy     ->cd(1);
  gPad        ->SetMargin(mar_l,mar_r,mar_b,mar_t);
  gStyle      ->SetOptStat(0);
  hEnergyQ[1] ->SetTitle("");
  hEnergyQ[1] ->SetStats(0);
  hEnergyQ[1]->SetMaximum( 1.1*std::max(GetHistMax(hEnergyQ[1]),GetHistMax(hEnergyQ[0]) ) );
  hEnergyQ[1]->SetMaximum( 1.1*std::max(float(hEnergyQ[1]->GetMaximum()),GetHistMax(hTrue_EnergyDep)) );
  hEnergyQ[1] ->DrawCopy("hist");
  hEnergyQ_mcerr->DrawCopy("same E2");
  hEnergyQ[0] ->DrawCopy("same P E X0");
  hEnergyQ[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDep->DrawCopy("hist same");
  
  TPaveText* hdr1 = MakeTextBox(leg_x1, leg_y2, textSize, 4);
  hdr1->AddText(Form("#bf{LArIAT %s}",runtag.c_str()));
  hdr1->AddText("Cosmic Michel e^{+/-}");
  //hdr1 ->AddText(Form("%s Dataset",runtag.c_str()));
  //hdr1 ->AddText(Form("%i Events",(int)hEnergyQ[0]->GetEntries()));
  hdr1 ->AddText("Q-only shower energy")->SetTextColor(kBlue);
  hdr1 ->AddText("w/recomb. correction")->SetTextColor(kBlue);
  hdr1 ->Draw();

  TPaveText* pt1 = MakeTextBox(leg_x1, leg_y2-4.*textSize - 0.02, textSize, 2);
  pt1 ->AddText(Form("Data-MC_{reco} #chi^{2}_{#nu} = %5.2f", GetChi2(hEnergyQ[0],hEnergyQ[1])));
  pt1 ->AddText(Form("Data-MC_{true} #chi^{2}_{#nu} = %5.2f", GetChi2(hEnergyQ[0],hTrue_EnergyDep)));
  //pt1 ->AddText("Cuts:");
  //pt1 ->AddText(Form("  #DeltaT > %3.1f #mus", fdTcut/1000.));
  //pt1 ->AddText(Form("  3D shower (Npts #geq %1d)",fMinNumPts3D ));
  pt1 ->Draw();
  
  TLegend* leg1  = MakeLegend(leg_x1+0.08, leg_y2 - 6.*textSize-0.07, textSize, 3);
  leg1 ->AddEntry(hTrue_EnergyDep, "MC true #it{E}_{#it{dep}}", "L");   
  leg1 ->AddEntry(hEnergyQ[1], "MC reco", "LF");
  leg1 ->AddEntry(hEnergyQ[0], "Data", "LPE");
  leg1 ->Draw();
 
 
  
  // ------------------------------------------
  // Q+L energy 
  cEnergy     ->cd(2);
  gPad        ->SetMargin(mar_l,mar_r,mar_b,mar_t);
  gStyle      ->SetOptStat(0);
  hEnergyQL[1] ->SetTitle("");
  hEnergyQL[1] ->SetStats(0);
  hEnergyQL[1]->SetMaximum( 1.1*std::max(GetHistMax(hEnergyQL[1]),GetHistMax(hEnergyQL[0]) ) );
  hEnergyQL[1]->SetMaximum( 1.1*std::max(float(hEnergyQL[1]->GetMaximum()),GetHistMax(hTrue_EnergyDep)) );
  hEnergyQL[1] ->DrawCopy("hist");
  hEnergyQL_mcerr->DrawCopy("same E2");
  hEnergyQL[0] ->DrawCopy("same P E X0");
  hEnergyQL[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDep->DrawCopy("hist same");
  
  TPaveText* hdr2 = MakeTextBox(leg_x1, leg_y2, textSize, 4);
  hdr2->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  hdr2->AddText("Cosmic Michel e^{+/-}");
  //hdr2 ->AddText(Form("%i Events",(int)hEnergyQL[0]->GetEntries()));
  hdr2 ->AddText("Combined Q+L")->SetTextColor(kMagenta+2);
  hdr2 ->AddText("shower energy")->SetTextColor(kMagenta+2);
  hdr2 ->Draw();

  TPaveText* pt2 = MakeTextBox(leg_x1, leg_y2-4.*textSize - 0.02, textSize, 2);
  pt2 ->AddText(Form("Data-MC_{reco} #chi^{2}_{#nu} = %5.2f", GetChi2(hEnergyQL[0],hEnergyQL[1])));
  pt2 ->AddText(Form("Data-MC_{true} #chi^{2}_{#nu} = %5.2f", GetChi2(hEnergyQL[0],hTrue_EnergyDep)));
//  pt2 ->AddText("Cuts:");
//  pt2 ->AddText(Form("  #DeltaT > %3.1f #mus", fdTcut/1000.));
//  pt2 ->AddText(Form("  3D shower (Npts #geq %1d)",fMinNumPts3D ));
  pt2 ->Draw();
  
  TLegend* leg2  = MakeLegend(leg_x1+0.08, leg_y2 - 6.*textSize-0.07, textSize, 3);
  leg2 ->AddEntry(hTrue_EnergyDep, "MC true #it{E}_{#it{dep}}", "L");   
  leg2 ->AddEntry(hEnergyQL[1], "MC reco", "LF");
  leg2 ->AddEntry(hEnergyQL[0], "Data", "LPE");
  leg2 ->Draw();
    
  


  // ===============================================================
  // 3b) Array of different energy spectra comparing to LogL method
  //      * Q-only
  //      * L-only
  //      * Q+L
  //      * Q+L (likelihood)
  std::cout
  <<"\n-----------------------------------------------\n"
  <<"Comparison of all energy spectra\n";
  TCanvas* cEnergyTrk2 = new TCanvas("energyTrk2","energyTrk2",600,600);
  TCanvas* cEnergy2a = new TCanvas("energy2a","energy2a",600,600);
  TCanvas* cEnergy2b = new TCanvas("energy2b","energy2b",600,600);
  TCanvas* cEnergy2c = new TCanvas("energy2c","energy2c",600,600);
  TCanvas* cEnergy2d = new TCanvas("energy2d","energy2d",600,600);

  hEnergyTrk[1]       ->SetTitleSize(16,"t");
  hEnergyQ[1]         ->SetTitleSize(16,"t");
  hEnergyL[1]         ->SetTitleSize(16,"t");
  hEnergyQL[1]        ->SetTitleSize(16,"t");
  hEnergyQL_LogL[1]   ->SetTitleSize(16,"t");
  FormatAxes(hEnergyTrk[1],     axisTitleSize, axisLabelSize, 1.0, 1.4);
  FormatAxes(hEnergyQ[1],       axisTitleSize, axisLabelSize, 1.0, 1.4);
  FormatAxes(hEnergyL[1],       axisTitleSize, axisLabelSize, 1.0, 1.4);
  FormatAxes(hEnergyQL[1],      axisTitleSize, axisLabelSize, 1.0, 1.4);
  FormatAxes(hEnergyQL_LogL[1], axisTitleSize, axisLabelSize, 1.0, 1.4);
  hEnergyTrk[0]     ->SetMarkerSize(0.8);
  hEnergyQ[0]       ->SetMarkerSize(0.8);
  hEnergyL[0]       ->SetMarkerSize(0.8);
  hEnergyQL[0]      ->SetMarkerSize(0.8);
  hEnergyQL_LogL[0] ->SetMarkerSize(1.2);

  // ------------------------------------------
  // Trk energy
  cEnergyTrk2   ->cd();
  gPad        ->SetMargin(mar_l,mar_r,0.12,0.08); gStyle  ->SetOptStat(0);
  hEnergyTrk[1]->SetTitle("");
  hEnergyTrk[1] ->DrawCopy("hist");
  hEnergyTrk_mcerr->DrawCopy("same E2");
  hEnergyTrk[0] ->DrawCopy("same P E X0");
  hEnergyTrk[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDepTrk->DrawCopy("hist same");

  TPaveText* b3_hd1trk = MakeTextBox(0.55, 0.90, 0.045,2);
  b3_hd1trk->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  b3_hd1trk->AddText("Cosmic Michel e^{+/-}");
  b3_hd1trk->Draw();

  TPaveText* b3_hd1btrk = MakeTextBox(0.55, 0.90-2*0.045-0.02, 0.045, 2);
  b3_hd1btrk->AddText("Q-only ion. trk energy")->SetTextColor(kBlue);
  b3_hd1btrk->AddText("w/recomb. correction")->SetTextColor(kBlue);
  b3_hd1btrk->Draw();
  
  TLegend* b3_leg1trk  = MakeLegend(0.6, 0.90-4*0.045-0.10, 0.045, 3);
  b3_leg1trk ->AddEntry(hTrue_EnergyDepTrk, "MC true E_{dep}^{ion}", "L");   
  b3_leg1trk ->AddEntry(hEnergyTrk[1], "MC reco", "LF");
  b3_leg1trk ->AddEntry(hEnergyTrk[0], "Data", "LPE");
  b3_leg1trk ->Draw();

  
  // ------------------------------------------
  // Q-only energy, uniform recomb (left)
  //cEnergy2     ->cd(1);
  cEnergy2a   ->cd();
  gPad        ->SetMargin(mar_l,mar_r,0.12,0.08); gStyle  ->SetOptStat(0);
  hEnergyQ[1]->SetTitle("");
  hEnergyQ[1] ->DrawCopy("hist");
  hEnergyQ_mcerr->DrawCopy("same E2");
  hEnergyQ[0] ->DrawCopy("same P E X0");
  hEnergyQ[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDep->DrawCopy("hist same");

  TPaveText* b3_hd1 = MakeTextBox(0.55, 0.90, 0.045,2);
  b3_hd1->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  b3_hd1->AddText("Cosmic Michel e^{+/-}");
  b3_hd1->Draw();

  TPaveText* b3_hd1b = MakeTextBox(0.55, 0.90-2*0.045-0.02, 0.045, 2);
  b3_hd1b->AddText("Q-only shower energy")->SetTextColor(kBlue);
  b3_hd1b->AddText("w/recomb. correction")->SetTextColor(kBlue);
  b3_hd1b->Draw();
  
  TLegend* b3_leg1  = MakeLegend(0.6, 0.90-4*0.045-0.10, 0.045, 3);
  b3_leg1 ->AddEntry(hTrue_EnergyDep, "MC true E_{dep}", "L");   
  b3_leg1 ->AddEntry(hEnergyQ[1], "MC reco", "LF");
  b3_leg1 ->AddEntry(hEnergyQ[0], "Data", "LPE");
  b3_leg1 ->Draw();
  
  // ------------------------------------------
  // L-only energy, uniform recomb (left)
  cEnergy2b   ->cd();
  gPad        ->SetMargin(mar_l,mar_r,0.12,0.08); gStyle  ->SetOptStat(0);
  hEnergyL[1]->SetTitle("");
  hEnergyL[1] ->DrawCopy("hist");
  hEnergyL_mcerr->DrawCopy("same E2");
  hEnergyL[0] ->DrawCopy("same P E X0");
  hEnergyL[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDep->DrawCopy("hist same");

  TPaveText* b3_hd2 = MakeTextBox(0.55, 0.90, 0.045,2);
  b3_hd2->AddText(Form("#bf{LArIAT: %s}",runtag.c_str()));
  b3_hd2->AddText("Cosmic Michel e^{+/-}");
  b3_hd2->Draw();

  TPaveText* b3_hd2b = MakeTextBox(0.55, 0.90-2*0.045-0.02, 0.045, 2);
  b3_hd2b->AddText("L-only shower energy")->SetTextColor(kRed+2);
  b3_hd2b->AddText("w/recomb. correction")->SetTextColor(kRed+2);
  b3_hd2b->Draw();
  
  TLegend* b3_leg2  = MakeLegend(0.6, 0.90-4*0.045-0.10, 0.045, 3);
  b3_leg2 ->AddEntry(hTrue_EnergyDep, "MC true E_{dep}", "L");   
  b3_leg2 ->AddEntry(hEnergyL[1], "MC reco", "LF");
  b3_leg2 ->AddEntry(hEnergyL[0], "Data", "LPE");
  b3_leg2 ->Draw();
  


  // ------------------------------------------
  // Q+L (middle)
  cEnergy2c   ->cd();
  ScaleHistoDataMC(hEnergyQL[1],hTrue_EnergyDep);
  gPad        ->SetMargin(mar_l,mar_r,0.12,0.08); gStyle  ->SetOptStat(0);
  hEnergyQL[1]  ->SetTitle("");
  hEnergyQL[1] ->SetStats(0);
  hEnergyQL[1]->SetMaximum( 1.1*std::max(GetHistMax(hEnergyQL[1]),GetHistMax(hEnergyQL[0]) ) );
  hEnergyQL[1]->SetMaximum( 1.1*std::max(float(hEnergyQL[1]->GetMaximum()),GetHistMax(hTrue_EnergyDep)) );
  hEnergyQL[1] ->DrawCopy("hist");
  hEnergyQL_mcerr->DrawCopy("same E2");
  hEnergyQL[0] ->DrawCopy("same P E X0");
  hEnergyQL[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDep->DrawCopy("hist same");

  TPaveText* b3_hd3 = MakeTextBox(0.55, 0.90, 0.045,2);
  b3_hd3->AddText(Form("#bf{LArIAT %s}",runtag.c_str()));
  b3_hd3->AddText("Cosmic Michel e^{+/-}");
  b3_hd3->Draw();

  TPaveText* b3_hd3b = MakeTextBox(0.55, 0.90-2*0.045-0.02, 0.045, 2);
  b3_hd3b->AddText("Combined Q+L")->SetTextColor(kMagenta+2);
  b3_hd3b->AddText("shower energy")->SetTextColor(kMagenta+2);
  b3_hd3b->Draw();
  
  TLegend* b3_leg3  = MakeLegend(0.6, 0.90-4*0.045-0.10, 0.045, 3);
  b3_leg3 ->AddEntry(hTrue_EnergyDep, "MC true E_{dep}", "L");   
  b3_leg3 ->AddEntry(hEnergyQL[1], "MC reco", "LF");
  b3_leg3 ->AddEntry(hEnergyQL[0], "Data", "LPE");
  b3_leg3 ->Draw();

  // ------------------------------------------
  // Q+L Log-Likelihood, uniform recomb
  //cEnergy2     ->cd(3);
  cEnergy2d->cd();
  gPad        ->SetMargin(mar_l,mar_r,0.12,0.08); gStyle  ->SetOptStat(0);
  ScaleHistoDataMC(hEnergyQL_LogL[1],hTrue_EnergyDep);
  hEnergyQL_LogL[1] ->SetTitle("");
  hEnergyQL_LogL[1] ->SetStats(0);
  hEnergyQL_LogL[1]->SetMaximum( 1.1*std::max(GetHistMax(hEnergyQL_LogL[1]),GetHistMax(hEnergyQL_LogL[0]) ) );
  hEnergyQL_LogL[1]->SetMaximum( 1.1*std::max(float(hEnergyQL_LogL[1]->GetMaximum()),GetHistMax(hTrue_EnergyDep)) );
  hEnergyQL_LogL[1] ->DrawCopy("hist");
  hEnergyQL_LogL_mcerr->DrawCopy("same E2");
  hEnergyQL_LogL[0] ->DrawCopy("same P E X0");
  hEnergyQL_LogL[1] ->DrawCopy("sameaxis"); // redraw axis
  hTrue_EnergyDep   ->DrawCopy("hist same");
  
  /*
  TPaveText* b3_hd3 = MakeTextBox(0.17, 0.88, 0.045, 3);
  b3_hd3->AddText(Form("#bf{LArIAT %s}",runtag.c_str()));
  b3_hd3->AddText("Cosmic Michel e^{+/-}");
  b3_hd3->AddText("Q+L (likelihood)")->SetTextColor(kGreen+2);
  b3_hd3 ->Draw();
  
  TLegend* b3_leg2  = MakeLegend(0.6, 0.88, 0.045, 3);
  b3_leg2 ->AddEntry(hTrue_EnergyDep, "MC true E_{dep}", "L");   
  b3_leg2 ->AddEntry(hEnergyQL_LogL[1], "MC reco", "LF");
  b3_leg2 ->AddEntry(hEnergyQL_LogL[0], "Data", "LPE");
  b3_leg2 ->Draw();
*/

  TPaveText* b3_hd4 = MakeTextBox(0.55, 0.90, 0.045,2);
  b3_hd4->AddText(Form("#bf{LArIAT %s}",runtag.c_str()));
  b3_hd4->AddText("Cosmic Michel e^{+/-}");
  b3_hd4->Draw();

  TPaveText* b3_hd4b = MakeTextBox(0.55, 0.90-2*0.045-0.02, 0.045, 3);
  b3_hd4b->AddText("Combined Q+L")->SetTextColor(kGreen+2);
  b3_hd4b->AddText("shower energy")->SetTextColor(kGreen+2);
  b3_hd4b->AddText("(likelihood method)")->SetTextColor(kGreen+2);
  b3_hd4b->Draw();
  
  TLegend* b3_leg4  = MakeLegend(0.6, 0.90-4*0.045-0.10, 0.045, 3);
  b3_leg4 ->AddEntry(hTrue_EnergyDep, "MC true E_{dep}", "L");   
  b3_leg4 ->AddEntry(hEnergyQL_LogL[1], "MC reco", "LF");
  b3_leg4 ->AddEntry(hEnergyQL_LogL[0], "Data", "LPE");
  b3_leg4 ->Draw();
  
  
  
  
  //**************************************
  // EQ vs EL

  TCanvas* cEnergy_QvsL[2];

  cEnergy_QvsL[0] = new TCanvas("cEnergy_QvsL_Data","Data",700,700);
  gPad->SetMargin(0.15,0.15,0.12,0.10);
  gPad->SetGrid();
  hEnergy_QvsL[0]->GetXaxis()->SetTitleOffset(1.3);
  hEnergy_QvsL[0]->SetTitle("Data");
  hEnergy_QvsL[0]->SetOption("colz");
  hEnergy_QvsL[0]->Draw();
 
  
  cEnergy_QvsL[1] = new TCanvas("cEnergy_QvsL_MC","MC",700,700);
  gPad->SetMargin(0.15,0.15,0.12,0.10);
  gPad->SetGrid();
  hEnergy_QvsL[1]->GetXaxis()->SetTitleOffset(1.3);
  hEnergy_QvsL[1]->SetTitle("MC");
  hEnergy_QvsL[1]->SetOption("colz");
  hEnergy_QvsL[1]->Draw();
  
  
  
  
  
  
  
  
  
  
  
  
  ///////////////////////////////
  // Write canvases to file
  outFile   ->cd();
  cQL       ->Write();
  cEnergyTrk->Write();
  cEnergy   ->Write();
  cEnergyTrk2->Write();
  cEnergy2a ->Write();
  cEnergy2b ->Write();
  cEnergy2c ->Write();
  cEnergy2d ->Write();
  ////////////////////////////////
 
 
  
  gStyle->SetOptFit(1); 
  
  if( doResolutionSlices ) { 
  
  // ===================================================================
  // RESOLUTION PLOTS
  //=====================================================================
  // 2.1) Analyze the Q and L resolution
  std::cout
  <<"===================================================\n"
  <<"Looking at Q and L resolutions for Michel events...\n";
  SetTrigEffParams(0);
  //fMaxMCEvts = -50000; //50000; //-1; //100000; //-1;
  //fMaxMCEvts = -1;
  RepMC();
    
  TF1* fResFit = new TF1("resFit","sqrt([0]^2/x + [1]^2)",1.,100.);
  fResFit->SetParameter(0,15.);
  fResFit->SetParameter(1,2.);
  fResFit->SetParLimits(0,0.,80.);
  fResFit->SetParLimits(1,0.,20.);
  fResFit->SetLineWidth(1);
  fResFit->SetNpx(1000);
  
  float Emin = 7.5; //2.5;
  float Emax = 52.5; //47.5;
  float y1 = 3;
  float y2 = 200; 
  //TCanvas* cQLres  = new TCanvas("QLres","QLres",900,450);
  //cQLres  ->Divide(2,1);  
  
  //auto gr_Qres = new TGraphAsymmErrors();  
  //auto gr_Lres = new TGraphAsymmErrors();  
  //auto gr_Qrms = new TGraphAsymmErrors();  
  //auto gr_Lrms = new TGraphAsymmErrors();  
  auto gr_Qsig = new TGraphAsymmErrors();  
  auto gr_Lsig = new TGraphAsymmErrors();  
  auto gr_Qrms = new TGraphAsymmErrors();  
  auto gr_Lrms = new TGraphAsymmErrors();  
  auto dummy = new TGraphAsymmErrors();
  
  std::vector<TF1> fvFitFuncQ;
  std::vector<TF1> fvFitFuncL;
  
  //gStyle->SetOptStat(1010);
  gStyle->SetOptFit(1111); 
    
    // Michel shower Q
    ResolutionSliceLoop(
      hEvs_Q, Emin, Emax, 9, 1.,
      true, "Q", "Michel Q",3,3,0.,2500e3,//-9,-9,
      1,15., false,
      -9, -9,
      gr_Qsig,
      gr_Qrms,
      fvFitFuncQ);
    
    // Michel shower L
    ResolutionSliceLoop(
      hEvs_L, Emin, Emax, 9, 1.,
      true, "L", "Michel L",3,3,0.,2500e3,
      0,-0.10, false,
      -9, -9,
      gr_Lsig,
      gr_Lrms,
      fvFitFuncL);
  
  
  // ------------------------------------------------------------
  // Plot overlays of the parameterized and actual distributions
  // for Q and L
  ParamDistOverlay(
    hEvs_Q, 7.5, 52.5, 9, 1., 0., 2500e3,
    "overlays_Q", "Michel Q", 3, 3, fFuncP_Q);
  ParamDistOverlay(
    hEvs_L, 7.5, 52.5, 9, 1., 0., 2500e3,
    "overlays_L", "Michel L", 3, 3, fFuncP_L);

  // -------------------------------------------------------------
  // Q-FITS:
  // Plot the parameters of the double-gaussian functions from the 
  // fits to Q at each energy slice
  std::cout<<"Beginning Q dist fit\n";
  
  int Nparams = 5;
  TGraphErrors* gr_Qparam[Nparams];
  for(size_t i=0; i<Nparams; i++)
    gr_Qparam[i] = new TGraphErrors();
  
  for(size_t i=0; i<fvFitFuncQ.size(); i++){
    double binEnergy, Qres, dE;
    gr_Qsig->GetPoint(i,binEnergy,Qres);
    dE = gr_Qsig->GetErrorX(i);
    std::cout<<"   "<<i<<"   E "<<binEnergy<<"\n";
    TF1 fitFunc = fvFitFuncQ[i];
    double* p = fitFunc.GetParameters();
    const double* ep = fitFunc.GetParErrors();
    double mu = p[1];
    double mu_err = ep[1];
    double sig = p[2]/p[1];
    double sig_err = sig*std::sqrt( std::pow(ep[2]/p[2],2) + std::pow(ep[1]/p[1],2) );
    double rBG = (p[3]*p[5])/(p[0]*p[2]);
    double rBG_err = rBG*std::sqrt( std::pow(ep[3]/p[3],2) + std::pow(ep[5]/p[5],2) + std::pow(ep[0]/p[0],2) + std::pow(ep[2]/p[2],2) );
//    double muDiff = p[4]-p[1];
//    double muDiff_err = std::sqrt( ep[4]*ep[4]+ ep[1]*ep[1] );
    double muBG = p[4]/p[1];
    double muBG_err = fabs(muBG)*std::sqrt( std::pow(ep[4]/p[4],2) + std::pow(ep[1]/p[1],2));
    double sigBG = p[5]/p[2];
    double sigBG_err = sigBG*std::sqrt( std::pow(ep[5]/p[5],2) + std::pow(ep[2]/p[2],2)  );
    std::cout
    <<"  Bin "<<i<<": Double-Gaus params \n"
    <<"    - sigma = "<<fitFunc.GetParameter(2)<<"\n"
    <<"    - rBG = "<<rBG<<" +/- "<<rBG_err<<"\n"
    <<"    - muBG = "<<muBG<<" +/- "<<muBG_err<<"\n"
    <<"    - sigBG = "<<sigBG<<" +/- "<<sigBG_err<<"\n";
    
    gr_Qparam[0]->SetPoint      ( gr_Qparam[0]->GetN(),   binEnergy, mu);
    gr_Qparam[0]->SetPointError ( gr_Qparam[0]->GetN()-1, dE, mu_err );
    gr_Qparam[1]->SetPoint      ( gr_Qparam[1]->GetN(),   binEnergy, sig);
    gr_Qparam[1]->SetPointError ( gr_Qparam[1]->GetN()-1, dE, sig_err );
    gr_Qparam[2]->SetPoint      ( gr_Qparam[2]->GetN(),   binEnergy, rBG);
    gr_Qparam[2]->SetPointError ( gr_Qparam[2]->GetN()-1, dE,rBG_err);
    gr_Qparam[3]->SetPoint      ( gr_Qparam[3]->GetN(),   binEnergy, muBG);
    gr_Qparam[3]->SetPointError ( gr_Qparam[3]->GetN()-1, dE,muBG_err);
    gr_Qparam[4]->SetPoint      ( gr_Qparam[4]->GetN(),   binEnergy, sigBG);
    gr_Qparam[4]->SetPointError ( gr_Qparam[4]->GetN()-1, dE,sigBG_err);
  }
  gStyle->SetOptFit(1);
  gr_Qparam[0]->SetMarkerStyle(20);
  gr_Qparam[0]->SetMarkerSize(0.8);
  gr_Qparam[0]->GetXaxis()->SetTitle("Deposited energy [MeV]");
  gr_Qparam[0]->GetXaxis()->SetTitleSize(0.06);
  gr_Qparam[0]->GetYaxis()->SetTitleSize(0.06);
  gr_Qparam[0]->GetYaxis()->SetTitleOffset(1.4);
  gr_Qparam[0]->GetXaxis()->SetLabelSize(0.05);
  gr_Qparam[0]->GetYaxis()->SetLabelSize(0.05);
  CopyTGraphFormat(gr_Qparam[0],gr_Qparam[1], true);
  CopyTGraphFormat(gr_Qparam[0],gr_Qparam[2], true);
  CopyTGraphFormat(gr_Qparam[0],gr_Qparam[3], true);
  CopyTGraphFormat(gr_Qparam[0],gr_Qparam[4], true);
  gr_Qparam[0]->GetYaxis()->SetTitle("Peak #mu [e]");
  gr_Qparam[1]->GetYaxis()->SetTitle("Peak #sigma/#mu");
  gr_Qparam[2]->GetYaxis()->SetTitle("A_{BG} / A_{peak}");
  gr_Qparam[3]->GetYaxis()->SetTitle("#mu_{BG} / #mu_{peak}");
  gr_Qparam[4]->GetYaxis()->SetTitle("#sigma_{BG} / #sigma_{peak}");
   
  gr_Qparam[0]->Fit(funcP_Q_fp[0]);
  gr_Qparam[1]->Fit(funcP_Q_fp[1]);
  gr_Qparam[2]->Fit(funcP_Q_fp[2]);
  gr_Qparam[3]->Fit(funcP_Q_fp[3]);
  gr_Qparam[4]->Fit(funcP_Q_fp[4]);

  TCanvas* c_Q_params = new TCanvas("Q_params","Q_params",1100,220);
  c_Q_params->Divide(5,1);
  gStyle->SetOptFit(1);
  for(size_t i=0; i<Nparams; i++){
    c_Q_params->cd(i+1);
    gPad->SetMargin(0.15,0.05,0.15,0.10);
    gr_Qparam[i]->GetYaxis()->SetRangeUser( 
      0.6*TMath::MinElement(gr_Qparam[i]->GetN(),gr_Qparam[i]->GetY()),
      1.3*TMath::MaxElement(gr_Qparam[i]->GetN(),gr_Qparam[i]->GetY()));
    gr_Qparam[i]->GetXaxis()->SetNdivisions(505);
    gr_Qparam[i]->GetYaxis()->SetNdivisions(505);
    gr_Qparam[i]->Draw("AP");
  }
  
  
  // -------------------------------------------------------------
  // L-FITS:
  // Plot the parameters of the single Guassian
  std::cout<<"Beginning L dist fit\n";
  
  Nparams = 2;
  TGraphErrors* gr_Lparam[Nparams];
  for(size_t i=0; i<Nparams; i++)
    gr_Lparam[i] = new TGraphErrors();
  
  for(size_t i=0; i<fvFitFuncL.size(); i++){
    double binEnergy, Qres, dE;
    gr_Lsig->GetPoint(i,binEnergy,Qres);
    dE = gr_Lsig->GetErrorX(i);
    std::cout<<"   "<<i<<"   E "<<binEnergy<<"\n";
    TF1 fitFunc = fvFitFuncL[i];
    double* p = fitFunc.GetParameters();
    const double* ep = fitFunc.GetParErrors();
    double mu = p[1];
    double mu_err = ep[1];
    double sig = p[2]/p[1];
    double sig_err = sig*std::sqrt( std::pow(ep[2]/p[2],2) + std::pow(ep[1]/p[1],2) );
    std::cout
    <<"  Bin "<<i<<": Double-Gaus params \n"
    <<"    - sigma = "<<fitFunc.GetParameter(2)<<"\n";
    
    gr_Lparam[0]->SetPoint      ( gr_Lparam[0]->GetN(),   binEnergy, mu);
    gr_Lparam[0]->SetPointError ( gr_Lparam[0]->GetN()-1, dE, mu_err );
    gr_Lparam[1]->SetPoint      ( gr_Lparam[1]->GetN(),   binEnergy, sig);
    gr_Lparam[1]->SetPointError ( gr_Lparam[1]->GetN()-1, dE, sig_err );
  }
  gStyle->SetOptFit(1);
  gr_Lparam[0]->SetMarkerStyle(20);
  gr_Lparam[0]->SetMarkerSize(0.8);
  gr_Lparam[0]->GetXaxis()->SetTitle("Deposited energy [MeV]");
  gr_Lparam[0]->GetXaxis()->SetTitleSize(0.06);
  gr_Lparam[0]->GetYaxis()->SetTitleSize(0.06);
  gr_Lparam[0]->GetYaxis()->SetTitleOffset(1.2);
  gr_Lparam[0]->GetXaxis()->SetLabelSize(0.05);
  gr_Lparam[0]->GetYaxis()->SetLabelSize(0.05);
  CopyTGraphFormat(gr_Lparam[0],gr_Lparam[1], true);
  gr_Lparam[0]->GetYaxis()->SetTitle("Peak #mu [#gamma]");
  gr_Lparam[1]->GetYaxis()->SetTitle("Peak #sigma/#mu");
   
  gr_Lparam[0]->Fit(funcP_L_fp[0]);
  gr_Lparam[1]->Fit(funcP_L_fp[1]);

  TCanvas* c_L_params = new TCanvas("L_params","L_params",600,300);
  c_L_params->Divide(2,1);
  gStyle->SetOptFit(1);
  for(size_t i=0; i<Nparams; i++){
    c_L_params->cd(i+1);
    gPad->SetMargin(0.15,0.05,0.15,0.10);
    gr_Lparam[i]->Draw("AP");
    gStyle->SetOptFit(1);
    gPad->Update();
  }


    




    // --------------------------------------------------------------------
    // 4.1) Michel electron resolution (including reconstruction effects).
    //      This will use the true deposited energy as reference since we
    //      are not making any requirements on containment.
    // --------------------------------------------------------------------
    std::cout
    <<"\n#########################################################\n"
    <<"Resolution slices: Michel electrons\n";
    float axts = 0.05;
    Emin = 7.5; //2.5;
    Emax = 52.5; //47.5;
    bool  useHybridRes = false;
    float fitRangeFac = 1.5;

    auto gr_sigma_Q = new TGraphAsymmErrors();
    auto gr_sigma_Q_Trk = new TGraphAsymmErrors();
    auto gr_sigma_QL = new TGraphAsymmErrors();
    auto gr_rms_Q = new TGraphAsymmErrors();
    auto gr_rms_Q_Trk = new TGraphAsymmErrors();
    auto gr_rms_QL = new TGraphAsymmErrors();
    auto gr_sigma_QL_LogL = new TGraphAsymmErrors();
    auto gr_rms_QL_LogL = new TGraphAsymmErrors();
    
    // Michel Trk
    ResolutionSliceLoop(
      hEvsRes_E_Trk, 7.5, 52.5, 9, 1.,
      false, "EQ_Trk", "Michel E_{Q}^{ion}",3,3,-1.2,1.2,
      0, 0.33, false,
      //1, 20, false,
      0, 47.5,
      gr_sigma_Q_Trk,
      gr_rms_Q_Trk);
    
    // Michel shower E_Q
    ResolutionSliceLoop(
      hEvsRes_E_Q, Emin, Emax, 9, 1.,
      true, "EQ", "Michel E_{Q}",3,3,-1.2,1.2,
      //1,-1, useHybridRes,
      0,0.33, useHybridRes,
      -9, -9,
      gr_sigma_Q,
      gr_rms_Q);
    
    // Michel shower E_QL
    ResolutionSliceLoop(
      hEvsRes_E_QL, Emin, Emax, 9, 1.,
      true, "EQL", "Michel E_{QL}",3,3,-1.2,1.2,
      //1,-1, useHybridRes,
      0,0.33, useHybridRes,
      -9, -9,
      gr_sigma_QL,
      gr_rms_QL);
    
    // Michel shower E_QL log-likelihood
    ResolutionSliceLoop(
      hEvsRes_E_QL_LogL, Emin, Emax, 9, 1.,
      true, "EQL_LogL", "Michel E_{QL}^{likelihood}",3,3,-1.2,1.2,
      0,0.33, useHybridRes,
      -9, -9,
      gr_sigma_QL_LogL,
      gr_rms_QL_LogL);
    
   
    // Plot formatting
    gr_sigma_Q    ->SetMarkerColor(kBlue);
    gr_sigma_Q    ->SetMarkerStyle(20);
    gr_sigma_Q    ->SetMarkerSize(1.);
    gr_sigma_Q    ->SetLineColor(kBlue);
    gr_sigma_Q    ->SetLineWidth(1);
    CopyTGraphFormat(gr_sigma_Q, gr_rms_Q);
    gr_rms_Q      ->SetLineStyle(2);
    gr_rms_Q      ->SetMarkerStyle(4);
    
    CopyTGraphFormat(gr_sigma_Q, gr_sigma_QL);
    gr_sigma_QL   ->SetMarkerStyle(21);
    gr_sigma_QL   ->SetMarkerColor(kViolet+2);
    gr_sigma_QL   ->SetLineColor(kViolet+2);
    CopyTGraphFormat(gr_sigma_QL, gr_rms_QL);
    gr_rms_QL      ->SetLineStyle(2);
    gr_rms_QL      ->SetMarkerStyle(25);
    
    CopyTGraphFormat(gr_sigma_QL, gr_sigma_QL_LogL);
    gr_sigma_QL_LogL->SetMarkerStyle(23);
    gr_sigma_QL_LogL->SetMarkerSize(1.5);
    gr_sigma_QL_LogL->SetMarkerColor(kGreen+2);
    gr_sigma_QL_LogL->SetLineColor(kGreen+2);
    CopyTGraphFormat(gr_sigma_QL_LogL, gr_rms_QL_LogL);
    gr_rms_QL_LogL      ->SetLineStyle(2);
    gr_rms_QL_LogL->SetMarkerStyle(32);
    gr_rms_QL_LogL->SetMarkerSize(1.5);
    
    CopyTGraphFormat(gr_sigma_Q, gr_sigma_Q_Trk);
    CopyTGraphFormat(gr_rms_Q, gr_rms_Q_Trk);

    
    // -----------------------------------------------------
    // PLOT
    std::cout
    <<"\n\n\n\n\n"
    <<"**********************************************************************\n"
    <<"Drawing shower resolution for Michel electrons\n";
    TCanvas* cEvsRes  = new TCanvas("EvsRes_michel","EvsRes_michel",600,600);
    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t);
    gPad->SetGrid();
    gStyle->SetGridColor(kGray);


    //gStyle->SetErrorX(0);
    auto mg = new TMultiGraph();
    mg->Add(gr_sigma_Q, "APL");
    mg->Add(gr_sigma_QL, "APL");
    mg->Add(gr_sigma_QL_LogL, "APL");
    mg->Add(gr_rms_Q, "APL");
    mg->Add(gr_rms_QL, "APL");
    mg->Add(gr_rms_QL_LogL, "APL");
    mg->Draw("a");
    mg->GetXaxis()->SetTitle("True energy deposited [MeV]");
    mg->GetYaxis()->SetTitle("Shower energy resolution [%]");
    mg->GetXaxis()->SetTitleSize(axisTitleSize);
    mg->GetXaxis()->SetTitleOffset(1.1);
    mg->GetYaxis()->SetTitleSize(axisTitleSize);
    mg->GetYaxis()->SetTitleOffset(1.3);
    mg->GetYaxis()->SetRangeUser(0,60);
    mg->GetXaxis()->SetRangeUser(0,60);
    //gPad->SetLogy(1);

    TLegend* legg = MakeLegend( 0.6, 1.-mar_t-0.02, textSize, 7,0.33); 
    legg  ->SetBorderSize(1);
    legg  ->SetFillStyle(1001);
    legg  ->AddEntry(gr_sigma_Q,"Q-only #sigma_{p}","PL");
    legg  ->AddEntry(gr_sigma_QL,"Q+L #sigma_{p}","PL");
    legg  ->AddEntry(gr_sigma_QL_LogL,"Q+L #sigma_{p} (likelihood)","PL");
    legg  ->AddEntry(gr_rms_Q,"Q-only RMS","PL");
    legg  ->AddEntry(gr_rms_QL,"Q+L RMS","PL");
    legg  ->AddEntry(gr_rms_QL_LogL,"Q+L RMS (likelihood)","PL");
    legg  ->Draw("same");
    TPaveText* headd = MakeTextBox(mar_l + 0.02, 1.-mar_t-0.02, textSize, 2);
    headd ->AddText("#bf{LArIAT Run IIB MC}");
    headd ->AddText("Cosmic Michel e^{+/-}");
    //headd ->AddText(Form("#tau_{e} = %4.2f ms",0.83));
    //headd ->AddText(Form("#sigma_{PE} = %4.1f%%",hPERes->GetRMS()*100.));
    headd ->Draw();
    
  
     



  ///////////////////////////////
  // Write canvases to file
  outFile   ->cd();
  c_Q_params       ->Write();
  c_L_params ->Write();
  cEvsRes   ->Write();
  ////////////////////////////////






    if( fDoBareElectrons ) {
      
    // --------------------------------------------------------------------
    // 4.2) Bare electron shower resolution plots. Here we will use the bare
    //      electron MC (no muon) for a simplified look at resolution of 
    //      electron showers for different scenarios.
    //
    //      We will assume a PE resolution of 0.1 pe (~ SiPM resolution)
    //      and no additional smearing beyond that.
    //    
    //      9x9 array:
    //          __________________
    //         |                  |
    //      ^  |                  |   
    //      |  |                  | 
    //      |  |                  |
    //     S/N |__________________|
    //          LY  --->
    std::cout
    <<"\n###############################################################\n"
    <<"Resolution slices: bare electrons, nominal\n";
  
    if( fMaxMCEvts_IsoEl > 0 ) fMaxMCEvts = fMaxMCEvts_IsoEl; 

    // Make E_vs_Q and L 2D histograms more finely binned along Y:
    //int E_bins = 50;
    int Eres_bins = 60;
    float Eres_x1   = 2.5; //2.5;
    float Eres_x2   = 62.5; //62.5;
    float res_bins = 320; // 1200
    float res_max = 0.8;
    hEvs_Q  = new TH2D("Evs_Q","Q vs. true Michel energy dep.;True energy deposited [MeV];Reco charge [e^{-}]",Eres_bins,Eres_x1,Eres_x2,200,0.,2000e3);
    hEvs_L  = new TH2D("Evs_L","L vs. true Michel energy dep.;True energy deposited [MeV];Reco light [#gamma]",Eres_bins,Eres_x1,Eres_x2,200,0.,2000e3);
    hEvsRes_E_Trk = new TH2D("EvsRes_E_Trk","Electron ion. track;True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
    hEvsRes_E_Q   = new TH2D("EvsRes_E_Q","Q-only shower energy (uniform R);True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
    hEvsRes_E_QL  = new TH2D("EvsRes_E_QL","Q+L shower energy;True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 
    hEvsRes_E_QL_LogL = new TH2D("EvsRes_E_QL_LogL","Q+L log-likelihood shower energy (uniform R);True energy deposited [MeV];( reco - true ) / true",Eres_bins, Eres_x1, Eres_x2,res_bins, -1.*res_max, res_max); 


    // Lower energy threshold
    Emin = 2.5;
  
    // Switches for using true vs. deposited energy
    // as the reference for resolution
    fUseTrueEnergy  = true;
    fRequireContainment = true;
    fChargeSignFilter = -1;
    
    // Turn off extra smearing (which was due to the 
    // Michel topology)
    /*
    float fSmearFactorDefaults[2];
    float fScaleFactorDefaults[2];
    fSmearFactorDefaults[0] = fSmearFactor[0];
    fSmearFactorDefaults[1] = fSmearFactor[1];
    fScaleFactorDefaults[0] = fScaleFactor[0];
    fScaleFactorDefaults[1] = fScaleFactor[1];
    fSmearFactor[0]=0.;
    fSmearFactor[1]=0.;
    */
   
    // All possible graphs
    //   LY:
    //    1 = 2 pe/MeV
    //    2 = 10 pe/MeV
    //    3 = 20 pe/MeV
    //    4 = 100 pe/MeV
    //   SN:
    //    1 = 10/1
    //    2 = 70/1 (LArIAT)
    //    3 = inf (no noise)
    
    const size_t N_ly = 4;
    const size_t N_sn = 3;

    // LY targets (these will determine what scalings 
    // we need to apply to QE)
    float lytarget[N_ly]={
      2.,
      10.,
      20.,
      100.};  
    
    // signal-to-noise for each sample (approximate)
    std::string sntag[N_sn]={
      "S/N #approx 7:1",
      "S/N #approx 10:1",
      "S/N #approx 70:1"};  
    
    // the files
    std::string mcfilenames[3]={
      "files/MichelAna_mc2_sn7to1.root",
      "files/MichelAna_mc2_sn10to1.root",
//      "files/MichelAna_mc2_electrons.root"};
      "files/MichelAna_mc2_nominal.root"};
    
    // Q,L (nominal)
    TGraphAsymmErrors* grEl_Qdist_mean_nom = new TGraphAsymmErrors;
    TGraphAsymmErrors* grEl_Qdist_rms_nom = new TGraphAsymmErrors;
    TGraphAsymmErrors* grEl_Qdist_sig_nom = new TGraphAsymmErrors;
    TGraphAsymmErrors* grEl_Ldist_mean_nom = new TGraphAsymmErrors;
    TGraphAsymmErrors* grEl_Ldist_rms_nom = new TGraphAsymmErrors;
    TGraphAsymmErrors* grEl_Ldist_sig_nom = new TGraphAsymmErrors;
    // Q,L (LY or SN scenario)
    TGraphAsymmErrors* grEl_Qdist_mean_sn[N_sn];
    TGraphAsymmErrors* grEl_Qdist_rms_sn[N_sn];
    TGraphAsymmErrors* grEl_Qdist_sig_sn[N_sn];
    TGraphAsymmErrors* grEl_Ldist_mean_ly[N_ly];
    TGraphAsymmErrors* grEl_Ldist_rms_ly[N_ly];
    TGraphAsymmErrors* grEl_Ldist_sig_ly[N_ly];

    // Resolution 
    TGraphAsymmErrors* grEl_sig_Q_nom  = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_sig_Q_Trk_nom  = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_sig_QL_nom = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_sig_QL_LogL_nom = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_sig_Q_sn_ly[N_sn][N_ly];
    TGraphAsymmErrors* grEl_sig_Q_Trk_sn[N_sn];
    TGraphAsymmErrors* grEl_sig_QL_sn_ly[N_sn][N_ly];
    TGraphAsymmErrors* grEl_sig_QL_LogL_sn_ly[N_sn][N_ly];
    TGraphAsymmErrors* grEl_rms_Q_nom  = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_Q_Trk_nom  = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_QL_nom = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_QL_LogL_nom = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_Q_sn_ly[N_sn][N_ly];
    TGraphAsymmErrors* grEl_rms_Q_Trk_sn[N_sn];
    TGraphAsymmErrors* grEl_rms_QL_sn_ly[N_sn][N_ly];
    TGraphAsymmErrors* grEl_rms_QL_LogL_sn_ly[N_sn][N_ly];
  
    /*
    TGraphAsymmErrors* grEl_sig_Q_Trk_constrecomb = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_Q_Trk_constrecomb = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_sig_Q_constrecomb = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_sig_QL_constrecomb = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_Q_constrecomb = new TGraphAsymmErrors();
    TGraphAsymmErrors* grEl_rms_QL_constrecomb = new TGraphAsymmErrors();
    */
    
    // ------------------------------------------------
    // Likelihood fitter parameterizations 
    // 1st layer: SN setting
    //  2nd layer: parameter
    //    3rd layer: paramter values
    std::vector< std::vector < float >> vParams_Q_nom = 
    {
      { -1.2e4, 2.77e4},
      { 0.408,   1.107,   0.0174}
    };
    std::vector< std::vector < float >> vParams_L_nom = 
    {
      { -1.2e4,    2.32e4},
      { 0.511, 0.9205, 0.0516 }
    };


    std::vector< std::vector< std::vector< float >>> vParams_Q = 
    { // SN = 7:1
      {
        { -1.52e4, 2.63e4 },        // A
        { 0.395,   0.699, 0.0260 } // B
      }, 
      // SN = 10:1
      {
        { -1.51e4, 2.70e4 },        // A
        { 0.388,   0.856, 0.0258 } // B
      },
      // SN = 70:1
      {
        { -1.21e4, 2.77e4 },        // A
        { 0.429,   1.13,  0.0159 } // B
      }
    };

    std::vector< std::vector< std::vector< float >>> vParams_L = 
    { // LY 0
      {
        { -2.82e4,  2.39e4   },        // A
        { 1.32,    0.706,   0.0199 } // B
      }, 
      // LY 1
      {
        { -1.26e4,   2.36e4  },        // A
        { 0.526,   0.691, 0.0147 } // B
      },
      // LY 2
      {
        { -1.21e4,   2.36e4  },        // A
        {  0.360,  0.631, 0.00764 } // B
      },
      // LY 3 (need to define)
      {
        { -1.25e4,    2.36e4  },        // A
        {  0.291,  0.814, 0.0140 } // B
      }
    };


    // ---------------------------------------
    // Switch to nominal Q and L distributions for 
    // the likelihood fitter
    fFuncP_Q = fFuncP_Q_Gaus;
    fFuncP_L = fFuncP_L_Gaus;
    for(size_t i=0; i<funcP_Q_Gaus_fp.size(); i++){
      for(size_t j=0; j<funcP_Q_Gaus_fp[i]->GetNpar(); j++)
        funcP_Q_Gaus_fp[i]->SetParameter(j,vParams_Q_nom[i][j]);
    }
    for(size_t i=0; i<funcP_L_Gaus_fp.size(); i++){
      for(size_t j=0; j<funcP_L_Gaus_fp[i]->GetNpar(); j++)
        funcP_L_Gaus_fp[i]->SetParameter(j,vParams_L_nom[i][j]);
    }
  

    
    // First find the default LY
    TFile* file = new TFile("files/MichelAna_mc2_nominal.root","read");
    TTree* tree = (TTree*)file->Get("michelana/anatree");
    setBranches(tree);
    SetTrigEffParams(0);
    Loop( tree, 1, true );
    float LY_nom = hTrue_LY->GetMean();
    std::cout<<"Mean LY for bare electron shower MC: "<<LY_nom<<"\n";
    
      ResolutionSliceLoop(
        hEvs_Q, Emin, Emax, 10, 1.,
        true, "Q_e_nom", "Electron Q",3,3,0.,2500e3,
        0, -0.10, false,
        -9,-9,
        grEl_Qdist_sig_nom,
        grEl_Qdist_rms_nom,
        grEl_Qdist_mean_nom);
      ResolutionSliceLoop(
        hEvs_L, Emin, Emax, 10, 1.,
        true, "L_e_nom", "Electron L",3,3,0.,2500e3,
        0, -0.10, false,
        -9,-9,
        grEl_Ldist_sig_nom,
        grEl_Ldist_rms_nom,
        grEl_Ldist_mean_nom);
   
    useHybridRes = true;
    float FitThresh = -0.20;


    // Isolated electrons trk nominal
    ResolutionSliceLoop(
      hEvsRes_E_Trk, 2.5, 42.5, 8, 1.,
      true, "EQ_Trk_e_nom", "Electron E_{Q}^{ion}",4,2,-0.8,0.8,
      0, 0.33, false,
        Emin, 42.5,
      grEl_sig_Q_Trk_nom,
      grEl_rms_Q_Trk_nom);
    
    // Isolated electrons shower E_Q nominal
    ResolutionSliceLoop(
      hEvsRes_E_Q, Emin, Emax, 10, 1.,
      true, "EQ_e_nom", "Electron E_{Q}",3,3,-0.3,0.3,
      0, FitThresh, useHybridRes,
        Emin, 42.5,
      grEl_sig_Q_nom,
      grEl_rms_Q_nom);
    
    // Isolated electrons shower E_QL nominal
    ResolutionSliceLoop(
      hEvsRes_E_QL, Emin, Emax, 10, 1.,
      true, "EQL_e_nom", "Electron E_{QL}",3,3,-0.3,0.3,
      0, FitThresh, useHybridRes,
        Emin, 42.5,
      grEl_sig_QL_nom,
      grEl_rms_QL_nom);
    
    // Isolated electrons shower E_QL LogL nominal
    ResolutionSliceLoop(
      hEvsRes_E_QL_LogL, Emin, Emax, 10, 1.,
      true, "EQL_LogL_e_nom", "Electron E_{QL}^{likelihood}",3,3,-0.3,0.3,
      0, FitThresh, useHybridRes,
        Emin, 42.5,
      grEl_sig_QL_LogL_nom,
      grEl_rms_QL_LogL_nom);
   
    /*
    // Constant recombination E_Q
    file = new TFile("files/MichelAna_mc2_constrecomb.root","read");
    tree = (TTree*)file->Get("michelana/anatree");
    setBranches(tree);
    Loop( tree, 1, true );
    ResolutionSliceLoop(
      hEvsRes_E_Trk, Emin, Emax,
      false, "EQ_Trk_e_constrecomb", "Q-only Track Energy",3,3,1.2,
      0, 0.33, false,
      grEl_sig_Q_Trk_constrecomb,
      grEl_rms_Q_Trk_constrecomb);
    ResolutionSliceLoop(
      hEvsRes_E_Q, Emin, Emax,
      false, "EQ_e_constrecomb", "Q-only Shower Energy",3,3,1.2,
      0, 0.33, false,
      grEl_sig_Q_constrecomb,
      grEl_rms_Q_constrecomb);
    ResolutionSliceLoop(
      hEvsRes_E_QL, Emin, Emax,
      false, "EQL_e_constrecomb", "Q+L Shower Energy",3,3,1.2,
      0, -0.33, true,
      grEl_sig_QL_constrecomb,
      grEl_rms_QL_constrecomb);
    */
    
    // ------------------------------------------------------------
    // Plot overlays of the parameterized and actual distributions
    // for Q and L
    ParamDistOverlay(
      hEvs_Q, 7.5, 52.5, 9, 1., 0., 2500e3,
      "overlays_e_Q_nom", "Electron shower Q", 3, 3, fFuncP_Q);
    ParamDistOverlay(
      hEvs_L, 7.5, 52.5, 9, 1., 0., 2500e3,
      "overlays_e_L_nom", "Electron shower L", 3, 3, fFuncP_L);
    
    // Turn off the LArIAT smearing and assume smearing
    // based on a light collection system capable of achieving
    // 0.1pe resolution (like SiPMs) 
    float fSmearFactorDefaults[2];
    float fScaleFactorDefaults[2];
    fSmearFactorDefaults[0] = fSmearFactor[0];
    fSmearFactorDefaults[1] = fSmearFactor[1];
    fScaleFactorDefaults[0] = fScaleFactor[0];
    fScaleFactorDefaults[1] = fScaleFactor[1];
    fSmearFactor[0] = 0.1;
    fSmearFactor[1] = 0.1;
    fSmearTruePE    = true;
    fSmearMode      = 1;    
    fCorrectQuenching = true;
   
   
    // Now do all hypothetical scenarios 
    for(size_t i_sn = 0; i_sn < N_sn; i_sn++) { 

      std::cout
      <<"\n\n\n"
      <<"Reading in file: "<<mcfilenames[i_sn]<<"\n";
      file = new TFile(mcfilenames[i_sn].c_str(),"read");
      tree = (TTree*)file->Get("michelana/anatree");
      setBranches(tree);
    
      grEl_sig_Q_Trk_sn[i_sn] = new TGraphAsymmErrors();
      grEl_rms_Q_Trk_sn[i_sn] = new TGraphAsymmErrors();
      
      grEl_Qdist_sig_sn[i_sn] = new TGraphAsymmErrors();
      grEl_Qdist_rms_sn[i_sn] = new TGraphAsymmErrors();
      grEl_Qdist_mean_sn[i_sn] = new TGraphAsymmErrors();
      

      // LY loop
      for(size_t i_ly = 0; i_ly < N_ly; i_ly++) { 
        
        // initialize the histograms for this S/N
        grEl_sig_Q_sn_ly[i_sn][i_ly] = new TGraphAsymmErrors();
        grEl_sig_QL_sn_ly[i_sn][i_ly] = new TGraphAsymmErrors();
        grEl_sig_QL_LogL_sn_ly[i_sn][i_ly] = new TGraphAsymmErrors();
        grEl_rms_Q_sn_ly[i_sn][i_ly] = new TGraphAsymmErrors();
        grEl_rms_QL_sn_ly[i_sn][i_ly] = new TGraphAsymmErrors();
        grEl_rms_QL_LogL_sn_ly[i_sn][i_ly] = new TGraphAsymmErrors();
        
        // scale QE factors 
        fScaleFactor[0] = fScaleFactorDefaults[0]*( lytarget[i_ly] / LY_nom);
        fScaleFactor[1] = fScaleFactorDefaults[1]*( lytarget[i_ly] / LY_nom);
       
        // tell us what's up! 
        std::cout
        <<"\n\n"
        <<"-------------------------------------------\n"
        <<sntag[i_sn]<<"\n"   
        <<"Target LY = "<<lytarget[i_ly]<<"\n"
        <<"Nominal LY = "<<LY_nom<<"\n"
        <<"Scale factors: "<<fScaleFactor[0]<<", "<<fScaleFactor[1]<<"\n";
        
        
        // Switch to appropriate Q and L distribution for this combination
        // of noise and LY
        for(size_t i=0; i<funcP_Q_Gaus_fp.size(); i++){
          for(size_t j=0; j<funcP_Q_Gaus_fp[i]->GetNpar(); j++)
            funcP_Q_Gaus_fp[i]->SetParameter(j,vParams_Q[i_sn][i][j]);
        }
        for(size_t i=0; i<funcP_L_Gaus_fp.size(); i++){
          for(size_t j=0; j<funcP_L_Gaus_fp[i]->GetNpar(); j++)
            funcP_L_Gaus_fp[i]->SetParameter(j,vParams_L[i_ly][i][j]);
        }
        
        //if( i_ly == 3 ) fDebugFlag = true;
        //else            fDebugFlag = false;
        
        // loop the tree
        Loop( tree, 1, true );
  
        if( i_ly == 2 ) {
         
          ResolutionSliceLoop(
            hEvsRes_E_Trk, 2.5, 42.5, 8, 1.,
            false, Form("Trk_e_sn%lu_ly%lu",i_sn, i_ly), "Electron E_{Q}^{ion}",4,2,-1.2,1.2,
            0, 0.33, false,
            Emin, 42.5,
            grEl_sig_Q_Trk_sn[i_sn],
            grEl_rms_Q_Trk_sn[i_sn]);
          
          ResolutionSliceLoop(
            hEvs_Q, Emin, Emax, 10, 1.,
            false, Form("Q_e_sn%lu",i_sn), "Electron Q",3,3,-9.,-9.,
            0, -0.10, false,
            -9,-9,
            grEl_Qdist_sig_sn[i_sn],
            grEl_Qdist_rms_sn[i_sn],
            grEl_Qdist_mean_sn[i_sn]);
        }

        if( i_sn == 2 ) {
          grEl_Ldist_sig_ly[i_ly]= new TGraphAsymmErrors();
          grEl_Ldist_rms_ly[i_ly] = new TGraphAsymmErrors();
          grEl_Ldist_mean_ly[i_ly] = new TGraphAsymmErrors();
          ResolutionSliceLoop(
            hEvs_L, Emin, Emax, 10, 1.,
            false, Form("L_e_ly%lu", i_ly), "Electron L",3,3,-9.,-9.,
            0, -0.10, false,
            -9,-9,
            grEl_Ldist_sig_ly[i_ly],
            grEl_Ldist_rms_ly[i_ly],
            grEl_Ldist_mean_ly[i_ly]);
        }
        
        ResolutionSliceLoop(
          hEvsRes_E_Q, Emin, Emax,10,1., 
          true, Form("EQ_e_sn%lu_ly%lu",i_sn, i_ly), "Electron E_{Q}",3,3,-0.8,0.8,
          0, FitThresh, useHybridRes,
            Emin, 42.5,
          grEl_sig_Q_sn_ly[i_sn][i_ly],
          grEl_rms_Q_sn_ly[i_sn][i_ly]);

        ResolutionSliceLoop(
          hEvsRes_E_QL,Emin, Emax,10,1.,
          true, Form("EQL_e_sn%lu_ly%lu",i_sn, i_ly), "Electron E_{QL}",3,3,-0.8,0.8,
          0, FitThresh, useHybridRes,
            Emin, 42.5,
          grEl_sig_QL_sn_ly[i_sn][i_ly],
          grEl_rms_QL_sn_ly[i_sn][i_ly]);
        
        ResolutionSliceLoop(
          hEvsRes_E_QL_LogL,Emin, Emax,10,1.,
          true, Form("EQL_LogL_e_sn%lu_ly%lu",i_sn, i_ly), "Electron E_{QL}^{likelihood}",3,3,-0.8,0.8,
          0, FitThresh, useHybridRes,
            Emin, 42.5,
          grEl_sig_QL_LogL_sn_ly[i_sn][i_ly],
          grEl_rms_QL_LogL_sn_ly[i_sn][i_ly]);

      } // Done with all 3 LY scenarios
      std::cout<<"Done with all "<<N_ly<<" LY scenarios for S/N # "<<i_sn<<"\n";

    }// Done with all 3 S/N scenarios
    std::cout<<"Done with all "<<N_sn<<" S/N scenarios!\n";






    
    // =======================================================
    // First, we need to parameterize Q and L as we did for the 
    // Michel sample
    
    // formatting: 
    // Q mean
    grEl_Qdist_mean_nom ->SetMarkerStyle(20);
    grEl_Qdist_mean_nom ->SetMarkerSize(0.6);
    grEl_Qdist_mean_nom ->SetMarkerColor(kBlack);
    grEl_Qdist_mean_nom ->SetLineColor(kBlack);
    grEl_Qdist_mean_nom ->GetXaxis()->SetNdivisions(505);
    grEl_Qdist_mean_nom ->GetYaxis()->SetNdivisions(505);
    FormatAxes(grEl_Qdist_mean_nom, 0.06,0.05,1.,1.5); 
    grEl_Qdist_mean_nom ->GetXaxis()->SetTitle("Electron energy [MeV]");
    grEl_Qdist_mean_nom ->GetYaxis()->SetTitle("#mu_{Q} [e]");
    // Q sigma
    CopyTGraphFormat( grEl_Qdist_mean_nom, grEl_Qdist_sig_nom,true);
    grEl_Qdist_sig_nom->GetYaxis()->SetTitle("#sigma_{Q}");
    // L mean
    CopyTGraphFormat(grEl_Qdist_mean_nom, grEl_Ldist_mean_nom,true);    
    grEl_Ldist_mean_nom->GetYaxis()->SetTitle("#mu_{L} [#gamma]");
    // L sigma
    CopyTGraphFormat( grEl_Qdist_sig_nom, grEl_Ldist_sig_nom,true);
    grEl_Ldist_sig_nom->GetYaxis()->SetTitle("#sigma_{L}");
   
    
    // Divide the widths by the means to get fractional sigma "rsig"
    
    // nominal cases:
    TGraphAsymmErrors* grEl_Qdist_rsig_nom = GetFracSig(grEl_Qdist_sig_nom, grEl_Qdist_mean_nom, "#sigma_{Q} / #mu_{Q}");
    TGraphAsymmErrors* grEl_Ldist_rsig_nom = GetFracSig(grEl_Ldist_sig_nom, grEl_Ldist_mean_nom, "#sigma_{L} / #mu_{L}");
    
    // Now, for each scenario   
    TGraphAsymmErrors* grEl_Qdist_rsig_sn[N_sn];
    TGraphAsymmErrors* grEl_Ldist_rsig_ly[N_ly];
    for(int i=0; i<N_sn; i++){
      CopyTGraphFormat(grEl_Qdist_mean_nom, grEl_Qdist_mean_sn[i], true);
      CopyTGraphFormat(grEl_Qdist_sig_nom, grEl_Qdist_sig_sn[i], true);
      grEl_Qdist_rsig_sn[i] = GetFracSig(grEl_Qdist_sig_sn[i], grEl_Qdist_mean_sn[i], "#sigma_{Q} / #mu_{Q}" );
    }
    for(int i=0; i<N_ly; i++){
      CopyTGraphFormat(grEl_Ldist_mean_nom, grEl_Ldist_mean_ly[i], true);
      CopyTGraphFormat(grEl_Ldist_sig_nom, grEl_Ldist_sig_ly[i], true);
      grEl_Ldist_rsig_ly[i] = GetFracSig(grEl_Ldist_sig_ly[i], grEl_Ldist_mean_ly[i], "#sigma_{L} / #mu_{L}" );
    }
      
   /* 
    for(size_t i=0; i<3; i++){
      grEl_Qdist_rsig_sn[i] = new TGraphAsymmErrors();
      grEl_Ldist_rsig_ly[i] = new TGraphAsymmErrors();
      //CopyTGraphFormat(grEl_Qdist_sig_sn[i],grEl_Qdist_rsig_sn[i],true);
      //CopyTGraphFormat(grEl_Ldist_sig_ly[i],grEl_Ldist_rsig_ly[i],true);
      for(size_t j=0; j<grEl_Qdist_sig_sn[i]->GetN(); j++){
        double E, dE, sig, sig_err, mean, mean_err, rsig, rsig_err;
        dE = grEl_Qdist_sig_sn[i]->GetErrorX(j);
        
        grEl_Qdist_sig_sn[i]->GetPoint(j,E,sig);
        grEl_Qdist_mean_sn[i]->GetPoint(j,E,mean);
        sig_err = grEl_Qdist_sig_sn[i]->GetErrorY(j);
        mean_err = grEl_Qdist_mean_sn[i]->GetErrorY(j);
        rsig = sig / mean / 100.;
        rsig_err = rsig*std::sqrt( std::pow(sig_err/sig,2) + std::pow(mean_err/mean,2));
        grEl_Qdist_rsig_sn[i]->SetPoint(  grEl_Qdist_rsig_sn[i]->GetN(),E,  rsig ); 
        grEl_Qdist_rsig_sn[i]->SetPointError(  grEl_Qdist_rsig_sn[i]->GetN()-1, dE,dE, rsig_err, rsig_err );

        grEl_Ldist_sig_ly[i]->GetPoint(j,E,sig);
        grEl_Ldist_mean_ly[i]->GetPoint(j,E,mean);
        sig_err = grEl_Ldist_sig_ly[i]->GetErrorY(j);
        mean_err = grEl_Ldist_mean_ly[i]->GetErrorY(j);
        rsig = sig / mean / 100.;
        rsig_err = rsig*std::sqrt( std::pow(sig_err/sig,2) + std::pow(mean_err/mean,2));
        grEl_Ldist_rsig_ly[i]->SetPoint(  grEl_Ldist_rsig_ly[i]->GetN(),E,  rsig ); 
        grEl_Ldist_rsig_ly[i]->SetPointError(  grEl_Ldist_rsig_ly[i]->GetN()-1, dE,dE, rsig_err, rsig_err );
        
        // nom cases
        if( i==0 ) {
          
          grEl_Qdist_sig_nom->GetPoint(j,E,sig);
          grEl_Qdist_mean_nom->GetPoint(j,E,mean);
          sig_err = grEl_Qdist_sig_nom->GetErrorY(j);
          mean_err = grEl_Qdist_mean_nom->GetErrorY(j);
          rsig = sig / mean / 100.;
          rsig_err = rsig*std::sqrt( std::pow(sig_err/sig,2) + std::pow(mean_err/mean,2));
          grEl_Qdist_rsig_nom->SetPoint(  grEl_Qdist_rsig_nom->GetN(),E,  rsig ); 
          grEl_Qdist_rsig_nom->SetPointError(  grEl_Qdist_rsig_nom->GetN()-1, dE,dE, rsig_err, rsig_err );
          
          grEl_Ldist_sig_nom->GetPoint(j,E,sig);
          grEl_Ldist_mean_nom->GetPoint(j,E,mean);
          sig_err = grEl_Ldist_sig_nom->GetErrorY(j);
          mean_err = grEl_Ldist_mean_nom->GetErrorY(j);
          rsig = sig / mean / 100.;
          rsig_err = rsig*std::sqrt( std::pow(sig_err/sig,2) + std::pow(mean_err/mean,2));
          grEl_Ldist_rsig_nom->SetPoint(  grEl_Ldist_rsig_nom->GetN(),E,  rsig ); 
          grEl_Ldist_rsig_nom->SetPointError(  grEl_Ldist_rsig_nom->GetN()-1, dE,dE, rsig_err, rsig_err );
      
          //CopyTGraphFormat(grEl_Qdist_rsig_sn[0],grEl_Qdist_rsig_nom,true);
          //CopyTGraphFormat(grEl_Ldist_rsig_ly[0],grEl_Ldist_rsig_nom,true);
        }
      }
    }
    */
    
    
    // -------------------------------
    // First, nominal case
    std::cout<<"Beginning Q and L dist fist for the nominal case (LArIAT LY and S/N)\n";
    TCanvas* c_el_params_nom = new TCanvas("c_el_params_nom","c_el_params_nom",600,600);
    c_el_params_nom->Divide(2,2);
    //gStyle->SetOptStat(1);
    gStyle->SetOptFit(1111);
      
      // Q
      c_el_params_nom->cd(1);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_mean_nom->Fit(funcP_Q_Gaus_fp[0]);
      grEl_Qdist_mean_nom->DrawClone("AP");
      c_el_params_nom->cd(2);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_rsig_nom->Fit(funcP_Q_Gaus_fp[1]);
      grEl_Qdist_rsig_nom->DrawClone("AP");
      
      // L
      c_el_params_nom->cd(3);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_mean_nom->Fit(funcP_L_Gaus_fp[0]);
      grEl_Ldist_mean_nom->DrawClone("AP");
      c_el_params_nom->cd(4);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_rsig_nom->Fit(funcP_L_Gaus_fp[1]);
      grEl_Ldist_rsig_nom->DrawClone("AP");
      gStyle->SetOptFit(1);
      gPad->Update();

      
    

    // -------------------------------
    // Now scenarios
    //
    // First, Q:
    std::cout<<"Beginning Q dist fit for isolated showers\n";
    TCanvas* c_el_paramsQ = new TCanvas("c_el_paramsQ","c_el_paramsQ",400,600);
    c_el_paramsQ->Divide(2,3);
    //gStyle->SetOptStat(1);
    gStyle->SetOptFit(1111);
    
      // Case 1:
      c_el_paramsQ->cd(1);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_mean_sn[0]->Fit(funcP_Q_Gaus_fp[0]);
      grEl_Qdist_mean_sn[0]->DrawClone("AP");
    
      c_el_paramsQ->cd(2);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_rsig_sn[0]->Fit(funcP_Q_Gaus_fp[1]);
      grEl_Qdist_rsig_sn[0]->DrawClone("AP");
    
      // Case 2:
      c_el_paramsQ->cd(3);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_mean_sn[1]->Fit(funcP_Q_Gaus_fp[0]);
      grEl_Qdist_mean_sn[1]->DrawClone("AP");
      
      c_el_paramsQ->cd(4);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_rsig_sn[1]->Fit(funcP_Q_Gaus_fp[1]);
      grEl_Qdist_rsig_sn[1]->DrawClone("AP");

      // Case 3
      c_el_paramsQ->cd(5);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_mean_sn[2]->Fit(funcP_Q_Gaus_fp[0]);
      grEl_Qdist_mean_sn[2]->DrawClone("AP");

      c_el_paramsQ->cd(6);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Qdist_rsig_sn[2]->Fit(funcP_Q_Gaus_fp[1]);
      grEl_Qdist_rsig_sn[2]->DrawClone("AP");
      gStyle->SetOptFit(1111);
   
   
    
    // -------------------------------
    // then L:
    std::cout<<"Beginning L dist fit for isolated showers\n";
    TCanvas* c_el_paramsL = new TCanvas("c_el_paramsL","c_el_paramsL",400,800);
    c_el_paramsL->Divide(2,4);
    gStyle->SetOptFit(1);
    
      // Case 1
      c_el_paramsL->cd(1);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_mean_ly[0]->Fit(funcP_L_Gaus_fp[0]);
      grEl_Ldist_mean_ly[0]->DrawClone("AP");
    
      c_el_paramsL->cd(2);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_rsig_ly[0]->Fit(funcP_L_Gaus_fp[1]);
      grEl_Ldist_rsig_ly[0]->DrawClone("AP");
    
      // Case 2
      c_el_paramsL->cd(3);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_mean_ly[1]->Fit(funcP_L_Gaus_fp[0]);
      grEl_Ldist_mean_ly[1]->DrawClone("AP");
    
      c_el_paramsL->cd(4);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      grEl_Ldist_rsig_ly[1]->Fit(funcP_L_Gaus_fp[1]);
      grEl_Ldist_rsig_ly[1]->DrawClone("AP");

      // Case 3
      c_el_paramsL->cd(5);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_mean_ly[2]->Fit(funcP_L_Gaus_fp[0]);
      grEl_Ldist_mean_ly[2]->DrawClone("AP");
    
      c_el_paramsL->cd(6);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_rsig_ly[2]->Fit(funcP_L_Gaus_fp[1]);
      grEl_Ldist_rsig_ly[2]->DrawClone("AP");
      gStyle->SetOptFit(1111);
      
      // Case 4
      c_el_paramsL->cd(7);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_mean_ly[3]->Fit(funcP_L_Gaus_fp[0]);
      grEl_Ldist_mean_ly[3]->DrawClone("AP");
    
      c_el_paramsL->cd(8);
      gPad->SetMargin(0.20, 0.05, mar_b, 0.10);
      //gStyle->SetOptStat(1);
      gStyle->SetOptFit(1111);
      grEl_Ldist_rsig_ly[3]->Fit(funcP_L_Gaus_fp[1]);
      grEl_Ldist_rsig_ly[3]->DrawClone("AP");
      gStyle->SetOptFit(1111);


  
    
    
    
    
    
    
    
    
    
    // -----------------------------------------------------
    // B) Track resolutions at all 3 S/N
    // -----------------------------------------------------
    std::cout
    <<"\n\n\n\n\n"
    <<"**********************************************************************\n"
    <<"Drawing electron ionization-only (track) resolutions for different S/N\n";
    TCanvas* cEvsRes_trk  = new TCanvas("EvsRes_trk","EvsRes_trk",500,500);
    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t);
    gStyle->SetOptFit(0);
  
    FormatTGraph(grEl_sig_Q_Trk_sn[0], kAzure+5, kAzure+5,  20, 1, 0.8); // 10-1 
    //FormatTGraph(grEl_sig_Q_Trk_sn[1], kBlue, kBlue, 22, 1, 0.8); // 70-1
    //FormatTGraph(grEl_sig_Q_Trk_sn[2], kBlue+2, kBlue+2, 21, 1, 0.8); // inf
    FormatTGraph(grEl_sig_Q_Trk_sn[1], kBlue, kBlue,        20, 1, 0.8); // 70-1
    FormatTGraph(grEl_sig_Q_Trk_sn[2], kBlue+3, kBlue+3,    20, 1, 0.8); // inf
   
    TF1* fResFitTrk[3];
    fResFitTrk[0] = (TF1*)fResFit->Clone("resFitTrk_0");
    fResFitTrk[1] = (TF1*)fResFit->Clone("resFitTrk_1");
    fResFitTrk[2] = (TF1*)fResFit->Clone("resFitTrk_2");
   
    fResFitTrk[0]->SetLineColor( grEl_sig_Q_Trk_sn[0]->GetLineColor()); //fResFitTrk[0]->DrawCopy("same"); 
    fResFitTrk[1]->SetLineColor( grEl_sig_Q_Trk_sn[1]->GetLineColor()); //fResFitTrk[1]->DrawCopy("same"); 
    fResFitTrk[2]->SetLineColor( grEl_sig_Q_Trk_sn[2]->GetLineColor()); //fResFitTrk[2]->DrawCopy("same");
    fResFitTrk[0]->SetLineWidth(2);
    fResFitTrk[1]->SetLineWidth(2);
    fResFitTrk[2]->SetLineWidth(2);
    fResFitTrk[0]->SetLineStyle(1);
    fResFitTrk[1]->SetLineStyle(2);
    fResFitTrk[2]->SetLineStyle(7);
    grEl_sig_Q_Trk_sn[0]->Fit(fResFitTrk[0]);
    grEl_sig_Q_Trk_sn[1]->Fit(fResFitTrk[1]);
    grEl_sig_Q_Trk_sn[2]->Fit(fResFitTrk[2]);
    
    //gStyle->SetErrorX(0);
    gPad->SetGrid();
    gStyle->SetGridColor(kGray);
    auto mg_trk = new TMultiGraph();
    mg_trk->Add(grEl_sig_Q_Trk_sn[0],"AP");
    mg_trk->Add(grEl_sig_Q_Trk_sn[1],"AP");
    mg_trk->Add(grEl_sig_Q_Trk_sn[2],"AP");
    mg_trk->Draw("a");
    mg_trk->GetXaxis()->SetTitle("True electron ionization energy deposited [MeV]");
    mg_trk->GetYaxis()->SetTitle("Ionization energy resolution [%]");
    FormatAxes(mg_trk, axisTitleSize, axisLabelSize,1.1,1.3);  
    mg_trk->GetYaxis()->SetRangeUser(0,15);
   
    TLegend* lg_trk = MakeLegend( 0.6, 1.-mar_t-0.02, textSize, 4, 0.33); 
    lg_trk  ->SetBorderSize(1);
    lg_trk  ->SetFillStyle(1001);
    lg_trk  ->SetNColumns(2);
    lg_trk  ->AddEntry(grEl_sig_Q_Trk_sn[0],"S/N ~ 7:1", "P");
    lg_trk  ->AddEntry(fResFitTrk[0],       "Fit",        "L");
    lg_trk  ->AddEntry(grEl_sig_Q_Trk_sn[1],"S/N ~ 10:1","P");
    lg_trk  ->AddEntry(fResFitTrk[1],       "Fit",        "L");
    lg_trk  ->AddEntry(grEl_sig_Q_Trk_sn[2],"S/N ~ 70:1","P");
    lg_trk  ->AddEntry(fResFitTrk[2],       "Fit",        "L");
    lg_trk  ->Draw("same");
    TPaveText* hd_trk = MakeTextBox(mar_l + 0.02, 1.-mar_t-0.02, textSize, 2);
    hd_trk ->AddText("#bf{LArIAT MC}");
    hd_trk ->AddText("Isolated Electrons");
    hd_trk ->Draw();
    
    // Display the fit results
    float y0 = 1.-mar_t-0.02-4*textSize-0.02;
    TPaveText* hd_fittrk =  MakeTextBox(0.68, y0, 0.03, 10, 0.25 );
    hd_fittrk -> AddText("S/N ~ 7:1 fit:")->SetTextColor( fResFitTrk[0]->GetLineColor() );
      hd_fittrk -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitTrk[0]->GetParameter(0), fResFitTrk[0]->GetParError(0)))->SetTextColor( fResFitTrk[0]->GetLineColor() );
      hd_fittrk -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitTrk[0]->GetParameter(1), fResFitTrk[0]->GetParError(1)))->SetTextColor( fResFitTrk[0]->GetLineColor() );
    hd_fittrk -> AddText("S/N ~ 10:1 fit:")->SetTextColor( fResFitTrk[1]->GetLineColor() );
      hd_fittrk -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitTrk[1]->GetParameter(0), fResFitTrk[1]->GetParError(0)))->SetTextColor( fResFitTrk[1]->GetLineColor() );
      hd_fittrk -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitTrk[1]->GetParameter(1), fResFitTrk[1]->GetParError(1)))->SetTextColor( fResFitTrk[1]->GetLineColor() );
    hd_fittrk -> AddText("S/N ~ 70:1 fit:")->SetTextColor( fResFitTrk[2]->GetLineColor() );
      hd_fittrk -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitTrk[2]->GetParameter(0), fResFitTrk[2]->GetParError(0)))->SetTextColor( fResFitTrk[2]->GetLineColor() );
      hd_fittrk -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitTrk[2]->GetParameter(1), fResFitTrk[2]->GetParError(1)))->SetTextColor( fResFitTrk[2]->GetLineColor() );
    hd_fittrk->Draw();    
    
    
    
    
    
    
    // -----------------------------------------------------
    // C) Nominal shower resolution 
    // -----------------------------------------------------
    std::cout
    <<"\n\n\n\n\n"
    <<"****************************************************************\n"
    <<"Drawing nominal shower resolution (bare electrons)\n";
    TCanvas* cEvsRes_shwr  = new TCanvas("EvsRes_shwr","EvsRes_shwr",500,500);
    gPad->SetMargin(mar_l, mar_r, mar_b, mar_t);
    gPad->SetGrid();
    gStyle->SetGridColor(kGray);
    gStyle->SetOptStat(0);
    
    TF1* fResFitQ = (TF1*)fResFit->Clone("resFitQ");
    fResFitQ->SetLineColor( kBlue );
    fResFitQ->SetLineWidth(2);
    TF1* fResFitQL = (TF1*)fResFit->Clone("resFitQL");
    fResFitQL->SetLineColor( kViolet+2 );
    fResFitQL->SetLineWidth(2);
    TF1* fResFitQL_LogL = (TF1*)fResFitQL->Clone("resFitQL_LogL");
    fResFitQL_LogL->SetLineColor( kGreen+2 );
    fResFitQL_LogL->SetLineWidth(2);
      
      // Fits
      grEl_sig_Q_nom->Fit(fResFitQ);
      //fResFitQ->DrawCopy("same");
      std::cout
      <<"   Q-only fit:\n"
      <<"     Chi2/ndf = "<<fResFitQ->GetChisquare()/fResFitQ->GetNDF()<<"\n"
      <<"     A = "<<fResFitQ->GetParameter(0)<<" +/- "<< fResFitQ->GetParError(0)<<"\n"
      <<"     B = "<<fResFitQ->GetParameter(1)<<" +/- "<< fResFitQ->GetParError(1)<<"\n";
      grEl_sig_QL_nom->Fit(fResFitQL);
      //fResFitQL->DrawCopy("same");
      std::cout
      <<"   Q+L fit:\n"
      <<"     Chi2/ndf = "<<fResFitQL->GetChisquare()/fResFitQL->GetNDF()<<"\n"
      <<"     A = "<<fResFitQL->GetParameter(0)<<" +/- "<< fResFitQL->GetParError(0)<<"\n"
      <<"     B = "<<fResFitQL->GetParameter(1)<<" +/- "<< fResFitQL->GetParError(1)<<"\n";
      grEl_sig_QL_LogL_nom->Fit(fResFitQL_LogL);
      //fResFitQL->DrawCopy("same");
      std::cout
      <<"   Q+L fit:\n"
      <<"     Chi2/ndf = "<<fResFitQL_LogL->GetChisquare()/fResFitQL_LogL->GetNDF()<<"\n"
      <<"     A = "<<fResFitQL_LogL->GetParameter(0)<<" +/- "<< fResFitQL_LogL->GetParError(0)<<"\n"
      <<"     B = "<<fResFitQL_LogL->GetParameter(1)<<" +/- "<< fResFitQL_LogL->GetParError(1)<<"\n";
    
    FormatTGraph(grEl_sig_Q_nom,        kBlue,      kBlue,      20, 1,  0.8, 2);
    FormatTGraph(grEl_sig_QL_nom,       kViolet+2,  kViolet+2,  21, 1,  0.8, 2);
    FormatTGraph(grEl_sig_QL_LogL_nom,  kGreen+2,   kGreen+2,   23, 1,  1.2, 2);
    FormatTGraph(grEl_rms_Q_nom,        kBlue,      kBlue,      20,  2,  0.8, 2);
    FormatTGraph(grEl_rms_QL_nom,       kViolet+2,  kViolet+2,  21,  2,  0.8, 2);
    FormatTGraph(grEl_rms_QL_LogL_nom,  kGreen+2,   kGreen+2,   23,  2,  1.2, 2);
   
    auto mg_shwr = new TMultiGraph();
    mg_shwr->Add(grEl_sig_Q_nom,"AP");
    mg_shwr->Add(grEl_sig_QL_nom,"AP");
    mg_shwr->Add(grEl_sig_QL_LogL_nom,"AP");
//    mg_shwr->Add(grEl_rms_Q_nom,"APL");
//    mg_shwr->Add(grEl_rms_QL_nom,"APL");
    mg_shwr->Draw("a");
    mg_shwr->GetXaxis()->SetTitle("True electron energy [MeV]");
    mg_shwr->GetYaxis()->SetTitle("Shower energy resolution [%]");
    FormatAxes(mg_shwr, axisTitleSize, axisLabelSize,1.1,1.3); 
    mg_shwr->GetYaxis()->SetRangeUser(0,8);
    //gPad->SetLogy();
    
    TPaveText* hd_shwr = MakeTextBox(mar_l + 0.02, 1.-mar_t-0.02, textSize, 2);
    hd_shwr ->AddText("#bf{LArIAT MC}");
    hd_shwr ->AddText("Isolated Electrons");
//    hd_shwr ->AddText("Electrons");
    hd_shwr ->Draw();
      
    
    TLegend* lg_shwr = MakeLegend( 0.55, 1.-mar_t-0.02, textSize, 4, 0.32 ); 
    lg_shwr  ->SetBorderSize(1);
    lg_shwr  ->SetFillStyle(1001);
    lg_shwr ->SetNColumns(2);
    lg_shwr  ->AddEntry(grEl_sig_Q_nom, "Q-only","P");
    lg_shwr  ->AddEntry(fResFitQ,          "Fit",          "L" );
    lg_shwr  ->AddEntry(grEl_sig_QL_nom, "Q+L","P");
    lg_shwr  ->AddEntry(fResFitQL,          "Fit",          "L" );
    lg_shwr  ->AddEntry(grEl_sig_QL_LogL_nom, "Q+L (likelihood)","P");
    lg_shwr  ->AddEntry(fResFitQL_LogL,          "Fit",          "L" );
  //  lg_shwr  ->AddEntry(grEl_rms_Q_nom, "Q-only RMS","PL");
  //  lg_shwr  ->AddEntry((TObject*)0, "","");
  //  lg_shwr  ->AddEntry(grEl_rms_QL_nom, "Q+L RMS","PL");
  //  lg_shwr  ->AddEntry((TObject*)0, "","");
    lg_shwr  ->Draw("same");

    //gStyle->SetErrorX(0);
      
    // Display the fit results
    y0 = 1.-mar_t-0.02-4*textSize-0.02;
    TPaveText* hd_fitnom =  MakeTextBox(0.68, y0, 0.03, 10, 0.28 );
    hd_fitnom -> AddText("Q-only fit:")->SetTextColor(kBlue);
      hd_fitnom -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitQ->GetParameter(0), fResFitQ->GetParError(0)))->SetTextColor(kBlue);
      hd_fitnom -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitQ->GetParameter(1), fResFitQ->GetParError(1)))->SetTextColor(kBlue);
    hd_fitnom -> AddText("Q+L fit:")->SetTextColor(kViolet+2);
      hd_fitnom -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitQL->GetParameter(0), fResFitQL->GetParError(0)))->SetTextColor(kViolet+2);
      hd_fitnom -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitQL->GetParameter(1), fResFitQL->GetParError(1)))->SetTextColor(kViolet+2);
    hd_fitnom -> AddText("Q+L likelihood fit:")->SetTextColor(kGreen+2);
      hd_fitnom -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitQL_LogL->GetParameter(0), fResFitQL_LogL->GetParError(0)))->SetTextColor(kGreen+2);
      hd_fitnom -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitQL_LogL->GetParameter(1), fResFitQL_LogL->GetParError(1)))->SetTextColor(kGreen+2);
    hd_fitnom -> Draw();
  
  
   
   
    // ----------------------------------------------------
    // C2)  Plot fractional improvement going from Q-only to
    //      Q+L likelihood for the 4 different light yields.
    std::cout
    <<"============================================\n"
    <<"Making frac improvement plot\n";
    TCanvas* cFracImprove = new TCanvas("FracImprove","FracImprove",600,600);
    TGraphAsymmErrors* g_fracimprov[N_ly];
    for(size_t i=0; i<N_ly; i++){
      g_fracimprov[i] = new TGraphAsymmErrors();
      g_fracimprov[i] ->SetTitle(Form("%4.0f pe/MeV",lytarget[i]));
      //g_fracimprov[i] ->GetXaxis()->SetTitle("Electron energy [MeV]");
      //g_fracimprov[i] ->GetYaxis()->SetTitle("Improvement in resolution [%]");
      for(size_t ipt=0; ipt<grEl_sig_Q_sn_ly[1][i]->GetN(); ipt++){
        double E, dE;
        double sig_ref, dsig_ref, sig, dsig;
        grEl_sig_Q_sn_ly[1][i]->GetPoint(ipt,E,sig_ref);
        dE        = grEl_sig_Q_sn_ly[1][i]->GetErrorX(ipt);
        dsig_ref  = grEl_sig_Q_sn_ly[1][i]->GetErrorY(ipt);
        grEl_sig_QL_LogL_sn_ly[1][i]->GetPoint(ipt,E,sig);
        dsig      = grEl_sig_QL_LogL_sn_ly[1][i]->GetErrorY(ipt);
        double f = -100.*( sig/sig_ref - 1. );
        double df = 100.*(sig/sig_ref)*std::sqrt( std::pow(dsig/sig,2) + std::pow(dsig_ref/sig_ref,2) );
        
        //double a = sig-sig_ref;
        //double da = std::sqrt( std::pow(dsig,2) + std::pow(dsig_ref,2) );
        //double b = sig_ref;
        //double db = dsig_ref;
        //double f = -100.*(a/b); // negative change = positive improvement
        //double df = fabs(f)*std::sqrt( std::pow( da/a,2 ) + std::pow( db/b,2 ) );
        int n = g_fracimprov[i]->GetN();
        g_fracimprov[i]   ->SetPoint( n, E, f);
        g_fracimprov[i]   ->SetPointError(n, dE,dE,df,df);
        std::cout
        <<"  pt "<<n<<"  E= "<<E<<"  dE= "<<dE<<"   f = "<<f<<"  df= "<<df<<"\n";
      }
    }

    // Formatting
    g_fracimprov[0] -> SetMarkerStyle(3);
    g_fracimprov[0] -> SetLineWidth(2);
    g_fracimprov[0] -> SetMarkerColor(kBlack);
    g_fracimprov[0] -> SetLineColor(kBlack);
    CopyTGraphFormat( g_fracimprov[0], g_fracimprov[1]);
    CopyTGraphFormat( g_fracimprov[0], g_fracimprov[2]);
    CopyTGraphFormat( g_fracimprov[0], g_fracimprov[3]);
    
    g_fracimprov[0] -> SetLineStyle(3);
    g_fracimprov[1] -> SetLineStyle(4);
    g_fracimprov[2] -> SetLineStyle(2);
    g_fracimprov[3] -> SetLineStyle(1);
   
    gPad->SetMargin(0.15,0.05,0.15,0.05); 
    gPad->SetGrid();
    gStyle->SetGridColor(kGray);

    TMultiGraph* mg_fracimprov = new TMultiGraph();
    mg_fracimprov->Add( g_fracimprov[0], "APL");
    mg_fracimprov->Add( g_fracimprov[1], "APL");
    mg_fracimprov->Add( g_fracimprov[2], "APL");
    mg_fracimprov->Add( g_fracimprov[3], "APL");
    mg_fracimprov->Draw("a");
    mg_fracimprov->GetXaxis()->SetTitle("True electron energy [MeV]");
    mg_fracimprov->GetYaxis()->SetTitle("Improvement in resolution [%]");
    FormatAxes(mg_fracimprov, 0.045, 0.04, 1.1, 1.3);

    TLegend* lg_fracimprov = MakeLegend( 0.5, 0.88, 0.04, 4, 0.38 ); 
    lg_fracimprov ->SetBorderSize(1);
    lg_fracimprov ->SetFillStyle(1001);
    lg_fracimprov ->AddEntry( g_fracimprov[0], g_fracimprov[0]->GetTitle(), "PL");
    lg_fracimprov ->AddEntry( g_fracimprov[1], g_fracimprov[1]->GetTitle(), "PL");
    lg_fracimprov ->AddEntry( g_fracimprov[2], g_fracimprov[2]->GetTitle(), "PL");
    lg_fracimprov ->AddEntry( g_fracimprov[3], g_fracimprov[3]->GetTitle(), "PL");
    lg_fracimprov ->Draw("same");
   
   
   
   
   
   
   
   
   
     
    
    // -----------------------------------------------------
    // D) The big honkin' array
    // -----------------------------------------------------
    std::cout
    <<"\n\n\n\n\n"
    <<"****************************************************************\n"
    <<"Drawing big honkin' array of resolution scenarios\n";
    TCanvas* cEvsRes_array = new TCanvas("EvsRes_array","EvsRes_array",800,800);
    cEvsRes_array->Divide(3,3,0.005,0.005,0);
   
    axisTitleSize = 0.05;
    axisLabelSize = 0.05;
    textSize      = 0.04;
    
    float legend_x1 = 0.58;
    float legend_width = 1.-mar_r-legend_x1-0.02;

    int sn_index[9]={2, 2, 2, 1, 1, 1, 0, 0, 0}; 
    int ly_index[9]={0, 1, 2, 0, 1, 2, 0, 1, 2}; 
    //gStyle->SetErrorX(0);
    TMultiGraph* mg_array[9];
    TPaveText* hd_array[9];
    TLegend* lg_array[9]; 
    TPaveText* hd_fit_q[9];
    TPaveText* hd_fit_ql[9];
    TPaveText* hd_fit_ql_logl[9];

    for(size_t i=0; i<9; i++){
      
      int ily = ly_index[i];
      int isn = sn_index[i];

      std::cout
      <<"\n"
      <<"-------------------------------------------------------\n"
      <<"Drawing resolution for LY = "<<lytarget[ily]<< ", "<<sntag[isn]<<"\n";
     
      if( grEl_sig_Q_sn_ly[isn][ily]->GetN() == 0 ) continue;
      
    
      TGraphAsymmErrors** grQ         = &grEl_sig_Q_sn_ly[isn][ily];
      TGraphAsymmErrors** grQL        = &grEl_sig_QL_sn_ly[isn][ily];
      TGraphAsymmErrors** grQL_LogL   = &grEl_sig_QL_LogL_sn_ly[isn][ily];
    
      // Initialize with guesses
      fResFitQ->SetParameter(0, 18.);
      fResFitQ->SetParameter(1, 3.);
      fResFitQL->SetParameter(0, 20.);
      fResFitQL->SetParameter(1, 1.);
      fResFitQL_LogL->SetParameter(0, 20.);
      fResFitQL_LogL->SetParameter(1, 1.);
     
      
      std::cout
      <<"\n"
      <<"   ..................................................\n"
      <<"   Fits to resolution (using Gaussian fit)\n"
      <<"   ...................................................\n";
      (*grQ)->Fit(fResFitQ);
      std::cout
      <<"   Q-only fit:\n"
      <<"     Chi2/ndf = "<<fResFitQ->GetChisquare()/fResFitQ->GetNDF()<<"\n"
      <<"     A = "<<fResFitQ->GetParameter(0)<<" +/- "<< fResFitQ->GetParError(0)<<"\n"
      <<"     B = "<<fResFitQ->GetParameter(1)<<" +/- "<< fResFitQ->GetParError(1)<<"\n";
      (*grQL)->Fit(fResFitQL);
      std::cout
      <<"   Q+L fit:\n"
      <<"     Chi2/ndf = "<<fResFitQL->GetChisquare()/fResFitQL->GetNDF()<<"\n"
      <<"     A = "<<fResFitQL->GetParameter(0)<<" +/- "<< fResFitQL->GetParError(0)<<"\n"
      <<"     B = "<<fResFitQL->GetParameter(1)<<" +/- "<< fResFitQL->GetParError(1)<<"\n";
      (*grQL_LogL)->Fit(fResFitQL_LogL);
      std::cout
      <<"   Q+L fit:\n"
      <<"     Chi2/ndf = "<<fResFitQL_LogL->GetChisquare()/fResFitQL_LogL->GetNDF()<<"\n"
      <<"     A = "<<fResFitQL_LogL->GetParameter(0)<<" +/- "<< fResFitQL_LogL->GetParError(0)<<"\n"
      <<"     B = "<<fResFitQL_LogL->GetParameter(1)<<" +/- "<< fResFitQL_LogL->GetParError(1)<<"\n"
      <<"   ...................................................\n\n";

      
      FormatTGraph((*grQ),        kBlue,      kBlue,      20, 1,  0.5, 2);
      FormatTGraph((*grQL),       kViolet+2,  kViolet+2,  21, 1,  0.5, 2);
      FormatTGraph((*grQL_LogL),  kGreen+2,   kGreen+2,   23, 1,  0.8, 2);
     
      float textSizeLeg = 0.03;
      float heightLeg = 3;
      float marginTopLeg = 0.02;
      float heightFit = 3.;
      float sepFit = 0.01;
      float fitX = legend_x1 + 0.09;
      
      // Plot fit curves
      cEvsRes_array->cd(i+1);
      gPad->SetMargin(mar_l, mar_r, mar_b, mar_t);
      //gStyle->SetErrorX(0);
      gStyle->SetOptFit(0);

      mg_array[i] = new TMultiGraph();
      mg_array[i]->Add((*grQ), "AP");
      //mg_array[i]->Add((*grQL),"AP");
      mg_array[i]->Add((*grQL_LogL),"AP");
      mg_array[i]->Draw("a");
      mg_array[i]->GetXaxis()->SetTitle("True electron energy [MeV]");
      mg_array[i]->GetYaxis()->SetTitle("Shower energy resolution [%]");
      mg_array[i]->GetYaxis()->SetRangeUser(0.,17.5);
      //mg_array[i]->GetYaxis()->SetRangeUser(1,80);
      //gPad->SetLogy();
      FormatAxes(mg_array[i], axisTitleSize, axisLabelSize,1.1,1.3);  
      
      gPad->SetGrid();
      gStyle->SetGridColor(kGray);
      gPad->Update();
      gPad->RedrawAxis();
      
      hd_array[i] = MakeTextBox(mar_l + 0.02, 1.-mar_t-0.02, textSize, 3.5, 0.35);
//      hd_array[i] ->AddText("#bf{LArIAT MC}");
      hd_array[i] ->AddText("Isolated Electrons");
      hd_array[i]  ->AddText(Form("LY = %2.0f pe/MeV",lytarget[ly_index[i]]));
      hd_array[i]  ->AddText(Form("%s",sntag[sn_index[i]].c_str()));
      //hd_array[i]  ->SetFillStyle(1001);
      //hd_array[i] ->SetFillColor(kWhite);
//      hd_array[i] ->SetBorderSize(1);
//      hd_array[i]  ->SetFillColorAlpha(kWhite,0.60);
      hd_array[i] ->Draw();

      lg_array[i] = MakeLegend( legend_x1, 1.-mar_t-marginTopLeg, textSizeLeg, heightLeg, legend_width); 
      lg_array[i]  ->SetBorderSize(1);
      lg_array[i]  ->SetFillStyle(1001);
      lg_array[i]  ->SetNColumns(2);
      lg_array[i]  ->AddEntry((*grQ),  "Q-only",  "P");
      lg_array[i]  ->AddEntry(fResFitQ,            "Fit",              "L");
      //lg_array[i]  ->AddEntry((*grQL),  "Q+L",  "P");
      //lg_array[i]  ->AddEntry(fResFitQL,           "Fit",              "L");
      lg_array[i]  ->AddEntry((*grQL_LogL),  "Q+L (likelihood)",  "P");
      lg_array[i]  ->AddEntry(fResFitQL_LogL,           "Fit",              "L");
      lg_array[i]  ->Draw("same");
      

      // Make 3 Text Boxes for the fit info

      hd_fit_q[i] =  MakeTextBox(fitX, 1.-mar_t-marginTopLeg-textSizeLeg*heightLeg-0.01, textSizeLeg, heightFit, legend_width);
      hd_fit_q[i] -> AddText("Q-only fit:")->SetTextColor(kBlue);
      hd_fit_q[i] -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitQ->GetParameter(0), fResFitQ->GetParError(0)))->SetTextColor(kBlue);
      hd_fit_q[i] -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitQ->GetParameter(1), fResFitQ->GetParError(1)))->SetTextColor(kBlue);
      hd_fit_q[i] -> Draw();
      
      /*
      hd_fit_ql[i] = MakeTextBox(fitX, 1.-mar_t-marginTopLeg-textSizeLeg*heightLeg-0.01-textSizeLeg*heightFit-sepFit, textSizeLeg, heightFit, legend_width);
      hd_fit_ql[i] -> AddText("Q+L fit:")->SetTextColor(kViolet+2);
      hd_fit_ql[i] -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitQL->GetParameter(0), fResFitQL->GetParError(0)))->SetTextColor(kViolet+2);
      hd_fit_ql[i] -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitQL->GetParameter(1), fResFitQL->GetParError(1)))->SetTextColor(kViolet+2);
      hd_fit_ql[i] -> Draw();
      */

      //hd_fit_ql_logl[i] = MakeTextBox(fitX, 1.-mar_t-marginTopLeg-textSizeLeg*heightLeg-0.01-2.*textSizeLeg*heightFit-2*sepFit, textSizeLeg, heightFit, legend_width);
      hd_fit_ql_logl[i] = MakeTextBox(fitX, 1.-mar_t-marginTopLeg-textSizeLeg*heightLeg-0.01-textSizeLeg*heightFit-sepFit, textSizeLeg, heightFit, legend_width);
      hd_fit_ql_logl[i] -> AddText("Q+L likelihood fit:")->SetTextColor(kGreen+2);
      hd_fit_ql_logl[i] -> AddText(Form("A = %5.2f #pm %5.2f %%", fResFitQL_LogL->GetParameter(0), fResFitQL_LogL->GetParError(0)))->SetTextColor(kGreen+2);
      hd_fit_ql_logl[i] -> AddText(Form("B = %5.2f #pm %5.2f %%", fResFitQL_LogL->GetParameter(1), fResFitQL_LogL->GetParError(1)))->SetTextColor(kGreen+2);
      hd_fit_ql_logl[i] -> Draw();
      
      // Lock the pad so formatting doesn't change
      gPad->SetEditable(false);
      
    }








    

  








  //////////////////////////////
  // Write to file 
  outFile       ->cd();
  c_el_params_nom ->Write();
  c_el_paramsQ ->Write();
  c_el_paramsL ->Write();
  cEvsRes_trk->Write();
  cEvsRes_shwr->Write();
  cEvsRes_array->Write();
  cFracImprove->Write();
  /////////////////////////////

  } // do isolated electron showers

  } // End resolution slices
 
  //gStyle->SetOptFit(1111);
  fSmearTruePE=false;

  outFile->Close();

  double duration = (time(0) - startTime)/60.;
  printf("\n\n\n");
  printf("======================================\n");
  printf("Done with making calorimetry plots\n");
  printf("Total time: %6.2f minutes\n",(time(0)-startTime)/60.);
  printf("======================================\n");
   
}












//##########################################################################
// Recombination plots
// (WORK IN PROGRESS)
void RecombPlots(){
  
  float mar_l  = 0.15;
  float mar_r  = 0.10;
  float mar_t  = 0.05;
  float leg_x1 = 0.55;
  float leg_x2 = 1.-mar_r-0.02;
  float leg_y1 = 0.7;
  float leg_y2 = 1.-mar_t-0.02;
  float axisTitleSize = 0.045;
  float textSize = 0.035;
  
  // First use MC to figure out what scale factors to apply to Q and L
  Loop(1);
  fQCorrFactor = 1./(1.+hQRes->GetMean());
  fLCorrFactor = 1./(1.+hLRes->GetMean());
  std::cout<<"QCorrFactor = "<<fQCorrFactor<<"\n";
  std::cout<<"LCorrFactor = "<<fLCorrFactor<<"\n";

  TCanvas* cr = new TCanvas("recomb","recomb",600,600);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);
  gStyle->SetOptStat(0);
  
  float dataMax = GetHistMax(hR[0]);
  float mcMax   = GetHistMax(hR[1]);
  float max     = std::max(dataMax, mcMax);
  
//  hR[1]->Scale(hR[0]->GetEntries() / hR[1]->Integral() );
  hR[1]->SetLineWidth(0);
  hR[1]->SetFillColor(kGray);
  hR[1]->GetXaxis()->SetTitleSize(axisTitleSize);
  hR[1]->GetYaxis()->SetTitleSize(axisTitleSize);
  hR[1]->SetMaximum( max*1.2);
  hR[1]->Draw("hist");
  hR[0]->SetMarkerStyle(20);
  hR[0]->SetMarkerColor(kBlack);
  hR[0]->SetLineColor(kBlack);
  hR[0]->SetMarkerSize(0.8);
  hR[0]->Draw("same P E X0");
  hR[1]->Draw("same axis");
  
  TPaveText* header = MakeTextBox(mar_l + 0.02, leg_y2, textSize, 3);
    header->AddText("#bf{LArIAT Preliminary}");
    header->AddText(Form("%s Dataset",runtag.c_str()));
    header->AddText(Form("%d Events",(int)hR[0]->GetEntries()));
    header->Draw();

  TPaveText* cutsr = MakeTextBox(leg_x1,leg_y2, textSize, 3);
    cutsr->AddText("Cuts: ");
    cutsr->AddText(Form("  3D shower (Npts #geq %d)",fMinNumPts3D));
    cutsr->AddText(Form("  #DeltaT > %3.1f #mus",fdTcut/1000.));
//    cutsr->AddText("  E_{Q} < 55 MeV");
    cutsr->Draw();
  
  TLegend* leg = MakeLegend( leg_x1+0., leg_y2-4.*textSize, textSize, 2); 
    leg  ->AddEntry(hR[1], Form("MC reco: #it{#bar{r}} = %5.3f",hR[1]->GetMean()), "F"); 
    leg ->AddEntry(hR[0],  Form("Data: #it{#bar{r}} = %5.3f",hR[0]->GetMean()),"EPL");
    leg->Draw();

  
}

//################################################################################
// Plot the decay time and perform fit
void TimePlots(){

  Loop(0);

  //if( hDecayTime[0]->GetEntries() == 0 || hDecayTime[1]->GetEntries() == 0 ) {
  //  std::cout<<"ERROR: one of the histograms is empty.\n";
  //  return;
  //}
  
  
  float mar_l  = 0.15;
  float mar_r  = 0.10;
  float mar_t  = 0.05;
  float leg_x1 = 0.45;
  float leg_x2 = 1.-mar_r-0.02;
  float leg_y1 = 0.7;
  float leg_y2 = 1.-mar_t-0.02;
  float axisTitleSize = 0.045;
  float textSize = 0.035;
  
  TLegend* lTime[2];
  TCanvas* cTime =  new TCanvas("dt_data","dt_data",650,600);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);

  float statH = 0.11;
  float statW = 0.20;
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(1112);
  gStyle->SetStatY(leg_y2);
  gStyle->SetStatX(leg_x2);
  gStyle->SetStatW(statW);
  gStyle->SetStatH(statH);

  std::cout<<"Defining function...\n";
  
  TF1*      fit[2];
  TF1*      fminus[2];
  TF1*      fplus[2];
  TF1*      fbg[2];

  fminus[0] = new TF1("fminus_data","[0]*exp(-x/[1])", fDecayFit_T1, fDecayFit_T2);
  fminus[1] = (TF1*)fminus[0]->Clone("fminus_mc");
  fplus[0] = (TF1*)fminus[0]->Clone("fplus_data");
  fplus[1] = (TF1*)fminus[0]->Clone("fplus_mc");
  fbg[0]    = new TF1("fbg_data","[0]", fDecayFit_T1, fDecayFit_T2 );
  fbg[1]    = (TF1*)fbg[0]->Clone("fbg_mc");

  fit[0]  = new TF1("DecayFit_Data","[0]+[1]*exp(-x/[3])+[2]*exp(-x/[4])",fDecayFit_T1, fDecayFit_T2);
  std::cout<<"Function was defined\n";

  fDecayFit_CPlus   = GetHistMax(hDecayTime[0])*0.6;
  fDecayFit_CMinus  = GetHistMax(hDecayTime[0])*0.4;
  fit[0]->SetParameter(0, fDecayFit_BG );
  fit[0]->SetParLimits(0, 0, 200);
  fit[0]->SetParameter(1, fDecayFit_CPlus );
  fit[0]->SetParLimits(1, 0, hDecayTime[0]->GetMaximum());
  fit[0]->SetParameter(2, fDecayFit_CMinus );
  fit[0]->SetParLimits(2, 0, hDecayTime[0]->GetMaximum());
  fit[0]->SetParameter(3, fDecayFit_TauPlus );
  fit[0]->SetParameter(4, fDecayFit_TauMinus );
  if( fDecayFit_FixFreeLifetime ) fit[0]->FixParameter(3, fDecayFit_TauPlus );
  if( fDecayFit_FixBG )           fit[0]->FixParameter(0, 0.);
  fit[1] = (TF1*)fit[0]->Clone("DecayFit_MC");
  std::cout<<"All parameters were set\n"; 

  //fit[0]->FixParameter(4, fDecayFit_TauMinus);
  //fit[1]->FixParameter(4, fDecayFit_TauMinus);

  int index = 0;
  double chi2[2];
  
  for(int i=0; i<2; i++){
  fit[i]->SetParName(0,"Constant BG");
  fit[i]->SetParName(1,"C_{#mu+}");
  fit[i]->SetParName(2,"C_{#mu-}");
  fit[i]->SetParName(3,"#mu^{+} lifetime (fixed) [ns]");
  fit[i]->SetParName(4,"#mu^{-} lifetime [ns]");
  }
  
  //=======================================================================
  // 1) Data fit
  //=======================================================================
  index = 0;
  //cTime->cd(0);
  gPad->SetLogy();
  TFitResultPtr rdata = hDecayTime[index]->Fit(fit[index],"RS");
  chi2[index] = rdata->Chi2();
  rdata->Print("V");
  
  hDecayTime[index]->SetMarkerStyle(20);
  hDecayTime[index]->SetMarkerSize(0.4);
  hDecayTime[index]->SetLineColor(kBlack);
  hDecayTime[index]->SetLineStyle(1);
  hDecayTime[index]->GetXaxis()->SetTitleSize(axisTitleSize);
  hDecayTime[index]->GetYaxis()->SetTitleSize(axisTitleSize);
  hDecayTime[index]->SetMaximum( hDecayTime[index]->GetBinContent(hDecayTime[index]->GetMaximumBin())*8 ); 
  hDecayTime[index]->Draw("P E X0");
  
  fminus[index]->SetParameter(0, fit[index]->GetParameter(2));
  fminus[index]->SetParameter(1, fit[index]->GetParameter(4));
  fplus[index]->SetParameter(0, fit[index]->GetParameter(1));
  fplus[index]->SetParameter(1, fit[index]->GetParameter(3));
  fbg[index]->SetParameter(0, fit[index]->GetParameter(0));

  fminus[index] ->SetLineColor(kBlue);
  fplus[index]  ->SetLineColor(kGreen+2);
  fbg[index]    ->SetLineColor(kBlack);
  fminus[index] ->SetLineStyle(2);
  fplus[index]  ->SetLineStyle(2);
  fbg[index]    ->SetLineStyle(2);
  fplus[index]    ->Draw("same");
  fminus[index]    ->Draw("same");
  //fbg[index]    ->Draw("same");
 
  lTime[index]  = MakeLegend( leg_x2 - 1.2*statW, leg_y2 - 2.*statH - 0.08, textSize, 4); 
  lTime[index]  ->AddEntry(hDecayTime[index], "Data", "LPE"); 
  lTime[index]  ->AddEntry(fit[index], "Overall fit", "L"); 
//  lTime[index]  ->AddEntry(fbg[index], "Background", "L"); 
  lTime[index]  ->AddEntry(fplus[index], "#mu^{+} (free)", "L"); 
  lTime[index]  ->AddEntry(fminus[index], "#mu^{-} (muonic Ar)", "L"); 
  lTime[index]  ->Draw();

  TPaveText* header = MakeTextBox(mar_l + 0.02, leg_y2, textSize, 2);
    header->AddText("#bf{LArIAT Preliminary}");
    sprintf(buffer,"%s Dataset",runtag.c_str());
    header->AddText(buffer);
    header->Draw();


  //TPaveText* cuts = MakeTextBox(mar_l + 0.02, leg_y2-2.*textSize-0.01, textSize, n);

  int n = 0;


  if( fDecayFit_ReqShower1 ) n += 2;
  if( fAmpCut[0] > 0 ) n++;
  if( fAmpCut[1] > 0 ) n++;
  TPaveText* cuts = MakeTextBox(mar_l + 0.02, leg_y2-2.*textSize-0.01, 0.03, 5);
    cuts->AddText("Cuts: ");
    cuts->AddText(" - 1 stopping track");
    cuts->AddText(" - Michel optical ID");
    if( fDecayFit_ReqShower1 ) 
    cuts->AddText(" - 2D shower reco'd"); //in either plane");
    if( fAmpCut[0] > 0 ) cuts->AddText(Form(" - HMM amp > %3.0f ADC",fAmpCut[0]));
    if( fAmpCut[1] > 0 ) cuts->AddText(Form(" - ETL amp > %3.0f ADC",fAmpCut[1]));
    cuts->Draw();
  
  //=======================================================================
  // 1) MC fit
  //=======================================================================
  index = 1;
  TCanvas* cTimeMC = new TCanvas("dt_mc","dt_mc",650,600);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);
  gPad->SetLogy();
  TFitResultPtr rmc = hDecayTime[index]->Fit(fit[index],"RS");
  chi2[index] = rmc->Chi2();
  rmc->Print("V");
  
  hDecayTime[index]->SetFillColor(kGray);
  hDecayTime[index]->SetLineStyle(1);
  hDecayTime[index]->SetMaximum( hDecayTime[index]->GetMaximum()*4 );
  hDecayTime[index]->Draw("hist"); // "hist"
  fit[1]->Draw("same");
  
  // TO-DO: draw components

}

//###############################################################################
void CalibrationPlots(){
    
    char xTitle[100] = "Time (days)";
    char yTitle[100] = "Single Photoelectron Response [ADC]";
    char yTitleTau[100] = "Average Muon Late-light Lifetime #tau [ns]";
    float leftMargin = 0.10;
    float rightMargin = 0.10;
    float topMargin = 0.05;
    float leg_x1 = 0.12;
    float leg_y2 = 0.88;
    float textSize=0.035;
    //float legX1 = 0.5;

    // ==================================================================
    // Run 1

    if( fRunMode == 1 ) { 
      
      fChargeRes  = 0.;
      fPhelRes    = 0.;
       
      auto c1 = new TCanvas("c1","c1",850,650);
      c1->cd();
      gPad->SetRightMargin(rightMargin);
      gPad->SetLeftMargin(leftMargin);
      //gPad->SetTopMargin(topMargin);
      auto gr  = new TGraphErrors();
    
      gr->SetPoint(gr->GetN(),    1.5, 60.70);  gr->SetPointError(gr->GetN()-1,  1.5, 0.43);
      gr->SetPoint(gr->GetN(),    4.5, 59.44);  gr->SetPointError(gr->GetN()-1,  1.5, 0.32);
      gr->SetPoint(gr->GetN(),    7.5, 57.45);  gr->SetPointError(gr->GetN()-1,  1.5, 0.54);
      gr->SetPoint(gr->GetN(),    10.5,56.99);  gr->SetPointError(gr->GetN()-1,  1.5, 0.65);
   
      gr  ->SetMinimum(62);
      gr  ->SetMaximum(72);
      gr  ->GetXaxis()->SetRangeUser(0.,12.);
      gr  ->SetMarkerStyle(20);
      gr  ->SetMarkerColor(kBlue);
      gr  ->SetLineColor(kBlue);
      gr  ->GetXaxis()->SetTitle(xTitle);
      gr  ->GetYaxis()->SetTitle(yTitle);
      gr  ->GetYaxis()->SetTitleOffset(1.2);
      gr  ->GetXaxis()->SetTitleOffset(1.2);
      gr  ->Draw("AP");

      TPaveText* pt = MakeTextBox(leg_x1, leg_y2, textSize, 2);
      pt  ->AddText("#bf{LArIAT Preliminary}");
      pt  ->AddText(Form("%s Cosmic Michel Electron Dataset",runtag.c_str()));
      pt  ->Draw();
    
      auto leg  = MakeLegend( leg_x1, leg_y2 - 2.5*textSize, textSize, 1);
      leg       ->AddEntry(gr, "ETL (2-inch) PMT",       "lep");
      leg       ->Draw();
      
      // -----------------------------------------------------
      // Tau plots
      auto c2 = new TCanvas("c2","c2",850,650);
      c2->cd();
      gPad->SetLeftMargin(leftMargin);
      gPad->SetRightMargin(rightMargin);
      //gPad->SetTopMargin(topMargin);
      auto grtau  = new TGraphErrors();

      grtau->SetPoint(grtau->GetN(),  1.5, 1265);  grtau->SetPointError(grtau->GetN()-1, 1.5,3);
      grtau->SetPoint(grtau->GetN(),  4.5, 1251);  grtau->SetPointError(grtau->GetN()-1, 1.5, 3);
      grtau->SetPoint(grtau->GetN(),  7.5, 1247);  grtau->SetPointError(grtau->GetN()-1,1.5, 3);
      grtau->SetPoint(grtau->GetN(),    10.5,1270);  grtau->SetPointError(grtau->GetN()-1,1.5,4);
      grtau  ->GetXaxis()->SetRangeUser(0.,12.);
      grtau  ->GetYaxis()->SetRangeUser(1200.,1300.);
      grtau ->GetYaxis()->SetTitleOffset(1.2);
      grtau  ->SetMarkerStyle(20);
      grtau  ->SetMarkerColor(kBlue+2);
      grtau  ->SetLineColor(kBlue+2);
      grtau  ->GetXaxis()->SetTitle(xTitle);
      grtau  ->GetYaxis()->SetTitle(yTitleTau);
      grtau  ->GetYaxis()->SetTitleOffset(1.3);
      grtau  ->GetXaxis()->SetTitleOffset(1.2);
      grtau  ->Draw("AP");
      
      TPaveText* pt2 = MakeTextBox(leg_x1, leg_y2, textSize, 2);
      pt2  ->AddText("#bf{LArIAT Preliminary}");
      pt2  ->AddText(Form("%s Cosmic Michel Electron Dataset",runtag.c_str()));
      pt2  ->Draw();
    
      auto leg2  = MakeLegend( leg_x1, leg_y2 - 2.5*textSize, textSize, 1);
      leg2       ->AddEntry(grtau, "ETL (2-inch) PMT",       "lep");
      leg2       ->Draw();
    
      

    }
   
    
    // ==================================================================
    // Run 2b

    if( fRunMode == 2 ) {
      auto c1 = new TCanvas("c1","c1",850,650);
      c1->cd();
      gPad->SetRightMargin(rightMargin);
      gPad->SetLeftMargin(leftMargin);
      //gPad->SetTopMargin(topMargin);

      TGraphErrors* gr[2];
      
      gr[0]  = new TGraphErrors();
      gr[0]->SetPoint(gr[0]->GetN(), 1.5,  35.07); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.03);
      gr[0]->SetPoint(gr[0]->GetN(), 4.5,  34.67); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.03);
      gr[0]->SetPoint(gr[0]->GetN(), 7.5,  34.69); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.04);
      gr[0]->SetPoint(gr[0]->GetN(), 10.5, 34.11); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.04);
      gr[0]->SetPoint(gr[0]->GetN(), 14.0, 34.04); gr[0]->SetPointError(gr[0]->GetN()-1,2.0,0.03);
      gr[0]->SetPoint(gr[0]->GetN(), 25.5, 43.42); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.04);
      gr[0]->SetPoint(gr[0]->GetN(), 28.5, 43.06); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.03);
      gr[0]->SetPoint(gr[0]->GetN(), 31.0, 43.10); gr[0]->SetPointError(gr[0]->GetN()-1,1.0,0.08);
//      gr[0]->SetPoint(gr[0]->GetN(), 38.5, 43.04); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.03);
      gr[0]->SetPoint(gr[0]->GetN(), 41.5, 42.32); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.04);
      gr[0]->SetPoint(gr[0]->GetN(), 44.5, 42.04); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.04);
      gr[0]->SetPoint(gr[0]->GetN(), 47.5, 41.68); gr[0]->SetPointError(gr[0]->GetN()-1,1.5,0.05);
      gr[0]  ->GetXaxis()->SetRangeUser(0.,49.);
      gr[0]  ->SetMarkerStyle(20);
      gr[0]  ->SetMarkerColor(kRed);
      gr[0]  ->SetLineColor(kRed);
      gr[0]  ->GetXaxis()->SetTitle(xTitle);
      gr[0]  ->GetYaxis()->SetTitle(yTitle);
      gr[0]  ->GetYaxis()->SetTitleOffset(1.2);
      gr[0]  ->GetXaxis()->SetTitleOffset(1.2);
      gr[0]  ->Draw("AP");
      
      gr[1]  = new TGraphErrors();
      gr[1]->SetPoint(gr[1]->GetN(), 1.5,  38.43); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.07);
      gr[1]->SetPoint(gr[1]->GetN(), 4.5,  38.28); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.06);
      gr[1]->SetPoint(gr[1]->GetN(), 7.5,  37.39); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.07);
      gr[1]->SetPoint(gr[1]->GetN(), 10.5, 36.82); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.07);
      gr[1]->SetPoint(gr[1]->GetN(), 14.0, 36.72); gr[1]->SetPointError(gr[1]->GetN()-1,2.0,0.06);
      gr[1]->SetPoint(gr[1]->GetN(), 25.5, 40.87); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.08);
      gr[1]->SetPoint(gr[1]->GetN(), 28.5, 40.60); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.07);
      gr[1]->SetPoint(gr[1]->GetN(), 31.0, 40.11); gr[1]->SetPointError(gr[1]->GetN()-1,1.0,0.20);
 //     gr[1]->SetPoint(gr[1]->GetN(), 38.5, 42.13); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.06);
      gr[1]->SetPoint(gr[1]->GetN(), 41.5, 39.66); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.08);
      gr[1]->SetPoint(gr[1]->GetN(), 44.5, 39.16); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.08);
      gr[1]->SetPoint(gr[1]->GetN(), 47.5, 39.02); gr[1]->SetPointError(gr[1]->GetN()-1,1.5,0.10);
      gr[1]  ->GetXaxis()->SetRangeUser(0.,49.);
      gr[1]  ->SetMarkerStyle(20);
      gr[1]  ->SetMarkerColor(kBlue);
      gr[1]  ->SetLineColor(kBlue);
      gr[1]  ->GetXaxis()->SetTitle(xTitle);
      gr[1]  ->GetYaxis()->SetTitle(yTitle);
      gr[1]  ->GetYaxis()->SetTitleOffset(1.2);
      gr[1]  ->GetXaxis()->SetTitleOffset(1.2);
      gr[1]  ->Draw("P");

      gr[0]->SetMinimum(32.);
      gr[0]->SetMaximum(48.);
      
      TPaveText* pt = MakeTextBox(leg_x1, leg_y2, textSize, 2);
      pt  ->AddText("#bf{LArIAT Preliminary}");
      pt  ->AddText(Form("%s Cosmic Michel Electron Dataset",runtag.c_str()));
      pt  ->Draw();
    
      auto leg  = MakeLegend( leg_x1, leg_y2 - 2.*textSize, textSize, 2);
      leg       ->AddEntry(gr[0], "HMM (3-inch) PMT",       "lep");
      leg       ->AddEntry(gr[1], "ETL (2-inch) PMT",       "lep");
      leg       ->Draw();
      
      // -----------------------------------------------------
      // Tau plots
      auto c2 = new TCanvas("c2","c2",850,650);
      c2->cd();
      gPad->SetRightMargin(rightMargin);
      gPad->SetLeftMargin(leftMargin);
      //gPad->SetTopMargin(topMargin);
      TGraphErrors* grtau[2];
      
      grtau[0] = new TGraphErrors();
    grtau[0]->SetPoint(grtau[0]->GetN(), 1.5, 1135.8); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,1.3);
    grtau[0]->SetPoint(grtau[0]->GetN(), 4.5,  1138.3); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,1.3);
    grtau[0]->SetPoint(grtau[0]->GetN(), 7.5,  1124.5); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,1.5);
    grtau[0]->SetPoint(grtau[0]->GetN(), 10.5, 1128.6); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,1.3);
    grtau[0]->SetPoint(grtau[0]->GetN(), 14.0, 1122.3); grtau[0]->SetPointError(grtau[0]->GetN()-1,2.0,1.2);
    grtau[0]->SetPoint(grtau[0]->GetN(), 25.5, 1189.0); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,2);
    grtau[0]->SetPoint(grtau[0]->GetN(), 28.5, 1186.5); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,1.5);
    grtau[0]->SetPoint(grtau[0]->GetN(), 31.0, 1189); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.0,3.);
//    grtau[0]->SetPoint(grtau[0]->GetN(), 38.5,1131 ); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,3.);
    grtau[0]->SetPoint(grtau[0]->GetN(), 41.5, 1125 ); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,2);
    grtau[0]->SetPoint(grtau[0]->GetN(), 44.5, 1107 ); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,2);
    grtau[0]->SetPoint(grtau[0]->GetN(), 47.5, 1081 ); grtau[0]->SetPointError(grtau[0]->GetN()-1,1.5,2);
    grtau[0]  ->GetXaxis()->SetRangeUser(0.,49.);
    grtau[0]  ->GetYaxis()->SetRangeUser(1000.,1350.);
    grtau[0]  ->SetMarkerStyle(20);
    grtau[0]  ->SetMarkerColor(kRed+2);
    grtau[0]  ->SetLineColor(kRed+2);
    grtau[0]  ->GetXaxis()->SetTitle(xTitle);
    grtau[0]  ->GetYaxis()->SetTitle(yTitleTau);
    grtau[0]  ->GetYaxis()->SetTitleOffset(1.3);
    grtau[0]  ->GetXaxis()->SetTitleOffset(1.2);
    grtau[0]  ->Draw("AP");
    
      grtau[1] = new TGraphErrors();
    grtau[1]->SetPoint(grtau[1]->GetN(), 1.5, 1221 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3.);
    grtau[1]->SetPoint(grtau[1]->GetN(), 4.5, 1203 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3.);
    grtau[1]->SetPoint(grtau[1]->GetN(), 7.5, 1233 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3);
    grtau[1]->SetPoint(grtau[1]->GetN(), 10.5,1214); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3);
    grtau[1]->SetPoint(grtau[1]->GetN(), 14.0,1209 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,2.0,2);
    grtau[1]->SetPoint(grtau[1]->GetN(), 25.5,1191 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3);
    grtau[1]->SetPoint(grtau[1]->GetN(), 28.5,1221 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3);
    grtau[1]->SetPoint(grtau[1]->GetN(), 31.0,1175 ); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.0,6);
//    grtau[1]->SetPoint(grtau[1]->GetN(), 38.5,1152); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,6);
    grtau[1]->SetPoint(grtau[1]->GetN(), 41.5, 1184); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3.);
    grtau[1]->SetPoint(grtau[1]->GetN(), 44.5, 1145); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,3.);
    grtau[1]->SetPoint(grtau[1]->GetN(), 47.5, 1079); grtau[1]->SetPointError(grtau[1]->GetN()-1,1.5,4.);
    grtau[1]  ->GetXaxis()->SetRangeUser(0.,49.);
    grtau[1]  ->SetMarkerStyle(20);
    grtau[1]  ->SetMarkerColor(kBlue+2);
    grtau[1]  ->SetLineColor(kBlue+2);
    grtau[1]  ->GetXaxis()->SetTitle(xTitle);
    grtau[1]  ->GetYaxis()->SetTitle(yTitleTau);
    grtau[1]  ->GetYaxis()->SetTitleOffset(1.3);
    grtau[1]  ->GetXaxis()->SetTitleOffset(1.2);
    grtau[1]  ->Draw("P");
      
      TPaveText* pt2 = MakeTextBox(leg_x1, leg_y2, textSize, 2);
      pt2  ->AddText("#bf{LArIAT Preliminary}");
      pt2  ->AddText(Form("%s Cosmic Michel Electron Dataset",runtag.c_str()));
      pt2  ->Draw();
    
      auto leg2  = MakeLegend( leg_x1, leg_y2 - 2.5*textSize, textSize, 2);
      leg2       ->AddEntry(grtau[0], "HMM (3-inch) PMT",       "lep");
      leg2       ->AddEntry(grtau[1], "ETL (2-inch) PMT",       "lep");
      leg2       ->Draw();
    
    
    
    }


 
}

//##############################################################################
void MakeCalibrationPlots(){
  hxLifetime = hTimeVsElectronLifetime->ProfileX();
  hxCharge = hTimeVsCharge->ProfileX();
  hxLight = hTimeVsLight->ProfileX();
  hxChargePerCm = hTimeVsChargePerCm->ProfileX();
  hxLightPerCm = hTimeVsLightPerCm->ProfileX();

  // "Reference" baseline level for the calibration plots
  // will be a truncated mean of the values inputted in the
  // time plots.
  std::vector<float> QPerCmvals;
  std::vector<float> LPerCmvals;
  for(size_t i=1; i<hxChargePerCm->GetNbinsX(); i++){
    if( hxChargePerCm->GetBinContent(i) == 0 ) continue;
    QPerCmvals.push_back( hxChargePerCm->GetBinContent(i) );
    LPerCmvals.push_back( hxLightPerCm->GetBinContent(i) );
  }

  refQPerCm = CalcTruncatedMean( QPerCmvals, 0.33, 0);
  refLPerCm = CalcTruncatedMean( LPerCmvals, 0.33, 0);

}


// ##############################################################################
// Plots comparing spatial distributions and visibilities of Michel showers
// in data vs. MC.
void SpatialPlots(){
  
  float textSize = 0.045;
  float mar_b  = 0.12;
  float mar_t  = 0.08;
  float mar_l  = 0.14;
  float mar_r  = 0.02;
  float axisTitleSize = 0.05;
  float axisTitleOffset = 1.2;
  float axisLabelSize = 0.045;
  float leg_x1 = mar_l + 0.03;
  float leg_y2 = 1.-mar_t-0.01; 
 
  std::string type[2]={"Data:","MC:"};
  
  // Scale
  hElShowerCentroid_X[1]->Scale( hElShowerCentroid_X[0]->GetEntries() / hElShowerCentroid_X[1]->Integral() );
  hElShowerCentroid_Y[1]->Scale( hElShowerCentroid_Y[0]->GetEntries() / hElShowerCentroid_Y[1]->Integral() );
  hElShowerCentroid_Z[1]->Scale( hElShowerCentroid_Z[0]->GetEntries() / hElShowerCentroid_Z[1]->Integral() );
  hElShowerVis[1]->Scale( hElShowerVis[0]->GetEntries() / hElShowerVis[1]->Integral() );
  
  // -----------------------------------
  // Plot X, Y, Z
  TCanvas*    c_centroid = new TCanvas("centroid","centroid",1050,400);
  c_centroid -> Divide(3,1);
  TLegend* leg[3];
  gStyle->SetOptStat(0);
  float mean[2];
  float rms[2];
  int index;
   
  index=1;
  c_centroid->cd(index);
    gPad->SetBottomMargin(mar_b); gPad->SetTopMargin(mar_t); gPad->SetRightMargin(mar_r); gPad->SetLeftMargin(mar_l);
    hElShowerCentroid_X[1]->SetLineColor(kRed-3);
    hElShowerCentroid_X[1]->SetLineWidth(2);
    hElShowerCentroid_X[1]->SetMaximum( std::max(hElShowerCentroid_X[1]->GetMaximum(),hElShowerCentroid_X[0]->GetMaximum())*1.3);
    hElShowerCentroid_X[1]->GetXaxis()->SetTitleSize(axisTitleSize);
    hElShowerCentroid_X[1]->GetXaxis()->SetLabelSize(axisLabelSize);
    hElShowerCentroid_X[1]->GetYaxis()->SetTitleSize(axisTitleSize);
    hElShowerCentroid_X[1]->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hElShowerCentroid_X[1]->GetYaxis()->SetLabelSize(axisLabelSize);
    hElShowerCentroid_X[1]->GetXaxis()->SetTitle("Shower centroid X [cm]");
    hElShowerCentroid_X[1]->GetYaxis()->SetTitle("Events");
    hElShowerCentroid_X[1]->Draw("hist E1 X0");
    hElShowerCentroid_X[0]->SetMarkerStyle(20);
    hElShowerCentroid_X[0]->SetMarkerColor(kBlack);
    hElShowerCentroid_X[0]->SetLineColor(kBlack);
    hElShowerCentroid_X[0]->SetLineWidth(2);
    hElShowerCentroid_X[0]->Draw("same P E X0");
    hElShowerCentroid_X[1]->Draw("sameaxis");
    mean[0] = hElShowerCentroid_X[0]->GetMean();
    rms[0]  = hElShowerCentroid_X[0]->GetRMS();
    mean[1] = hElShowerCentroid_X[1]->GetMean();
    rms[1]  = hElShowerCentroid_X[1]->GetRMS();
    leg[index-1]  = MakeLegend(leg_x1, leg_y2, textSize, 2);
    leg[index-1]  ->AddEntry(hElShowerCentroid_X[0],Form("%-5s: #mu = %4.1f cm, #sigma = %4.1f cm",type[0].c_str(),mean[0],rms[0]),"PL");
    leg[index-1]  ->AddEntry(hElShowerCentroid_X[1],Form("%-5s: #mu = %4.1f cm, #sigma = %4.1f cm",type[1].c_str(),mean[1],rms[1]),"FL");
    leg[index-1]  ->Draw();
  
  index=2;
  c_centroid->cd(index);
    gPad->SetBottomMargin(mar_b); gPad->SetTopMargin(mar_t); gPad->SetRightMargin(0.02); gPad->SetLeftMargin(mar_l);
    hElShowerCentroid_Y[1]->SetLineColor(kRed-3);
    hElShowerCentroid_Y[1]->SetLineWidth(2);
    hElShowerCentroid_Y[1]->SetMaximum( std::max(hElShowerCentroid_Y[1]->GetMaximum(),hElShowerCentroid_Y[0]->GetMaximum())*1.3);
    hElShowerCentroid_Y[1]->GetXaxis()->SetTitleSize(axisTitleSize);
    hElShowerCentroid_Y[1]->GetXaxis()->SetLabelSize(axisLabelSize);
    hElShowerCentroid_Y[1]->GetYaxis()->SetTitleSize(axisTitleSize);
    hElShowerCentroid_Y[1]->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hElShowerCentroid_Y[1]->GetYaxis()->SetLabelSize(axisLabelSize);
    hElShowerCentroid_Y[1]->GetXaxis()->SetTitle("Shower centroid Y [cm]");
    hElShowerCentroid_Y[1]->GetYaxis()->SetTitle("Events");
    hElShowerCentroid_Y[1]->Draw("hist E1 X0");
    hElShowerCentroid_Y[0]->SetMarkerStyle(20);
    hElShowerCentroid_Y[0]->SetMarkerColor(kBlack);
    hElShowerCentroid_Y[0]->SetLineColor(kBlack);
    hElShowerCentroid_Y[0]->SetLineWidth(2);
    hElShowerCentroid_Y[0]->Draw("same P E X0");
    hElShowerCentroid_Y[1]->Draw("sameaxis");
    mean[0] = hElShowerCentroid_Y[0]->GetMean();
    rms[0]  = hElShowerCentroid_Y[0]->GetRMS();
    mean[1] = hElShowerCentroid_Y[1]->GetMean();
    rms[1]  = hElShowerCentroid_Y[1]->GetRMS();
    leg[index-1]  = MakeLegend(leg_x1, leg_y2, textSize, 2);
    leg[index-1]  ->AddEntry(hElShowerCentroid_Y[0],Form("Data: #mu = %4.1f cm, #sigma = %4.1f cm",mean[0],rms[0]),"PL");
    leg[index-1]  ->AddEntry(hElShowerCentroid_Y[1],Form("MC: #mu = %4.1f cm, #sigma = %4.1f cm",mean[1],rms[1]),"FL");
    leg[index-1]  ->Draw();
  
  index=3;
  c_centroid->cd(index);
    gPad->SetBottomMargin(mar_b); gPad->SetTopMargin(mar_t); gPad->SetRightMargin(0.02); gPad->SetLeftMargin(mar_l);
    hElShowerCentroid_Z[1]->SetMaximum( std::max(hElShowerCentroid_Z[1]->GetMaximum(),hElShowerCentroid_Z[0]->GetMaximum())*1.3);
    hElShowerCentroid_Z[1]->SetLineColor(kRed-3);
    hElShowerCentroid_Z[1]->SetLineWidth(2);
    hElShowerCentroid_Z[1]->GetXaxis()->SetTitleSize(axisTitleSize);
    hElShowerCentroid_Z[1]->GetXaxis()->SetLabelSize(axisLabelSize);
    hElShowerCentroid_Z[1]->GetYaxis()->SetTitleSize(axisTitleSize);
    hElShowerCentroid_Z[1]->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hElShowerCentroid_Z[1]->GetYaxis()->SetLabelSize(axisLabelSize);
    hElShowerCentroid_Z[1]->GetXaxis()->SetTitle("Shower centroid Z [cm]");
    hElShowerCentroid_Z[1]->GetYaxis()->SetTitle("Events");
    hElShowerCentroid_Z[1]->Draw("hist E1 X0");
    hElShowerCentroid_Z[0]->SetMarkerStyle(20);
    hElShowerCentroid_Z[0]->SetMarkerColor(kBlack);
    hElShowerCentroid_Z[0]->SetLineColor(kBlack);
    hElShowerCentroid_Z[0]->SetLineWidth(2);
    hElShowerCentroid_Z[0]->Draw("same P E X0");
    hElShowerCentroid_Z[1]->Draw("sameaxis");
    mean[0] = hElShowerCentroid_Z[0]->GetMean();
    rms[0]  = hElShowerCentroid_Z[0]->GetRMS();
    mean[1] = hElShowerCentroid_Z[1]->GetMean();
    rms[1]  = hElShowerCentroid_Z[1]->GetRMS();
    leg[index-1]  = MakeLegend(leg_x1, leg_y2, textSize, 2);
    leg[index-1]  ->AddEntry(hElShowerCentroid_Z[0],Form("Data: #mu = %4.1f cm, #sigma = %4.1f cm",mean[0],rms[0]),"PL");
    leg[index-1]  ->AddEntry(hElShowerCentroid_Z[1],Form("MC: #mu = %4.1f cm, #sigma = %4.1f cm",mean[1],rms[1]),"FL");
    leg[index-1]  ->Draw();
  
  // -----------------------------------
  // Visibility
  TCanvas*    c_showervis = new TCanvas("showervis","showervis",600,600);
    gPad->SetBottomMargin(mar_b); gPad->SetTopMargin(mar_t); gPad->SetRightMargin(0.06); gPad->SetLeftMargin(mar_l);
    hElShowerVis[1]->SetMaximum( std::max(hElShowerVis[1]->GetMaximum(),hElShowerVis[0]->GetMaximum())*1.3);
    hElShowerVis[1]->SetLineColor(kRed-3);
    hElShowerVis[1]->SetLineWidth(2);
    hElShowerVis[1]->GetXaxis()->SetTitleSize(axisTitleSize);
    hElShowerVis[1]->GetXaxis()->SetLabelSize(axisLabelSize);
    hElShowerVis[1]->GetYaxis()->SetTitleSize(axisTitleSize);
    hElShowerVis[1]->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hElShowerVis[1]->GetYaxis()->SetLabelSize(axisLabelSize);
    hElShowerVis[1]->GetXaxis()->SetTitle("Fractional photon visibility (#times 10^{-3})");
    hElShowerVis[1]->GetYaxis()->SetTitle("Events");
//    hElShowerVis[1]->GetXaxis()->SetNdivisions(5,5,0);
    hElShowerVis[1]->Draw("hist E1 X0");
    hElShowerVis[0]->SetMarkerStyle(20);
    hElShowerVis[0]->SetMarkerSize(0.6);
    hElShowerVis[0]->SetMarkerColor(kBlack);
    hElShowerVis[0]->SetLineColor(kBlack);
    hElShowerVis[0]->SetLineWidth(2);
    hElShowerVis[0]->Draw("same P E X0");
    hElShowerVis[1]->Draw("sameaxis");
    //gPad->SetGridx();
    //gPad->SetGridy();
    mean[0] = hElShowerVis[0]->GetMean();
    rms[0]  = hElShowerVis[0]->GetRMS();
    mean[1] = hElShowerVis[1]->GetMean();
    rms[1]  = hElShowerVis[1]->GetRMS();
    TLegend* lg = MakeLegend(leg_x1, leg_y2, textSize, 2);
    lg  ->AddEntry(hElShowerVis[0],Form("Data: #mu = %5.3f #times 10^{-3}, #sigma = %5.3f #times 10^{-3}",mean[0],rms[0]),"PL");
    lg  ->AddEntry(hElShowerVis[1],Form("MC: #mu = %5.3f #times 10^{-3}, #sigma = %5.3f #times 10^{-3}",mean[1],rms[1]),"FL");
    lg  ->Draw();

}

// ##############################################################################
// Make MC weight histogram for smoothing out the spatial distribution in Z-X
void MakeWeightMap(){
  
  // ---------------------------------------
  // 2D map
  /*
  float sum = 0;
  float n = 0;
  for(size_t i=1; i<hTrueMuEndpt_ZX->GetXaxis()->GetNbins()+1; i++ ) {
    for(size_t j=1; j<hTrueMuEndpt_ZX->GetYaxis()->GetNbins()+1; j++ ) {
      sum += hTrueMuEndpt_ZX->GetBinContent(i,j);
      n ++; 
    }
  }
  std::cout<<"N = "<<n<<"\n";
  float mean = sum/n;
  for(size_t i=1; i<hTrueMuEndpt_ZX->GetXaxis()->GetNbins()+1; i++ ) {
    for(size_t j=1; j<hTrueMuEndpt_ZX->GetYaxis()->GetNbins()+1; j++ ) {
      float w = mean / hTrueMuEndpt_ZX->GetBinContent(i,j);
      hMCSpatialWgt_ZX->SetBinContent(i,j, w);
    }
  }
  */
  
  // --------------------------------------
  // 3D map
  float sum = 0;
  float n = 0; 
  for(size_t i=1; i<=hTrueMuEndpt->GetXaxis()->GetNbins(); i++ ){
    for(size_t j=1; j<=hTrueMuEndpt->GetYaxis()->GetNbins(); j++ ){
      for(size_t k=1; k<=hTrueMuEndpt->GetZaxis()->GetNbins(); k++ ){
        sum += hTrueMuEndpt->GetBinContent(i,j,k);
        n++;
      }
    }
  }
  float mean = sum/n;
  for(size_t i=1; i<=hTrueMuEndpt->GetXaxis()->GetNbins(); i++ ){
    for(size_t j=1; j<=hTrueMuEndpt->GetYaxis()->GetNbins(); j++ ){
      for(size_t k=1; k<=hTrueMuEndpt->GetZaxis()->GetNbins(); k++ ){
        float w = mean / hTrueMuEndpt->GetBinContent(i,j,k);
        hMCSpatialWgt->SetBinContent(i,j,k,w);
      }
    }
  }

 
  /* 
  // get mean
  float sum = 0;
  float n = 0;
  float sum_zx = 0;
  float n_zx = 0;
  for(size_t i=1; i<hTrueMuEndpt->GetXaxis()->GetNbins()+1; i++ ) {
    for(size_t j=1; j<hTrueMuEndpt->GetYaxis()->GetNbins()+1; j++ ) {
      sum_zx += hTrueMuEndpt_ZX->GetBinContent(i,j);
      n_zx++;
      for(size_t k=1; k<hTrueMuEndpt->GetZaxis()->GetNbins()+1; k++ ) {
        sum += hTrueMuEndpt->GetBinContent(i,j,k);
        n++; 
      }
    }
  }
  std::cout<<"N = "<<n<<"\n";
  float mean = sum/n;
  float mean_zx = sum_zx/n_zx;
  float w = 1.;
  for(size_t i=1; i<hTrueMuEndpt->GetXaxis()->GetNbins()+1; i++ ) {
    for(size_t j=1; j<hTrueMuEndpt->GetYaxis()->GetNbins()+1; j++ ) {
      w = mean_zx / hTrueMuEndpt_ZX->GetBinContent(i,j);
      hMCSpatialWgt_ZX->SetBinContent(i,j, w);
      for(size_t k=1; k<hTrueMuEndpt->GetZaxis()->GetNbins()+1; k++ ) {
        w = mean / hTrueMuEndpt->GetBinContent(i,j,k);
        hMCSpatialWgt->SetBinContent(i,j,k, w);
      }
    }
  }
  */

 
  // normalize 3D map of reco'd centroid
  if( hElShowerCentroid[0]->GetEntries() && hElShowerCentroid[1]->GetEntries() ) {
    hElShowerCentroid[1]->Scale( hElShowerCentroid[0]->GetEntries()/hElShowerCentroid[1]->Integral() );
    
    for(size_t i=1; i<hElShowerCentroid[0]->GetXaxis()->GetNbins()+1; i++ ){
      for(size_t j=1; j<hElShowerCentroid[0]->GetYaxis()->GetNbins()+1; j++ ){
        for(size_t k=1; k<hElShowerCentroid[0]->GetZaxis()->GetNbins()+1; k++ ){
          float w = hElShowerCentroid[0]->GetBinContent(i,j,k) / hElShowerCentroid[1]->GetBinContent(i,j,k);
          hMCSpatialWgt->SetBinContent(i,j,k,w);
        }
      }
    }

  }
  
  //hMCSpatialWgt_ZX->Write(); 
  //wgtFile->Close();




  
}


//###############################################################################
void StabilityPlots(){
  
  MakeCalibrationPlots();
  
  float textSize = 0.035;
  float mar_b  = 0.15;
  float axisTitleSize = 0.07;
  float axisTitleOffset = 0.5;
  float axisLabelSize = 0.07;

  
  TCanvas*    cStability = new TCanvas("stab","stab",1000,700);
  cStability -> Divide(1,3);
  gStyle->SetOptStat(0);

  cStability->cd(1);
    gPad  ->SetBottomMargin(mar_b);
    gPad  ->SetTopMargin(0.08);
    gPad  ->SetRightMargin(0.05);
    gPad  ->SetGridx();
    gPad  ->SetGridy();
    hxLifetime->GetXaxis()->SetTitleSize(axisTitleSize);
    hxLifetime->GetXaxis()->SetLabelSize(axisLabelSize);
    hxLifetime->GetYaxis()->SetTitleSize(axisTitleSize);
    hxLifetime->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hxLifetime->GetYaxis()->SetLabelSize(axisLabelSize);
    hxLifetime->SetTitle("");
    hxLifetime->GetYaxis()->SetTitle("Electron lifetime [#mus]");
    hxLifetime->SetMarkerStyle(4);
    hxLifetime->Draw();

  cStability->cd(2);
    gPad  ->SetBottomMargin(mar_b);
    gPad  ->SetTopMargin(0.08);
    gPad  ->SetRightMargin(0.05);
    gPad  ->SetGridx();
    gPad  ->SetGridy();
    hxCharge->SetTitle("");
    hxCharge->GetXaxis()->SetTitleSize(axisTitleSize);
    hxCharge->GetXaxis()->SetLabelSize(axisLabelSize);
    hxCharge->GetYaxis()->SetTitleSize(axisTitleSize);
    hxCharge->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hxCharge->GetYaxis()->SetLabelSize(axisLabelSize);
    hxCharge->SetLineColor(kBlue+2);
    hxCharge->SetMarkerStyle(21);
    hxCharge->SetMarkerSize(0.4);
    hxCharge->SetMarkerColor(kBlue+2);
    hxCharge->GetYaxis()->SetTitle("Michel charge [e-]");
    hxCharge->Draw();
    TLine* mc_aveQ = new TLine(gPad->GetUxmin(), hQ[1]->GetMean(), gPad->GetUxmax(), hQ[1]->GetMean() );
    mc_aveQ->SetX1( 0);
    mc_aveQ->SetX2( fDayEnd );
    mc_aveQ->SetY1( hQ[1]->GetMean() );
    mc_aveQ->SetY2( hQ[1]->GetMean() );
    mc_aveQ->SetLineColor(kBlue);
    mc_aveQ->SetLineWidth(1);
    mc_aveQ->Draw();


  cStability->cd(3);
    gPad  ->SetBottomMargin(mar_b);
    gPad  ->SetTopMargin(0.08);
    gPad  ->SetRightMargin(0.05);
    gPad  ->SetGridx();
    gPad  ->SetGridy();
    hxLight->SetTitle("");
    hxLight->GetXaxis()->SetTitleSize(axisTitleSize);
    hxLight->GetXaxis()->SetLabelSize(axisLabelSize);
    hxLight->GetYaxis()->SetTitleSize(axisTitleSize);
    hxLight->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hxLight->GetYaxis()->SetLabelSize(axisLabelSize);
    hxLight->SetLineColor(kRed+2);
    hxLight->SetMarkerStyle(21);
    hxLight->SetMarkerSize(0.4);
    hxLight->SetMarkerColor(kRed+2);
    hxLight->GetYaxis()->SetTitle("Michel light [#gamma]");
    hxLight->Draw();
    TLine* mc_aveL = new TLine(gPad->GetUxmin(), hL[1]->GetMean(), gPad->GetUxmax(), hL[1]->GetMean() );
    mc_aveL->SetX1( 0);
    mc_aveL->SetX2( fDayEnd );
    mc_aveL->SetY1( hL[1]->GetMean() );
    mc_aveL->SetY2( hL[1]->GetMean() );
    mc_aveL->SetLineColor(kRed);
    mc_aveL->SetLineWidth(1);
    mc_aveL->Draw();

    cStability->Update();
  

    // ====================================================================
    TCanvas*    cStabilityCal = new TCanvas("stabcal","stabcal",1000,700);
    cStabilityCal -> Divide(1,3);
    gStyle->SetOptStat(0);
  
  cStabilityCal->cd(1);
    gPad  ->SetBottomMargin(mar_b);
    gPad  ->SetTopMargin(0.08);
    gPad  ->SetRightMargin(0.05);
    gPad  ->SetGridx();
    gPad  ->SetGridy();
    hxLifetime->GetXaxis()->SetTitleSize(axisTitleSize);
    hxLifetime->GetXaxis()->SetLabelSize(axisLabelSize);
    hxLifetime->GetYaxis()->SetTitleSize(axisTitleSize);
    hxLifetime->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hxLifetime->GetYaxis()->SetLabelSize(axisLabelSize);
    hxLifetime->SetTitle("");
    hxLifetime->GetYaxis()->SetTitle("Electron lifetime [#mus]");
    hxLifetime->SetMarkerStyle(4);
    hxLifetime->Draw();

  cStabilityCal->cd(2);
    gPad  ->SetBottomMargin(mar_b);
    gPad  ->SetTopMargin(0.08);
    gPad  ->SetRightMargin(0.05);
    gPad  ->SetGridx();
    gPad  ->SetGridy();
    hxChargePerCm->SetTitle("");
    hxChargePerCm->GetXaxis()->SetTitleSize(axisTitleSize);
    hxChargePerCm->GetXaxis()->SetLabelSize(axisLabelSize);
    hxChargePerCm->GetYaxis()->SetTitleSize(axisTitleSize);
    hxChargePerCm->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hxChargePerCm->GetYaxis()->SetLabelSize(axisLabelSize);
    hxChargePerCm->SetLineColor(kBlack);
    hxChargePerCm->SetMarkerStyle(21);
    hxChargePerCm->SetMarkerSize(0.4);
    hxChargePerCm->SetMarkerColor(kBlack);
    hxChargePerCm->GetYaxis()->SetTitle("Crossing #mu Q density [e-/cm]");
    hxChargePerCm->Draw();
    TLine* refQ = new TLine(0., refQPerCm, fDayEnd, refQPerCm);
    refQ->SetLineColor(kRed);
    refQ->SetLineWidth(1);
    refQ->Draw();

  cStabilityCal->cd(3);
    gPad  ->SetBottomMargin(mar_b);
    gPad  ->SetTopMargin(0.08);
    gPad  ->SetRightMargin(0.05);
    gPad  ->SetGridx();
    gPad  ->SetGridy();
    hxLightPerCm->SetTitle("");
    hxLightPerCm->GetXaxis()->SetTitleSize(axisTitleSize);
    hxLightPerCm->GetXaxis()->SetLabelSize(axisLabelSize);
    hxLightPerCm->GetYaxis()->SetTitleSize(axisTitleSize);
    hxLightPerCm->GetYaxis()->SetTitleOffset(axisTitleOffset);
    hxLightPerCm->GetYaxis()->SetLabelSize(axisLabelSize);
    hxLightPerCm->SetLineColor(kBlack);
    hxLightPerCm->SetMarkerStyle(21);
    hxLightPerCm->SetMarkerSize(0.4);
    hxLightPerCm->SetMarkerColor(kBlack);
    hxLightPerCm->GetYaxis()->SetTitle("Crossing #mu L density [#gamma/cm]");
    hxLightPerCm->Draw();
    TLine* refL = new TLine(0., refLPerCm, fDayEnd, refLPerCm);
    refL->SetLineColor(kRed);
    refL->SetLineWidth(1);
    refL->Draw();

    cStabilityCal->Update();
  


}


//################################################################################
void EventCuts(){

  fMaxMCEvts=-1;

  TH1D* hOpEventCuts;
  TH1D* hEventCuts;// = (TH1D*)fileData
//  fTreeMC_ns   = (TTree*)fFileMC_ns->Get("michelana/anatree");
  
    
  printf("======================================================================\n");
  printf("Event cut parameters\n\n");
  printf("  Min mu cluster size     : %d \n", fMinMuClusterSize);
  printf("  Min-max el cluster size : %d - %d \n", fMinElClusterSize, fMaxElClusterSize);
  printf("  Min Bragg slope         : %5.0f\n",fMinBraggSlope); 
  printf("  Min mu linearity/frac   : %4.2f / %4.2f\n",fMinMuLinearity,fMinFracMuHitsLinear);
  printf("  Min mu hits for dir fit : %d\n",fMinMuClusterHitsEndFit); 
  printf("  Min-max decay angle     : %.0f - %.0f\n",fMinDecayAngle2D, fMaxDecayAngle2D);
  printf("  Min-max el shower size  : %d - %d \n",fMinElShowerSize, fMaxElShowerSize);
  printf("  Min shower frac         : %4.2f\n",fMinShowerFrac);
  printf("\n");
  printf("  Min num 3D pts          : %d\n",fMinNumPts3D);
  printf("  Frac hits 3D            : %.2f\n",fMinFracHits3D);
  printf("  Min proj. dist to wires : %.0f\n",fMinProjDist);
  printf("  Decay time cut          : %.0f\n",fdTcut);
  printf("\n\n");


  for(size_t i=0; i<2; i++){
    nEvt_TotalEvents=0;
    nEvt_1StpTrk=0;
    nEvt_OpID_timeRange=0;
    nEvt_OpID=0;
    nEvt_OpID_muPromptCut=0;
      int nEvt_OpID_RMS = 0;
      int nEvt_OpID_2hits = 0;
      int nEvt_OpID_dTmatch = 0;
      int nEvt_OpID_hitWidths = 0;
      int nEvt_OpID_hitTime   = 0;
      int nEvt_OpID_hitSat    = 0;
    nEvt_Shwr2D_BndFound=0;
    nEvt_Shwr2D_ClusterSize=0;
    nEvt_Shwr2D_ShowerSize=0;
    nEvt_Shwr2D_ExtraHits=0;
    nEvt_Shwr2D_Slope=0;
    nEvt_Shwr2D_MuLinearity=0;
    nEvt_Shwr2D_MuDir=0;
    nEvt_Shwr2D_DecayAngle=0;
    nEvt_Shwr3D_ShowersOnBothPlanes=0;
    nEvt_Shwr3D_Npts=0;
    nEvt_Shwr3D_FracHits3D=0;
    nEvt_Shwr3D_Direction=0;
    nEvt_Shwr3D_Centroid=0;
    nEvt_3DMuEnd=0;
    nEvt_dTcut=0;
    nEvt_DecayFit_TimeRange=0;
    nEvt_DecayFit_2DShower=0;
    nEvt_DecayFit_AmpCut=0;
    nEvt_Shwr2D_ShowerFrac=0;
  
    if( i==0 ) {
      hEventCuts = (TH1D*)fFile->Get("michelana/EventCuts");
      hOpEventCuts = (TH1D*)fFile->Get("michelana/EventCutsLight");
    } else
    if( i==1 ) {
      hEventCuts = (TH1D*)fFileMC->Get("michelana/EventCuts");
      hOpEventCuts = (TH1D*)fFileMC->Get("michelana/EventCutsLight");
    }

    nEvt_TotalEvents  = hEventCuts->GetBinContent(1);
    nEvt_1StpTrk      = hEventCuts->GetBinContent(2);
    nEvt_OpID_RMS         = hOpEventCuts->GetBinContent(3);
    nEvt_OpID_2hits       = hOpEventCuts->GetBinContent(4);
    nEvt_OpID_dTmatch       = hOpEventCuts->GetBinContent(5);
    nEvt_OpID_hitWidths       = hOpEventCuts->GetBinContent(6);
    nEvt_OpID_hitTime       = hOpEventCuts->GetBinContent(7);
    nEvt_OpID_hitSat       = hOpEventCuts->GetBinContent(8);
    //nEvt_OpID         = hOpEventCuts->GetBinContent(9);
    //nEvt_Shwr2D_BndFound = hEventCuts->GetBinContent(4);
   
    Loop(i);
    
    nEvt_OpID = nEvt_OpID_muPromptCut;
    
    printf("======================================================================\n");
    if( i==0 ) printf("%s Data \n",runtag.c_str());
    if( i==1 ) printf("%s MC \n",runtag.c_str());
    printf("\n");
    printf("Calorimetry event selection\n");
    printf("------------------------------------------------------------\n");
    printf("  Total events*                       : %d\n",nEvt_TotalEvents);
    printf("------------------------------------------------------------\n");
    printf("  1 Stopping Track                    : %d (%4.1f %%)\n",nEvt_1StpTrk,          100.*CalcFracChange(nEvt_TotalEvents,            nEvt_1StpTrk));
    printf("------------------------------------------------------------\n");
    printf("  OpID: RMS cut                       : %d (%4.1f %%)\n",nEvt_OpID_RMS,             100.*CalcFracChange(nEvt_1StpTrk,         nEvt_OpID_RMS));
    printf("  OpID: 2 op hits                     : %d (%4.1f %%)\n",nEvt_OpID_2hits,             100.*CalcFracChange(nEvt_OpID_RMS,         nEvt_OpID_2hits));
    printf("  OpID: dT match                      : %d (%4.1f %%)\n",nEvt_OpID_dTmatch,             100.*CalcFracChange(nEvt_OpID_2hits,         nEvt_OpID_dTmatch));
    printf("  OpID: hit width cuts                : %d (%4.1f %%)\n",nEvt_OpID_hitWidths,             100.*CalcFracChange(nEvt_OpID_dTmatch,         nEvt_OpID_hitWidths));
 //   printf("  OpID: hit timing cut                : %d (%4.1f %%)\n",nEvt_OpID_hitTime,             100.*CalcFracChange(nEvt_OpID_hitWidths,         nEvt_OpID_hitTime));
    printf("  OpID: hit saturation cut            : %d (%4.1f %%)\n",nEvt_OpID_hitSat,             100.*CalcFracChange(nEvt_OpID_hitTime,         nEvt_OpID_hitSat));
    printf("  OpID: decay time range              : %d (%4.1f %%)\n",nEvt_OpID_timeRange,             100.*CalcFracChange(nEvt_OpID_hitSat,         nEvt_OpID_timeRange));
//    printf("  OpID: muon prompt light cut         : %d (%4.1f %%)\n",nEvt_OpID_muPromptCut,             100.*CalcFracChange(nEvt_OpID_timeRange,         nEvt_OpID_muPromptCut));
    nEvt_OpID = nEvt_OpID_timeRange;
    printf("------------------------------------------------------------\n");
    printf("  2D shower: cluster bnd found        : %d (%4.1f %%)\n",nEvt_Shwr2D_BndFound,  100.*CalcFracChange(nEvt_OpID,         nEvt_Shwr2D_BndFound));
    printf("  2D shower: cluster size             : %d (%4.1f %%)\n",nEvt_Shwr2D_ClusterSize,  100.*CalcFracChange(nEvt_Shwr2D_BndFound, nEvt_Shwr2D_ClusterSize));
//    printf("  2D shower: no extra hits            : %d (%4.1f %%)\n",nEvt_Shwr2D_ExtraHits, 100.*CalcFracChange(nEvt_Shwr2D_ClusterSize, nEvt_Shwr2D_ExtraHits));
//    printf("  2D shower: Bragg slope              : %d (%4.1f %%)\n",nEvt_Shwr2D_Slope,     100.*CalcFracChange(nEvt_Shwr2D_ClusterSize, nEvt_Shwr2D_Slope));
    printf("  2D shower: mu linearity             : %d (%4.1f %%)\n",nEvt_Shwr2D_MuLinearity,100.*CalcFracChange(nEvt_Shwr2D_Slope, nEvt_Shwr2D_MuLinearity));
    printf("  2D shower: mu dir                   : %d (%4.1f %%)\n",nEvt_Shwr2D_MuDir,     100.*CalcFracChange(nEvt_Shwr2D_MuLinearity, nEvt_Shwr2D_MuDir));
    printf("  2D shower: decay angle              : %d (%4.1f %%)\n",nEvt_Shwr2D_DecayAngle,100.*CalcFracChange(nEvt_Shwr2D_MuDir,    nEvt_Shwr2D_DecayAngle));
    printf("  2D shower: shower size              : %d (%4.1f %%)\n",nEvt_Shwr2D_ShowerSize,  100.*CalcFracChange(nEvt_Shwr2D_DecayAngle, nEvt_Shwr2D_ShowerSize));
//    printf("  2D shower: shower completeness      : %d (%4.1f %%)\n",nEvt_Shwr2D_ShowerFrac,100.*CalcFracChange(nEvt_Shwr2D_ShowerSize,    nEvt_Shwr2D_ShowerFrac));
    printf("------------------------------------------------------------\n");
    if( fRequire3DMuEndpt ) 
      printf(" (3D shower: 3D mu endpoint           : %d (%4.1f %%) )\n",nEvt_3DMuEnd, 100.*CalcFracChange(nEvt_Shwr2D_ShowerSize,  nEvt_3DMuEnd));
    else 
      nEvt_3DMuEnd = nEvt_Shwr2D_ShowerFrac;
    printf("  3D shower: showers on both planes   : %d (%4.1f %%)\n",nEvt_Shwr3D_ShowersOnBothPlanes, 100.*CalcFracChange(nEvt_Shwr2D_ShowerSize,  nEvt_Shwr3D_ShowersOnBothPlanes));
    printf("  3D shower: num 3D pts               : %d (%4.1f %%)\n",nEvt_Shwr3D_Npts,                100.*CalcFracChange(nEvt_Shwr3D_ShowersOnBothPlanes,  nEvt_Shwr3D_Npts));
    printf("  3D shower: frac 3D pts              : %d (%4.1f %%)\n",nEvt_Shwr3D_FracHits3D,            100.*CalcFracChange(nEvt_Shwr3D_Npts, nEvt_Shwr3D_FracHits3D));
    printf("  3D shower: shower centroid cut      : %d (%4.1f %%)\n",nEvt_Shwr3D_Centroid,            100.*CalcFracChange(nEvt_Shwr3D_FracHits3D, nEvt_Shwr3D_Centroid));
    //printf("  3D shower: proj. shower direction   : %d (%4.1f %%)\n",nEvt_Shwr3D_Direction,            100.*CalcFracChange(nEvt_Shwr3D_Centroid, nEvt_Shwr3D_Direction));
    printf("------------------------------------------------------------\n");
    printf("  dT cut                              : %d (%4.1f %%)\n",nEvt_dTcut,              100.*CalcFracChange(nEvt_Shwr3D_Centroid, nEvt_dTcut));
    printf("------------------------------------------------------------\n");
    printf("Overall selection                     : %d (%4.1f %%)\n",nEvt_dTcut,100.*(1.+CalcFracChange(nEvt_TotalEvents, nEvt_dTcut)));
    printf("\n");
    printf("Decay fit event selection\n");
    printf("------------------------------------------------------------\n");
    printf("  Total events*                       : %d\n",nEvt_TotalEvents);
    printf("  1 Stp trk                           : %d (%4.1f %%)\n",nEvt_1StpTrk,              100.*CalcFracChange(nEvt_TotalEvents, nEvt_1StpTrk));
    printf("  Optical ID                          : %d (%4.1f %%)\n",nEvt_OpID,                 100.*CalcFracChange(nEvt_1StpTrk,     nEvt_OpID));
    printf("  Time range                          : %d (%4.1f %%)\n",nEvt_DecayFit_TimeRange,   100.*CalcFracChange(nEvt_OpID, nEvt_DecayFit_TimeRange));
    printf("  amplitude cut                       : %d (%4.1f %%)\n",nEvt_DecayFit_AmpCut,      100.*CalcFracChange(nEvt_DecayFit_TimeRange, nEvt_DecayFit_AmpCut));
    printf("\n\n");
  }

  float axts = 0.06;
  float axls = 0.05;
  
  TH1D* h_ptr[2];
  TLine* cutLine[9];
  TLine* cutLine2[9];

  TCanvas* c_cuts = new TCanvas("c_cuts","c_cuts",1200,600);
  c_cuts->Divide(4,2);

  int index;
  float cutx1, cutx2;
  float scaleFac;

  index = 0;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_MuClusterSize[0], h_evtCut_MuClusterSize[1], fMinMuClusterSize, -1., axts, axls);
  index ++;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_ElClusterSize[0], h_evtCut_ElClusterSize[1], fMinElClusterSize, -1., axts, axls);
  index ++;
  //MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_BraggSlope[0], h_evtCut_BraggSlope[1], fMinBraggSlope, -1., axts, axls);
  //index = 3;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_MuLinearity[0], h_evtCut_MuLinearity[1], fMinMuLinearity, -1., axts, axls);
  gPad->SetLogy(1);
  h_evtCut_MuLinearity[0]->GetYaxis()->SetRangeUser(0.1,10000);
  index ++;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_MuDir[0], h_evtCut_MuDir[1], fMinMuClusterHitsEndFit, -1., axts, axls);
  index ++;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_DecayAngle[0], h_evtCut_DecayAngle[1], fMinDecayAngle2D, fMaxDecayAngle2D, axts, axls);
  index ++;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_ShowerSize[0], h_evtCut_ShowerSize[1], fMinElShowerSize, -1., axts, axls);
  index ++;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_NumPts3D[0], h_evtCut_NumPts3D[1], fMinNumPts3D, -1., axts, axls);
  index ++;
  MakeEventCutSubplot(c_cuts, cutLine[index], cutLine2[index], index, h_evtCut_FracPts3D[0], h_evtCut_FracPts3D[1], fMinFracHits3D, -1., axts, axls);
  
  
  // ------------------------------------
  TCanvas* c_mupe = new TCanvas("c_mupe","c_mupe",800,400);
  c_mupe->Divide(2,1);
  TLine* lineMuCut[2];
  MakeEventCutSubplot(c_mupe, lineMuCut[0], lineMuCut[0], 0, h_evtCut_MuPE_prompt[0][0], h_evtCut_MuPE_prompt[1][0], fMinMuPromptPE[0], -1., axts, axls);
  MakeEventCutSubplot(c_mupe, lineMuCut[1], lineMuCut[1], 1, h_evtCut_MuPE_prompt[0][1], h_evtCut_MuPE_prompt[1][1], fMinMuPromptPE[1], -1., axts, axls);


  


}

void MakeEventCutSubplot( TCanvas* c, TLine* line1, TLine* line2, int index, TH1D* h_data, TH1D* h_mc, float cutx1, float cutx2, float axts, float axls){
  float scaleFac = float(h_data->GetEntries()) / h_mc->GetEntries();
  h_mc->Scale( scaleFac );
  c->cd(index+1);
  gPad->SetMargin(0.15,0.05,0.15,0.05);
  FormatAxes(h_data,axts,axls, 1., 1.4);
  h_data->GetXaxis()->SetNdivisions(505);
  h_data->SetMaximum( std::max( GetHistMax(h_data), GetHistMax(h_mc) )*1.2 );
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.5);
  h_data->SetMarkerColor(kBlack);
  h_data->SetLineColor(kBlack);
  h_data->SetLineWidth(1);
  h_data->GetYaxis()->SetRangeUser(0, h_data->GetMaximum() );
  h_data->Draw("P E X0");
  h_mc->SetLineColor(kRed);
  h_mc->SetLineWidth(1);
  h_mc->Draw("same hist");
  h_data->Draw("same P E X0");
  gStyle->SetOptStat(0);
  line1 = new TLine();
  line1->SetX1( cutx1 );
  line1->SetX2( cutx1 );
  line1  ->SetLineStyle(2);
  line1  ->SetLineColor(kBlue);
  line1  ->SetLineWidth(2);
  line1->SetY1( gPad->GetUymin() );
  line1->SetY1( gPad->GetUymax() );
  if( cutx1 > 0 ) line1->Draw();
  if( cutx2 > 0 ) {
    line2 = (TLine*)line1->Clone("cutline2");
    line2->SetX1( cutx2 );
    line2->SetX2( cutx2 );
    line2->Draw();
  }
}


//################################################################################
float CalcFracChange(int n1, int n2){
  return float(n2-n1)/float(n1);
}


//##################################################################################
void RepMC(){
  std::cout<<"RepMC()...\n";
  Loop(fTreeMC_ns,true,true);
}



//#################################################################################
void SetTrigEffParams(int s=1){
  if( s == 0 ) {
    fTrigEff_P[0]=0.;
    fTrigEff_P[1]=0.;
  } else {
    for(size_t i=0; i<2; i++){
      trigEffCut[i] ->SetParameter(0, fTrigEff_P[i]);
      trigEffCut[i] ->SetParameter(1, fTrigEff_K[i]);
      trigEffCut[i] ->SetParameter(2, 1.);
      trigEffCut[i] ->SetParameter(3, 1.);
      if( fTrigEff_P[i] <= 0 ) trigEffCut[i]->SetParameter(2,0.);
    }
  }
}


// ====================================================================================
// Calculates the chi-squared goodness of match between two histograms. It is assumed
// that the 1st histograms (h1) is data and the 2nd (h2) is MC. 
float GetChi2(TH1D* h1, TH1D* h2, bool verbose){
  return GetChi2(h1,h2,-9.e9,9.e9, false, verbose);
}
float GetChi2Weighted(TH1D* h1, TH1D* h2, bool verbose){
  return GetChi2(h1, h2, -9.e9, 9.e9, true, verbose);
}
float GetChi2(TH1D* h1, TH1D* h2, float x1, float x2, bool useWeight=false, bool verbose=false){
  
  size_t n1 = h1->GetNbinsX();
  size_t n2 = h2->GetNbinsX();
  if( n1 != n2 ) {
    std::cout<<"!!!! GetChi2(): number of bins does not match !!!\n";
    return 999.; 
  }

  float sum = 0.;
  float sumW = 0.;
  int N = 0;
  for(size_t i=1; i<=n1; i++){
    float x_l = float(h1->GetXaxis()->GetBinLowEdge(i));
    float x_u = float(h1->GetXaxis()->GetBinLowEdge(i+1));
    if( x_l < x1 || x_u > x2 ) continue;
    float b1 = float(h1->GetBinContent(i));
    float b2 = float(h2->GetBinContent(i));
    float e1 = h1->GetBinError(i);
    float e2 = h2->GetBinError(i);
    if( b1 == 0 ) continue;
    if( b1 < 3 ) e1 = std::sqrt(2); // Data
    float sig = std::sqrt( float(std::pow(e1,2) + std::pow(e2,2)) );
    if( sig > 0 ) {
      float W = 1.;
      if( useWeight ) W = b1; 
      sum += W*std::pow( (b1-b2) / sig, 2);
      N++;
      sumW += W;
    }
  }
  
  float out = 999.;
  if( sumW > 0. ) out = sum / sumW;
  
  if( verbose) std::cout<<"GetChi2(w="<<useWeight<<"): "<<h1->GetTitle()<<":  --> Chi2/NDF= "<<out<<"\n";

  return out;
}










//###########################################################################
void Optimize(int ch){
  
  cout
  <<"Preparing to optimize optical parameters... target = "<<fOptimizeTarget<<"\n";
  
  // Refill data histograms
  Loop(0);
  
  // Initialize TMinuit function to be used in minimization
  TMinuit *ptMinuit = new TMinuit(6); 
  ptMinuit->SetPrintLevel(1);
  ptMinuit->SetFCN(fcn_trigeff);

  // Set minimization strategy (default 1)
  double arglist[6];
  int ierflg = 0;
  arglist[0] = fSetStrategy;
  ptMinuit->mnexcm("SET STR", arglist, 1, ierflg);  
  
  // Set channel and argument variables
  fCh = ch;
  
  double p0   = fTrigEff_P[ch];
  double k0   = fTrigEff_K[ch];
  double sc0  = fScaleFactor[ch];
  double smp0 = 0.01; //fSmearFactorPr[ch]; // should be fixed
  double smt0 = fSmearFactor[ch];
  double fr0  = fFastRatio;         // should be fixed

  ptMinuit->DefineParameter(0, "P",         p0,     0.1,   0.2*p0,   2.0*p0  );
  ptMinuit->DefineParameter(1, "K",         k0,     0.5,  0.2*k0,   2.0*k0  );
  ptMinuit->DefineParameter(2, "scale",     sc0,    0.005, 0.1*sc0,  1.5*sc0 );
  ptMinuit->DefineParameter(3, "smear_pr",  smp0,   0.005,  0.,      0.5     );
  ptMinuit->DefineParameter(4, "smear_tot", smt0,   0.005,  0.,      0.5     );
  ptMinuit->DefineParameter(5, "fratio",    fr0,    0.01, 0.3,      0.6     );
  if( fFixP[ch]     )       ptMinuit->FixParameter(0);
  if( fFixK[ch]     )       ptMinuit->FixParameter(1);
  if( fFixScale )       ptMinuit->FixParameter(2);
  if( fFixSmearPrompt ) ptMinuit->FixParameter(3);
  if( fFixSmearTotal )  ptMinuit->FixParameter(4);
  if( fFixFastRatio )   ptMinuit->FixParameter(5);

  // Fix params based on optimization target
  // 
  // Use prompt light to optimize trigg eff params (P,K)
  // and the overall scale factor (fix smearing params)
  if( fOptimizeTarget == "prompt") {
    ptMinuit->FixParameter(4);
  } else 
  // Use total light to set the smearing (fix all others)
  if( fOptimizeTarget == "total") {
    ptMinuit->FixParameter(0);
    ptMinuit->FixParameter(1);
    ptMinuit->FixParameter(2);
    ptMinuit->FixParameter(3);
  } else 
  // Optimize both prompt and total simultaneously
  // (no params fixed unless explicitly fixed using
  // the "fFixX" params above)
  if( fOptimizeTarget == "both" ) {
    // *shrug*
  }
 
  // Run the minimization!
  ptMinuit->SetMaxIterations(500);
  if( fSetMinimizer == "simplex" ) ptMinuit->mnsimp();
  if( fSetMinimizer == "migrad" ) ptMinuit->Migrad();
  
  // Get the results and update the parameters
  double v,ev;
  ptMinuit->GetParameter(0, v, ev);
  fTrigEff_P[ch]=float(v);
  ptMinuit->GetParameter(1, v, ev);
  fTrigEff_K[ch]=float(v);
  ptMinuit->GetParameter(2, v, ev);
  fScaleFactor[ch]=float(v);
  ptMinuit->GetParameter(3, v, ev);
  fSmearFactorPr[ch]=float(v);
  ptMinuit->GetParameter(4, v, ev);
  fSmearFactor[ch]=float(v);
  ptMinuit->GetParameter(5, v, ev);
  fFastRatio=float(v);
  
  delete ptMinuit; 
  
  RepMC();

}

void OptimizePrompt(int ch, int n = 1){
  fOptimizeTarget = "prompt";
  for(int i=0; i<n; i++) {
    std::cout<<"\n\nOPTIMIZATION -- TARGET: "<<fOptimizeTarget<<", PASS: "<<i+1<<" / "<<n<<"\n\n";
    Optimize(ch);
  }
}
void OptimizePromptBothPMTs(int n = 1){
  fOptimizeTarget = "prompt";
  for(int i=0; i<n; i++) {
    std::cout<<"\n\nOPTIMIZATION PMT 0 -- TARGET: "<<fOptimizeTarget<<", PASS: "<<i+1<<" / "<<n<<"\n\n";
    Optimize(0);
    std::cout<<"\n\nOPTIMIZATION PMT 0 -- TARGET: "<<fOptimizeTarget<<", PASS: "<<i+1<<" / "<<n<<"\n\n";
    Optimize(1);
  }
}
void OptimizeTotal(int ch, int n = 1){
  fOptimizeTarget = "total";
  float def = fMinimizeBothPMTs;
  fMinimizeBothPMTs = false;
  for(int i=0; i<n; i++) {
    std::cout<<"\n\nOPTIMIZATION -- TARGET: "<<fOptimizeTarget<<", PASS: "<<i+1<<" / "<<n<<"\n\n";
    Optimize(ch);
  }
  fMinimizeBothPMTs = def;
}
void OptimizeBoth(int ch, int n = 1){
  fOptimizeTarget = "both";
  for(int i=0; i<n; i++) {
    std::cout<<"\n\nOPTIMIZATION -- TARGET: "<<fOptimizeTarget<<", PASS: "<<i+1<<" / "<<n<<"\n\n";
    Optimize(ch);
  }
}
void OptimizeTotalOnly(int ch, int n = 1){
  fOptimizeTarget = "totalonly";
  for(int i=0; i<n; i++) {
    std::cout<<"\n\nOPTIMIZATION -- TARGET: "<<fOptimizeTarget<<", PASS: "<<i+1<<" / "<<n<<"\n\n";
    Optimize(ch);
  }
}


// Function that is to be minimized. 
//   p0 = trig eff. P
//   p1 = trig eff. K
//   p2 = scale factor
//   p3 = prompt smearing
//   p4 = total smearing
//   p5 = fast ratio
static void fcn_trigeff(int& nDim, double* gout, double& result, double *par, int flg){
  result = trigeff(par[0],par[1],par[2],par[3],par[4],par[5]); 
}

double trigeff(double p, double k, double sc, double sm, double smtot, double fratio){
  
  // Set variables
  fTrigEff_P[fCh] = p;
  fTrigEff_K[fCh] = k;
  fSmearFactorPr[fCh] = sm;
  fScaleFactor[fCh]=sc;
  fSmearFactor[fCh] = smtot;
  fFastRatio = fratio;

  // Remake the MC histograms
  // (w/ smearing and scaling)
  RepMC();





  // Define output chi2 metric. 
  double out = 0; 
  int nn = 0;

  for(size_t i=0; i<2; i++){
    
    if( !fMinimizeBothPMTs && i != fCh ) continue;
    if( !fUsePmt[i] ) continue;

    // Get integrals of PE histograms
    double n_prompt_data  = hPE_prompt[0][i]->Integral();
    double n_total_data   = hPE_total_qc[0][i]->Integral();
    
    if( n_prompt_data > 0 && (fOptimizeTarget == "prompt" || fOptimizeTarget == "both") ) {
      out += GetChi2Weighted( hPE_prompt[0][i], hPE_prompt[1][i], true ); // useWeight = true 
      nn++;
    }
    if( n_total_data > 0 && (fOptimizeTarget == "total" || fOptimizeTarget == "both" || fOptimizeTarget == "totalonly") ) {
      //out += GetChi2Weighted( hPE_total_qc[0][fCh], hPE_total_qc[1][fCh], true); // useWeight = true
      out += GetChi2Weighted( hPE_total[0][i], hPE_total[1][i], true); // useWeight = true
      nn++;
    }
  
  }
  
  out /= nn;
  
  /*
  // Define output chi2 metric. 
  double out = 0; 
  int nn = 0;
  if( n_prompt_data > 0 && (fOptimizeTarget == "prompt" || fOptimizeTarget == "both") ) {
    out += GetChi2Weighted( hPE_prompt[0][fCh], hPE_prompt[1][fCh], true ); // useWeight = true 
    nn++;
  }
  if( n_total_data > 0 && (fOptimizeTarget == "total" || fOptimizeTarget == "both" || fOptimizeTarget == "totalonly") ) {
    //out += GetChi2Weighted( hPE_total_qc[0][fCh], hPE_total_qc[1][fCh], true); // useWeight = true
    out += GetChi2Weighted( hPE_total[0][fCh], hPE_total[1][fCh], true); // useWeight = true
    nn++;
  }
  out /= nn;
  */

  std::cout<<"P= "<< p << ", K= "<<k<<"  Scale="<<sc<<"    Smear(tot)= "<<smtot<<"     chi2= "<<out<<"\n";
  return out;
}




//#############################################################################
void MapChi2Space(int ch, double sc1, double sc2, double sig1, double sig2){
  
  fCh = ch;

  int bins_x = 10;
  int bins_y = 10;
  TH2D * h_map = new TH2D("h_map","#chi^{2} surface for Monte Carlo fit to Michel PE spectrum;Scaling (#epsilon);Nominal smearing resolution (#sigma_{0})",bins_x,sc1,sc2, bins_y,sig1,sig2);
      
  double dx = (sc2-sc1)/(double)bins_x;
  double dy = (sig2-sig1)/(double)bins_y;
  
  for(int i=0; i<bins_x; i++){
    for(int j=0; j<bins_y; j++){
      double x = sc1 + (i+0.5)*dx;
      double y = sig1 + (j+0.5)*dy;
      fScaleFactor[fCh] = x;
      fSmearFactor[fCh] = y;
      RepMC();
      ScaleMC();
      double chi2 = GetChi2Weighted( hPE_total[0][ch], hPE_total[1][ch]);
        // method 2: optimize total + prompt
        //double out1 = GetChi2( hPE_prompt[0][fCh], hPE_prompt[1][fCh]);
        //double out2 = GetChi2( hPE_total[0][fCh], hPE_total[1][fCh]);
        //double chi2 = out1 + out2;
      std::cout<<"(i,j) = "<<i<<"  "<<j<<"   x/y = "<<x<<"/"<<y<<"    chi2 = "<<chi2<<"\n";
      h_map->Fill(x,y,chi2);
    }
  }

  h_map->Draw("colz");
  gStyle->SetOptStat(0);

}


//##########################################################################################
TH1D* Make1DSlice(const TH2D* h2d, int bin1, int bin2){
  return Make1DSlice(h2d,bin1,bin2,"h");
}

TH1D* Make1DSlice(const TH2D* h2d, int bin1, int bin2, std::string name)
{
  float nbins = h2d->GetYaxis()->GetNbins();
  float xmin = h2d->GetYaxis()->GetXmin();
  float xmax = h2d->GetYaxis()->GetXmax();
  float mean = (h2d->GetXaxis()->GetBinCenter(bin1) + h2d->GetXaxis()->GetBinCenter(bin2))/2.;
  name = Form("%s_%d_%d",name.c_str(),bin1,bin2);
  TH1D* h = new TH1D(name.c_str(),Form("1D slice: bins %d-%d (x= %f)",bin1,bin2,mean),nbins,xmin,xmax);
  for(int i=bin1; i<= bin2; i++){
    for(int j=1; j<nbins; j++){
      int bincontent = h2d->GetBinContent(i,j);
      for(int jj=0; jj<bincontent; jj++) h->Fill(h2d->GetYaxis()->GetBinCenter(j)); //h2d->GetBinContent(i,j));
    }
  }
  
  return h;

}


//##########################################################################################
void Get1DSliceMeanRMS(TH2D* h2d, int bin1, int bin2, float &mean, float &rms)
{
  mean = 0.;
  rms = 0.;
  float nbins = h2d->GetYaxis()->GetNbins();
  float xmin = h2d->GetYaxis()->GetXmin();
  float xmax = h2d->GetYaxis()->GetXmax();
  float sum = 0.;
  int N = 0;
  for(int i=bin1; i<= bin2; i++){
    for(int j=1; j<nbins; j++){
      N += h2d->GetBinContent(i,j);
      sum += h2d->GetBinContent(i,j)*h2d->GetYaxis()->GetBinCenter(j);
    }
  }

  if( N > 0 ) mean = sum/N;
  
  float sumsq = 0.;
  for(int i=bin1; i<= bin2; i++){
    for(int j=1; j<nbins; j++){
      for(int jj=0; jj<h2d->GetBinContent(i,j); jj++) sumsq += std::pow(h2d->GetBinContent(i,j) - mean,2);
    }
  }
  rms = std::sqrt( sumsq / N);

}

//#########################################################################################
void FitToLandau(TH1D* h,float& mpv, float& mpv_err, float& sigma, float& sigma_err){
  FitToLandau(h, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), mpv, mpv_err, sigma, sigma_err);
}
void FitToLandau(TH1D* h, float r1, float r2, float& mpv, float& mpv_err, float& sigma, float& sigma_err){
  std::cout<<"Fitting to Landau, range "<<r1<<"-"<<r2<<"\n";
  /*
  TF1 landau("landau","landau",r1,r2);
  landau.SetParameter(0,h->GetMaximum());
  landau.SetParameter(1,h->GetMean());
  landau.SetParameter(2,500);
  h->Fit("landau","R");
  
  // refit
  landau.SetRange(landau.GetParameter(1)-3.*landau.GetParameter(2), landau.GetParameter(1)+3.*landau.GetParameter(2));
  h->Fit("landau","R");

  sigma     = landau.GetParameter(2);
  sigma_err = landau.GetParError(2);
  mpv     = landau.GetParameter(1)- 0.2228*sigma;
  mpv_err = sqrt( pow(landau.GetParError(1),2) + 0.2228*pow(sigma_err,2) );
  */
  
  TF1 landau("LandauGaus",LandauGaus,r1,r2,4);
  landau.SetParameter(0,h->GetRMS()/10.); // width scale of landau
  landau.SetParLimits(0,20.,h->GetRMS());
  landau.SetParameter(1,h->GetMean()*0.90);  // mpv
  landau.SetParameter(2,h->GetEntries()*100); // total integral
  landau.SetParameter(3,h->GetRMS()); // convoluted gaussian width
  landau.SetParLimits(3,h->GetRMS()*0.01,h->GetRMS()); // convoluted gaussian width
  landau.SetParName(0, "Landau width");
  landau.SetParName(1, "Landau MPV");
  landau.SetParName(2, "Landau intgrl");
  landau.SetParName(3, "Gaus width");
  h->Fit("LandauGaus","R");
  sigma     = landau.GetParameter(0);
  sigma_err = landau.GetParError(0);
  //refit
  float r = std::sqrt(std::pow(landau.GetParameter(0),2) + std::pow(landau.GetParameter(3),2) );
  landau.SetRange(landau.GetParameter(1)-3*r, landau.GetParameter(1)+5*r);
  h->Fit("LandauGaus","R");
  sigma     = landau.GetParameter(0);
  sigma_err = landau.GetParError(0);
  mpv     = landau.GetParameter(1); // MPV already corrected in function
  mpv_err = landau.GetParError(1);
}

//#########################################################################################
void FitToGaus(TH1D* h,float& mpv, float& mpv_err, float& sigma, float& sigma_err){
  FitToGaus(h, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), mpv, mpv_err, sigma, sigma_err);
}
void FitToGaus(TH1D* h, float r1, float r2, float& mpv, float& mpv_err, float& sigma, float& sigma_err){
  TF1 gaus("gaus","gaus",r1,r2);
  gaus.SetParameter(0,h->GetMaximum());
  gaus.SetParameter(1,h->GetMean());
  gaus.SetParameter(2,500);
  h->Fit("gaus","R");
  sigma     = gaus.GetParameter(2);
  sigma_err = gaus.GetParError(2);
  mpv     = gaus.GetParameter(1);
  mpv_err = gaus.GetParError(1);
}

//#########################################################################################
void FitToLandauGaus(TH1D* h,float& mpv, float& mpv_err, float& sigma, float& sigma_err){
  FitToLandauGaus(h, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), mpv, mpv_err, sigma, sigma_err);
}
void FitToLandauGaus(TH1D* h, float r1, float r2, float& mpv, float& mpv_err, float& sigma, float& sigma_err){
  TF1 landau("landau","landau(0)+gaus(3)",r1,r2);
  landau.SetParameter(0,h->GetMaximum());
  landau.SetParameter(1,h->GetMean());
  landau.SetParameter(2,500);
  landau.SetParameter(3,h->GetMaximum());
  landau.SetParameter(4,h->GetMean());
  landau.SetParameter(5,500);
  h->Fit("landau","R");
  sigma     = landau.GetParameter(2);
  sigma_err = landau.GetParError(2);
  mpv     = landau.GetParameter(1)- 0.2228*sigma;
  mpv_err = std::sqrt( std::pow(landau.GetParError(1),2) + 0.2228*std::pow(sigma_err,2) );
  mpv     = landau.GetParameter(4);
  mpv_err = landau.GetParError(4);
  /*
  TF1 landau("LandauGaus",LandauGaus,r1,r2,4);
  landau.SetParameter(0,h->GetRMS()/10.); // width scale of landau
  landau.SetParLimits(0,20.,h->GetRMS());
  landau.SetParameter(1,h->GetMean()*0.90);  // mpv
  landau.SetParameter(2,h->GetEntries()*100); // total integral
  landau.SetParameter(3,h->GetRMS()); // convoluted gaussian width
  landau.SetParLimits(3,20.,h->GetRMS()); // convoluted gaussian width
  h->Fit("LandauGaus","R");
  sigma     = landau.GetParameter(0);
  sigma_err = landau.GetParError(0);
  mpv     = landau.GetParameter(1);
  mpv_err = landau.GetParError(1);
  */
}

// #########################################################################################
void SetBoxModelParamNames(bool flag){
  if( flag ) {
    f_ModBoxRecomb->SetParName(0, "Box Model #alpha"); 
    f_ModBoxRecomb->SetParName(1, "Box Model #beta"); //[(kV/cm)(g/cm^{2})/MeV]");
    f_ModBoxRecomb->SetParName(2, "E-field"); // [kV/cm]");
  } else {
    f_ModBoxRecomb->SetParName(0, "p0"); 
    f_ModBoxRecomb->SetParName(1, "p1"); //[(kV/cm)(g/cm^{2})/MeV]");
    f_ModBoxRecomb->SetParName(2, "p2"); // [kV/cm]");
  }
}


//########################################################################
// TODO: 
void MuRangeCalibration(int mc, int ipl){
  if( mc==0 ) MuRangeCalibration( fFile, mc, ipl );
  if( mc==1 ) MuRangeCalibration( fFileMC, mc, ipl );
}

void MuRangeCalibration(TFile* file, int mc, int ipl){
  
  float rr_offset = -0.16; 
 
  std::cout<<"Plane "<<ipl<<" (isMC="<<mc<<")\n"; 
  std::cout<<"TFile  size "<<file->GetSize()<<"\n";
  
  float a = fBoxModelA;
  float b = fBoxModelB;
  if( !fUseModBoxRecomb ) {
    a = fBirksModelA;
    b = fBirksModelk;
  }
  if( fCalFixB ) 
    f_dADCdx_from_dEdx->FixParameter(2,b);  
  else
    f_dADCdx_from_dEdx->SetParameter(2,b);  
  f_dADCdx_from_dEdx->FixParameter(1,a);
  f_dADCdx_from_dEdx->FixParameter(3,fEField); // E-field
  f_dQdx_from_dEdx->FixParameter(0,a);
  f_dQdx_from_dEdx->FixParameter(1,b);
  f_dQdx_from_dEdx->FixParameter(2,fEField); // E-field

  TH2D* h_dQdx;
  TH2D* h_dADCdx;

  if( ipl == 1 ) {
    h_dQdx = (TH2D*)file->Get("michelana/dEdx_profiles/MuResRangeVsdQdx");
    h_dADCdx = (TH2D*)file->Get("michelana/dEdx_profiles/MuResRangeVsdADCdx");
  } else 
  if( ipl == 0 ) {
    h_dQdx = (TH2D*)file->Get("michelana/dEdx_profiles/MuResRangeVsdQdx_Pl0");
    h_dADCdx = (TH2D*)file->Get("michelana/dEdx_profiles/MuResRangeVsdADCdx_Pl0");
  }

  //std::string typeName[2]={"Data","MC"};
  //std::string planeName[2]={"Induction","Collection"};
  TGraphErrors* gr = new TGraphErrors();
  
  MakeResRangeGraphs(h_dADCdx, gr);
  
  gStyle->SetOptFit(0);
  
  TCanvas* cfit = new TCanvas(Form("c_pl%d",ipl),Form("c_pl%d",ipl),650,650);
  f_dADCdx_from_dEdx->SetParameter(0, fCalAreaConstants[mc][ipl]);
  gr->Fit("dADCdx_from_dEdx");
  gPad->SetMargin(0.15,0.05,0.15,0.05); // left, right, bottom, top
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.7);
  gr->GetXaxis()->SetTitle("(dE/dx)_{hyp} [MeV/cm]");
  gr->GetYaxis()->SetTitle("dQ/dx [ADC/cm]");
  gr->GetXaxis()->SetTitleSize(0.045);
  gr->GetYaxis()->SetTitleSize(0.045);
  gr->GetYaxis()->SetTitleOffset(1.4);
  gr->Draw("AP");

  std::string tag;
  if( mc == 0 ) tag = Form("%s Data",runtag.c_str());
  if( mc == 1 ) tag = Form("%s MC",runtag.c_str());
//  if( mc == 1 ) tag = "MC";
   

  std::string planeTag;
  if( ipl == 0 ) planeTag = "Induction plane";
  if( ipl == 1 ) planeTag = "Collection plane";

  TPaveText* head = MakeTextBox(0.17,0.93,0.035,3);
  head  ->AddText(Form("#bf{LArIAT %s}",tag.c_str()));
  head  ->AddText(    "Stopping cosmic #mu" );
  head  ->AddText(Form("%s",planeTag.c_str()));
  head  ->Draw();
  
  TLegend* tl  = MakeLegend(0.17, 0.82, 0.035, 5);
  tl          ->AddEntry(gr,"Reco","PL");
//  tl          ->AddEntry(f_dADCdx_from_dEdx,"#splitline{Fit w/ fixed Mod. Box}{params from ArgoNeuT}","L");
//  tl          ->AddEntry(f_dADCdx_from_dEdx,"params from ArgoNeuT","");
  tl          ->AddEntry(f_dADCdx_from_dEdx,"Fit w/ fixed Mod. Box","L")->SetTextColor(kRed);
  tl          ->AddEntry(f_dADCdx_from_dEdx,"params from ArgoNeuT","")->SetTextColor(kRed);
  tl          ->AddEntry(f_dADCdx_from_dEdx,Form("#it{#chi}^{#it{2}} / ndf = %4.2f",f_dADCdx_from_dEdx->GetChisquare()/(f_dADCdx_from_dEdx->GetNDF()-1)),"")->SetTextColor(kRed);
  tl          ->AddEntry(f_dADCdx_from_dEdx,Form("#it{C}_{#it{cal}} = %7.5f(%2.0f) ADC/e^{-}",f_dADCdx_from_dEdx->GetParameter(0),1e5*f_dADCdx_from_dEdx->GetParError(0)),"")->SetTextColor(kRed);
  tl          ->Draw();

  

}

//###########################################################################################
void MakeResRangeGraphs(TH2D* h2d, TGraphErrors* gr ) {

  float rr_offset = -0.16; 
  float rr_err0 =0.35; // 0.26
  float rr1 = fCalRR1;
  float rr2 = fCalRR2;
  size_t binwin = 2;

  int N = h2d->GetXaxis()->GetNbins();
  float drr = (0.5+float(binwin))*h2d->GetXaxis()->GetBinWidth(1);
  
  TCanvas* cgr[100];
  TH1D* h[100];
  int index = 0;

  // Determine first bin
  size_t binStart = binwin+1;
  for(size_t i=binStart; i<N; i++){
    if( h2d->GetXaxis()->GetBinCenter(i-binwin) >= rr1 ) {
      binStart = i;
      break;
    }
  }

  // Loop over bins
  for(size_t i=binStart; i<N; i += binwin*2+1){
    
    if( i+binwin*2+1 > N ) continue;
   
    // Get res range
    float rr = h2d->GetXaxis()->GetBinCenter(i);
    
    std::cout<<"  bin "<<i<<"   rr "<<rr<<"\n";

    if( rr >= rr1 && rr <= rr2 ) {
     
      // Apply offset
      rr -= rr_offset; 
      float mpv, mpv_err, sigma, sigma_err, mean, rms;
      float llfac = 2.; 
      float ulfac = 10.; 

      float rr_err = rr_err0;
      std::cout<<"rr_err0 = "<<rr_err0<<"    drr= "<<drr<<"    total err = "<<rr_err<<"\n";
      float ra = std::max((float)0.1, rr - rr_err);
      float rb = (rr + rr_err);
      float dEdxHyp     = 0.5*(MuMPV_from_rr(rr+drr) + MuMPV_from_rr(rr-drr));
      float dEdxHyp_err = 0.5*fabs(MuMPV_from_rr(ra) - MuMPV_from_rr(rb));
      std::cout<<"    dEdxHyp "<<dEdxHyp<<"   err: "<<dEdxHyp_err<<"\n";
        
//      h[index] = h2d->ProjectionY( "h", i-binwin, i+binwin ); 
      h[index] = Make1DSlice( h2d, i-binwin, i+binwin );
      mean = h[index]->GetMean();
      rms = h[index]->GetRMS();
      float max_x= h[index]->GetBinCenter(h[index]->GetMaximumBin());
      //FitToLandau(h[index], mean - llfac*rms, mean + ulfac*rms, mpv, mpv_err, sigma, sigma_err);
      FitToLandau(h[index], max_x - llfac*rms, max_x + ulfac*rms, mpv, mpv_err, sigma, sigma_err);
      float dQdx = mpv;
      float dQdx_err = mpv_err;
      if( dQdx > 0. ) {
          std::cout<<"ADDING POINT "<<dEdxHyp<<", "<<dQdx<<"\n";
          gr->SetPoint     (gr->GetN(),   dEdxHyp,      dQdx); 
          gr->SetPointError(gr->GetN()-1, dEdxHyp_err,  dQdx_err); 
      }
      
      if( index < 99 ) {
        std::cout<<"Drawing canvas\n"; 
        cgr[index]  = new TCanvas(Form("cgr_%d",index),Form("dADC/dx, rr = %5.2f",rr),600,600);
        gStyle->SetOptStat(1110);
        h[index]->Draw();
      }
      index++;
      
    } // endif rr is within range
  
  
  } // done looping over bins

}


//###########################################################################################
float CalcTruncatedMean( std::vector<float>& v, float p, int skew = 0){

 // Sort the vector
 if( p > 0 ) std::sort( v.begin(), v.end() );

  int N = v.size();

  // Calculate truncated mean
  int ntrim = p*N;
  if( N < 3 ) ntrim = 0;
  float sum = 0.;
  int NN = 0;
  int low = ntrim;
  int high = N - ntrim;
  if( skew == -1 ) low  = 0;
  if( skew ==  1 ) high = N;
  for(int i=0; i < N; i++){
    if( i >= low && i < high ){
      sum += v.at(i);
      NN++;
    }
  }

  if(NN > 0){
    return sum / NN;
  } else {
    return 0.;
  }

}


// =========================================================
// Modified Box Model recombination survival frac
//  x       = dE/dx
//  par[0]  = alpha (=0.93)
//  par[1]  = beta (=0.212)
//  par[2]  = Efield
double _ModBoxRecomb(double* x, double* par){
  double dEdx   = x[0];
  double alpha  = par[0];
  double beta   = par[1];
  double Efield = par[2];
  return ModBoxRecomb(dEdx, alpha, beta, Efield );
}

// =========================================================
// Modified Box Model recombination survival frac
//  x       = E-field
//  par[0]  = alpha (=0.93)
//  par[1]  = beta (=0.212)
//  par[2]  = dE/dx
double _ModBoxRecomb_Efield(double* x, double* par){
  double dEdx   = par[2];
  double alpha  = par[0];
  double beta   = par[1];
  double Efield = x[0]/1000.;
  return ModBoxRecomb(dEdx, alpha, beta, Efield );
}

double ModBoxRecomb(double dEdx, double alpha, double beta, double Efield){
  double xi = (beta*dEdx)/(fLArDensity*Efield);
  if( xi > 0. )
    return log(alpha+xi)/xi;
  else 
    return 0.;
}

// =========================================================
// Birks model recombination survival frac
//  x       = dE/dx
//  par[0]  = A (=0.80)
//  par[1]  = k (=0.0486)
//  par[2]  = Efield
double _BirksRecomb(double* x, double* par){
  double dEdx   = x[0];
  double A      = par[0];
  double k      = par[1];
  double Efield = par[2];
  return BirksRecomb(dEdx,A,k,Efield);
}

double _BirksRecomb_Efield(double* x, double* par){
  double dEdx   = par[2];
  double A      = par[0];
  double k      = par[1];
  double Efield = x[0]/1000.;
  return BirksRecomb(dEdx,A,k,Efield);
}

double BirksRecomb(double dEdx, double A, double k, double Efield){
  if( Efield > 0. ) 
    return A/(1. + (dEdx*k)/(fLArDensity*Efield));
  else 
    return 0.;
}




// =========================================================
// Log liklihood for optimized Q+L energy computation
double _logLEnergyQL(double* x, double* par){
  
  // Object: given a Q and an L, determine the most 
  // probable energy that produced that combination.
  // This incorporates resolution effects to use the
  // given information in the most informed way.

  // Input energy [MeV]
  double E = x[0];

  // Q, L
  double Qion = par[0];
  double Qph  = par[1];
  double L    = par[2];
  double Q    = Qion;
 
  if( Q == 0 ) {
    std::cout<<" Bad Q \n";
    return -9.;
  }
  
  // parameters controlling division of 
  // energy into light and charge
  double alpha  = par[3]; // excitation ratio (=0.21)
  double R      = par[4]; // Assumed recombination
  double W      = par[5]; // ionization energy of LAr

  // parameters affecting Q resolution
  double aveT   = par[6]; // average drift time of given event.
  double tau    = par[7]; // electron lifetime.
  double elRes  = par[8]; // expected jitter due to electronics
                          // (unit: frac % per e-).

  // parameters affecting L resolution
  double vis    = par[9];   // ave photon visibility for event.
  double peRes  = par[10];   // expected jitter due to electronics
                            // or in other words, the phoeoelectron
                            // resolution (unit: frac % per pe).
  
  // expected mean values
  R = (Qion*fRecombIon + Qph*fRecombPh)/Q;
  double q_mean = (E/W)*R;
  double l_mean = q_mean*( (1.+alpha)/R - 1. ); 

  // probabilities for given Q,L
  // (initialize to zero, meaning
  // both are at infinity)
  double probQ, probL;
  double minProb = 1e-300;
  double LogPq = TMath::Log(minProb);
  double LogPl = TMath::Log(minProb);
  
  // --------------------------------------------
  // Now for the modeling of Q and L. By default,
  // each is treated as a simple Gaussian with 
  // widths set based on given information about
  // drift time, tau, visibility, and detector
  // resolutions. Alternatively, we can use custom
  // functions.
  double useResFuncQ = par[11];
  double useResFuncL = par[12];

  if( useResFuncQ == 1 ) {

    //std::cout<<"Calling energy: "<<E<<"\n";
    //fFuncP_resQ->SetParameter("energy",E);
    //probQ = fFuncP_resQ->Eval( (Q-q_mean)/q_mean );
    fFuncP_Q->SetParameter("energy",E);
    probQ = fFuncP_Q->Eval( Q );
  
  } else {
  
    // statistical
    double q_width_stat = sqrt( q_mean );
    // from drift attenuation
    double rsurv = exp(-aveT/tau);
    double q_width_drift = sqrt( q_mean*rsurv*(1.- rsurv)) / rsurv;
    // from variation of R within event
    double q_width_recomb = q_mean*std::sqrt(   Qion*std::pow(fRecombIon-R,2)/Q 
                                              + Qph*std::pow(fRecombPh-R,2)/Q    );
    // from electronic jitter
    double q_width_elec = elRes*sqrt( q_mean );
    // Add up the contributing factors 
    double q_width = sqrt(  std::pow(q_width_stat,2) 
                          + std::pow(q_width_recomb,2)
                          + std::pow(q_width_drift,2) 
                          + std::pow(q_width_elec,2) 
                          );
  
    probQ = TMath::Gaus(Q, q_mean, q_width ); 

  }

  if( useResFuncL == 1 ) {
    
    // Set parameters 
    //fFuncP_resL->SetParameter("energy",E);
    
    //probL = fFuncP_resL->Eval( (L-l_mean)/l_mean );
    fFuncP_L->SetParameter("energy",E);
    probL = fFuncP_L->Eval( L );

  } else {
  
    // statistical
    double l_width_stat = sqrt( l_mean );
    // from propagation 
    double l_width_prop = sqrt( l_mean*vis*(1.-vis) ) / vis; 
    // recombination variation (will be as for charge
    // due to complementarity)
    double l_width_recomb = q_mean*std::sqrt(   Qion*std::pow(fRecombIon-R,2)/Q 
                                              + Qph*std::pow(fRecombPh-R,2)/Q    );
    // from quenching
    //double l_width_quench = 0.;
    // from electronic jitter
    double l_width_elec = peRes*sqrt( l_mean );
    // add up all the contributing factors
    double l_width = sqrt(  std::pow(l_width_stat,2) 
                        + std::pow(l_width_prop,2) 
                        + std::pow(l_width_recomb,2) 
                        + std::pow(l_width_elec,2) 
                        );
  
    //if( l_width == 0 || l_width/l_mean < 0.01 ) l_width = 0.01*l_mean;
    probL = TMath::Gaus( L, l_mean, l_width );
    
  }
   
  // Safeguard for cases where P was 0
//  if( LogPq > 0 ) LogPq = TMath::Log(0);
//  if( LogPl > 0 ) LogPl = TMath::Log(0);
//  if( !std::isfinite( LogPq )) LogPq = TMath::Log(minProb);
//  if( !std::isfinite( LogPl )) LogPl = TMath::Log(minProb);
  if( probQ > minProb )  LogPq = TMath::Log( probQ );
  if( probL > minProb )  LogPl = TMath::Log( probL ); 
 
  if( L <= 0 ) LogPl = 0.;

  if( fMaxMCEvts > 0 && fMaxMCEvts < 200 ) {
  std::cout
  <<"---------------------------\n"
  <<"  Call to LogLEnergyQL(): Q = "<<Q<<"    L = "<<L<<"\n"
  <<"  E = "<<E<<", R: "<<R<<"\n"
  <<"  fp_Q integral = "<<fFuncP_Q->Integral(-5000e3,5000e3)<<"\n"
  <<"  fp_L integral = "<<fFuncP_L->Integral(-5000e3,5000e3)<<"\n"
  //<<"  Q_pred = "<<q_mean<<" +/- "<<q_width<<"\n"
  //<<"  L_pred = "<<l_mean<<" +/- "<<l_width<<"\n"
  //<<"  Width contributions to Q: "<<q_width_stat<<"   "<<q_width_recomb<<"  "<<q_width_drift<<"  "<<q_width_elec<<"\n"
  //<<"  Width contributions to L: "<<l_width_stat<<"   "<<l_width_recomb<<"  "<<l_width_prop<<"  "<<q_width_elec<<"\n"
  //<<"  P(Q) = "<<TMath::Gaus(Q, q_mean, q_width)<<"\n"
  //<<"  P(L) = "<<TMath::Gaus(L, l_mean, l_width)<<"\n"
  <<"  -LogP(Q): "<<-1.*LogPq<<"\n"
  <<"  -LogP(L): "<<-1.*LogPl<<"\n"
  <<"  -------------------------\n"; 
  }

  // return the negative log-liklihood
  return -1.*LogPq -1.*LogPl;
 
}


double _funcP_Q(double* x, double* par){
  
  double Q = x[0];

  // hypothesized energy
  double E = par[0];
  
  // scale factor (default set to 1)
  double scale = par[1];
 
  // the parameters
  double muPeak = funcP_Q_fp[0]->Eval(E);
  double sigPeak = muPeak*funcP_Q_fp[1]->Eval(E);
  double rBG    = funcP_Q_fp[2]->Eval(E);
  double muBG  = muPeak*funcP_Q_fp[3]->Eval(E);
  double sigBG = sigPeak*funcP_Q_fp[4]->Eval(E);
 
  // limit BG normalization
  if( rBG <= 0 ) {
    rBG = 0;
    muBG = 0;
    sigBG = 1; 
  } else 
  if ( E > 50 ) {
    rBG = funcP_Q_fp[2]->Eval(50);
  } else 
  if( E < 10 ) {
    rBG = funcP_Q_fp[2]->Eval(10);
  }
  
  funcP_Q->SetParameter(0,muPeak);
  funcP_Q->SetParameter(1,sigPeak);
  funcP_Q->SetParameter(2,rBG);
  funcP_Q->SetParameter(3,muBG);
  funcP_Q->SetParameter(4,sigBG);
  return scale*funcP_Q->Eval( x[0] );  

}


double _funcP_L(double* x, double* par){
  
  // hypothesized energy
  double E = par[0];
  
  // scale factor (default set to 1)
  double scale = par[1];
  
  // the parameters
  double muPeak = funcP_L_fp[0]->Eval(E);
  double sigPeak = muPeak*funcP_L_fp[1]->Eval(E);
  
  funcP_L->SetParameter(0,muPeak);
  funcP_L->SetParameter(1,sigPeak);
  return scale*funcP_L->Eval( x[0] );  

}

double _funcP_Q_Gaus(double* x, double* par){
  double E = par[0];
  double scale = par[1];
  double muPeak = funcP_Q_Gaus_fp[0]->Eval(E);
  double sigPeak = muPeak*funcP_Q_Gaus_fp[1]->Eval(E);
  funcP_Q_Gaus->SetParameter(0,muPeak);
  funcP_Q_Gaus->SetParameter(1,sigPeak);
  return scale*funcP_Q_Gaus->Eval( x[0] );  
}

double _funcP_L_Gaus(double* x, double* par){
  double E = par[0];
  double scale = par[1];
  double muPeak = funcP_L_Gaus_fp[0]->Eval(E);
  double sigPeak = muPeak*funcP_L_Gaus_fp[1]->Eval(E);
  funcP_L_Gaus->SetParameter(0,muPeak);
  funcP_L_Gaus->SetParameter(1,sigPeak);
  return scale*funcP_L_Gaus->Eval( x[0] );  
}

/* 
double _funcP_resQ(double* x, double* par){

  double E = par[0];
  
  // Set parameters
  double sigPeak = funcP_resQ_fp[0]->Eval(E);
  double rBG = funcP_resQ_fp[1]->Eval(E);
  double muBG = funcP_resQ_fp[2]->Eval(E);
  double sigBG = funcP_resQ_fp[3]->Eval(E);

  // limit BG normalization
  if( rBG <= 0 ) {
    rBG = 0;
    muBG = 0;
    sigBG = 1; 
  } else 
  if ( E > 50 ) {
    rBG = funcP_resQ_fp[1]->Eval(50);
  }
  
  // BG can't be narrower than the peak!
  if( sigBG < sigPeak ) {
    sigBG = sigPeak;
  }
  
  funcP_resQ->SetParameter(0,sigPeak);
  funcP_resQ->SetParameter(1,rBG);
  funcP_resQ->SetParameter(2,muBG);
  funcP_resQ->SetParameter(3,sigBG);
  return funcP_resQ->Eval( x[0] );  

}

double _funcP_resL(double* x, double* par){
  double E = par[0];
  // Set parameters 
  double param[funcP_resL_fp.size()];
  for(size_t i=0; i<funcP_resL_fp.size(); i++)
    param[i] = funcP_resL_fp[i]->Eval(E);
  // name them
  double sig  = param[0];
  funcP_resL->SetParameter(0,sig);
  return funcP_resL->Eval( x[0] );  
}
*/  



// ========================================================
// Calculate dADC/dx expected at given muon res. range
//  x       = rr
//  par[0]  = ADCs per electron
double dADCdx_from_rr(double* x, double *par) {
  float rr = x[0];
  float AdcPerElectron = par[0];
  float dEdx = f_muStopPower->Eval(rr);
  float R = f_ModBoxRecomb->Eval(dEdx);
  if( !fUseModBoxRecomb ) R = f_BirksRecomb->Eval(dEdx);
  float dQdx = (dEdx/fWion)*R;
  return AdcPerElectron * dQdx;
}

// =======================================================
// Calculate dADC/dx expected from given dE/dx
//  x       = dE/dx
//  par[0]  = ADCs per electron
//  par[1]  = alpha or A
//  par[2]  = beta or k
//  par[3]  = EField 
double dADCdx_from_dEdx(double* x, double *par) {
  double AdcPerElectron  = par[0];
  double params[3];
  params[0]             = par[1];
  params[1]             = par[2];
  params[2]             = par[3];
  return AdcPerElectron * dQdx_from_dEdx(x,params);
}

// =======================================================
// Calculate dQ/dx expected from given dE/dx
//  x       = dE/dx
//  par[0]  = alpha or A
//  par[1]  = beta or k
//  par[2]  = EField 
double dQdx_from_dEdx(double* x, double *par) {
  double dEdx            = x[0];
  double params[3];
  params[0]             = par[0];
  params[1]             = par[1];
  params[2]             = par[2];
  double R = 0.;
  if( fUseModBoxRecomb ){
    R = f_ModBoxRecomb->EvalPar( x, params);
  } else {
    R = f_BirksRecomb->EvalPar( x, params);
  }
  return (dEdx/fWion) * R;
}

// ##############################################################################
// Values from: https://www.sciencedirect.com/science/article/pii/S0092640X01908617,
//              pg 237
void ParameterizeMuStoppingPower(){
  
  // Values directly from paper:
  //   ionization density in MeV cm^2/g (must multiple by density of 1.396
  //double id[10]= { 5.696,    4.466,    3.505,    2.732,    2.341,    1.771,    1.670,    1.570,    1.519 };
  //   res. range in units g/cm^2 (must divide by density)
  //double rr[10]= { 9.937e-1, 1.795,    3.329,    6.605,    1.058e1,  3.084e1,  4.250e1,  6.732e1,  1.063e2};

  // Different values seen here:
  // http://pdg.lbl.gov/2017/AtomicNuclearProperties/MUE/muE_liquid_argon.txt
  double p[18]  = { 10.,      12.,    14.,    17.,    20.,    25.,    30.,    35.,    40.,      45.,      50.,      55.,      60.,      70.,    80.,    90.,    100.,   120.  };  
  double id[18] = { 5.687,    4.979,  4.461,  3.901,  3.502,  3.042,  2.731,  2.508,  2.34,     2.21,     2.107,    2.023,    1.954,    1.848,  1.771,  1.713,  1.669,  1.608};
  double rr[18] = { 9.833e-1, 1.360,  1.786,  2.507,  3.321,  4.859,  6.598,  8.512,  1.058e1,  1.278e1,  1.510e1,  1.752e1,  2.004e1,  2.531e1,3.084e1,3.659e1,4.250e1,5.473e1}; 
  for(size_t i=0; i<18; i++){
    id[i] *= fLArDensity;
    rr[i] /= fLArDensity;
  }

  TCanvas* cmustop = new TCanvas("cmustop","cmustop",500,500);  
  TGraphErrors* gr = new TGraphErrors();
  for(size_t i=0; i<18; i++){
    gr->SetPoint(gr->GetN(), rr[i], id[i] );
    gr->SetPointError(gr->GetN()-1, rr[i]*0.001, id[i]*0.001 );
  }
  f_muStopPower->SetRange(1.,30.);
  gr->Fit("muStopPwr","R");
  gStyle->SetOptFit(1);
  gr->Draw("AP");
  gr->SetMarkerStyle(22);
  gr->GetXaxis()->SetRangeUser(0,30);
  gr->GetXaxis()->SetTitle("Muon residual range, #it{r} [cm]");
  gr->GetYaxis()->SetTitle("Mean dE/dx [MeV/cm]");
  //f_muStopPower->Draw("same");
 
  TLegend* tleg  = MakeLegend(0.5, 0.9, 0.045, 3);
  tleg         ->AddEntry(gr,"PDG values","P");
  tleg        ->AddEntry(f_muStopPower,"Parameterized fit","PL");
  tleg        ->AddEntry(f_muStopPower,"Ar^{b}+C","");
  tleg         ->Draw();
  
  
  TCanvas* cmustop2 = new TCanvas("cmustop2","cmustop2",500,500);  
  TGraphErrors* gr2 = new TGraphErrors();
  for(size_t i=0; i<18; i++){
    gr2->SetPoint(gr2->GetN(), rr[i], p[i] );
    gr2->SetPointError(gr2->GetN()-1, rr[i]*0.001, p[i]*0.001);
  }
  gr2->Draw("AP");

}





// #################################################################################
// Function from Limetime_module
// #################################################################################

double LandauGaus(Double_t *x,Double_t *par) {
    
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.


    // Numeric constants
    const Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    const Double_t mpshift  = -0.22278298;       // Landau maximum location
    
    // Control constants
    //      Double_t np = 18000.0;      // number of convolution steps
    //      Double_t sc =    8.0;      // convolution extends to +-sc Gaussian sigmas
    
    const Double_t np =   1000.0;      // number of convolution steps
    const Double_t sc =    8.0;      // convolution extends to +-sc Gaussian sigmas
    
    // Variables
    double xx;
    double mpc;
    double fland;
    double sum = 0.0;
    double xlow,xupp;
    double step;
    double i;
    
    // MP shift correction
    mpc = par[1] - mpshift * par[0];
    
    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    
    step = (xupp-xlow) / np;
    
    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
    
    return (par[2] * step * sum * invsq2pi / par[3]);
}

//#######################################################################################
// Crossing muon calibration test:
void CrsMuTest(){
  gStyle->SetOptFit(1);
  
  TFile* crsMuFile = new TFile("files/MichelAna_mc2b_crsmu.root","read");
  hdEdx_crsmu = (TH1D*)crsMuFile->Get("michelana/crsmu/CrsMu_dEdx")->Clone("crsmu_dedx");
  TH1D* hdEdx_true = (TH1D*)crsMuFile->Get("michelana/crsmu/True_dEdx_CrsMu");
  
  std::cout<<"True dE/dx = "<<hdEdx_true->GetMean()<<" +/- "<<hdEdx_true->GetRMS()/std::sqrt(hdEdx_true->GetEntries())<<" MeV/cm\n";
  
  TCanvas* cCrsMuTest = new TCanvas("crsMuTest","crsMuTest",500,500);
  cCrsMuTest->cd();

  float mpv = hdEdx_true->GetMean();
  float mpv_err;
  float sigma = 0.2;
  float sigma_err;
  FitToLandau(hdEdx_crsmu, 0., 5., mpv, mpv_err, sigma, sigma_err);
  hdEdx_crsmu->SetLineColor(kBlack);
  hdEdx_crsmu->SetMarkerColor(kBlack);
  hdEdx_crsmu->Draw("axis hist");
  
  //crsMuFile->Close();
  

}

//#######################################################################################
float CorrectForQuenching(float L_tot, float L_prompt, float tau){
  float quenchCorrFactor = 1.0;
  float estimatedTau    = TauConversion((*fMuContamCorr_EffTau)[1]);
  // scale up light to compensate for quenching
  if( estimatedTau < fTauT ) quenchCorrFactor = (fTauT*exp(-100./fTauT)) / (estimatedTau*exp(-100./estimatedTau));
  float lateLight = (L_tot - L_prompt)*quenchCorrFactor;
  return L_prompt + lateLight;
}

//#####################################
float TauConversion( float tau_eff ) {
  return 0.965*tau_eff - 41.7;
}

//######################################################################################
double MuMPV_from_rr( float rr ) {
  float A = f_muStopPower->GetParameter(0);
  float b = f_muStopPower->GetParameter(1);
  float C = f_muStopPower->GetParameter(2);
  float m = 105.65839;
  float T = (A/(b+1)) * std::pow(rr,b+1) + C*rr;
  float P = m*std::sqrt( std::pow(T/m + 1.,2)-1. );
  float mpv = f_MPV_mu->Eval(P);
  return mpv;
}



//######################################################################################
void FieldScanAna() {
 
  float total_scaleFac = 0.7;
  
  float textSize = 0.035;
  float mar_l  = 0.15;
  float mar_r  = 0.05;
  float mar_t  = 0.05;
  float mar_b  = 0.15;
  float leg_x1 = 0.50;
  float leg_x2 = 1.-mar_r-0.02;
  float leg_y1 = 0.7;
  float leg_y2 = 1.-mar_t-0.02;
  float text_x1 =  mar_l+0.02;
  float text_x2 =  0.4;
  float text_y1 = 0.75;
  float text_y2 = 1.-mar_t-0.03;
  float axisTitleSize = 0.045;
  
  std::cout<<"FIELD SCAN ANA\n";

  // First use MC to figure out what scale factors to apply to Q and L
  Loop(1);
  fQCorrFactor = 1./(1.+hQRes->GetMean());
  fLCorrFactor = 1./(1.+hLRes->GetMean());
  std::cout<<"QCorrFactor = "<<fQCorrFactor<<"\n";
  std::cout<<"LCorrFactor = "<<fLCorrFactor<<"\n";
  


  // No longer limit runs 
  int maxRunOrig = fMaxRun; 
  if( fMaxRun > 0 ) fMaxRun *= -1;


  const int Nfields = 6; 
  TFile* fieldScanFile[Nfields];
  TTree* fieldScanTree[Nfields];
  float fieldNom[Nfields]={ 270, 300, 350, 400, 450, 500  };
  float fieldEff[Nfields]={ 260, 288, 337, 386, 434, 487  };

  TGraphErrors* gr_murange[Nfields];
  auto gr_Q = new TGraphErrors();
  auto gr_L = new TGraphErrors();
  auto gr_T = new TGraphErrors();
  auto gr_R = new TGraphErrors();

  // Read in each fieldscan file
  for(size_t i=0; i<Nfields; i++){

    std::cout<<i<<"\n";

    fieldScanFile[i]    = new TFile(Form("fscan/MichelAna_%3.0fVcm.root",fieldNom[i]), "read");
    fieldScanTree[i]    = (TTree*)fieldScanFile[i]->Get("michelana/anatree");
    setBranches(fieldScanTree[i]);
    Loop(fieldScanTree[i],0,0);
    
    // -----------------------------------------------------------------
    // Make muon range plot (coll. plane)
    TCanvas* c_muRange = new TCanvas(
      Form("muRange_fieldScan_%3.0fVcm",fieldNom[i]),
      Form("muRange_fieldScan_%3.0fVcm",fieldNom[i]),
      500,500 );
    
    TH2D* h_dQdx = (TH2D*)fieldScanFile[i]->Get("michelana/dEdx_profiles/MuResRangeVsdQdx");
    gr_murange[i] = new TGraphErrors();
    MakeResRangeGraphs(h_dQdx, gr_murange[i]);
    c_muRange->cd();
    gPad->SetMargin(0.15,0.05,0.15,0.05);
    f_dQdx_from_dEdx->FixParameter(2,fieldEff[i]/1000.); 
    f_dQdx_from_dEdx->FixParameter(0,f_dQdx_from_dEdx->GetParameter(0));
    f_dQdx_from_dEdx->FixParameter(1,f_dQdx_from_dEdx->GetParameter(1));
    gr_murange[i]->Fit("dQdx_from_dEdx");
    gStyle->SetOptFit(1112);
    gr_murange[i]->SetMarkerStyle(20);
    gr_murange[i]->SetMarkerSize(0.7);
    gr_murange[i]->SetTitle(Form("%3.0f V/cm",fieldEff[i]));
    gr_murange[i]->GetXaxis()->SetTitle("(dEdx)_{hyp} [MeV/cm]");
    gr_murange[i]->GetYaxis()->SetTitle("dQ/dx [e-/cm]");
    gr_murange[i]->GetXaxis()->SetTitleSize(0.045);
    gr_murange[i]->GetYaxis()->SetTitleSize(0.045);
    gr_murange[i]->GetYaxis()->SetTitleOffset(1.4);
    gr_murange[i]->Draw("AP");
    

    // -----------------------------------------------------------------
    float Q = hQ[0]->GetMean();
    float dQ = hQ[0]->GetRMS()/std::sqrt(hQ[0]->GetEntries());
    float L = hL[0]->GetMean();
    float dL = hL[0]->GetRMS()/std::sqrt(hL[0]->GetEntries());
    float T = Q+L;
    float dT = std::sqrt(dQ*dQ + dL*dL);
    float R = hR[0]->GetMean();
    float dR = hR[0]->GetRMS()/std::sqrt( hR[0]->GetEntries());
    gr_Q->SetPoint(       gr_Q->GetN(),   fieldEff[i],  Q );
    gr_Q->SetPointError(  gr_Q->GetN()-1, fieldEff[i]*0.01,          dQ );
    gr_L->SetPoint(       gr_L->GetN(),   fieldEff[i],  L   );
    gr_L->SetPointError(  gr_L->GetN()-1, fieldEff[i]*0.01,          dL );
    gr_T->SetPoint(       gr_T->GetN(),   fieldEff[i],  T*total_scaleFac   );
    gr_T->SetPointError(  gr_T->GetN()-1, fieldEff[i]*0.01,          dT*total_scaleFac );
    gr_R->SetPoint(       gr_R->GetN(),   fieldEff[i],  R   );
    gr_R->SetPointError(  gr_R->GetN()-1, fieldEff[i]*0.01,          dR );
    
    std::cout
    <<"--------------------------------\n"
    <<"Field              : "<<fieldNom[i]<<" V/cm\n"
    <<"Entries in Q histo : "<<hQ[0]->GetEntries()<<"\n"
    <<"Entreis in L histo : "<<hL[0]->GetEntries()<<"\n"
    <<"Q                  : "<<Q<<"\n"
    <<"L                  : "<<L<<"\n"
    <<"T                  : "<<T<<"\n";
   
    fieldScanTree[i]->Delete(); 
    fieldScanFile[i]->Close();
  }

  // Format each graph
  gr_Q->SetTitle("Michel shower charge (Q)");
  gr_L->SetTitle("Michel shower light (L)");
  gr_T->SetTitle("Total (Q+L)");
  gr_Q->SetMarkerStyle(20);
  gr_Q->SetMarkerColor(kBlue);
  gr_Q->SetLineColor(kBlue);
  gr_L->SetMarkerStyle(20);
  gr_L->SetMarkerColor(kRed+2);
  gr_L->SetLineColor(kRed+2);
  gr_T->SetMarkerStyle(20);
  gr_T->SetMarkerColor(kViolet+2);
  gr_T->SetLineColor(kViolet+2);
  gr_R->SetMarkerStyle(20);
  gr_R->SetMarkerColor(kBlack);
  gr_R->SetLineColor(kBlack);
  gr_R->SetFillColor(kWhite);
 
  // -------------------------------------------------------------
  // Plot Q, L, and Q+L
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetLeftMargin(mar_l);
  gPad->SetRightMargin(mar_r);
  gPad->SetTopMargin(mar_t);
  gPad->SetBottomMargin(mar_b);
  auto mg = new TMultiGraph();
  mg->Add(gr_Q, "APL");
  mg->Add(gr_L, "APL");
  mg->Add(gr_T, "APL");
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("Electric field [V/cm]");
  mg->GetYaxis()->SetTitle("Reconstructed quanta");
  mg->GetYaxis()->SetRangeUser(750e3,1650e3);
  mg->GetXaxis()->SetTitleSize(axisTitleSize);
  mg->GetYaxis()->SetTitleSize(axisTitleSize);
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetYaxis()->SetTitleOffset(1.5);
  TPaveText* t1a = MakeTextBox(text_x1, leg_y2, textSize, 2);
  t1a            ->AddText("#bf{LArIAT: Run IIB}");
  t1a            ->AddText("Cosmic Michel Electron Dataset");
  t1a           ->Draw();
  TLegend* l1a  = MakeLegend(text_x1, leg_y2 - (2.+0.5)*textSize, textSize, 3);
  l1a         ->AddEntry(gr_Q,"Charge (Q)","PL");
  l1a         ->AddEntry(gr_L,"Light (L)","PL");
  l1a         ->AddEntry(gr_T,"Total (Q+L)#times0.7","PL");
  l1a         ->Draw();
  
  // -------------------------------------------------------------
  // Plot recombination
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  gPad->SetLeftMargin(mar_l);
  gPad->SetRightMargin(mar_r);
  gPad->SetTopMargin(mar_t);
  gPad->SetBottomMargin(mar_b);
  gr_R->GetXaxis()->SetTitle("Electric field [V/cm]");
  gr_R->GetYaxis()->SetTitle("Recombination survival fraction");
  gr_R->GetXaxis()->SetTitleSize(axisTitleSize);
  gr_R->GetYaxis()->SetTitleSize(axisTitleSize);
  gr_R->GetXaxis()->SetTitleOffset(1.1);
  gr_R->GetYaxis()->SetTitleOffset(1.5);
  gr_R->GetYaxis()->SetRangeUser(0.59,0.72);
//  gr_R->Draw("AP");

  // unfix params?
//  f_ModBoxRecomb_Efield->ReleaseParameter(0);
//  f_ModBoxRecomb_Efield->ReleaseParameter(2);
//  f_ModBoxRecomb_Efield->ReleaseParameter(1);
//  gStyle->SetOptFit(1101); // prob chi2 errors names

  /*
  TGraphErrors* gr_Rbox=(TGraphErrors*)gr_R->Clone("Rbox");
  gr_Rbox->Fit(f_ModBoxRecomb_Efield);
  TGraphErrors* gr_Rbirks=(TGraphErrors*)gr_R->Clone("Rbirks");
  gr_Rbirks->Fit(f_BirksRecomb_Efield,"+");
  gr_Rbirks->GetFunction("_BirksRecomb_Efield")->SetLineColor(kBlue);
  gPad->Update();
  TPaveStats *p2 = (TPaveStats*)gr_Rbirks->FindObject("stats");
  p2->SetY1NDC(0.1);
  p2->SetY2NDC(0.3);
  p2->SetTextColor(kBlue);
  gPad->Modified();
  */

  TGraphErrors* gr_R2 = (TGraphErrors*)gr_R->Clone("R2");
  gr_R->Fit(f_ModBoxRecomb_Efield);
  gr_R->Fit(f_BirksRecomb_Efield,"+");
  gr_R->Fit(f_BirksRecombArgoNeuT_Efield,"+");
  gr_R->Draw("AP");
  gr_R2->Draw("P same");
 
  std::cout
  <<"==============================================\n"
  <<"Mod Box (ArgoNeuT) Chi2  = "<<f_ModBoxRecomb_Efield->GetChisquare()/(f_ModBoxRecomb_Efield->GetNDF()-1)<<"\n"
  <<"Birks (ArgoNeuT) Chi2    = "<<f_BirksRecombArgoNeuT_Efield->GetChisquare()/(f_BirksRecombArgoNeuT_Efield->GetNDF()-1)<<"\n"
  <<"Birks (ICARUS) Chi2      = "<<f_BirksRecomb_Efield->GetChisquare()/(f_BirksRecomb_Efield->GetNDF()-1)<<"\n"
  <<"==============================================\n";

  float statH = 0.1;
  float statW = 0.25;
  gStyle->SetStatX(1-mar_r-0.02);
  gStyle->SetStatY(mar_b + 2.*statH);
  gStyle->SetStatW(statW);
  gStyle->SetStatH(statH);
  gStyle->SetOptStat(0);
 
  TPaveText* t1b = MakeTextBox(text_x1, leg_y2, textSize, 2);
  t1b            ->AddText("#bf{LArIAT: Run IIB}");
  t1b            ->AddText("Cosmic Michel Electron Dataset");
  t1b           ->Draw();
  
  textSize = 0.03; 
  TLegend* l1b  = MakeLegend(text_x1, leg_y2 - (2.+0.5)*textSize, textSize, 5);
  l1b         ->AddEntry(gr_R,"Data","PL");
  l1b         ->AddEntry(f_BirksRecomb_Efield,"Birks (ICARUS params)","L");
  l1b         ->AddEntry(f_BirksRecombArgoNeuT_Efield,"Birks (ArgoNeuT params)","L");
  l1b         ->AddEntry(f_ModBoxRecomb_Efield,"Mod. Box (ArgoNeuT params)","L");
  l1b         ->AddEntry((TObject*)0,Form("dE/dx = %3.1f MeV/cm",f_ModBoxRecomb_Efield->GetParameter(2)),"");
  l1b         ->Draw();


  fMaxRun = maxRunOrig;

}

// ###############################################################################3
// Plot the SPEs and average waveform effective lifetime
void PlotSpeAndTau(){
  
  /*  
  // From calibration file:
  // HMM PMT
  Run1  Run2	  SER     SER_err   width   eff. tau
  9438  9474        38.16   0.08      12.50  1188.     
  9476  9493        37.74   0.08      12.40  1177.
  9494  9514        37.71   0.08      12.10  1182.
  9515  9553        37.16   0.08      12.20  1177.
  9554  9631        37.07   0.08      12.00  1175.

  // ETL PMT
  Run1  Run2	  SER     SER_err   width   eff. tau
	6245  6262        62.90   0.50      32.50   1260.
	6263  6284        63.60   0.50      31.60   1278.
	6285  6326        61.60   0.70      30.80   1261.
	6327  6345        59.20   0.80      30.40   1270.
  9438  9474        42.10   0.20      17.80  1263.  
  9476  9493        41.90   0.20      17.40  1273.
  9494  9514        41.00   0.20      17.00  1251.
  9515  9553        40.60   0.20      17.00  1231.
  9554  9631        40.30   0.20      16.80  1259.
  */


  // =============================================================
  // Run 1
  const int N1=4;
  double days1[N1]          ={  1.5,  4.5,  7.5,  10.5  };
  double days1_err[N1]      ={  1.5,  1.5,  1.5,  1.5  };
  double spe1etl[N1]        ={  62.9,	63.6,	61.6,	59.2 	};
  double spe1etl_err[N1]    ={  0.50,	0.5,	0.70,	0.80	};
  double tau1etl[N1]        ={  1260.,1278.,1261.,1270. };
  double tau1etl_err[N1]    ={  4.,   3.,   4.,   4.    }; 

  // =========================================
  // Run 2 (abridged dataset)
  const int N2=5;
  double days2[N2]          ={  1.5,  4.5,  7.5,  10.5, 14.0  };
  double days2_err[N2]      ={  1.5,  1.5,  1.5,  1.5,  2.    };
  double spe2etl[N2]        ={  42.10,41.90,41.00,40.60,40.30 };
  double spe2etl_err[N2]    ={  0.2,  0.2,  0.2,  0.2,  0.2   };
  double spe2hmm[N2]        ={  38.16,37.74,37.71,37.16,37.07  };
  double spe2hmm_err[N2]    ={  0.08, 0.08, 0.08, 0.08, 0.08  };
  double tau2etl[N2]        ={  1263.,1273.,1251.,1231.,1259. };
  double tau2etl_err[N2]    ={  4.,   4.,   4.,   4.,   4.    };
  double tau2hmm[N2]        ={  1188.,1177.,1182.,1177.,1175. };
  double tau2hmm_err[N2]    ={  2.,   3.,   3.,   2.,   3.    };
  

  float textSize = 0.035;
  float mar_l  = 0.15;
  float mar_r  = 0.05;
  float mar_t  = 0.05;
  float leg_x1 = 0.50;
  float leg_x2 = 1.-mar_r-0.02;
  float leg_y1 = 0.7;
  float leg_y2 = 1.-mar_t-0.02;
  float text_x1 =  mar_l+0.02;
  float text_x2 =  0.4;
  float text_y1 = 0.75;
  float text_y2 = 1.-mar_t-0.03;
  float axisTitleSize = 0.045;

  // ..........................................................................
  // Plot Run 1
  TCanvas* c1 = new TCanvas("c1","",1000,500);
  c1->Divide(2,1);
  TGraphErrors* g_spe1etl = new TGraphErrors();
  TGraphErrors* g_tau1etl = new TGraphErrors();
  for(size_t i=0; i<N1; i++ ) {
    g_spe1etl->SetPoint(        g_spe1etl->GetN(),    days1[i],     spe1etl[i] );
    g_spe1etl->SetPointError(   g_spe1etl->GetN()-1,  0.99*days1_err[i], spe1etl_err[i] );
    g_tau1etl->SetPoint(        g_tau1etl->GetN(),    days1[i],     tau1etl[i] );
    g_tau1etl->SetPointError(   g_tau1etl->GetN()-1,  0.99*days1_err[i], tau1etl_err[i] );
  }
  
  c1->cd(1);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);
  g_spe1etl->GetXaxis()->SetTitle("Days elapsed");
  g_spe1etl->GetYaxis()->SetTitle("Single photoelectron response [ADC]");
  g_spe1etl->SetMarkerStyle(20);
  g_spe1etl->SetMarkerColor(kBlue);  
  g_spe1etl->SetLineColor(kBlue); 
  g_spe1etl->Draw("AP");
  g_spe1etl->GetXaxis()->SetRangeUser(0.,12.);
  g_spe1etl->GetYaxis()->SetRangeUser(60,90);
  TPaveText* t1a = MakeTextBox(text_x1, leg_y2, textSize, 2);
  t1a            ->AddText("#bf{LArIAT: Run I}");
  t1a            ->AddText("Cosmic Michel Electron Dataset");
  t1a           ->Draw();
  TLegend* l1a  = MakeLegend(text_x1, leg_y2 - 2.*textSize, textSize, 1);
  l1a         ->AddEntry(g_spe1etl,"ETL (2-inch) PMT","PL");
  l1a         ->Draw();

  c1->cd(2);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);
  g_tau1etl->GetXaxis()->SetTitle("Days elapsed");
  g_tau1etl->GetYaxis()->SetTitle("Measured muon late-light lifetime #tau_{meas} [ns]");
  g_tau1etl->SetMarkerStyle(20);
  g_tau1etl->SetMarkerColor(kBlue+2);  
  g_tau1etl->SetLineColor(kBlue+2);  
  g_tau1etl->Draw("AP");
  g_tau1etl->GetXaxis()->SetRangeUser(0.,12.);
  g_tau1etl->GetYaxis()->SetRangeUser(1200.,1360.);
  TPaveText* t1b = MakeTextBox(text_x1, leg_y2, textSize, 2);
  t1b            ->AddText("#bf{LArIAT: Run I}");
  t1b            ->AddText("Cosmic Michel Electron Dataset");
  t1b           ->Draw();
  TLegend* l1b  = MakeLegend(text_x1, leg_y2 - 2.*textSize, textSize, 1);
  l1b         ->AddEntry(g_tau1etl,"ETL (2-inch) PMT","PL");
  l1b         ->Draw();
  
  
  
  // ..........................................................................
  // Plot Run 2b
  TCanvas* c2 = new TCanvas("c2","",1000,500);
  c2->Divide(2,1);
  TGraphErrors* g_spe2etl = new TGraphErrors();
  TGraphErrors* g_tau2etl = new TGraphErrors();
  TGraphErrors* g_spe2hmm = new TGraphErrors();
  TGraphErrors* g_tau2hmm = new TGraphErrors();
  for(size_t i=0; i<N2; i++ ) {
    g_spe2etl->SetPoint(        g_spe2etl->GetN(),    days2[i],     spe2etl[i] );
    g_spe2etl->SetPointError(   g_spe2etl->GetN()-1,  0.99*days2_err[i], spe2etl_err[i] );
    g_tau2etl->SetPoint(        g_tau2etl->GetN(),    days2[i],     tau2etl[i] );
    g_tau2etl->SetPointError(   g_tau2etl->GetN()-1,  0.99*days2_err[i],  tau2etl_err[i] );
    g_spe2hmm->SetPoint(        g_spe2hmm->GetN(),    days2[i],     spe2hmm[i] );
    g_spe2hmm->SetPointError(   g_spe2hmm->GetN()-1,  0.99*days2_err[i], spe2hmm_err[i] );
    g_tau2hmm->SetPoint(        g_tau2hmm->GetN(),    days2[i],     tau2hmm[i] );
    g_tau2hmm->SetPointError(   g_tau2hmm->GetN()-1,  0.99*days2_err[i],  tau2hmm_err[i] );
  }
  
  c2->cd(1);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);
  g_spe2etl->GetXaxis()->SetTitle("Days elapsed");
  g_spe2etl->GetYaxis()->SetTitle("Single photoelectron response [ADC]");
  g_spe2etl->SetMarkerStyle(20);
  g_spe2etl->SetMarkerColor(kBlue);  
  g_spe2etl->SetLineColor(kBlue); 
  g_spe2etl->Draw("AP");
  g_spe2etl->GetXaxis()->SetRangeUser(0.,16.);
  g_spe2etl->GetYaxis()->SetRangeUser(32,46);
  g_spe2hmm->SetMarkerStyle(20);
  g_spe2hmm->SetMarkerColor(kRed);  
  g_spe2hmm->SetLineColor(kRed); 
  g_spe2hmm->Draw("P");
  TPaveText* t2a = MakeTextBox(text_x1, leg_y2, textSize, 2);
  t2a            ->AddText("#bf{LArIAT: Run IIB}");
  t2a            ->AddText("Cosmic Michel Electron Dataset");
  t2a           ->Draw();
  TLegend* l2a  = MakeLegend(text_x1, leg_y2 - 2.*textSize, textSize, 2);
  l2a         ->AddEntry(g_spe2hmm,"HMM (e-inch) PMT","PL");
  l2a         ->AddEntry(g_spe2etl,"ETL (2-inch) PMT","PL");
  l2a         ->Draw();

  c2->cd(2);
  gPad->SetLeftMargin(mar_l); gPad->SetRightMargin(mar_r); gPad->SetTopMargin(mar_t);
  g_tau2etl->GetXaxis()->SetTitle("Days elapsed");
  g_tau2etl->GetYaxis()->SetTitle("Measured muon late-light lifetime #tau_{meas} [ns]");
  g_tau2etl->SetMarkerStyle(20);
  g_tau2etl->SetMarkerColor(kBlue+2);  
  g_tau2etl->SetLineColor(kBlue+2);  
  g_tau2etl->Draw("AP");
  g_tau2etl->GetXaxis()->SetRangeUser(0.,16.);
  g_tau2etl->GetYaxis()->SetRangeUser(1100.,1400.);
  g_tau2hmm->SetMarkerStyle(20);
  g_tau2hmm->SetMarkerColor(kRed+2);  
  g_tau2hmm->SetLineColor(kRed+2); 
  g_tau2hmm->Draw("P");
  TPaveText* t2b = MakeTextBox(text_x1, leg_y2, textSize, 2);
  t2b            ->AddText("#bf{LArIAT: Run IIB}");
  t2b            ->AddText("Cosmic Michel Electron Dataset");
  t2b           ->Draw();
  TLegend* l2b  = MakeLegend(text_x1, leg_y2 - 2.*textSize, textSize, 2);
  l2b         ->AddEntry(g_tau2hmm,"HMM (3-inch) PMT","PL");
  l2b         ->AddEntry(g_tau2etl,"ETL (2-inch) PMT","PL");
  l2b         ->Draw();

   
  

}

//###################################################################################
void PeakExclusionFit(TH1D* h, TF1* f, float x1, float x2 ){
  std::cout
  <<"---------------------------------------------------\n"
  <<"PeakExclusionFit(): initial params "<<f->GetParameter(0)<<", "<<f->GetParameter(1)<<", "<<f->GetParameter(2)<<"\n"
  <<"Excluding region between "<<x1<<" and "<<x2<<"\n";
  float sum = 0; 
  TGraph gr;
  for(size_t i=1; i<=h->GetXaxis()->GetNbins(); i++){
    float x = h->GetXaxis()->GetBinCenter(i);
    if( (x < x1 || x > x2 ) && h->GetBinContent(i) > 0 ) {
      gr.SetPoint(gr.GetN(), x, h->GetBinContent(i));
      sum += h->GetBinContent(i);
    }
  }
  gr.Fit(f);
  std::cout
  <<"Fit included "<<gr.GetN()<<" points ("<<sum<<" events)\n"
  <<"final params "<<f->GetParameter(0)<<", "<<f->GetParameter(1)<<"\n"
  <<"---------------------------------------------------\n";
}



//####################################################################################
void MultiGausFit(TH1D* h, TF1* f2g, TF1* fbg, float param){
  MultiGausFit(h,f2g,fbg,1,param);
}

void MultiGausFit(TH1D* h, TF1* f2g, TF1* fbg, int mode, float param){
  
  // ---------------------------------------------------
  // Mode:
  //            0 == Simple fit to single Gaussian using max bin +/- thresh-based range
  //  (default) 1 == Initial fit to 2-Gaus, then fit signal peak on top of fixed BG
  //            2 == Initial fit to 2-Gaus, then zero-out BG and fit signal peak
  //            3 == Use initial 2-Gaus fit w/ no re-fitting 
  
  
  fbg ->ReleaseParameter(0);
  fbg ->ReleaseParameter(1);
  fbg ->ReleaseParameter(2);
  f2g ->ReleaseParameter(3);
  f2g ->ReleaseParameter(4);
  f2g ->ReleaseParameter(5);
  //fbg ->ReleaseParameter(0);
  //fbg ->ReleaseParameter(1);
  //fbg ->ReleaseParameter(2);

  // Initialize parameters  
  float max_bin = h->GetMaximumBin();
  float peak_res= h->GetBinCenter(max_bin);
  float max     = GetHistMax(h);
  f2g->SetParLimits(0,  0,1.2*max);
  f2g->SetParameter(0,  0.95*max );
  f2g->SetParameter(1,  h->GetBinCenter(max_bin));
  f2g->SetParLimits(2,  0.01, 1.5*h->GetRMS() );
  f2g->SetParameter(2,  std::max( float(h->GetRMS()/2.), float(0.01)) );
  f2g->SetParLimits(3,  0,0.5*max );
  f2g->SetParameter(3,  0.05*max );
  f2g->SetParameter(4,  h->GetMean() );
  f2g->SetParLimits(5,  h->GetRMS(), h->GetRMS()*5. );
  f2g->SetParameter(5,  h->GetRMS()*2. );

  float r1 = h->GetXaxis()->GetXmin();
  float r2 = h->GetXaxis()->GetXmax();

 
  // 0)
  // Simple case of single gaus fit.
  if( mode == 0 ) {

    f2g->FixParameter(3, 0);
    f2g->FixParameter(4, 0);
    f2g->FixParameter(5, 0);
    //f2g->SetParameter(2,h->GetRMS()/10.);
    fbg ->SetParameter(0, f2g->GetParameter(3));
    fbg ->SetParameter(1, f2g->GetParameter(4));
    fbg ->SetParameter(2, f2g->GetParameter(5));
    
    // Find where histogram drops below given threshold
    // of the peak height (for use in defining fit range)
    size_t minBinsForFit = 8;
    if( param >= 0 ) {
      minBinsForFit /= 2;
      for(size_t ii=max_bin+minBinsForFit; ii<h->GetXaxis()->GetNbins(); ii++){
        if( h->GetBinContent(ii) < max*param && h->GetBinContent(ii+1) < max*param ) {
          //|| (h->GetBinLowEdge(ii) >= peak_res + limit ) ) {
          r2 = h->GetBinLowEdge(ii+1);
          break;
        }
      }
      for(size_t ii=max_bin-minBinsForFit; ii>0; ii--){
        if( h->GetBinContent(ii) < max*param && h->GetBinContent(ii-1) < max*param ) {
          //|| (h->GetBinLowEdge(ii) <= peak_res - limit ) ) {
          r1 = h->GetBinLowEdge(ii);
          break;
        }
      }
    }
    
    /*
    float limit = 0.30;
    if( param >= 0 ) {
      for(size_t ii=max_bin; ii<h->GetXaxis()->GetNbins(); ii++){
        if( (h->GetBinContent(ii) < max*param && h->GetBinContent(ii+1) < max*param )
          || (h->GetBinLowEdge(ii) >= peak_res + limit ) ) {
          r2 = h->GetBinLowEdge(ii+1);
          break;
        }
      }
      for(size_t ii=max_bin; ii>0; ii--){
        if( (h->GetBinContent(ii) < max*param && h->GetBinContent(ii-1) < max*param )
          || (h->GetBinLowEdge(ii) <= peak_res - limit ) ) {
          r1 = h->GetBinLowEdge(ii);
          break;
        }
      }
    }
    */
    f2g->SetRange(r1,r2);
    h->Fit( f2g, "QR");
    // In rare instances the error on the Gaus width is
    // humungous for now discernable reason. If this happens, 
    // refit with different weighting scheme.
    if( f2g->GetParError(2)/f2g->GetParameter(2) > 0.5 ) 
      h->Fit( f2g, "QRW");


  } else {
    // Initial 2-Gaus fit to find where the peak is
    std::cout
    <<"Using initial parameters to 2-Gauss: \n"
    <<"  p0: "<<f2g->GetParameter(0)<<"\n"
    <<"  p1: "<<f2g->GetParameter(1)<<"\n"
    <<"  p2: "<<f2g->GetParameter(2)<<"\n"
    <<"  p3: "<<f2g->GetParameter(3)<<"\n"
    <<"  p4: "<<f2g->GetParameter(4)<<"\n"
    <<"  p5: "<<f2g->GetParameter(5)<<"\n";
    //f2g ->SetRange(-1.2,1.2);
    h   ->Fit( f2g, "Q" );
  
    if( mode == 1 ) {

      // Fit peak on top of BG 
      //fbg ->SetParLimits(0,0,0.5*max);
      //fbg ->SetParLimits(2,h->GetRMS(),h->GetRMS()*5.);
      //if( excludePeak ) {
      //  PeakExclusionFit(h,fbg,
      //    f2g->GetParameter(1)-3*f2g->GetParameter(2),
      //    f2g->GetParameter(1)+3*f2g->GetParameter(2));
      //}
      
      // Fix the Gaussian parameters and refit
      //if2g->FixParameter(3, fbg->GetParameter(0));
      //f2g->FixParameter(4, fbg->GetParameter(1));
      //f2g->FixParameter(5, fbg->GetParameter(2));


      if( param > 0 ) {
        r1 = f2g->GetParameter(1)-param*f2g->GetParameter(2);
        r2 = f2g->GetParameter(1)+param*f2g->GetParameter(2);
      }
      f2g->SetRange(r1,r2);
      h->Fit( f2g, "QR"); 
      if( f2g->GetParError(5)/f2g->GetParameter(5) > 0.5 ) 
        h->Fit( f2g, "QRW");
      
      fbg ->SetParameter(0, f2g->GetParameter(3));
      fbg ->SetParameter(1, f2g->GetParameter(4));
      fbg ->SetParameter(2, f2g->GetParameter(5));
      
    } else
    if( mode == 2 ) {
      // Zero-out BG and fit to peak
      fbg ->SetParameter(0, 0.);
      fbg ->SetParameter(1, 0.);
      fbg ->SetParameter(2, 0.);
      f2g->FixParameter(3, 0.);
      f2g->FixParameter(4, 0.);
      f2g->FixParameter(5, 0.);
      f2g->SetRange(
        f2g->GetParameter(1)-param*f2g->GetParameter(2),
        f2g->GetParameter(1)+param*f2g->GetParameter(2));
      h->Fit( f2g, "QRW"); 
      if( f2g->GetParError(5)/f2g->GetParameter(5) > 0.3 ) 
        h->Fit( f2g, "QRW");
    }
  
  }

}





// ####################################################################################
// Helper function for the resolution slice methods below. This function takes in a 
// 2D histogram and slices it up into some number of 1D projections.
void MakeSlices( const TH2D* h2, 
                std::vector<TH1D*>& vec_h1, 
                std::vector<float>& vec_x,
                std::vector<int>& vec_nbins,
                std::string name, 
                float xmin, float xmax, int nslices, float dx){
  
  // Find first bin within the range and count up *all*
  // bins in this range. If the provided "nslices" is negative,
  // then make a slice out of every bin.
  float         bwidth    = h2->GetXaxis()->GetBinWidth(1);
  const size_t  nbins     = h2->GetXaxis()->GetNbins();
  
  // By default, dx is a single bin width.
  if( dx <= 0 ) dx = bwidth;
  
  std::cout
  <<"\n"
  <<"################################################################################\n"
  <<"Slicing up histogram: "<<h2->GetTitle()<<" ("<<nbins<<" bins total )\n"
  //<<"Looping "<<N<<" bins between "<<xmin<<" and "<<xmax<<"\n"
  <<"Creating "<<nslices<<" slices of width "<<dx<<"\n";
  
  // Define the centerpoint of the first slice region
  float intervalwidth = (xmax-xmin) / float (nslices);
  float xcen          = xmin + intervalwidth/2.;
  float x1            = xcen - dx/2.;
  float x2            = xcen + dx/2.; 
  int   i_slice       = 0;
  int   startbin      = -1;
  bool  inBinGroup    = false;
  float xsum          = 0.;

  for(size_t i=1; i<nbins; i++){

    float binx      = h2->GetXaxis()->GetBinCenter(i);
    float lowEdge   = h2->GetXaxis()->GetBinLowEdge(i);
    float highEdge  = lowEdge + bwidth;
    bool  flag       = (binx >= x1 && binx < x2);
    TH1D* h_tmp;
    
    if( flag )  std::cout<<"   "<<i<<"   binx = "<<binx<<"    xcen = "<<xcen<<"    x1-x2 "<<x1<<"-"<<x2<<"\n";
      
    // if we're already in a bin group and 
    // we reach a bin that doesn't satisfy
    // the range, then pack things up.
    if( inBinGroup && !flag ) {
      
      //h_tmp     = h2->ProjectionY( name.c_str(), startbin, i-1 ); 
      h_tmp     = Make1DSlice( h2, startbin, i-1, name.c_str() );
      vec_h1    .push_back( h_tmp );
      vec_x     .push_back( xsum / (i-startbin) );
      vec_nbins .push_back( i-startbin );
      
      i_slice++;
      xsum = 0.;
      startbin = -1;
      x1 += intervalwidth;
      x2 += intervalwidth;
      flag = (binx >= x1 && binx < x2);
      inBinGroup = false;
    }
    
    if( i_slice >= nslices ) break;

    if( flag ) {
      xsum += binx;
      if( startbin < 0 ) startbin = i;
      inBinGroup = true;  
    } 
    
  } // done looping over all bins

}


//######################################################################################
// RESLOOP
void ResolutionSliceLoop(
  const TH2D* h2d, float Emin, float Emax,int N, float dE, 
  bool makeFitArray, std::string canvasName, std::string tag, int nx, int ny, float xmin, float xmax,
  int mode, float param, bool useHybridRes,
  float Emin_gr, float Emax_gr,
  TGraphAsymmErrors* gr_sig,
  TGraphAsymmErrors* gr_rms){
  std::vector<TF1> fin;
  ResolutionSliceLoop(h2d,Emin,Emax,N, dE,
                      makeFitArray,canvasName,tag,nx,ny,xmin,xmax,
                      mode,param,useHybridRes,
                      Emin_gr, Emax_gr,
                      gr_sig,
                      gr_rms,
                      fin);
}

void ResolutionSliceLoop(
  const TH2D* h2d, float Emin, float Emax, int nslice, float dE, 
  bool makeFitArray, std::string canvasName, std::string tag, int nx, int ny, float xmin, float xmax,
  int mode, float param, bool useHybridRes,
  float Emin_gr, float Emax_gr,
  TGraphAsymmErrors* gr,
  TGraphAsymmErrors* gr_rms,
  std::vector<TF1> &fvFunc )
{
  // Axis title and label size
  float axls = 0.05;
  float axts = 0.05;

  // Make vector of 1D histograms corresponding to each slice,
  // as well as the slice central energy and num bins in slice.
  std::vector<TH1D*> vec_h1d;
  std::vector<float> vec_E;
  std::vector<int> vec_nbins;
  
  // Make a default canvas
  TCanvas* cdefault = new TCanvas("default","default",200,200);
  cdefault->cd();
  
  //gStyle->SetOptStat(1110);
  //gStyle->SetOptFit(1111);

  // Slice up the 2D histo
  MakeSlices( h2d, vec_h1d, vec_E, vec_nbins, canvasName, Emin, Emax, nslice, dE);

  // Make the multi-gauss functions which we'll be fitting to
  // (and a separate single-gaussian BG function)
  TF1 f2g("f2g","[0]*exp( - pow( (x-[1])/( sqrt(2)*[2]), 2)) + [3]*exp( - pow( (x-[4])/( sqrt(2)*[5]), 2)) ");
  TF1 fbg("fbg","[0]*exp( -pow((x-[1])/(sqrt(2)*[2]),2))");
  fbg .SetLineColor(kOrange-1);
  fbg .SetLineStyle(2);
   
  // Define the main canvas array
  TCanvas* cArray;
  if( makeFitArray ){
    int npx_x = 250*nx;
    int npx_y = 250*ny;
    //gStyle->SetOptStat(1110);
    //gStyle->SetOptFit(1111);
    cArray = new TCanvas(Form("arr_%s",canvasName.c_str()),Form("arr_%s",canvasName.c_str()),npx_x,npx_y);
    cArray->Divide(nx,ny);
    cArray->cd(1);
  }
  
  // Loop through the slices and do the fits
  for(size_t i=0; i<vec_h1d.size(); i++){

    // by default, cd into some throw-away canvas to prevent weird
    // formatting effects on the main canvas array (ROOT sucks)
    cdefault->cd();
    
    TH1D* h_tmp     = vec_h1d.at(i);
    float histMax   = GetHistMax(h_tmp);
//    if( histMax < 80 ) h_tmp->Rebin(2); 
   
    // get info about this 1D projection
    float binwidth  = vec_h1d.at(i)->GetXaxis()->GetBinWidth(1);
    float rms       = vec_h1d.at(i)->GetRMS();
    float nentries  = vec_h1d.at(i)->GetEntries();
    float E         = vec_E.at(i);
    float nbins     = vec_nbins.at(i);
    float Ebinwidth = h2d->GetXaxis()->GetBinWidth(1);
    float Eerr      = 0.5*nbins*Ebinwidth/std::sqrt(12);

    std::cout
    <<"\n"
    <<"-----------------------------------------------------------------------------------\n"
    <<"Fitting energy slice from "<<nbins<<" bins, <E> = "<<E<<" MeV ("<<nentries<<" entries, RMS = "<<rms<<")\n";
 
    if( nentries > 10 ) {
     
      // send projection to all-purpose gaussian fitting routine
      MultiGausFit( h_tmp, &f2g, &fbg, mode, param); 

      // chi2
      float chi2      = f2g.GetChisquare()/(f2g.GetNDF()-1);
    
      // Gaussian fit sigma
      float sig       = f2g.GetParameter(2)*100.;
      float sig_err   = f2g.GetParError(2)*100.;
      float sig_err_l = sig_err;
      float sig_err_u = sig_err;
      float sigRMS    = rms*100.;

      // incorporate the difference between RMS and sigma 
      // into the sigma error ("hybrid" approach)
      if( useHybridRes ) {
      
        float sig0    = sig;

        // set resolution to be the midpoint between the 
        // fitted Gaussian at the peak and the overall RMS.
        // The new error on this resolution will be a quadrature
        // sum of the original sigma error and half the difference
        // between sigma and RMS
        /*
        sig     = 0.5*(sig + sigRMS);
        sig_err = std::sqrt( std::pow( sig_err,2 ) + std::pow( 0.5*(sigRMS-sig0), 2 ) );
        sig_err_l = sig_err;
        sig_err_u = sig_err;    
        */
        
        // weighting based on chi^2 of fit
        float w_sig = 1.;
        float w_rms = 1.;
        if( fHybridRes_UseChisquareWeighting ) 
          w_rms = std::sqrt(chi2);
        
        float W = w_sig + w_rms;
        float min = std::min(sig,sigRMS);
        float max = std::max(sig,sigRMS);
        
        sig      = (w_sig*sig + w_rms*sigRMS)/(w_sig+w_rms);
        sig_err_l = std::sqrt( std::pow((w_sig/W)*sig_err,2) + std::pow(sig-min,2) );
        sig_err_u = std::sqrt( std::pow((w_sig/W)*sig_err,2) + std::pow(sig-max,2) );

        // .................................
        // Temporary 2/26/19
        sig = sigRMS;
        sig_err_l = 0;
        sig_err_u = 0;
      
      } 
    

      // add points to all the TGraphs
      if( ( Emin_gr < 0 || E >= Emin_gr ) && ( Emax_gr < 0 || E <= Emax_gr )) {
        std::cout<<"  Setting point on sigma graph:   "<<sig<<" (+ "<<sig_err_u<<" - "<<sig_err_l<<")\n";
        std::cout<<"  Setting point on RMS graph:     "<<sigRMS<<"\n";
        std::cout<<"  slice energy            "<<E<<" +/- "<<Eerr<<"\n";
        gr      ->SetPoint      (gr->GetN(),    E, sig);
        gr      ->SetPointError (gr->GetN()-1,  Eerr, Eerr, sig_err_l, sig_err_u);
        gr_rms  ->SetPoint      (gr_rms->GetN(),E, sigRMS);
        gr_rms  ->SetPointError (gr_rms->GetN()-1, Eerr, Eerr, 0., 0.);
     
        std::cout<<"  saving function\n"; 
        // save the fit function   
        fvFunc.push_back(f2g);
      }
 
  
      // plot it in the array
      if( makeFitArray && i < nx*ny ) {
        std::cout<<"  Preparing to plot\n";
        float x1 = h_tmp->GetXaxis()->GetXmin();
        float x2 = h_tmp->GetXaxis()->GetXmax();
        if( xmin != xmax && xmax > xmin ) {
          x1 = xmin;
          x2 = xmax;
        }
        cArray->cd(i+1);
        gPad->SetMargin(0.10,0.10,0.15,0.10);
        FormatAxes(h_tmp, axts, axls, 1.0, 1.0);
        h_tmp->GetXaxis()->SetNdivisions(505);
        h_tmp->SetTitle(Form("%s : %4.1f +/- %4.1f MeV",tag.c_str(),E,dE/2.));
        h_tmp->SetTitleSize(0.06);
        h_tmp->GetXaxis()->SetTitle(h2d->GetYaxis()->GetTitle());
        h_tmp->GetXaxis()->SetRangeUser(x1,x2);
        h_tmp->DrawCopy();
        if( fbg.GetParameter(0) ) {
          std::cout<<"  Background gaus has N = "<<fbg.GetParameter(0)<<" --> plotting!\n";
          fbg.SetRange( x1, x2);
          fbg.DrawCopy("LSAME");
        }
        gStyle->SetOptStat(1110);
        gPad->Update();
        gPad->SetEditable(false);
      
        cdefault->cd();

      }// done plotting in array

    }// end check for nentries
  
  }// done looping through all slices 
  
  cdefault->cd();

  if( fSavePlots && makeFitArray ) {
    outFile->cd();
    cArray->Write();
  }

}

void ResolutionSliceLoop(
  const TH2D* h2d, float Emin, float Emax, int nslice, float dE, 
  bool makeFitArray, std::string canvasName, std::string tag, int nx, int ny, float xmin, float xmax,
  int mode, float param, bool useHybridRes,
  float Emin_gr, float Emax_gr,
  TGraphAsymmErrors* gr,
  TGraphAsymmErrors* gr_rms,
  TGraphAsymmErrors* gr_mean ){
  
  std::vector<TF1> funcs;
  ResolutionSliceLoop( h2d, Emin, Emax, nslice, dE,
  makeFitArray, canvasName, tag, nx, ny, xmin, xmax,
  mode, param, useHybridRes, Emin_gr, Emax_gr,
  gr, gr_rms, funcs);

  for(size_t i=0; i<funcs.size(); i++){
    double E, dE, dummy;
    gr->GetPoint(i,E,dummy);
    if( ( Emin_gr < 0 || E >= Emin_gr ) && ( Emax_gr < 0 || E <= Emax_gr )) {
      dE = gr->GetErrorX(i);
      gr_mean->SetPoint( gr_mean->GetN(), E, funcs[i].GetParameter(1));
      gr_mean->SetPointError( gr_mean->GetN()-1, dE, dE, funcs[i].GetParError(1), funcs[i].GetParError(1) ); 
    } 
  
  }
}



void ParamDistOverlay(
  const TH2D* h2d, float Emin, float Emax, int nslice, float dE, float min, float max, 
  std::string canvasName, std::string tag, int nx, int ny,
  TF1* func_in)
{
  
  // Axis title and label size
  float axls = 0.05;
  float axts = 0.05;

  // Make vector of 1D histograms corresponding to each slice,
  // as well as the slice central energy and num bins in slice.
  std::vector<TH1D*> vec_h1d;
  std::vector<float> vec_E;
  std::vector<int> vec_nbins;
  
  // Make a default canvas
  TCanvas* cdefault = new TCanvas("default","default",200,200);

  // Slice up the 2D histo
  MakeSlices( h2d, vec_h1d, vec_E, vec_nbins, canvasName, Emin, Emax, nslice, dE);
  
  // Define the main canvas array
  int npx_x = 250*nx;
  int npx_y = 250*ny;
  TCanvas* cArray = new TCanvas(Form("arr_%s",canvasName.c_str()),Form("arr_%s",canvasName.c_str()),npx_x,npx_y);
  cArray->Divide(nx,ny);
  gStyle->SetOptStat(0);
  cArray->cd(1);
  //gStyle->SetOptStat(1110);
  
  // Loop through the slices and do the fits
  for(size_t i=0; i<vec_h1d.size(); i++){

    // by default, cd into some throw-away canvas to prevent weird
    // formatting effects on the main canvas array (ROOT sucks)
    cdefault->cd();
   
    // get info about this 1D projection
    TH1D* h_tmp     = vec_h1d.at(i);
    float binwidth  = vec_h1d.at(i)->GetXaxis()->GetBinWidth(1);
    float rms       = vec_h1d.at(i)->GetRMS();
    float nentries  = vec_h1d.at(i)->GetEntries();
    float E         = vec_E.at(i);
    float nbins     = vec_nbins.at(i);
    float Ebinwidth = h2d->GetXaxis()->GetBinWidth(1);
    float Eerr      = 0.5*nbins*Ebinwidth/std::sqrt(12);
   
    // format and scale scale the function
    func_in -> SetLineColor(kGreen+2);
    func_in ->SetParameter("scale",h_tmp->Integral()*binwidth);
    func_in ->SetParameter("energy",E);
    
    // plot it
    if( i < nx*ny ) {
      cArray->cd(i+1);
      //gStyle->SetTitleFont(42,"t");
      //gStyle->SetOptStat(0);
      gPad->SetMargin(0.10,0.10,0.15,0.10);
      h_tmp->SetTitle(Form("%s: %4.1f +/- %4.1f MeV",tag.c_str(),E,dE/2.));
      h_tmp->SetTitleSize(0.06);
      h_tmp->GetXaxis()->SetTitle(h2d->GetYaxis()->GetTitle());
      h_tmp->GetXaxis()->SetRangeUser(min,max);
      h_tmp->SetLineWidth(2);
      h_tmp->GetXaxis()->SetNdivisions(505);
      FormatAxes(h_tmp, axts, axls, 1.0, 1.0);
      h_tmp->DrawCopy();
      func_in->DrawCopy("LSAME"); 
      gStyle->SetOptStat(0);
      gPad->SetEditable(false);
      cdefault->cd();
    }
      
  
  
 }// done looping through all slices 
  
  if( fSavePlots ) {
  outFile->cd();
  cArray->Write();
  }

}







void PlotRecombCurves () {
  
  TF1* fmodbox = new TF1("modbox","log([0]+([1]*x)/(1.370*[2]))/( ([1]*x)/(1.370*[2]) )",0.,100); 
  fmodbox  ->SetParameter(0,fBoxModelA);
  fmodbox  ->SetParameter(1,fBoxModelB);
  fmodbox  ->SetParameter(2,fEField);
  fmodbox   ->SetNpx(2000);
  fmodbox   ->SetLineColor(kBlack);
  TF1* fbirks = new TF1("birks","[0]/(1+(x*[1])/(1.370*[2]) )",0.,100);
  fbirks  ->SetParameter(0,fBirksModelA);
  fbirks  ->SetParameter(1,fBirksModelk);
  fbirks  ->SetParameter(2,fEField);
  fbirks   ->SetNpx(100);
  fbirks    ->SetLineStyle(2);
  fbirks   ->SetLineColor(kBlack);

  TFile* file = new TFile("files/MichelAna_mc2_nonoise.root","read");
  TH1D* h_dEdx = (TH1D*)file->Get("largeant/energyPerLength");
  
  TCanvas* c_recomb = new TCanvas("recomb","recomb",600,600);
  gPad->SetMargin(0.15,0.05,0.15,0.05);
  gStyle->SetOptStat(0);

  TF1* birks_lariat = (TF1*)fbirks->Clone("birks_lariat");
  TF1* modbox_lariat = (TF1*)fmodbox->Clone("modbox_lariat");
  TF1* birks_ub = (TF1*)fbirks->Clone("birks_ub");
  TF1* modbox_ub = (TF1*)fmodbox->Clone("modbox_ub");
  birks_ub->FixParameter(2,0.277);
  modbox_ub->FixParameter(2,0.277);
  
  modbox_lariat->SetLineColor(kBlue+2);
  birks_lariat->SetLineColor(kBlue+2);
  birks_lariat->SetLineStyle(2);
  modbox_ub->SetLineColor(kOrange+2);
  birks_ub->SetLineColor(kOrange+2);
  birks_ub->SetLineStyle(2);

  h_dEdx->Scale( 0.95 / h_dEdx->GetMaximum() );
  h_dEdx->SetFillColor(kGreen-8);
  h_dEdx->SetFillStyle(1001);
  h_dEdx->SetLineColor(kBlack);

  h_dEdx->SetTitle("");
  h_dEdx->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
  h_dEdx->GetYaxis()->SetTitle("Recombination survival fraction");
  h_dEdx->GetXaxis()->SetTitleSize(0.05);
  h_dEdx->GetYaxis()->SetTitleSize(0.05);
  h_dEdx->GetXaxis()->SetLabelSize(0.04);
  h_dEdx->GetYaxis()->SetLabelSize(0.04);
  h_dEdx->GetXaxis()->SetTitleOffset(1.2);
  h_dEdx->GetYaxis()->SetTitleOffset(1.2);
  h_dEdx->GetXaxis()->SetRangeUser(0,6);
  h_dEdx->GetYaxis()->SetRangeUser(0,1);
  h_dEdx->Draw("hist");

  modbox_lariat->Draw("same");
  birks_lariat->Draw("same");
  //modbox_ub->Draw("same");
  //birks_ub->Draw("same");
 
  TLegend* ll = MakeLegend( 0.5, 0.93, 0.035, 5 ); 
  ll->SetBorderSize(1);
  ll->SetFillStyle(1001);
  ll->AddEntry(fbirks,"Birks model at 484 V/cm", "L");
  ll->AddEntry(fbirks,"(ICARUS parameters)", "");
  ll->AddEntry(fmodbox,"Mod. Box model at 484 V/cm","L");
  ll->AddEntry(fmodbox,"(ArgoNeuT parameters)", "");
  ll->AddEntry(h_dEdx,"True dE/dx of Geant4 steps","LF");
//  ll->AddEntry(h_dEdx,"(arb. normalized)","");
  ll->Draw();
  
  
  
  TCanvas* c_recomb2 = new TCanvas("recomb2","recomb2",600,600);
  modbox_lariat->Draw();
  birks_lariat->Draw("same");

  
  
  
   
}



// ###################################################################################
void MuWfmPlots(){
  
  TCanvas* c_muwfm[2][2];
  TF1* f_latefit;
  TH1D* h_dat[2];
  TH1D* h_sim[2];
  TPaveText* hd[2][2];
  TLegend* leg[2][2];

  float taumeas = 0.;
  float err_taumeas = 0.;

  f_latefit = new TF1("latefit","[0]*exp(-x/[1])+[2]*exp(-x/3550)");
  f_latefit->SetParName(0,"A");
  f_latefit->SetParName(1,"Quenched #tau (fixed)");
  f_latefit->SetParName(2,"B");

  // Loop over PMTs
  for(size_t i=0; i<2; i++){
    
    c_muwfm[i][0] = new TCanvas(Form("c_muwfm_pmt%lu_dat",i),Form("c_muwfm_pmt%lu_dat",i),600,600);
    c_muwfm[i][1] = new TCanvas(Form("c_muwfm_pmt%lu_sim",i),Form("c_muwfm_pmt%lu_sim",i),600,600);

    if( !fUsePmt[i] ) continue; 
    h_dat[i] = (TH1D*)fFile->Get(Form("michelana/diagnostics/optical/%lu_AveWfm_crsmu",i));
    h_sim[i] = (TH1D*)fFileMC->Get(Form("michelana/diagnostics/optical/%lu_AveWfm",i));
    h_dat[i] -> Scale( 100./GetHistMax(h_dat[i]) );
    h_sim[i] -> Scale( 100./GetHistMax(h_sim[i]) );
    
    // --------------------------------------
    // Data
    std::cout
    <<"------------------------\n"
    <<"Data PMT "<<i<<"\n";
    c_muwfm[i][0]->cd();
    gStyle->SetOptFit(1112);
    gPad->SetMargin(0.15,0.10,0.10,0.10);
    gPad->SetLogy(1);
    h_dat[i]->GetYaxis()->SetTitleOffset(1.5);
    h_dat[i]->GetYaxis()->SetTitle("Normalized Average Waveform [A.U.]");
    h_dat[i]->GetYaxis()->SetRangeUser(1e-2,200);
    // Find tau_meas
    f_latefit->FixParameter(2,0.);
    f_latefit->SetParameter(0,3.);
    f_latefit->SetParLimits(0,0.,80.);
    f_latefit->ReleaseParameter(1);
    f_latefit->SetParameter(1,1200);
    f_latefit->SetRange(400,2000);
    h_dat[i]->Fit("latefit","RNQ");
    taumeas       = f_latefit->GetParameter(1);
    err_taumeas   = f_latefit->GetParError(1);
    std::cout<<"tau_meas (400-2000ns) = "<<taumeas<<" +/- "<<err_taumeas<<"\n";
    // Now fit to full 7us
    f_latefit->ReleaseParameter(2);
    f_latefit->SetParameter(2,0.001);
    f_latefit->FixParameter(1,TauConversion(taumeas));
    f_latefit->SetRange(400,7000);
    h_dat[i]->Fit(f_latefit,"R");
    h_dat[i]->Draw();    
    f_latefit->DrawCopy("same"); 
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    hd[i][0] = MakeTextBox( 0.5, 0.90-0.02, 0.03, 4);
    hd[i][0] -> AddText("#bf{LArIAT Preliminary}"); 
    hd[i][0] -> AddText(Form("#bf{%s Data: %s}",runtag.c_str(),pmttag[i].c_str()));
    hd[i][0] -> AddText("Crossing Cosmic Muons");
    hd[i][0] -> AddText(Form("#tau_{meas} (0.2-4#mus) = %4.0f #pm %1.0f ns",taumeas,err_taumeas));
    hd[i][0] -> Draw();
    
    leg[i][0] = MakeLegend( 0.5, 0.90-4*0.03-0.02, 0.03, 2); 
    leg[i][0] ->AddEntry(h_dat[i],"PMT signal");
    leg[i][0] ->AddEntry(f_latefit,"Fit to S_{late}(t)");
    leg[i][0] ->Draw();
    
    
    
    
    // --------------------------------------
    // Simulated
    std::cout
    <<"------------------------\n"
    <<"MC PMT "<<i<<"\n";
    c_muwfm[i][1]->cd();
    gStyle->SetOptFit(1112);
    gPad->SetMargin(0.15,0.10,0.10,0.10);
    gPad->SetLogy(1);
    h_sim[i]->GetYaxis()->SetTitleOffset(1.5);
    h_sim[i]->GetYaxis()->SetTitle("Normalized Average Waveform [A.U.]");
    h_sim[i]->GetYaxis()->SetRangeUser(1e-2,200);
    // Find tau_meas
    f_latefit->FixParameter(2,0.);
    f_latefit->SetParameter(0,3.);
    f_latefit->SetParLimits(0,0.,80.);
    f_latefit->ReleaseParameter(1);
    f_latefit->SetParameter(1,1200);
    f_latefit->SetRange(400,2000);
    h_sim[i]->Fit("latefit","RNQ");
    taumeas       = f_latefit->GetParameter(1);
    err_taumeas   = f_latefit->GetParError(1);
    std::cout<<"tau_meas (400-2000ns) = "<<taumeas<<" +/- "<<err_taumeas<<"\n";
    // Now fit to full 7us
    f_latefit->ReleaseParameter(2);
    f_latefit->SetParameter(2,0.001);
    f_latefit->FixParameter(1,TauConversion(taumeas));
    f_latefit->SetRange(400,7000);
    h_sim[i]->Fit(f_latefit,"R");
    h_sim[i]->Draw();    
    f_latefit->DrawCopy("same"); 
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    hd[i][1] = MakeTextBox( 0.5, 0.90-0.02, 0.03, 3);
    hd[i][1] -> AddText("#bf{LArIAT Preliminary}"); 
    hd[i][1] -> AddText(Form("#bf{%s MC}",runtag.c_str()));
    hd[i][1] -> AddText(Form("#tau_{meas} (0.2-4#mus) = %4.0f #pm %1.0f ns",taumeas,err_taumeas));
    hd[i][1] -> Draw();
    
    leg[i][1] = MakeLegend( 0.5, 0.90-3*0.03-0.02, 0.03, 2); 
    leg[i][1] ->AddEntry(h_sim[i],"PMT signal");
    leg[i][1] ->AddEntry(f_latefit,"Fit to S_{late}(t)");
    leg[i][1] ->Draw();


   
   
    /* 
    c_muwfm[i][1]->cd();
    gPad->SetLogy(1);
    h_sim[i]->Draw();    
    */

  }



}






