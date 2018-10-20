#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TMath.h"
#include <TMinuit.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitter.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TLine.h"
#include <time.h>
#include <vector>
#include <math.h>


float kPi = TMath::Pi();

//##########################################################################
// Physical constants used in analysis macro.
float fLArDensity     = 1.370;  // [g/cm^3] (used in LArSoft (90.7K). ArgoNeuT used 1.383)
float fExcRatio       = 0.21;
float fTauT           = 1300.;  // triplet lifetime [ns]:w
float fWph = 19.5e-6; // eV
float fWion = 23.6e-6;
float fBoxModelB      = 0.212;       // ArgoNeuT
float fBoxModelA      = 0.93;        // ArgoNeuT
float fBirksModelA    = 0.80;      // ICARUS
float fBirksModelk    = 0.0486;    // ICARUS
//##########################################################################

float fRecomb = 0.65; 
float fRecombPh = 0.59;
float fRecombIon = 0.67; // 0.69 from Mod Box

TFile* fFile;
TFile* fFileMC;
TFile* fFileMC_ns;
TTree* fTree[2];
TTree* fTreeMC_ns;
TH1D* hDecayTime[2];
TH1D* hdT_mc;
TH1D* hPE_total[2][2];
TH1D* hPE_prompt[2][2];
TH1D* hPE_prompt_pretrigMC[2];
TH1D* hPE_prompt_posttrigMC[2];
TH1D* hPE_prompt_nosmearMC[2];
TH1D* hPE_total_nosmearMC[2];
TH1D* hPE_total_qc[2][2];
TH1D* hCrsMuPE_total[2][2];
TH1D* hCrsMuPE_prompt[2][2];
TH1D* hPE[2];
TH1D* hQ[2];
TH1D* hL[2];
TH1D* hQRes;
TH1D* hLRes;

// Event cut histograms
TH1D* h_evtCut_MuPE_prompt[2][2];
TH1D* h_evtCut_MuClusterSize[2];
TH1D* h_evtCut_ElClusterSize[2];
TH1D* h_evtCut_BraggSlope[2];
TH1D* h_evtCut_MuLinearity[2];
TH1D* h_evtCut_MuDir[2];
TH1D* h_evtCut_DecayAngle[2];
TH1D* h_evtCut_ShowerSize[2];
TH1D* h_evtCut_NumPts3D[2];
TH1D* h_evtCut_FracPts3D[2];
TH1D* h_evtCut_ProjDist[2];

TH1D* hElShowerCentroid_X[2];
TH1D* hElShowerCentroid_Y[2];
TH1D* hElShowerCentroid_Z[2];
TH1D* hElShowerVis[2];
TH3D* hElShowerCentroid[2];

TH2D* hTrue_Dist_vs_Containment;

void CrsMuTest();

TH3D* hTrueMuEndpt;
TH3D* hMCSpatialWgt;
TH2D* hTrueMuEndpt_ZX;
TH2D* hMCSpatialWgt_ZX;
float GetMCSpatialWgt(float x, float y, float z);

TH1D* hdEdx_crsmu;

TH1D* hTrue_NContainShwr;

// Energy spectra
TH1D* hEnergyTrk[2];
TH1D* hEnergyQ[2];
TH1D* hEnergyQ_2R[2];
TH1D* hEnergyQL[2];
TH1D* hEnergyQL_LogL[2];
TH1D* hEnergyQL_LogL_Fail[2];

// 2D resolution histograms (vs. energy)
TH2D* hEvs_Qtrue;
TH2D* hEvs_Q;
TH2D* hEvs_L;
TH2D* hEvs_Ltrue;
TH2D* hEvsRes_E_Trk;
TH2D* hEvsRes_E_Q;
TH2D* hEvsRes_E_Q_2R;
TH2D* hEvsRes_E_QL;
TH2D* hEvsRes_E_QL_LogL;

TH2D* hEvsRes_Q;
TH2D* hEvsRes_L;

//TH2D* h_Qcol_vs_QcolRes;
//TH2D* h_Pe_vs_PeRes;

TH2D* hEnergy_QvsQL[2];
TH1D* hPERes;
TH1D* hEnergyDepResQ;
TH1D* hEnergyDepResQL;
TH1D* hEnergyResQ;
TH1D* hEnergyResQL;
TH2D* hEnergyVsResPE;
TH2D* hEnergyVsResQ;
TH2D* hEnergyVsResL;
//TH2D* hEnergyVsResEQ;
//TH2D* hEnergyVsResEQ_Trk;
//TH2D* hEnergyVsResEQL;
TH1D* hTrue_EnergyDepTrk;
TH1D* hTrue_EnergyDep;
TH1D* hTrue_EnergyDep_PreCuts;
TH1D* hTrue_LY;
TH1D* hTrue_Energy;
TH2D* hTrue_EnergyVsContainment;
TH2D* hTrue_EnergyVsEnergyDepTrk;
TH2D* hTrue_ZX_Vs_Energy;
TH2D* hTrue_ZX_Vs_Containment;
TH1D* hdEdx[2];
//TH1D* hBeta[2];
//TH1D* hTrue_Beta;
TH1D* hR[2];
TH1D* hTrue_R;
TH1D* hTrue_RTrack;
TH1D* hTrue_RShower;
TH1D* hPromptFrac[2];
TH1D* hPromptFracRaw[2];
TH1D* hLightYield[2];
TH2D* hLY_entries_ZX[2];
TH2D* hLY_entries_ZY[2];
TH2D* hLY_ZX[2];
TH2D* hE_ZX[2];
TH2D* hLY_ZY[2];
TH2D* hE_ZY[2];

TH2D* hTrue_Energy_vs_RecoEnergy;
TH2D* hTrue_EnergyDep_vs_RecoEnergy;

TH1D* hShowerAngleMean[2];
TH1D* hShowerAngleRMS[2];
TH2D* hShowerAngle[2];

TH2D* hRunVsEnergy;
TH2D* hRunVsCharge;
TH2D* hRunVsLight;
TH2D* hRunVsElectronLifetime;
TH2D* hDecayTimeVs;
TH2D* hQvsL[2];
TH2D* hTimeVsCharge;
TH2D* hTimeVsLight;
TH2D* hTimeVsElectronLifetime;
TH1D* hChargePerCm;
TH1D* hLightPerCm;
TH2D* hTimeVsChargePerCm;
TH2D* hTimeVsLightPerCm;
TProfile* hxLifetime;
TProfile* hxCharge;
TProfile* hxLight;
TProfile* hxChargePerCm;
TProfile* hxLightPerCm;
TProfile* hx_dQdx;
TProfile* hx_dADCdx;
TProfile* hx_dEdx;

TH2D* hPE_total_vs_prompt[2][2];

TH2D* hPE_total_compare;


TFile* outFile;
TF1* trigEffCut[2];
std::string filenameData;
std::string filenameMC;
std::string filenameMC_ns;
std::string runtag;
std::string pmttag[2];
std::string diamtag[2];

char buffer[100];

std::string fOptimizeTarget;
int fSetStrategy;
int fSetChi2Mode;
int fFixScale;
int fFixFastRatio;
int fFixK[2];
int fFixP[2];
int fFixSmearPrompt;
int fFixSmearTotal;
int fMinimizeBothPMTs;
std::string fSetMinimizer;

float TauConversion(float);

bool fApplyCalibCorr;

int fCh;
int fMaxMCEvents;
int fRunMode;
int fMinRun;
int fMaxRun;
int fRunStart;
int fRunEnd;
int fTimeStart;
int fTimeEnd;
float fDayEnd;

bool fRequire3DMuEndpt;

float fdTcut = 2000;
float fGateDelay = 300;
float fMaxDecayTime = 7200.;
float fTrigEff_P[2];
float fTrigEff_K[2];
float fMinMuPromptPE[2];

bool  fPE_Require3DShower;

float fTrigEff_X1[2];
float fTrigEff_X2[2];

float fDecayFit_T1;
float fDecayFit_T2;
float fDecayFit_MinDecayTime;
float fDecayFit_MaxDecayTime;
float fDecayFit_BG;
float fDecayFit_CPlus;
float fDecayFit_CMinus;
float fDecayFit_TauPlus;
float fDecayFit_TauMinus;
bool  fDecayFit_FixFreeLifetime;
bool  fDecayFit_FixBG;
bool  fDecayFit_ReqShower1;
bool  fDecayFit_ReqShower0;


float fFastRatio;

int   fSmearMode; 
float fSmearFactor[2];
float fSmearFactorPr[2];
float fScaleFactor[2];
float fAmpCut[2];

void Init();
void Init(int mode);
void Init(std::string, std::string, std::string);

float GetHistMax(TH1D*);
void  rebin2(TH1 *h, int ngx, int ngy);

double ModBoxRecomb(double, double, double, double);
double _ModBoxRecomb(double*, double*);
double _ModBoxRecomb_Efield(double*, double*);
double BirksRecomb(double, double, double, double);
double _BirksRecomb(double*, double*);
double _BirksRecomb_Efield(double*, double*);
double dADCdx_from_rr(double*, double*);
double dADCdx_from_dEdx(double*, double*);
double dQdx_from_dEdx(double*, double*);
double MuMPV_from_rr(float rr);
double LandauGaus(Double_t *x,Double_t *par);

TF1* f_LogLEnergyQL;
double _logLEnergyQL(double*,double*);
double _funcP_resQ(double*,double*);
double _funcP_resL(double*,double*);
double _funcP_Q(double*,double*);
double _funcP_L(double*,double*);
double _funcP_Q_Gaus(double*,double*);
double _funcP_L_Gaus(double*,double*);

void ResolutionSliceLoop(
  const TH2D* h2d, float Emin, float Emax,int n, float dE, 
  bool makeFitArray, std::string canvasName, std::string tag, int nx, int ny, float xmin, float xmax,
  int mode, float param, bool hybridRes,
  TGraphAsymmErrors* gr,
  TGraphAsymmErrors* gr_rms,
  std::vector<TF1> &f);

void ResolutionSliceLoop(
  const TH2D* h2d, float Emin, float Emax, int n, float dE,
  bool makeFitArray, std::string canvasName, std::string tag, int nx, int ny, float xmin, float xmax,
  int mode, float param, bool hybridRes,
  TGraphAsymmErrors* gr,
  TGraphAsymmErrors* gr_rms);

void ResolutionSliceLoop(
  const TH2D* h2d, float Emin, float Emax, int n, float dE,
  bool makeFitArray, std::string canvasName, std::string tag, int nx, int ny, float xmin, float xmax,
  int mode, float param, bool hybridRes,
  TGraphAsymmErrors* gr,
  TGraphAsymmErrors* gr_rms,
  TGraphAsymmErrors* gr_mean);

void ParamDistOverlay(
  const TH2D* h2d, float Emin, float Emax, int nslice, float dE, 
  std::string canvasName, std::string tag, int nx, int ny,
  TF1* func_in );

float CorrectForQuenching(float, float, float);

// ===========================

void RepMC();
void LightPlots();
void EnergyPlots(bool in=false);
void dEdxPlots();
void TimePlots();
void AddTextLine();
void AddTextLine(TLatex* t, float x, float y, int lineNum, std::string line);
TH1D* Make1DSlice(const TH2D* h2d, int, int);
TH1D* Make1DSlice(const TH2D* h2d, int, int, std::string);
void Get1DSliceMeanRMS(TH2D* h2d, int, int, float&, float&);
void FitToLandau(TH1D* h,float& mpv, float& mpv_err, float& sigma, float& sigma_err);
void FitToLandau(TH1D* h, float r1, float r2, float& mpv, float& mpv_err, float& sigma, float& sigma_err);
void FitToGaus(TH1D* h,float& mpv, float& mpv_err, float& sigma, float& sigma_err);
void FitToGaus(TH1D* h, float r1, float r2, float& mpv, float& mpv_err, float& sigma, float& sigma_err);
void FitToLandauGaus(TH1D* h,float& mpv, float& mpv_err, float& sigma, float& sigma_err);
void FitToLandauGaus(TH1D* h, float r1, float r2, float& mpv, float& mpv_err, float& sigma, float& sigma_err);

void MakeResRangeGraphs(TH2D*, TGraphErrors*);
void MuRangeCalibration(TFile*, int mc, int ipl);
void MuRangeCalibration(int mc, int ipl);

void PeakExclusionFit(TH1D* h, TF1* f, float x1, float x2 );
void MultiGausFit(TH1D* h, TF1* f2g, TF1* fbg, float param);
void MultiGausFit(TH1D* h, TF1* f2g, TF1* fbg, int mode, float param);


void MakeWeightMap();

void DrawPE();
void ScaleMC();

void MapChi2Space(int ch, double sc1, double sc2, double sig1, double sig2);

void MakeCalibrationPlots();
void StabilityPlots();
float GetChargeCorrFactor(float);
float GetLightCorrFactor(float);

void MakeEventCutSubplot( TCanvas* c, TLine* line1, TLine* line2, int index, TH1D* h_data, TH1D* h_mc, float cutx1, float cutx2, float axts, float axls);

void AnaMichel();

void Loop(int);
void Loop(TTree*, bool, bool);
 
float GetChi2Weighted(TH1D* h1, TH1D* h2, bool verbose=false);
float GetChi2(TH1D* h1, TH1D* h2, bool verbose=false);
float GetChi2(TH1D* h1, TH1D* h2, float x1, float x2, bool useWeight=false, bool verbose=false);
 
float CalcFracChange(int, int);
void SetTrigEffParams(int s=1);

float CalcTruncatedMean( std::vector<float>& v, float p, int skew = 0);

void Optimize(int ch);
void OptimizePrompt(int ch, int n=1);
void OptimizePromptBothPMTs(int n=1);
void OptimizeTotal(int ch, int n=1);
void OptimizeBoth(int ch, int n=1);
void OptimizeTotalOnly(int ch, int n=1);
static void   fcn_trigeff(int& nDim, double* gout, double& result, double *par, int flg);
double        trigeff(double, double,double,double,double, double);


// Event cut params
float fMaxTotalCharge;
int fMinClusterSize;
int fMinElClusterSize;
int fMaxElClusterSize;
int fMinElShowerSize;
int fMaxElShowerSize;
float fMinElTrackFrac;
int fMinMuClusterSize;
int fMaxExtraHits;
float fMinMuLinearity;
float fMinFracMuHitsLinear;
int fMinMuClusterHitsEndFit;
float fMinDecayAngle2D;
float fMaxDecayAngle2D;  
int   fMinNumPts3D;
float fMinFracHits3D;
float fMinShowerFrac;
float fMaxCovAtBnd;
float fMinProjDist;
float fPMTExclRadius;
float fMinBraggSlope;
 
  //========================================
  // Event identifiers 
  int                   fRunNumber;
  int                   fSubRunNumber;
  int                   fEventNumber;
  int                   fEventTime;
  bool                  fMichelOpticalID;
  bool                  fTriggered;
  bool                  fIsRealData;

  // Optical information
//  float                 fWfmBaseline[2];
  float                 fWfmRMS[2];
  float                 fSPE[2];
  float                 fSPE_err[2];
  int                   fNumOpHits[2];
  int                   fNumOpHits0[2];
  float                 fdT[2];
  float                 fAmplitude[2];
  float                 fMuWidth[2];
  float                 fWidth[2];
  
  // Hit times, trigger condition, saturation flag
  std::vector<short>    fvHitTimes0[2];
  std::vector<short>    fvHitTimes[2];
  std::vector<bool>     fvIsHitAtTrigger[2];
  std::vector<bool>     fvIsHitSaturated[2]; 

  // Raw hit integrals
  std::vector<float>    fvHitADC_100ns[2];
  std::vector<float>    fvHitADC_200ns[2];
  std::vector<float>    fvHitADC_300ns[2];
  std::vector<float>    fvHitADC_400ns[2];
  std::vector<float>    fvHitADC_500ns[2];
  std::vector<float>    fvHitADC_600ns[2];
  std::vector<float>    fvHitADC_700ns[2];
  std::vector<float>    fvHitADC_900ns[2];
  std::vector<float>    fvHitADC_1200ns[2];
  std::vector<float>    fvHitADC_1500ns[2];
  std::vector<float>    fvHitADC_1800ns[2];
  std::vector<float>    fvHitADC_7000ns[2];

  // Prefit hit integrals, hit amplitude, width
  std::vector<float>    fvHitAmplitudes[2];
  std::vector<float>    fvHitWidth[2];
  std::vector<float>    fvHitADCpf_100ns[2];
  std::vector<float>    fvHitADCpf_7000ns[2];
  std::vector<short>    fvPrepulseX1[2];
  std::vector<float>    fvPrepulseZeroPoint[2];
  std::vector<float>    fvPrepulseBaseline[2];
  std::vector<float>    fvPrepulseRMS[2];
  std::vector<float>    fvPrepulseSlowNorm[2];
  std::vector<float>    fvPrepulseSlowTau[2];

  // Derived quantities when dT is defined (2 hits)
  float                 fPE_promptRaw[2];
  float                 fPE_totalRaw[2];
  float                 fPE_prompt[2];
  float                 fPE_total[2];
  float                 fPE_total_qc[2];

  float                 fPromptFracRaw[2];
  float                 fPromptFrac[2];
  float                 fMuAmplitude[2];
  float                 fMuPE_prompt[2];
  bool                  fMuPulseSaturated[2];
  bool                  fElPulseSaturated[2];
  float                 fElShowerPhel;
  float                 fElShowerPhel_prompt;
  float                 fElShowerPhel_qc;
  float                 fPEPrompt;
  float                 fDecayTime;
  float                 fTotalPhelCh[2];
  float                 fTotalPhel;
    
  // Track information
  int                   fNumTrack;
  std::vector<TVector3> fvTrackVertex;
  std::vector<TVector3> fvTrackEnd;
  std::vector<float>    fvTrackVertex_X;
  std::vector<float>    fvTrackVertex_Y;
  std::vector<float>    fvTrackVertex_Z;
  std::vector<float>    fvTrackEnd_X;
  std::vector<float>    fvTrackEnd_Y;
  std::vector<float>    fvTrackEnd_Z;
  std::vector<bool>     fvIsTrackStopping;
  std::vector<bool>     fvIsTrackPassing;
  std::vector<bool>     fvIsTrackContained;
  std::vector<float>    fvTrackEnergy;
  std::vector<float>    fvTrackLength;
  std::vector<int>      fvTrackID;
  int                   fNumTrackStopping;
  int                   fNumTrackCrossing;
  int                   fNumTrackContained;
  std::vector<int>      fMuTrackHits;
  int                   fMuTrackIndex;
  int                   fMuTrackID;
  float                 fMuTrackLength;
  float                 fMuTrackEnergy;
  float                 fMuTrackZenithAngle;
  float                 fMuTrackVertex_X;
  float                 fMuTrackVertex_Y;
  float                 fMuTrackVertex_Z;
  float                 fMuTrackEnd_X;
  float                 fMuTrackEnd_Y;
  float                 fMuTrackEnd_Z;
  std::vector<double>    fvMuTrackdEdx;
  std::vector<double>    fvMuTrackResidualRange;
  bool                  fMuTrackIsBraggPeaked;
  bool                  fMuTrackIsCaloOrdered;
  TVector3              fMuTrackVertex;
  TVector3              fMuTrackEnd;

  // 2D clustering data
  int                   fNumPlaneHits[2];
  float                 fAveX;
  float                 fAveDriftTime;
  
  float                 fBraggSlope;

  int                   fClusterSize;
  int                   fMuClusterSize;
  int                   fElClusterSize;
  int                   fExtraHits;
  float                 fMuAveLinearity;
  float                 fCovAtBnd;
  float                 fFracMuHitsLinear;
  int                   fMuClusterHitsEndFit;
  float                 fDecayAngle2D;
  int                   fElShowerSize;
  float                 fElShowerFrac; 
  float                 fTotalCharge[2];
  float                 fTotalChargeCol[2];
  float                 fElShowerCharge;
  float                 fElShowerChargeCol;
  float                 fElTrackCharge;
  float                 fElShowerEnergy;
  float                 fElTrackEnergy;                 
  
  float                 fElDir2D_W;
  float                 fElDir2D_X;
  
  float                 fBraggSlope_Pl0;
  int                   fClusterSize_Pl0;
  int                   fMuClusterSize_Pl0;
  int                   fElClusterSize_Pl0;
  int                   fExtraHits_Pl0;
  float                 fMuAveLinearity_Pl0;
  float                 fFracMuHitsLinear_Pl0;
  int                   fMuClusterHitsEndFit_Pl0;
  float                 fDecayAngle2D_Pl0;
  int                   fElShowerSize_Pl0;
  float                 fElShowerFrac_Pl0; 
  float                 fElShowerCharge_Pl0;
  float                 fElTrackCharge_Pl0;
  float                 fElShowerEnergy_Pl0;
  float                 fElTrackEnergy_Pl0;                 
  
  // 3D:
  float                 fFracHits3D;
  int                   fNumPts3D;
  float                 fElShowerVis;
  float                 fElShowerVisCh[2];
  TVector3              fElShowerCentroid;
  float                 fElShowerCentroid_X;
  float                 fElShowerCentroid_Y;
  float                 fElShowerCentroid_Z;
  float                 fTrue_ElShowerCentroid_X;
  float                 fTrue_ElShowerCentroid_Y;
  float                 fTrue_ElShowerCentroid_Z;
  float                 fElShowerPhotons;
  float                 fElShowerEnergyQL;
  TVector3              fMuEnd3D;
  float                 fMuEnd3D_X;
  float                 fMuEnd3D_Y;
  float                 fMuEnd3D_Z;
  float                 fMuEndHitTimeDiff;


  // Truth information
  int                   fTrue_MuCharge;
  float                 fTrue_MuStopTime;
  TVector3              fTrue_MuTrackVertex;
  TVector3              fTrue_MuTrackEnd;
  float                 fTrue_MuTrackVertex_X;
  float                 fTrue_MuTrackVertex_Y;
  float                 fTrue_MuTrackVertex_Z;
  float                 fTrue_MuTrackEnd_X;
  float                 fTrue_MuTrackEnd_Y;
  float                 fTrue_MuTrackEnd_Z;
  float                 fTrue_MuEnergyDep;
  int                   fTrue_NumBremmPhotons;
  float                 fTrue_ElEnergy;
  float                 fTrue_ElShowerPhotons;
  float                 fTrue_ElShowerVisCh[2];
  float                 fTrue_ElTrackPhotons;
  float                 fTrue_ElShowerPhotonsPrompt;
  float                 fTrue_ElShowerPhotonsLate;
  float                 fTrue_ElTrackEnergyDep;
  float                 fTrue_ElTrackCharge;
  float                 fTrue_ElShowerEnergyDep;
  float                 fTrue_ElShowerCharge;
  float                 fTrue_ElShowerChargeCol;
  TVector3              fTrue_ElMomentum;
  float                 fTrue_ElMomentum_X;
  float                 fTrue_ElMomentum_Y;
  float                 fTrue_ElMomentum_Z;
  float                 fTrue_ElAngle;
  bool                  fTrue_IsBareElShower;
  bool                  fTrue_IsElTrkContained;
  bool                  fTrue_IsElShwrContained;
  float                 fTrue_TotalCharge;
  float                 fTrue_TotalChargeCol;
  float                 fTrue_TotalEnergyDep;
  float                 fTrue_dT;
  float                 fTrue_Amplitude[2];
  float                 fTrue_MuPE_prompt[2];
  float                 fTrue_MuContamination_prompt[2];
  float                 fTrue_MuContamination_total[2];
  float                 fTrue_PE[2];
  float                 fTrue_PE_prompt[2];
  float                 fTrue_PE_total[2];
  float                 fTrue_ElShowerPhel;

  float                 fMuContam_prompt[2];
  float                 fMuContam_total[2];

  float                 fCrsMuLength;
  float                 fCrsMuCharge;
  float                 fCrsMuPhotons; 
  
  float                 fProjDist3DToWires;
  float                 fProjPtOnWires_Y;
  float                 fProjPtOnWires_Z;
  
  float                 fElectronLifetime;
  
  
  std::vector<float>*   fMuContamCorr_EffTau=0;
  std::vector<float>*   fQE_ScaleFactor=0;
//  std::vector<float>*   fMuContamCorr_EffTau;
//  std::vector<float>*   fQE_ScaleFactor;
   

// #####################################################
// Event cut variables

  int nEvt_TotalEvents;
  int nEvt_OpID;
  int nEvt_OpID_timeRange;
  int nEvt_OpID_muPromptCut;
  int nEvt_1StpTrk;
  int nEvt_Shwr2D_BndFound;
  int nEvt_Shwr2D_ClusterSize;
  int nEvt_Shwr2D_ExtraHits;
  int nEvt_Shwr2D_Slope;
  int nEvt_Shwr2D_MuLinearity;
  int nEvt_Shwr2D_MuDir;
  int nEvt_Shwr2D_DecayAngle;
  int nEvt_Shwr2D_ShowerSize;
  int nEvt_Shwr2D_ShowerFrac;
  int nEvt_Shwr3D_ShowersOnBothPlanes;
  int nEvt_Shwr3D_Npts;
  int nEvt_Shwr3D_FracHits3D;
  int nEvt_Shwr3D_Centroid;
  int nEvt_Shwr3D_Direction;
  int nEvt_3DMuEnd;
  int nEvt_dTcut;
  int nEvt_DecayFit_2DShower;
  int nEvt_DecayFit_TimeRange;
  int nEvt_DecayFit_AmpCut;





// ======================================================================================
void setBranches(TTree *tree){
 
  tree  ->SetBranchAddress("IsRealData",&fIsRealData); 
  tree   ->SetBranchAddress("MuContamCorr_EffTau",&fMuContamCorr_EffTau);
  tree   ->SetBranchAddress("QE_ScaleFactor",&fQE_ScaleFactor);
  tree   ->SetBranchAddress("ElectronLifetime",&fElectronLifetime);
  tree   ->SetBranchAddress("EventTime",&fEventTime);
  tree   ->SetBranchAddress("RunNumber",     &fRunNumber);
//  tree   ->SetBranchAddress("SubRunNumber",  &fSubRunNumber);
//  tree   ->SetBranchAddress("EventNumber",   &fEventNumber);
  tree   ->SetBranchAddress("MichelOpticalID",&fMichelOpticalID);
//  tree   ->SetBranchAddress("Triggered",&fTriggered);
 // tree   ->SetBranchAddress("NumOpHits",     &fNumOpHits);
//  tree   ->SetBranchAddress("NumOpHits0",     &fNumOpHits0);
  //tree   ->SetBranchAddress("SPE",           &fSPE);
//  tree   ->SetBranchAddress("SPE_err",       &fSPE_err);
  //tree   ->SetBranchAddress("dT",            &fdT);
  tree   ->SetBranchAddress("Amplitude",     &fAmplitude);
  tree   ->SetBranchAddress("Width",         &fWidth);
  //tree   ->SetBranchAddress("MuWidth",       &fMuWidth);
  tree   ->SetBranchAddress("PE_prompt",     &fPE_prompt);
  tree   ->SetBranchAddress("PE_total",      &fPE_total);
  tree   ->SetBranchAddress("PE_total_qc",        &fPE_total_qc);
  tree   ->SetBranchAddress("PE_promptRaw",     &fPE_promptRaw);
  tree   ->SetBranchAddress("PE_totalRaw",      &fPE_totalRaw);
  tree   ->SetBranchAddress("PromptFrac",    &fPromptFrac);
  tree   ->SetBranchAddress("MuAmplitude",   &fMuAmplitude);
  tree   ->SetBranchAddress("MuPE_prompt",   &fMuPE_prompt);
  tree   ->SetBranchAddress("MuPulseSaturated",&fMuPulseSaturated);
  tree   ->SetBranchAddress("ElPulseSaturated",&fElPulseSaturated);
  tree   ->SetBranchAddress("ElShowerPhel",      &fElShowerPhel);
  tree   ->SetBranchAddress("ElShowerPhel_prompt",      &fElShowerPhel_prompt);
  tree   ->SetBranchAddress("ElShowerPhel_qc",      &fElShowerPhel_qc);
  tree   ->SetBranchAddress("MuContam_prompt",  &fMuContam_prompt);
  tree    ->SetBranchAddress("TotalPhelCh",     &fTotalPhelCh);
  tree   ->SetBranchAddress("MuContam_total",  &fMuContam_total);
//  tree   ->SetBranchAddress("NumTrack", &fNumTrack);
//  tree   ->SetBranchAddress("NumTrackStopping", &fNumTrackStopping);
//  tree   ->SetBranchAddress("NumTrackCrossing", &fNumTrackCrossing);
  tree   ->SetBranchAddress("DecayTime",     &fDecayTime);
  tree   ->SetBranchAddress("True_ElEnergy",        &fTrue_ElEnergy);
 
//  tree   ->SetBranchAddress("NumPlaneHits",&fNumPlaneHits);
  
  tree   ->SetBranchAddress("ClusterSize",      &fClusterSize);
  tree   ->SetBranchAddress("MuClusterSize",      &fMuClusterSize);
  tree   ->SetBranchAddress("ElClusterSize",      &fElClusterSize);
//  tree   ->SetBranchAddress("ExtraHits",&fExtraHits);
  tree   ->SetBranchAddress("MuAveLinearity",&fMuAveLinearity);
  tree   ->SetBranchAddress("FracMuHitsLinear",&fFracMuHitsLinear);
  tree   ->SetBranchAddress("MuClusterHitsEndFit",&fMuClusterHitsEndFit);
  tree   ->SetBranchAddress("DecayAngle2D",&fDecayAngle2D);
  tree   ->SetBranchAddress("ElShowerSize",&fElShowerSize);
  tree   ->SetBranchAddress("ElShowerFrac",      &fElShowerFrac);  
  tree   ->SetBranchAddress("TotalCharge", &fTotalCharge);
  tree   ->SetBranchAddress("TotalChargeCol", &fTotalChargeCol);
  tree   ->SetBranchAddress("ElShowerCharge",        &fElShowerCharge);
  tree   ->SetBranchAddress("ElShowerChargeCol",        &fElShowerChargeCol);
  tree   ->SetBranchAddress("ElShowerEnergy",        &fElShowerEnergy);
  tree   ->SetBranchAddress("ElTrackCharge",        &fElTrackCharge);
  tree   ->SetBranchAddress("ElTrackEnergy",        &fElTrackEnergy);
  
  tree   ->SetBranchAddress("ElShowerSize_Pl0",&fElShowerSize_Pl0);
  /*
  tree   ->SetBranchAddress("ClusterSize_Pl0",      &fClusterSize_Pl0);
  tree   ->SetBranchAddress("MuClusterSize_Pl0",      &fMuClusterSize_Pl0);
  tree   ->SetBranchAddress("ElClusterSize_Pl0",      &fElClusterSize_Pl0);
  tree   ->SetBranchAddress("ExtraHits_Pl0",&fExtraHits_Pl0);
  tree   ->SetBranchAddress("MuAveLinearity_Pl0",&fMuAveLinearity_Pl0);
  tree   ->SetBranchAddress("FracMuHitsLinear_Pl0",&fFracMuHitsLinear_Pl0);
  tree   ->SetBranchAddress("MuClusterHitsEndFit_Pl0",&fMuClusterHitsEndFit_Pl0);
  tree   ->SetBranchAddress("DecayAngle2D_Pl0",&fDecayAngle2D_Pl0);
  tree   ->SetBranchAddress("ElShowerAngleMean_Pl0",&fElShowerAngleMean_Pl0);
  tree   ->SetBranchAddress("ElShowerAngleRMS_Pl0",&fElShowerAngleRMS_Pl0);
  tree   ->SetBranchAddress("ElShowerFrac_Pl0",      &fElShowerFrac_Pl0);  
  tree   ->SetBranchAddress("ElShowerCharge_Pl0",        &fElShowerCharge_Pl0);
  tree   ->SetBranchAddress("ElShowerEnergy_Pl0",        &fElShowerEnergy_Pl0);
  tree   ->SetBranchAddress("ElTrackCharge_Pl0",        &fElTrackCharge_Pl0);
  tree   ->SetBranchAddress("ElTrackEnergy_Pl0",        &fElTrackEnergy_Pl0);
  */
  tree   ->SetBranchAddress("MuEndHitTimeDiff",          &fMuEndHitTimeDiff);

  tree    ->SetBranchAddress("CovAtBnd", &fCovAtBnd);

  tree    ->SetBranchAddress("BraggSlope",&fBraggSlope);
  //tree    ->SetBranchAddress("BraggSlope_Pl0",&fBraggSlope_Pl0);
  
  tree    ->SetBranchAddress("MuTrackEnd_X",&fMuTrackEnd_X);
  tree    ->SetBranchAddress("MuTrackEnd_Y",&fMuTrackEnd_Y);
  tree    ->SetBranchAddress("MuTrackEnd_Z",&fMuTrackEnd_Z);

//  tree   ->SetBranchAddress("ElShowerEnergyQL",        &fElShowerEnergyQL);
  tree  ->SetBranchAddress("True_dT",                 &fTrue_dT);
  tree   ->SetBranchAddress("ElShowerPhotons",        &fElShowerPhotons);
  tree   ->SetBranchAddress("True_MuTrackEnd_X",     &fTrue_MuTrackEnd_X);
  tree   ->SetBranchAddress("True_MuTrackEnd_Y",     &fTrue_MuTrackEnd_Y);
  tree   ->SetBranchAddress("True_MuTrackEnd_Z",     &fTrue_MuTrackEnd_Z);
  tree   ->SetBranchAddress("True_ElShowerEnergyDep",&fTrue_ElShowerEnergyDep);
  tree   ->SetBranchAddress("MuEnd3D_X",             &fMuEnd3D_X);
  tree   ->SetBranchAddress("MuEnd3D_Y",             &fMuEnd3D_Y);
  tree   ->SetBranchAddress("MuEnd3D_Z",             &fMuEnd3D_Z);
//  tree  ->SetBranchAddress("ElDir2D_W",             &fElDir2D_W);
//  tree  ->SetBranchAddress("ElDir2D_X",             &fElDir2D_X);
  tree   ->SetBranchAddress("FracHits3D",             &fFracHits3D);
  tree   ->SetBranchAddress("NumPts3D",             &fNumPts3D);
  tree   ->SetBranchAddress("ElShowerVis",             &fElShowerVis);
  tree   ->SetBranchAddress("ElShowerVisCh",             &fElShowerVisCh);
  tree   ->SetBranchAddress("True_ElShowerVisCh",            &fTrue_ElShowerVisCh);
  tree   ->SetBranchAddress("True_PE_total",     &fTrue_PE_total);
  tree   ->SetBranchAddress("True_PE_prompt",     &fTrue_PE_prompt);
  tree   ->SetBranchAddress("True_ElShowerPhotons", &fTrue_ElShowerPhotons);
  tree   ->SetBranchAddress("True_ElShowerPhotonsPrompt", &fTrue_ElShowerPhotonsPrompt);
  tree   ->SetBranchAddress("True_ElShowerPhotonsLate", &fTrue_ElShowerPhotonsLate);
  tree   ->SetBranchAddress("True_ElShowerCharge", &fTrue_ElShowerCharge);
  tree   ->SetBranchAddress("True_ElShowerChargeCol", &fTrue_ElShowerChargeCol);
  tree  ->SetBranchAddress("True_TotalChargeCol",&fTrue_TotalChargeCol);
  tree   ->SetBranchAddress("ElShowerCentroid_X",&fElShowerCentroid_X);
  tree   ->SetBranchAddress("ElShowerCentroid_Y",&fElShowerCentroid_Y);
  tree   ->SetBranchAddress("ElShowerCentroid_Z",&fElShowerCentroid_Z);
  tree   ->SetBranchAddress("True_ElShowerCentroid_X",&fTrue_ElShowerCentroid_X);
  tree   ->SetBranchAddress("True_ElShowerCentroid_Y",&fTrue_ElShowerCentroid_Y);
  tree   ->SetBranchAddress("True_ElShowerCentroid_Z",&fTrue_ElShowerCentroid_Z);
//  tree  ->SetBranchAddress("ProjDist3DToWires",&fProjDist3DToWires);
//  tree  ->SetBranchAddress("ProjPtOnWires_Y",&fProjPtOnWires_Y);
//  tree  ->SetBranchAddress("ProjPtOnWires_Z",&fProjPtOnWires_Z);
//  tree  ->SetBranchAddress("AveX",&fAveX);
  tree  ->SetBranchAddress("AveDriftTime",&fAveDriftTime);
  tree  ->SetBranchAddress("True_IsElTrkContained",&fTrue_IsElTrkContained);
  tree  ->SetBranchAddress("True_IsElShwrContained",&fTrue_IsElShwrContained);
  tree  ->SetBranchAddress("True_ElTrackEnergyDep",&fTrue_ElTrackEnergyDep);
  //tree   ->SetBranchAddress("CrsMuLength",&fCrsMuLength);
  //tree   ->SetBranchAddress("CrsMuCharge",&fCrsMuCharge);
  //tree   ->SetBranchAddress("CrsMuPhotons",&fCrsMuPhotons);
}



// ######################################################################
// Various tools and algorithms to make life easier

float GetHistMax(TH1D* h) {
  return h->GetBinContent(h->GetMaximumBin() );
}

void AddTextLine(TLatex* t, float x, float y, int lineNum, std::string line){
  t->DrawLatex( x, y-(lineNum-1)*t->GetTextSize(), line.c_str());
}


TPaveText* MakeTextBox(float x, float y, float textSize, float numLines, float width){
  TPaveText *pt = new TPaveText(x, y-numLines*textSize, x+width, y, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextSize(textSize);
  pt->SetTextFont(42);
  pt->SetTextAlign(12);
  return pt;
}

TPaveText* MakeTextBox(float x, float y, float textSize, float numLines ){
  return MakeTextBox(x,y,textSize,numLines,0.4);
}


TLegend* MakeLegend(float x1, float y2, float textSize, float numLines, float width){
    TLegend *leg = new TLegend(x1, y2-numLines*textSize, x1+width, y2 );
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(textSize);
    return leg;
}

TLegend* MakeLegend(float x1, float y2, float textSize, float numLines){
  return MakeLegend(x1, y2, textSize, numLines, 0.3);
}

void FormatAxes(TH1D* h, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  h->GetXaxis()->SetLabelSize(axisLabelSize);
  h->GetYaxis()->SetLabelSize(axisLabelSize);
  h->GetXaxis()->SetTitleSize(axisTitleSize);
  h->GetYaxis()->SetTitleSize(axisTitleSize);
  h->GetXaxis()->SetTitleOffset(xOffset);
  h->GetYaxis()->SetTitleOffset(yOffset);
}
void FormatAxes(TGraphErrors* g, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  g->GetXaxis()->SetLabelSize(axisLabelSize);
  g->GetYaxis()->SetLabelSize(axisLabelSize);
  g->GetXaxis()->SetTitleSize(axisTitleSize);
  g->GetYaxis()->SetTitleSize(axisTitleSize);
  g->GetXaxis()->SetTitleOffset(xOffset);
  g->GetYaxis()->SetTitleOffset(yOffset);
}
void FormatAxes(TMultiGraph* g, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  g->GetXaxis()->SetLabelSize(axisLabelSize);
  g->GetYaxis()->SetLabelSize(axisLabelSize);
  g->GetXaxis()->SetTitleSize(axisTitleSize);
  g->GetYaxis()->SetTitleSize(axisTitleSize);
  g->GetXaxis()->SetTitleOffset(xOffset);
  g->GetYaxis()->SetTitleOffset(yOffset);
}
void FormatAxes(TGraphAsymmErrors* g, float axisTitleSize, float axisLabelSize, float xOffset, float yOffset ){
  g->GetXaxis()->SetLabelSize(axisLabelSize);
  g->GetYaxis()->SetLabelSize(axisLabelSize);
  g->GetXaxis()->SetTitleSize(axisTitleSize);
  g->GetYaxis()->SetTitleSize(axisTitleSize);
  g->GetXaxis()->SetTitleOffset(xOffset);
  g->GetYaxis()->SetTitleOffset(yOffset);
}

void FormatTGraph(TGraphAsymmErrors* g, Color_t mc, Color_t lc, int ms, int ls, float msize){
  g->SetMarkerColor(mc);
  g->SetLineColor(lc);
  g->SetMarkerStyle(ms);
  g->SetLineStyle(ls);
  g->SetMarkerSize(msize);
}

void FormatTGraph(TGraphErrors* g, Color_t mc, Color_t lc, int ms, int ls, float msize){
  g->SetMarkerColor(mc);
  g->SetLineColor(lc);
  g->SetMarkerStyle(ms);
  g->SetLineStyle(ls);
  g->SetMarkerSize(msize);
}

void CopyTGraphFormat(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, bool copyTitles = false){
  g2->SetMarkerColor(g1->GetMarkerColor());
  g2->SetMarkerStyle(g1->GetMarkerStyle());
  g2->SetMarkerSize(g1->GetMarkerSize());
  g2->SetLineColor(g1->GetLineColor());
  g2->SetLineStyle(g1->GetLineStyle());
  g2->SetLineWidth(g1->GetLineWidth());
  FormatAxes(g2, 
    g1->GetXaxis()->GetTitleSize(),
    g1->GetXaxis()->GetLabelSize(), 
    g1->GetXaxis()->GetTitleOffset(),
    g1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    g2->GetXaxis()->SetTitle( g1->GetXaxis()->GetTitle() );
    g2->GetXaxis()->SetTitleOffset( g1->GetXaxis()->GetTitleOffset() );
    g2->GetYaxis()->SetTitle( g1->GetYaxis()->GetTitle() );
    g2->GetYaxis()->SetTitleOffset( g1->GetYaxis()->GetTitleOffset() );
  }
}

void CopyTGraphFormat(TGraphErrors* g1, TGraphErrors* g2, bool copyTitles = false){
  g2->SetMarkerColor(g1->GetMarkerColor());
  g2->SetMarkerStyle(g1->GetMarkerStyle());
  g2->SetMarkerSize(g1->GetMarkerSize());
  g2->SetLineColor(g1->GetLineColor());
  g2->SetLineStyle(g1->GetLineStyle());
  g2->SetLineWidth(g1->GetLineWidth());
  FormatAxes(g2, 
    g1->GetXaxis()->GetTitleSize(),
    g1->GetXaxis()->GetLabelSize(), 
    g1->GetXaxis()->GetTitleOffset(),
    g1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    g2->GetXaxis()->SetTitle( g1->GetXaxis()->GetTitle() );
    g2->GetXaxis()->SetTitleOffset( g1->GetXaxis()->GetTitleOffset() );
    g2->GetYaxis()->SetTitle( g1->GetYaxis()->GetTitle() );
    g2->GetYaxis()->SetTitleOffset( g1->GetYaxis()->GetTitleOffset() );
  }
}

void CopyHistoFormat(TH1D* h1, TH1D* h2, bool copyTitles = false){
  h2->SetFillColor( h1->GetFillColor() );
  h2->SetFillStyle( h1->GetFillStyle() );
  h2->SetMarkerColor(h1->GetMarkerColor());
  h2->SetMarkerStyle(h1->GetMarkerStyle());
  h2->SetMarkerSize(h1->GetMarkerSize());
  h2->SetLineColor(h1->GetLineColor());
  h2->SetLineStyle(h1->GetLineStyle());
  h2->SetLineWidth(h1->GetLineWidth());
  h2->GetXaxis()->SetTitleSize( h1->GetXaxis()->GetTitleOffset() );
  h2->SetTitleSize( h1->GetTitleSize() );
  FormatAxes(h2, 
    h1->GetXaxis()->GetTitleSize(),
    h1->GetXaxis()->GetLabelSize(), 
    h1->GetXaxis()->GetTitleOffset(),
    h1->GetYaxis()->GetTitleOffset());
  if( copyTitles ) {
    h2->SetTitle( h1->GetTitle() );
    h2->GetXaxis()->SetTitle( h1->GetXaxis()->GetTitle() );
    h2->GetXaxis()->SetTitleOffset( h1->GetXaxis()->GetTitleOffset() );
    h2->GetYaxis()->SetTitle( h1->GetYaxis()->GetTitle() );
    h2->GetYaxis()->SetTitleOffset( h1->GetYaxis()->GetTitleOffset() );
  }
}

void ScaleHistoDataMC(TH1D* h1, TH1D* h2) {
    double nd = h1->Integral();
    double nmc = h2->Integral();
    if( nd >= 0 && nmc >= 0 )
      h2->Scale( nd / nmc );
}

void ScaleHistoDataMC(TH1D** h) {
  ScaleHistoDataMC(h[0],h[1]);
}

