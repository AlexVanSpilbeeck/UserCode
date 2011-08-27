#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>


#include "TKey.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TLegend.h"
#include "TList.h"
#include "TStyle.h"
#include "TObject.h"
#include "TMean.cc"

#include "PixelNameTranslator.cc"

using namespace std;

TString RunNumber = "_Run";
bool PrintImages = true;
bool PrintAllFormats = false;
bool printSummary = true;
bool printTexSummary = true;

//******** Declaration of fonctions ********
vector<TString> makedirlist(TFile*,TString);
void write(TH1F* , Bool_t = 0 , Bool_t = 0 , TString = "rmou"); // logY , largePad , optstatoption
void write(TH2F* , Bool_t = 0);                                 // logZ
void writeAllLevels(TH1F*,TH1F*,TH1F*,TString);
int rocId(int, int);
void getBarrelPos(TString& , int& , int& , int& );
int  getDetId(TString);

//Making PixNameTranslator object available from everywhere.
PixelNameTranslator pix;
  
  
class map_module_level{
  public:
  vector<TH2F> barrel;
  vector<TH2F> endcap;
  
  TH1F SUM_barrel;
  TH1F SUM_endcap;
  TH1F SUM_whole;
  TH1F TOT_barrel;
  TH1F TOT_endcap;
  
  TH1F SUB_whole;
  vector<TMean> m_whole;
  
  map_module_level(TString name , int nbin = 10 , double xmin = 0 , double xmax = 10){
    barrel.assign(3,TH2F());
    barrel[0] = TH2F("MAP_"+name+"_layer1" , name+"_layer1;module;ladder" , 9 , -4.5 , 4.5 , 21 , -10.5 , 10.5);
    barrel[1] = TH2F("MAP_"+name+"_layer2" , name+"_layer2;module;ladder" , 9 , -4.5 , 4.5 , 33 , -16.5 , 16.5);
    barrel[2] = TH2F("MAP_"+name+"_layer3" , name+"_layer3;module;ladder" , 9 , -4.5 , 4.5 , 45 , -22.5 , 22.5);
    
    endcap.assign(4,TH2F());
    endcap[0] = TH2F("MAP_"+name+"_disk-2" , name+"_disk-2;2*blade+panel-1;plaquette" , 50, -24., 26., 4, 1., 5.);
    endcap[1] = TH2F("MAP_"+name+"_disk-1" , name+"_disk-1;2*blade+panel-1;plaquette" , 50, -24., 26., 4, 1., 5.);
    endcap[2] = TH2F("MAP_"+name+"_disk+1" , name+"_disk+1;2*blade+panel-1;plaquette" , 50, -24., 26., 4, 1., 5.);
    endcap[3] = TH2F("MAP_"+name+"_disk+2" , name+"_disk+2;2*blade+panel-1;plaquette" , 50, -24., 26., 4, 1., 5.);
    
    SUM_barrel = TH1F("SUM_"+name+"_barrel" , "SUM_"+name+"_barrel;;mean per module VS ordered detid" , 900 , 0 , 900);
    SUM_endcap = TH1F("SUM_"+name+"_endcap" , "SUM_"+name+"_endcap;;mean per module VS ordered detid" , 900 , 0 , 900);
    SUM_whole  = TH1F("SUM_"+name+"_whole"  , "SUM_"+name+"_whole ;;mean per module VS ordered detid" , 1450, 0 , 1450);
    
    TOT_barrel = TH1F("TOT_"+name+"_barrel" , "TOT_"+name+"_barrel;module means;# modules" , nbin , xmin , xmax);
    TOT_endcap = TH1F("TOT_"+name+"_endcap" , "TOT_"+name+"_endcap;module means;# modules" , nbin , xmin , xmax);
    
    SUB_whole  = TH1F("SUB_"+name+"_whole" , "SUB_"+name+"_whole;;mean per subdetectors" , 7 , 0 , 7);
    m_whole.assign(7 , TMean());
    
  }
  
  void fill(int detid , TString& dirname , double val){
    int layer = 0, ladder = 0, module = 0;
    int disk = 0, blade = 0, panel = 0 , plaquette = 0 ;
    ostringstream detid_str(""); detid_str << detid;
    
    SUM_whole.Fill(detid_str.str().c_str() , val);
      
    if(pix.detID2Bpix(detid , layer , ladder , module)){
      this->barrel[layer-1].Fill(module , ladder , val);
      SUM_barrel.Fill(detid_str.str().c_str() , val);
      TOT_barrel.Fill(val);
      m_whole[layer-1].Add(val);
      //cout << detid << "  " << layer << "  " << ladder << "  " << module << "  " << endl;
    }
    else if(pix.detID2Fpix(detid, disk, blade, panel, plaquette)){
      if(disk>0) disk=disk-1;
      this->endcap[disk+2].Fill(2*blade + panel -1 , plaquette , val);
      SUM_endcap.Fill(detid_str.str().c_str() , val);
      TOT_endcap.Fill(val);
      m_whole[3+disk+2].Add(val);
    }
  }
  
  void Write(){
    for(unsigned i = 0 ; i < barrel.size() ; ++i)
      write(&(barrel[i]));
    for(unsigned i = 0 ; i < endcap.size() ; ++i)
      write(&(endcap[i]));
    
    
    SUM_barrel.LabelsDeflate();
    SUM_barrel.LabelsOption("a");
    SUM_barrel.GetXaxis()->SetNdivisions(-1,0);
    SUM_endcap.LabelsDeflate();
    SUM_endcap.LabelsOption("a");
    SUM_endcap.GetXaxis()->SetNdivisions(-1,0);
    SUM_whole.LabelsDeflate();
    SUM_whole.LabelsOption("a");
    SUM_whole.GetXaxis()->SetNdivisions(-1,0);
    SUM_barrel.GetXaxis()->SetTickLength(0);
    SUM_endcap.GetXaxis()->SetTickLength(0);
    SUM_whole.GetXaxis()->SetTickLength(0);
    write(&SUM_barrel,0,0,"n");
    write(&SUM_endcap,0,0,"n");
    write(&SUM_whole,0,1,"n");
    
    write(&TOT_barrel , 1);
    write(&TOT_endcap , 1);
    
    SUB_whole.Fill("layer 1" , m_whole[0].GetMean());
    SUB_whole.Fill("layer 2" , m_whole[1].GetMean());
    SUB_whole.Fill("layer 3" , m_whole[2].GetMean());
    SUB_whole.Fill("disk -2" , m_whole[3].GetMean());
    SUB_whole.Fill("disk -1" , m_whole[4].GetMean());
    SUB_whole.Fill("disk +1" , m_whole[5].GetMean());
    SUB_whole.Fill("disk +2" , m_whole[6].GetMean());
    for(int i=0 ; i< 7 ; ++i)
      SUB_whole.SetBinError( i + 1 , m_whole[i].GetRMS());
    write(&SUB_whole,0,0,"n");
    
    
  }
  
  void Divide(map_module_level& m){
    for(int i = 0 ; i < this->barrel.size() ; ++i)
      this->barrel[i].Divide(&(m.barrel[i]));
    for(int i = 0 ; i < this->endcap.size() ; ++i)
      this->endcap[i].Divide(&(m.endcap[i]));
  
    //not finished !!!
  
  }
  
  private:

};


class pix_val{
  public:
  
    Double_t gain , pedestal , err_gain , err_pedestal , plateau , dynamic_range , low_point , high_point , 
                 chi2NDF , chi2_prob;
    Int_t fit_result , npoints ;
    
};




//******** MAIN FUNCTION ********
void make_ComparisonPlots(TString fname="" , TString run="" , TString fname2="" , TString run2="" ){

  Bool_t diff = 0;
  
  
  RunNumber+=run;

  pix.init();

  gStyle->SetPalette(1);
  
  TFile* file2 = TFile::Open(fname2,"READ");
  if(fname2 != "" && run2 != "") diff = 1;
  if(diff) RunNumber+="-"+run2;
  
  
  double frac_flag_high_gain = 5; //flag pixels when gain is higher than this
  if(diff) frac_flag_high_gain = 1;
 
  TFile* file = TFile::Open(fname,"READ");
  file->cd();
  
  // Getting list of modules to loop on
  std::vector<TString> dirlist;
  dirlist = makedirlist(file,"Module");
  TDirectory* dir;
  TList* list;
  
  Double_t bmin = 999999999999 , bmax = 0 , fmin = 999999999999 , fmax = 0;
  
  
  //*******************************
  //Declaration of histograms
  
  //gStyle->SetOptStat(0);
  
  //----------------------------------------
  //From File
  TH2F* Gain2d = new TH2F();
  TH2F* ErrorGain2d = new TH2F();
  TH2F* Pedestal2d = new TH2F();
  TH2F* ErrorPedestal2d = new TH2F();
  TH2F* GainSaturate2d = new TH2F();
  TH2F* GainDynamicRange2d = new TH2F();
  TH2F* GainFitResult2d = new TH2F();
  TH2F* GainChi2NDF2d = new TH2F();
  TH2F* GainChi2Prob2d = new TH2F();
  TH2F* GainHighPoint2d = new TH2F();
  TH2F* GainLowPoint2d = new TH2F();
  TH1F* GainNPoints1d = new TH1F();
  
  //From Second File
  TH2F* Gain2d_2 		= new TH2F();
  TH2F* ErrorGain2d_2 		= new TH2F();
  TH2F* Pedestal2d_2 		= new TH2F();
  TH2F* ErrorPedestal2d_2 	= new TH2F();
  TH2F* GainSaturate2d_2 	= new TH2F();
  TH2F* GainDynamicRange2d_2	= new TH2F();
  TH2F* GainFitResult2d_2 	= new TH2F();
  TH2F* GainChi2NDF2d_2 	= new TH2F();
  TH2F* GainChi2Prob2d_2 	= new TH2F();
  TH2F* GainHighPoint2d_2 	= new TH2F();
  TH2F* GainLowPoint2d_2 	= new TH2F();
  TH1F* GainNPoints1d_2 	= new TH1F();
  
  
  //----------------------------------------
  //Summaries
  double nmod = 1416;
  TH1F* SUM_Gain2d = new TH1F("SUM_Gain2d","SUM_Gain2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_ErrorGain2d = new TH1F("SUM_ErrorGain2d","SUM_ErrorGain2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_Pedestal2d = new TH1F("SUM_Pedestal2d","SUM_Pedestal2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_ErrorPedestal2d = new TH1F("SUM_ErrorPedestal2d","SUM_ErrorPedestal2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainSaturate2d = new TH1F("SUM_GainSaturate2d","SUM_GainSaturate2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainDynamicRange2d = new TH1F("SUM_GainDynamicRange2d","SUM_GainDynamicRange2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainFitResult2d = new TH1F("SUM_GainFitResult2d","SUM_GainFitResult2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainChi2NDF2d = new TH1F("SUM_GainChi2NDF2d","SUM_GainChi2NDF2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainChi2Prob2d = new TH1F("SUM_GainChi2Prob2d","SUM_GainChi2Prob2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainHighPoint2d = new TH1F("SUM_GainHighPoint2d","SUM_GainHighPoint2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainLowPoint2d = new TH1F("SUM_GainLowPoint2d","SUM_GainLowPoint2d",(int)nmod,0.5,nmod+0.5);
  TH1F* SUM_GainNPoints1d = new TH1F("SUM_GainNPoints1d","SUM_GainNPoints1d",(int)nmod,0.5,nmod+0.5);
  
  
  //----------------------------------------
  //Binning
  
  int nbin_gain = 400 ; 	  double xmin_gain = 0 ,	  xmax_gain = 10;
  int nbin_errgain = 400 ;	  double xmin_errgain = 0 ,	  xmax_errgain = 0.01;
  int nbin_ped = 1000 ; 	  double xmin_ped = -100 ,	  xmax_ped = 150;
  int nbin_errped = 1000 ;	  double xmin_errped = -1 ,	  xmax_errped = 1;
  
  if(diff){
    nbin_gain = 400 ; 	  xmin_gain = -6 ;	xmax_gain = 6;
    nbin_errgain = 400 ;  xmin_errgain = -0.1 ; xmax_errgain = 0.1;
    nbin_ped = 1000 ; 	  xmin_ped = -50 ;	xmax_ped = 50;
    nbin_errped = 1000 ;  xmin_errped = -10 ;	xmax_errped = 10;
  }
  
  //----------------------------------------
  //Maps
  map_module_level MAP_Gain("Gain" , 50 , xmin_gain , xmax_gain);
  map_module_level MAP_ErrorGain("ErrorGain" , 50 , xmin_errgain , xmax_errgain);
  map_module_level MAP_Pedestal("Pedestal" , 50 , xmin_ped , xmax_ped);
  map_module_level MAP_ErrorPedestal("ErrorPedestal" , 50 , xmin_errped , xmax_errped);
  map_module_level MAP_GainSaturate("GainSaturate" , 50 , -10 , 10);
  map_module_level MAP_GainDynamicRange("GainDynamicRange" , 50 , -25 , 25);
  map_module_level MAP_GainFitResult("GainFitResult" , 50 , -1 , 1);
  map_module_level MAP_GainChi2NDF("GainChi2NDF" , 50 , -10 , 10);
  map_module_level MAP_GainChi2Prob("GainChi2Prob" , 50 , -0.1 , 0.1);
  map_module_level MAP_GainHighPoint("GainHighPoint" , 50 , -100 , 100);
  map_module_level MAP_GainLowPoint("GainLowPoint" , 50 , -10 , 10);
  map_module_level MAP_GainNPoints("GainNPoints" , 50 , -2 , 2);
    
  map_module_level MAP_frac_gain_high("frac_gain_high" , 100 , 0 , 1);
  
  //----------------------------------------
  //Pixel Level
  TH1F* TOT_Gain = new TH1F("TOT_Gain","TOT_Gain_PixelLevel;Gain;# pixels",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGain = new TH1F("TOT_ErrorGain","TOT_ErrorGain_PixelLevel;Error on 1/Gain;#pixels",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_Pedestal = new TH1F("TOT_Pedestal","TOT_Pedestal_PixelLevel;Pedestal;# pixels",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestal = new TH1F("TOT_ErrorPedestal","TOT_ErrorPedestal_PixelLevel;Error on Pedestal;# pixels",nbin_errped,xmin_errped,xmax_errped);
  
  TH1F* TOT_GainBPix = new TH1F("TOT_GainBPix","TOT_GainBPix_PixelLevel;Gain;# pixels",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGainBPix = new TH1F("TOT_ErrorGainBPix","TOT_ErrorGainBPix_PixelLevel;Error on 1/Gain;# pixels",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_PedestalBPix = new TH1F("TOT_PedestalBPix","TOT_PedestalBPix_PixelLevel;Pedestal;# pixels",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestalBPix = new TH1F("TOT_ErrorPedestalBPix","TOT_ErrorPedestalBPix_PixelLevel;Error on Pedestal;# pixels",nbin_errped,xmin_errped,xmax_errped);
  TH1F* TOT_GainFPix = new TH1F("TOT_GainFPix","TOT_GainFPix_PixelLevel;Gain;# pixels",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGainFPix = new TH1F("TOT_ErrorGainFPix","TOT_ErrorGainFPix_PixelLevel;Error on 1/Gain;# pixels",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_PedestalFPix = new TH1F("TOT_PedestalFPix","TOT_PedestalFPix_PixelLevel;Pedestal;# pixels",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestalFPix = new TH1F("TOT_ErrorPedestalFPix","TOT_ErrorPedestalFPix_PixelLevel;Error on Pedestal;# pixels",nbin_errped,xmin_errped,xmax_errped);
  
  //----------------------------------------
  //Column Level
  TH1F* TOT_GainPerCol = new TH1F("TOT_GainPerCol","TOT_Gain_ColumnLevel;Gain;# columns",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_PedestalPerCol = new TH1F("TOT_PedestalPerCol","TOT_Pedestal_ColumnLevel;Gain;# columns",nbin_ped,xmin_ped,xmax_ped);
  
  TH1F* TOT_GainPerColBPix = new TH1F("TOT_GainPerColBPix","TOT_GainBPix_ColumnLevel;Gain;# columns",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGainPerColBPix = new TH1F("TOT_ErrorGainPerColBPix","TOT_ErrorGainBPix_ColumnLevel;Error on 1/Gain;# columns",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_PedestalPerColBPix = new TH1F("TOT_PedestalPerColBPix","TOT_PedestalBPix_ColumnLevel;Pedestal;# columns",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestalPerColBPix = new TH1F("TOT_ErrorPedestalPerColBPix","TOT_ErrorPedestalBPix_ColumnLevel;Error on Pedestal;# columns",nbin_errped,xmin_errped,xmax_errped);
  TH1F* TOT_GainPerColFPix = new TH1F("TOT_GainPerColFPix","TOT_GainFPix_ColumnLevel;Gain;# columns",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGainPerColFPix = new TH1F("TOT_ErrorGainPerColFPix","TOT_ErrorGainFPix_ColumnLevel;Error on 1/Gain;# columns",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_PedestalPerColFPix = new TH1F("TOT_PedestalPerColFPix","TOT_PedestalFPix_ColumnLevel;Pedestal;# columns",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestalPerColFPix = new TH1F("TOT_ErrorPedestalPerColFPix","TOT_ErrorPedestalFPix_ColumnLevel;Error on Pedestal;# columns",nbin_errped,xmin_errped,xmax_errped);
  
  //----------------------------------------
  //ROC Level
  TH1F* TOT_GainPerROC = new TH1F("TOT_GainPerROC","TOT_Gain_ROCLevel;Gain;# ROC",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_PedestalPerROC = new TH1F("TOT_PedestalPerROC","TOT_Pedestal_ROCLevel;Gain;# ROC",nbin_ped,xmin_ped,xmax_ped);
  
  TH1F* TOT_GainPerROCBPix = new TH1F("TOT_GainPerROCBPix","TOT_GainBPix_ROCLevel;Gain;# ROC",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGainPerROCBPix = new TH1F("TOT_ErrorGainPerROCBPix","TOT_ErrorGainBPix_ROCLevel;Error on 1/Gain;# ROC",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_PedestalPerROCBPix = new TH1F("TOT_PedestalPerROCBPix","TOT_PedestalBPix_ROCLevel;Pedestal;# ROC",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestalPerROCBPix = new TH1F("TOT_ErrorPedestalPerROCBPix","TOT_ErrorPedestalBPix_ROCLevel;Error on Pedestal;# ROC",nbin_errped,xmin_errped,xmax_errped);
  TH1F* TOT_GainPerROCFPix = new TH1F("TOT_GainPerROCFPix","TOT_GainFPix_ROCLevel;Gain;# ROC",nbin_gain,xmin_gain,xmax_gain);
  TH1F* TOT_ErrorGainPerROCFPix = new TH1F("TOT_ErrorGainPerROCFPix","TOT_ErrorGainFPix_ROCLevel;Error on 1/Gain;# ROC",nbin_errgain,xmin_errgain,xmax_errgain);
  TH1F* TOT_PedestalPerROCFPix = new TH1F("TOT_PedestalPerROCFPix","TOT_PedestalFPix_ROCLevel;Pedestal;# ROC",nbin_ped,xmin_ped,xmax_ped);
  TH1F* TOT_ErrorPedestalPerROCFPix = new TH1F("TOT_ErrorPedestalPerROCFPix","TOT_ErrorPedestalFPix_ROCLevel;Error on Pedestal;# ROC",nbin_errped,xmin_errped,xmax_errped);
  
  TH2F* CorrelationGainPed = new TH2F("CorrelationGainPed","CorrelationGainPed;gain;pedestal",nbin_gain,xmin_gain,xmax_gain,nbin_ped,xmin_ped,xmax_ped);
  TH2F* CorrelationGainPedBPix = new TH2F("CorrelationGainPedBPix","CorrelationGainPedBPix;gain;pedestal",nbin_gain,xmin_gain,xmax_gain,nbin_ped,xmin_ped,xmax_ped);
  TH2F* CorrelationGainPedFPix = new TH2F("CorrelationGainPedFPix","CorrelationGainPedFPix;gain;pedestal",nbin_gain,xmin_gain,xmax_gain,nbin_ped,xmin_ped,xmax_ped);


  //******** ERRORS **********
  
  //Pixel Level 
  TH2F* CorrelationError = new TH2F("CorrelationError","CorrelationError;Error on 1/gain;Error on pedestal" , nbin_errgain , xmin_errgain , xmax_errgain , nbin_errped , xmin_errped , xmax_errped);
  TH2F* ErrorVsGain = new TH2F("ErrorVsGain","ErrorVsGain;Gain;Error on 1/gain" , nbin_gain , xmin_gain , xmax_gain , nbin_errgain , xmin_errgain , xmax_errgain);
  TH2F* ErrorVsPedestal = new TH2F("ErrorVsPedestal","ErrorVsPedestal;Pedestal;Error on pedestal", nbin_ped , xmin_ped , xmax_ped , nbin_errped , xmin_errped , xmax_errped);
  TH2F* CorrelationErrorBPix = new TH2F("CorrelationErrorBPix","CorrelationErrorBPix;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainBPix = new TH2F("ErrorVsGainBPix","ErrorVsGainBPix;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalBPix = new TH2F("ErrorVsPedestalBPix","ErrorVsPedestalBPix;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  TH2F* CorrelationErrorFPix = new TH2F("CorrelationErrorFPix","CorrelationErrorFPix;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainFPix = new TH2F("ErrorVsGainFPix","ErrorVsGainFPix;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalFPix = new TH2F("ErrorVsPedestalFPix","ErrorVsPedestalFPix;Pedestal;Error on pedestal",350,-100,250,200,0,100);
 
  //Column Level
  TH2F* CorrelationGainPedPerCol = new TH2F("CorrelationGainPedPerCol","CorrelationGainPedPerCol;gain;pedestal",400,0,10,350,-100,250);
  TH2F* CorrelationErrorPerCol = new TH2F("CorrelationErrorPerCol","CorrelationErrorPerColPerCol;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainPerCol = new TH2F("ErrorVsGainPerCol","ErrorVsGainPerCol;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalPerCol = new TH2F("ErrorVsPedestalPerCol","ErrorVsPedestalPerCol;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  TH2F* CorrelationGainPedPerColBPix = new TH2F("CorrelationGainPedPerColBPix","CorrelationGainPedPerColBPix;gain;pedestal",400,0,10,350,-100,250);
  TH2F* CorrelationErrorPerColBPix = new TH2F("CorrelationErrorPerColBPix","CorrelationErrorPerColBPix;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainPerColBPix = new TH2F("ErrorVsGainPerColBPix","ErrorVsGainPerColBPix;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalPerColBPix = new TH2F("ErrorVsPedestalPerColBPix","ErrorVsPedestalPerColBPix;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  TH2F* CorrelationGainPedPerColFPix = new TH2F("CorrelationGainPedPerColFPix","CorrelationGainPedPerColFPix;gain;pedestal",400,0,10,350,-100,250);
  TH2F* CorrelationErrorPerColFPix = new TH2F("CorrelationErrorPerColFPix","CorrelationErrorPerColFPix;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainPerColFPix = new TH2F("ErrorVsGainPerColFPix","ErrorVsGainPerColFPix;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalPerColFPix = new TH2F("ErrorVsPedestalPerColFPix","ErrorVsPedestalPerColFPix;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  
  //ROC Level
  TH2F* CorrelationGainPedPerROC = new TH2F("CorrelationGainPedPerROC","CorrelationGainPedPerROC;gain;pedestal",400,0,10,350,-100,250);
  TH2F* CorrelationErrorPerROC = new TH2F("CorrelationErrorPerROC","CorrelationErrorPerROC;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainPerROC = new TH2F("ErrorVsGainPerROC","ErrorVsGainPerROC;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalPerROC = new TH2F("ErrorVsPedestalPerROC","ErrorVsPedestalPerROC;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  TH2F* CorrelationGainPedPerROCBPix = new TH2F("CorrelationGainPedPerROCBPix","CorrelationGainPedPerROCBPix;gain;pedestal",400,0,10,350,-100,250);
  TH2F* CorrelationErrorPerROCBPix = new TH2F("CorrelationErrorPerROCBPix","CorrelationErrorPerROCBPix;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainPerROCBPix = new TH2F("ErrorVsGainPerROCBPix","ErrorVsGainPerROCBPix;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalPerROCBPix = new TH2F("ErrorVsPedestalPerROCBPix","ErrorVsPedestalPerROCBPix;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  TH2F* CorrelationGainPedPerROCFPix = new TH2F("CorrelationGainPedPerROCFPix","CorrelationGainPedPerROCFPix;gain;pedestal",400,0,10,350,-100,250);
  TH2F* CorrelationErrorPerROCFPix = new TH2F("CorrelationErrorPerROCFPix","CorrelationErrorPerROCFPix;Error on 1/gain;Error on pedestal",100,0,1,200,0,100);
  TH2F* ErrorVsGainPerROCFPix = new TH2F("ErrorVsGainPerROCFPix","ErrorVsGainPerROCFPix;Gain;Error on 1/gain",400,0,10,100,0,1);
  TH2F* ErrorVsPedestalPerROCFPix = new TH2F("ErrorVsPedestalPerROCFPix","ErrorVsPedestalPerROCFPix;Pedestal;Error on pedestal",350,-100,250,200,0,100);
  
  
  //Counting overflows for mean/RMS
  TOT_GainBPix->StatOverflows(kTRUE);
  TOT_ErrorGainBPix->StatOverflows(kTRUE);
  TOT_PedestalBPix->StatOverflows(kTRUE);
  TOT_ErrorPedestalBPix->StatOverflows(kTRUE);
  TOT_GainFPix->StatOverflows(kTRUE);
  TOT_ErrorGainFPix->StatOverflows(kTRUE);
  TOT_PedestalFPix->StatOverflows(kTRUE);
  TOT_ErrorPedestalFPix->StatOverflows(kTRUE);
  
  
   
  //TH2F *temp;
  //std::vector< TH2F* > *histo2save;
  
  //TH1F* test = new TH1F("tt","tt",100,0,100);
  
  int NModules = 0;
  int Npix = 0;
  int NBpix = 0;
  int NFpix = 0;
  int NpixGoodFit = 0;
  int NBpixGoodFit = 0;
  int NFpixGoodFit = 0;
  
  //********* Loop on Modules *******
  //cout<<"starting loop on good dirs"<<endl;
  for(unsigned int i=0; i<dirlist.size() ; ++i)
  { 
    NModules++;
   //if(NModules==647  || NModules==648) continue;
    
   //cout<<NModules<<"  "<<dirlist[i]<<endl;
   
     //if(i>10) break;
    
    //****************************************************************************************
    dir = file->GetDirectory(dirlist[i]); 
    list= dir->GetListOfKeys();
    for(int ikey=0;ikey<list->GetEntries();ikey++){
      TKey *thekey = (TKey*)list->At(ikey);
      if(thekey==0) continue;
      TString keyname=thekey->GetName();
      keyname.ReplaceAll(" ","");
      TString keytype=thekey->GetClassName();
      keytype.ReplaceAll(" ","");
      
      //Getting histo from file
      if(keytype=="TH2F"){
	if(keyname.Contains("Gain2d")) Gain2d                         = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("ErrorGain2d")) ErrorGain2d               = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("Pedestal2d")) Pedestal2d                 = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("ErrorPedestal2d")) ErrorPedestal2d       = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainSaturate2d")) GainSaturate2d         = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainDynamicRange2d")) GainDynamicRange2d = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainFitResult2d")) GainFitResult2d       = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainChi2NDF2d")) GainChi2NDF2d           = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainChi2Prob2d")) GainChi2Prob2d         = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainHighPoint2d")) GainHighPoint2d       = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainLowPoint2d")) GainLowPoint2d         = (TH2F*) dir->Get(keyname);

      
	/*if(keyname.Contains("344081672")){
	  temp = new TH2F();
	  temp = (TH2F*) dir->Get(keyname);
	  histo2save->push_back(temp);
	}*/
      }
      
      if(keytype=="TH1F"){
	if(keyname.Contains("GainNPoints1d")) GainNPoints1d           = (TH1F*) dir->Get(keyname);
      }
    }
    
    //****************************************************************************************
    TString dirname = dirlist[i];
    
    if(diff){
      dirname.ReplaceAll("Run "+run,"Run "+run2);
      dir = file2->GetDirectory(dirname);
      if(dir==0){
        cout << dirname << " was not found ..." << endl;
        continue;
      }
      list= dir->GetListOfKeys();
      for(int ikey=0;ikey<list->GetEntries();ikey++){
        TKey *thekey = (TKey*)list->At(ikey);
        if(thekey==0) continue;
        TString keyname=thekey->GetName();
        keyname.ReplaceAll(" ","");
        TString keytype=thekey->GetClassName();
        keytype.ReplaceAll(" ","");
        
        //Getting histo from file
        if(keytype=="TH2F"){
	if(keyname.Contains("Gain2d")) 			Gain2d_2             = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("ErrorGain2d")) 		ErrorGain2d_2        = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("Pedestal2d")) 		Pedestal2d_2         = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("ErrorPedestal2d"))		ErrorPedestal2d_2    = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainSaturate2d")) 		GainSaturate2d_2     = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainDynamicRange2d")) 	GainDynamicRange2d_2 = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainFitResult2d")) 	GainFitResult2d_2    = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainChi2NDF2d")) 		GainChi2NDF2d_2      = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainChi2Prob2d")) 		GainChi2Prob2d_2     = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainHighPoint2d")) 	GainHighPoint2d_2    = (TH2F*) dir->Get(keyname);
	if(keyname.Contains("GainLowPoint2d")) 		GainLowPoint2d_2     = (TH2F*) dir->Get(keyname);

        
	/*if(keyname.Contains("344081672")){
	  temp = new TH2F();
	  temp = (TH2F*) dir->Get(keyname);
	  histo2save->push_back(temp);
	}*/
        }
        
        if(keytype=="TH1F"){
	if(keyname.Contains("GainNPoints1d")) 		GainNPoints1d_2      = (TH1F*) dir->Get(keyname);
        }
      }
    }
    
    
    int NX = Gain2d->GetNbinsX();
    int NY = Gain2d->GetNbinsY();
    
    int detid = getDetId(Gain2d->GetName());
    cout << "Module: " << detid << endl;
    
    if(dirname.Contains("Barrel")){
      if(detid>bmax) bmax=detid;
      if(detid<bmin) bmin=detid;
    }
    if(dirname.Contains("Endcap")){
      if(detid>fmax) fmax=detid;
      if(detid<fmin) fmin=detid;
    }
    
    /*int layer = 0, ladder = 0, module = 0;
    pix.detID2Bpix(detid , layer , ladder , module);
    cout << "    " << layer << "  " << ladder << "  " << module << endl;
    getBarrelPos(dirlist[i] , layer , ladder , module);
    cout << "    " << layer << "  " << ladder << "  " << module << endl;
    */
    
    pix_val pix;
    
    
   /* if(NModules==1 || NModules==100 || NModules==1000){
      output->cd();
      write(Gain2d);
      write(ErrorGain2d);
      write(Pedestal2d);
      write(ErrorPedestal2d);
      write(GainSaturate2d);
      write(GainDynamicRange2d);
      write(GainFitResult2d);
      write(GainChi2NDF2d);
      write(GainChi2Prob2d);
      file->cd();
    }
    */
    TMean* MGain2d             = new TMean();
    TMean* MErrorGain2d        = new TMean();
    TMean* MPedestal2d         = new TMean();
    TMean* MErrorPedestal2d    = new TMean();
    TMean* MGainSaturate2d     = new TMean();
    TMean* MGainDynamicRange2d = new TMean();
    TMean* MGainFitResult2d    = new TMean();
    TMean* MGainChi2NDF2d      = new TMean();
    TMean* MGainChi2Prob2d     = new TMean();
    TMean* MGainHighPoint2d    = new TMean();
    TMean* MGainLowPoint2d     = new TMean();
    
    TMean* M_frac_gain_high    = new TMean();    
    
    double gainsperROC[16]={0};
    double errorgainsperROC[16]={0};
    double pedestalsperROC[16]={0};
    double errorpedestalsperROC[16]={0};
    int nperROC[16]={0};
    
    //******** Looping over pixels ********
    for(int xpix=1;xpix<=NX;xpix++)
    {
      double gainspercol[2]={0};
      double errorgainspercol[2]={0};
      double pedestalspercol[2]={0};
      double errorpedestalspercol[2]={0};
      double npercol[2]={0};
      int changecol=0;
      
      for(int ypix=1;ypix<=NY;ypix++){
      
	      /*if(i==0 && Gain2d->GetBinContent(xpix,ypix)==0){
	        test->Fill(Gain2d->GetBinContent(xpix,ypix));
	        cout<<xpix<<"  "<<ypix<<"  "<<Gain2d->GetBinContent(xpix,ypix)<<endl;
	      }*/
	      
	 if(ypix%80==0) changecol=1;
	   
	 Npix++;
	 if(dirlist[i].Contains("Barrel")) NBpix++;
	 if (dirlist[i].Contains("Endcap")) NFpix++;
	 	
	 pix.gain	   = Gain2d->GetBinContent(xpix,ypix)		  ;
	 pix.err_gain	   = ErrorGain2d->GetBinContent(xpix,ypix)	  ; 
	 pix.pedestal	   = Pedestal2d->GetBinContent(xpix,ypix)	  ; 
	 pix.err_pedestal  = ErrorPedestal2d->GetBinContent(xpix,ypix)	  ;
	 pix.plateau	   = GainSaturate2d->GetBinContent(xpix,ypix)	  ;
	 pix.dynamic_range = GainDynamicRange2d->GetBinContent(xpix,ypix) ;
	 pix.fit_result    = GainFitResult2d->GetBinContent(xpix,ypix)	  ;
	 pix.chi2NDF	   = GainChi2NDF2d->GetBinContent(xpix,ypix)	  ;
	 pix.chi2_prob	   = GainChi2Prob2d->GetBinContent(xpix,ypix)	  ;
	 pix.high_point    = GainHighPoint2d->GetBinContent(xpix,ypix)	  ;
	 pix.low_point	   = GainLowPoint2d->GetBinContent(xpix,ypix)	  ;
	 
	 if(diff){
	   pix.gain	     -= Gain2d_2->GetBinContent(xpix,ypix);
	   pix.err_gain	     -= ErrorGain2d_2->GetBinContent(xpix,ypix);
	   pix.pedestal	     -= Pedestal2d_2->GetBinContent(xpix,ypix);
	   pix.err_pedestal  -= ErrorPedestal2d_2->GetBinContent(xpix,ypix);
	   pix.plateau       -= GainSaturate2d_2->GetBinContent(xpix,ypix);
	   pix.dynamic_range -= GainDynamicRange2d_2->GetBinContent(xpix,ypix);
	   pix.fit_result    -= GainFitResult2d_2->GetBinContent(xpix,ypix);
	   pix.chi2NDF       -= GainChi2NDF2d_2->GetBinContent(xpix,ypix);
	   pix.chi2_prob     -= GainChi2Prob2d_2->GetBinContent(xpix,ypix);
	   pix.high_point    -= GainHighPoint2d_2->GetBinContent(xpix,ypix);
	   pix.low_point     -= GainLowPoint2d_2->GetBinContent(xpix,ypix);
         }
	      
	 if(GainFitResult2d->GetBinContent(xpix,ypix)>0){
	   
	   if(diff && GainFitResult2d_2->GetBinContent(xpix,ypix)<=0) continue;
	   
	   
           if( fabs(pix.gain) > frac_flag_high_gain ) M_frac_gain_high->Add(1);
	   else                                       M_frac_gain_high->Add(0);
	   
	   if(TMath::Abs(pix.gain)>30)
	     cout << "      col,row " << xpix << "," << ypix <<
	     " gain: "<< pix.gain << " ped: "<< pix.pedestal << endl;
	   if(fabs(pix.gain) > 30)      continue;
	   if(fabs(pix.pedestal) > 200) continue;
	   
	   NpixGoodFit++;
	   if(dirlist[i].Contains("Barrel")) NBpixGoodFit++;
	   if (dirlist[i].Contains("Endcap")) NFpixGoodFit++;
	   	   	   
	   npercol[changecol]++;
	   gainspercol[changecol]		+=pix.gain;
	   errorgainspercol[changecol]		+=pix.err_gain;
	   pedestalspercol[changecol]		+=pix.pedestal;
	   errorpedestalspercol[changecol]	+=pix.err_pedestal;
	   
	   nperROC[rocId(xpix-1,ypix-1)]++;
	   gainsperROC[rocId(xpix-1,ypix-1)]		+=pix.gain;
	   errorgainsperROC[rocId(xpix-1,ypix-1)]	+=pix.err_gain;
	   pedestalsperROC[rocId(xpix-1,ypix-1)]	+=pix.pedestal;
	   errorpedestalsperROC[rocId(xpix-1,ypix-1)]	+=pix.err_pedestal;
	   
	   MGain2d		->Add(pix.gain);
	   MErrorGain2d		->Add(pix.err_gain);
	   MPedestal2d		->Add(pix.pedestal);
	   MErrorPedestal2d	->Add(pix.err_pedestal);
	   MGainSaturate2d	->Add(pix.plateau);
	   MGainDynamicRange2d	->Add(pix.dynamic_range);
	   MGainFitResult2d	->Add(pix.fit_result);
	   MGainChi2NDF2d	->Add(pix.chi2NDF);
	   MGainChi2Prob2d	->Add(pix.chi2_prob);
	   MGainHighPoint2d	->Add(pix.high_point);
	   MGainLowPoint2d	->Add(pix.low_point);


	   
	   TOT_Gain->Fill(pix.gain);
	   TOT_ErrorGain->Fill(pix.err_gain);
	   TOT_Pedestal->Fill(pix.pedestal);
	   TOT_ErrorPedestal->Fill(pix.err_pedestal);
	   
	   CorrelationError->Fill(pix.err_gain , pix.err_pedestal);
	   ErrorVsGain->Fill(pix.gain , pix.err_gain);
	   ErrorVsPedestal->Fill(pix.pedestal , pix.err_pedestal);
	   CorrelationGainPed->Fill(pix.gain , pix.pedestal);
	   
	   if(dirlist[i].Contains("Barrel")){
	     TOT_GainBPix->Fill(pix.gain);
	     TOT_ErrorGainBPix->Fill(pix.err_gain);
	     TOT_PedestalBPix->Fill(pix.pedestal);
	     TOT_ErrorPedestalBPix->Fill(pix.err_pedestal);
	     
	     
	     CorrelationErrorBPix->Fill(pix.err_gain , pix.err_pedestal);
	     ErrorVsGainBPix->Fill(pix.gain , pix.err_gain);
	     ErrorVsPedestalBPix->Fill(pix.pedestal , pix.err_pedestal);
	     CorrelationGainPedBPix->Fill(pix.gain , pix.pedestal);
	     
	   }
	   else if (dirlist[i].Contains("Endcap")) {
	     TOT_GainFPix->Fill(pix.gain);
	     TOT_ErrorGainFPix->Fill(pix.err_gain);
	     TOT_PedestalFPix->Fill(pix.pedestal);
	     TOT_ErrorPedestalFPix->Fill(pix.err_pedestal);
	     
	     CorrelationErrorFPix->Fill(pix.err_gain , pix.err_pedestal);
	     ErrorVsGainFPix->Fill(pix.gain , pix.err_gain);
	     ErrorVsPedestalFPix->Fill(pix.pedestal , pix.err_pedestal);
	     CorrelationGainPedFPix->Fill(pix.gain , pix.pedestal);
	     
	   }
         }
       }//end of loop over Y pixels
       
       //Looping over column
       for(int nCOL=0;nCOL<2;nCOL++){
         if(npercol[nCOL]!=0){
	   
	   gainspercol[nCOL]/=double(npercol[nCOL]);
	   pedestalspercol[nCOL]/=double(npercol[nCOL]);
	   errorgainspercol[nCOL]/=double(npercol[nCOL]);
	   errorpedestalspercol[nCOL]/=double(npercol[nCOL]);
	   
	   
           TOT_GainPerCol->Fill(gainspercol[nCOL]);
           TOT_PedestalPerCol->Fill(pedestalspercol[nCOL]);
	   
	   CorrelationGainPedPerCol->Fill(gainspercol[nCOL],pedestalspercol[nCOL]);
	   CorrelationErrorPerCol->Fill(errorgainspercol[nCOL],errorpedestalspercol[nCOL]);
	   ErrorVsGainPerCol->Fill(gainspercol[nCOL],errorgainspercol[nCOL]);
	   ErrorVsPedestalPerCol->Fill(pedestalspercol[nCOL],errorpedestalspercol[nCOL]);
	   
	   
	   if(dirlist[i].Contains("Barrel")){
             TOT_GainPerColBPix->Fill(gainspercol[nCOL]);
             TOT_PedestalPerColBPix->Fill(pedestalspercol[nCOL]);
	     TOT_ErrorGainPerColBPix->Fill(errorgainspercol[nCOL]);
             TOT_ErrorPedestalPerColBPix->Fill(errorpedestalspercol[nCOL]);
	     
	     CorrelationGainPedPerColBPix->Fill(gainspercol[nCOL],pedestalspercol[nCOL]);
	     CorrelationErrorPerColBPix->Fill(errorgainspercol[nCOL],errorpedestalspercol[nCOL]);
	     ErrorVsGainPerColBPix->Fill(gainspercol[nCOL],errorgainspercol[nCOL]);
	     ErrorVsPedestalPerColBPix->Fill(pedestalspercol[nCOL],errorpedestalspercol[nCOL]);
	   }
	   if(dirlist[i].Contains("Endcap")){
             TOT_GainPerColFPix->Fill(gainspercol[nCOL]);
             TOT_PedestalPerColFPix->Fill(pedestalspercol[nCOL]);
	     TOT_ErrorGainPerColFPix->Fill(errorgainspercol[nCOL]);
             TOT_ErrorPedestalPerColFPix->Fill(errorpedestalspercol[nCOL]);
	     
	     CorrelationGainPedPerColFPix->Fill(gainspercol[nCOL],pedestalspercol[nCOL]);
	     CorrelationErrorPerColFPix->Fill(errorgainspercol[nCOL],errorpedestalspercol[nCOL]);
	     ErrorVsGainPerColFPix->Fill(gainspercol[nCOL],errorgainspercol[nCOL]);
	     ErrorVsPedestalPerColFPix->Fill(pedestalspercol[nCOL],errorpedestalspercol[nCOL]);
	   }
         }
       }//end of loop over col
       
     }//end of loop over X for pixels
     
     //Looping over ROC
     for(int nROC=0;nROC<16;nROC++){
       if(nperROC[nROC]!=0){
	   
	 gainsperROC[nROC]/=double(nperROC[nROC]);
	 pedestalsperROC[nROC]/=double(nperROC[nROC]);
	 errorgainsperROC[nROC]/=double(nperROC[nROC]);
	 errorpedestalsperROC[nROC]/=double(nperROC[nROC]);
	   
	   
         TOT_GainPerROC->Fill(gainsperROC[nROC]);
         TOT_PedestalPerROC->Fill(pedestalsperROC[nROC]);
	   
	 CorrelationGainPedPerROC->Fill(gainsperROC[nROC],pedestalsperROC[nROC]);
	 CorrelationErrorPerROC->Fill(errorgainsperROC[nROC],errorpedestalsperROC[nROC]);
	 ErrorVsGainPerROC->Fill(gainsperROC[nROC],errorgainsperROC[nROC]);
	 ErrorVsPedestalPerROC->Fill(pedestalsperROC[nROC],errorpedestalsperROC[nROC]);
	 	   
	 if(dirlist[i].Contains("Barrel")){
           TOT_GainPerROCBPix->Fill(gainsperROC[nROC]);
           TOT_PedestalPerROCBPix->Fill(pedestalsperROC[nROC]);
	   TOT_ErrorGainPerROCBPix->Fill(errorgainsperROC[nROC]);
           TOT_ErrorPedestalPerROCBPix->Fill(errorpedestalsperROC[nROC]);
	     
	   CorrelationGainPedPerROCBPix->Fill(gainsperROC[nROC],pedestalsperROC[nROC]);
	   CorrelationErrorPerROCBPix->Fill(errorgainsperROC[nROC],errorpedestalsperROC[nROC]);
	   ErrorVsGainPerROCBPix->Fill(gainsperROC[nROC],errorgainsperROC[nROC]);
	   ErrorVsPedestalPerROCBPix->Fill(pedestalsperROC[nROC],errorpedestalsperROC[nROC]);
	 }
	 if(dirlist[i].Contains("Endcap")){
           TOT_GainPerROCFPix->Fill(gainsperROC[nROC]);
           TOT_PedestalPerROCFPix->Fill(pedestalsperROC[nROC]);
	   TOT_ErrorGainPerROCFPix->Fill(errorgainsperROC[nROC]);
           TOT_ErrorPedestalPerROCFPix->Fill(errorpedestalsperROC[nROC]);
	     
	   CorrelationGainPedPerROCFPix->Fill(gainsperROC[nROC],pedestalsperROC[nROC]);
	   CorrelationErrorPerROCFPix->Fill(errorgainsperROC[nROC],errorpedestalsperROC[nROC]);
	   ErrorVsGainPerROCFPix->Fill(gainsperROC[nROC],errorgainsperROC[nROC]);
	   ErrorVsPedestalPerROCFPix->Fill(pedestalsperROC[nROC],errorpedestalsperROC[nROC]);
	 }
       }
     }
     
     SUM_Gain2d->SetBinContent(NModules,MGain2d->GetMean());
     SUM_ErrorGain2d->SetBinContent(NModules,MErrorGain2d->GetMean());
     SUM_Pedestal2d->SetBinContent(NModules,MPedestal2d->GetMean());
     SUM_ErrorPedestal2d->SetBinContent(NModules,MErrorPedestal2d->GetMean());
     SUM_GainSaturate2d->SetBinContent(NModules,MGainSaturate2d->GetMean());
     SUM_GainDynamicRange2d->SetBinContent(NModules,MGainDynamicRange2d->GetMean());
     SUM_GainFitResult2d->SetBinContent(NModules,MGainFitResult2d->GetMean());
     SUM_GainChi2NDF2d->SetBinContent(NModules,MGainChi2NDF2d->GetMean());
     SUM_GainChi2Prob2d->SetBinContent(NModules,MGainChi2Prob2d->GetMean());
     SUM_GainHighPoint2d->SetBinContent(NModules,MGainHighPoint2d->GetMean());
     SUM_GainLowPoint2d->SetBinContent(NModules,MGainLowPoint2d->GetMean());
     if(diff == 0) SUM_GainNPoints1d->SetBinContent(NModules,GainNPoints1d->GetMean());
     else          SUM_GainNPoints1d->SetBinContent(NModules,GainNPoints1d->GetMean() - GainNPoints1d_2->GetMean());
          
     MAP_Gain			.fill(detid , dirlist[i] , MGain2d->GetMean());    
     MAP_ErrorGain		.fill(detid , dirlist[i] , MErrorGain2d->GetMean());		    
     MAP_Pedestal		.fill(detid , dirlist[i] , MPedestal2d->GetMean());		     
     MAP_ErrorPedestal		.fill(detid , dirlist[i] , MErrorPedestal2d->GetMean());   
     MAP_GainSaturate		.fill(detid , dirlist[i] , MGainSaturate2d->GetMean());      
     MAP_GainDynamicRange	.fill(detid , dirlist[i] , MGainDynamicRange2d->GetMean());
     MAP_GainFitResult		.fill(detid , dirlist[i] , MGainFitResult2d->GetMean());    
     MAP_GainChi2NDF		.fill(detid , dirlist[i] , MGainChi2NDF2d->GetMean());     
     MAP_GainChi2Prob		.fill(detid , dirlist[i] , MGainChi2Prob2d->GetMean());      
     MAP_GainHighPoint		.fill(detid , dirlist[i] , MGainHighPoint2d->GetMean());    
     MAP_GainLowPoint		.fill(detid , dirlist[i] , MGainLowPoint2d->GetMean());
     if(diff==0) MAP_GainNPoints.fill(detid , dirlist[i] , GainNPoints1d->GetMean());
     else        MAP_GainNPoints.fill(detid , dirlist[i] , GainNPoints1d->GetMean() - GainNPoints1d_2->GetMean());	  
   
     MAP_frac_gain_high         .fill(detid , dirlist[i] , M_frac_gain_high->GetMean());
   
     if(TMath::Abs(MGain2d->GetMean())>20) cout<<"************  "<<MGain2d->GetMean()<<" mod "<<NModules<<"  "<<dirlist[i]<<endl;
     if(MPedestal2d->GetMean()==0) cout<<"=============="<<" mod "<<NModules<<dirlist[i]<<endl;
     
     cout << "   gain: " << MGain2d->GetMean() << "    pedestal: " << MPedestal2d->GetMean() << endl;
     
     delete MGain2d; 	    
     delete MErrorGain2d;	 
     delete MPedestal2d;	 
     delete MErrorPedestal2d;   
     delete MGainSaturate2d;     
     delete MGainDynamicRange2d;
     delete MGainFitResult2d;    
     delete MGainChi2NDF2d;     
     delete MGainChi2Prob2d;     
     delete MGainHighPoint2d;    
     delete MGainLowPoint2d;
     delete M_frac_gain_high;
	      
   }//End of Loop over modules
   
   
  
  
  //********* WRITING IN FILE ***************
  cout << "Writing to file ..." << endl;
  
  TFile* output = new TFile("Comp"+RunNumber+".root","RECREATE");
  //output->cd();
  
 
  output->cd();
  
  //test->Write();
  
  Bool_t   log = 0;
  if(diff) log = 1;
  
  //Pixel Level
  write(TOT_Gain , log);
  write(TOT_ErrorGain , log);
  write(TOT_Pedestal , log);
  write(TOT_ErrorPedestal , log);
  write(TOT_GainBPix , log);
  write(TOT_ErrorGainBPix , log);
  write(TOT_PedestalBPix , log);
  write(TOT_ErrorPedestalBPix , log);
  write(TOT_GainFPix , log);
  write(TOT_ErrorGainFPix , log);
  write(TOT_PedestalFPix , log);
  write(TOT_ErrorPedestalFPix , log);
  
  //Column Level
  write(TOT_GainPerCol , log);
  write(TOT_PedestalPerCol , log);
  write(TOT_GainPerColBPix , log);
  write(TOT_ErrorGainPerColBPix , log);
  write(TOT_PedestalPerColBPix , log);
  write(TOT_ErrorPedestalPerColBPix , log);
  write(TOT_GainPerColFPix , log);
  write(TOT_ErrorGainPerColFPix , log);
  write(TOT_PedestalPerColFPix , log);
  write(TOT_ErrorPedestalPerColFPix , log);
  
  //ROC Level
  write(TOT_GainPerROC , log);
  write(TOT_PedestalPerROC , log);
  write(TOT_GainPerROCBPix , log);
  write(TOT_ErrorGainPerROCBPix , log);
  write(TOT_PedestalPerROCBPix , log);
  write(TOT_ErrorPedestalPerROCBPix , log);
  write(TOT_GainPerROCFPix , log);
  write(TOT_ErrorGainPerROCFPix , log);
  write(TOT_PedestalPerROCFPix , log);
  write(TOT_ErrorPedestalPerROCFPix , log);
  
  //All levels
  writeAllLevels(TOT_GainBPix,TOT_GainPerColBPix,TOT_GainPerROCBPix,"TOT_GainAllLevelBPix");
  writeAllLevels(TOT_ErrorGainBPix,TOT_ErrorGainPerColBPix,TOT_ErrorGainPerROCBPix,"TOT_ErrorGainAllLevelBPix");
  writeAllLevels(TOT_PedestalBPix,TOT_PedestalPerColBPix,TOT_PedestalPerROCBPix,"TOT_PedestalAllLevelBPix");
  writeAllLevels(TOT_ErrorPedestalBPix,TOT_ErrorPedestalPerColBPix,TOT_ErrorPedestalPerROCBPix,"TOT_ErrorPedestalAllLevelBPix");
  writeAllLevels(TOT_GainFPix,TOT_GainPerColFPix,TOT_GainPerROCFPix,"TOT_GainAllLevelFPix");
  writeAllLevels(TOT_ErrorGainFPix,TOT_ErrorGainPerColFPix,TOT_ErrorGainPerROCFPix,"TOT_ErrorGainAllLevelFPix");
  writeAllLevels(TOT_PedestalFPix,TOT_PedestalPerColFPix,TOT_PedestalPerROCFPix,"TOT_PedestalAllLevelFPix");
  writeAllLevels(TOT_ErrorPedestalFPix,TOT_ErrorPedestalPerColFPix,TOT_ErrorPedestalPerROCFPix,"TOT_ErrorPedestalAllLevelFPix");
    
  //Pixel level
  write(CorrelationGainPed);
  write(CorrelationGainPedBPix);
  write(CorrelationGainPedFPix);
  write(CorrelationError);
  write(ErrorVsGain);
  write(ErrorVsPedestal);
  write(CorrelationErrorBPix);
  write(ErrorVsGainBPix);
  write(ErrorVsPedestalBPix);
  write(CorrelationErrorFPix);
  write(ErrorVsGainFPix);
  write(ErrorVsPedestalFPix);
  	
  //Column level
  write(CorrelationGainPedPerCol);
  write(CorrelationGainPedPerColBPix);
  write(CorrelationGainPedPerColFPix);
  write(CorrelationErrorPerCol);
  write(ErrorVsGainPerCol);
  write(ErrorVsPedestalPerCol);
  write(CorrelationErrorPerColBPix);
  write(ErrorVsGainPerColBPix);
  write(ErrorVsPedestalPerColBPix);
  write(CorrelationErrorPerColFPix);
  write(ErrorVsGainPerColFPix);
  write(ErrorVsPedestalPerColFPix);
  	
  //ROC level
  write(CorrelationGainPedPerROC);
  write(CorrelationGainPedPerROCBPix);
  write(CorrelationGainPedPerROCFPix);
  write(CorrelationErrorPerROC);
  write(ErrorVsGainPerROC);
  write(ErrorVsPedestalPerROC);
  write(CorrelationErrorPerROCBPix);
  write(ErrorVsGainPerROCBPix);
  write(ErrorVsPedestalPerROCBPix);
  write(CorrelationErrorPerROCFPix);
  write(ErrorVsGainPerROCFPix);
  write(ErrorVsPedestalPerROCFPix);
  
  //Summary plots
  gStyle->SetOptStat(0);
  write(SUM_Gain2d);
  write(SUM_ErrorGain2d);
  write(SUM_Pedestal2d);
  write(SUM_ErrorPedestal2d);
  write(SUM_GainSaturate2d);
  write(SUM_GainDynamicRange2d);
  write(SUM_GainFitResult2d);
  write(SUM_GainChi2NDF2d);
  write(SUM_GainChi2Prob2d);
  write(SUM_GainHighPoint2d);
  write(SUM_GainLowPoint2d);
  write(SUM_GainNPoints1d);
  
  MAP_Gain.Write();
  MAP_ErrorGain.Write();
  MAP_Pedestal.Write();
  MAP_ErrorPedestal.Write();
  MAP_GainSaturate.Write();
  MAP_GainDynamicRange.Write();
  MAP_GainFitResult.Write();
  MAP_GainChi2NDF.Write();
  MAP_GainChi2Prob.Write();
  MAP_GainHighPoint.Write();
  MAP_GainLowPoint.Write();
  MAP_GainNPoints.Write();
  
  MAP_frac_gain_high.Write();
  
  //for(int i=0;i<histo2save->size();i++){;}
    //(histo2save[i]).Write();
    
    
    

  //********* SUMMARY WRITING ***************
  cout<<"Looked @ "<<NModules<<" modules"<<endl;
  cout<<"Total number of pixels with good fit : "<<NpixGoodFit<<" over "<<Npix<<endl;
  cout<<"In BPIX : "<<NBpixGoodFit<<" over "<<NBpix<<endl;
  cout<<"In FPIX : "<<NFpixGoodFit<<" over "<<NFpix<<endl<<endl;
  cout<<"Gain Mean : "<<endl;
  cout<<"BPIX : "<<TOT_GainBPix->GetMean()<<endl;
  cout<<"FPIX : "<<TOT_GainFPix->GetMean()<<endl;
  cout<<"Gain Error Mean : "<<endl;
  cout<<"BPIX : "<<TOT_ErrorGainBPix->GetMean()<<endl;
  cout<<"FPIX : "<<TOT_ErrorGainFPix->GetMean()<<endl<<endl;
  cout<<"Pedestal Mean : "<<endl;
  cout<<"BPIX : "<<TOT_PedestalBPix->GetMean()<<endl;
  cout<<"FPIX : "<<TOT_PedestalFPix->GetMean()<<endl;
  cout<<"Pedestal Error Mean : "<<endl;
  cout<<"BPIX : "<<TOT_ErrorPedestalBPix->GetMean()<<endl;
  cout<<"FPIX : "<<TOT_ErrorPedestalFPix->GetMean()<<endl<<endl;
  
  ofstream summary;
  ofstream texsummary;
  
  if(printSummary)
    summary.open("Summary"+RunNumber+".txt");
  if(printTexSummary)
    texsummary.open("texSummary"+RunNumber+".tex");
    
  if(summary.is_open()){
    cout << "Writing summary ... " << endl;
    summary<<"Looked @ "<<NModules<<" modules"<<endl;
    summary<<"Total number of pixels with good fit : "<<NpixGoodFit<<" over "<<Npix<<" -> "<<(double)NpixGoodFit/(double)Npix<<" %"<<endl;
    summary<<"In BPIX : "<<NBpixGoodFit<<" over "<<NBpix<<" -> "<<(double)NBpixGoodFit/(double)NBpix<<" %"<<endl;
    summary<<"In FPIX : "<<NFpixGoodFit<<" over "<<NFpix<<" -> "<<(double)NFpixGoodFit/(double)NFpix<<" %"<<endl<<endl;
    summary<<"Gain Mean : "<<endl;
    summary<<"BPIX : "<<TOT_GainBPix->GetMean()<<" +- "<<TOT_GainBPix->GetRMS()<<endl;
    summary<<"FPIX : "<<TOT_GainFPix->GetMean()<<" +- "<<TOT_GainFPix->GetRMS()<<endl;
    summary<<"Gain Error Mean : "<<endl;
    summary<<"BPIX : "<<TOT_ErrorGainBPix->GetMean()<<" +- "<<TOT_ErrorGainBPix->GetRMS()<<endl;
    summary<<"FPIX : "<<TOT_ErrorGainFPix->GetMean()<<" +- "<<TOT_ErrorGainFPix->GetRMS()<<endl<<endl;
    summary<<"Pedestal Mean : "<<endl;
    summary<<"BPIX : "<<TOT_PedestalBPix->GetMean()<<" +- "<<TOT_PedestalBPix->GetRMS()<<endl;
    summary<<"FPIX : "<<TOT_PedestalFPix->GetMean()<<" +- "<<TOT_PedestalFPix->GetRMS()<<endl;
    summary<<"Pedestal Error Mean : "<<endl;
    summary<<"BPIX : "<<TOT_ErrorPedestalBPix->GetMean()<<" +- "<<TOT_ErrorPedestalBPix->GetRMS()<<endl;
    summary<<"FPIX : "<<TOT_ErrorPedestalFPix->GetMean()<<" +- "<<TOT_ErrorPedestalFPix->GetRMS()<<endl<<endl;
    summary.close();
  }
  
  if(texsummary.is_open()){
    cout << "Writing tex summary ... " << endl;
    texsummary<<"Looked @ "<<NModules<<" modules"<< " \\\\" << endl;
    texsummary<<"Total number of pixels with good fit : "<<NpixGoodFit<<" over "<<Npix<<" $->$  "<<(double)NpixGoodFit/(double)Npix<<" \\%"<< " \\\\" << endl;
    texsummary<<"In BPIX : "<<NBpixGoodFit<<" over "<<NBpix<<" $->$ "<<(double)NBpixGoodFit/(double)NBpix<<" \\%"<< " \\\\" << endl;
    texsummary<<"In FPIX : "<<NFpixGoodFit<<" over "<<NFpix<<" $->$ "<<(double)NFpixGoodFit/(double)NFpix<<" \\%"<< " \\\\" << endl<< " \\\\" << endl;
    texsummary<<"Gain Mean : "<< " \\\\" << endl;
    texsummary<<"BPIX : "<<TOT_GainBPix->GetMean()<<" $\\pm$ "<<TOT_GainBPix->GetRMS()<< " \\\\" << endl;
    texsummary<<"FPIX : "<<TOT_GainFPix->GetMean()<<" $\\pm$ "<<TOT_GainFPix->GetRMS()<< " \\\\" << endl;
    texsummary<<"Gain Error Mean : "<< " \\\\" << endl;
    texsummary<<"BPIX : "<<TOT_ErrorGainBPix->GetMean()<<" $\\pm$ "<<TOT_ErrorGainBPix->GetRMS()<< " \\\\" << endl;
    texsummary<<"FPIX : "<<TOT_ErrorGainFPix->GetMean()<<" $\\pm$ "<<TOT_ErrorGainFPix->GetRMS()<< " \\\\" << endl<< " \\\\" << endl;
    texsummary<<"Pedestal Mean : "<< " \\\\" << endl;
    texsummary<<"BPIX : "<<TOT_PedestalBPix->GetMean()<<" $\\pm$ "<<TOT_PedestalBPix->GetRMS()<< " \\\\" << endl;
    texsummary<<"FPIX : "<<TOT_PedestalFPix->GetMean()<<" $\\pm$ "<<TOT_PedestalFPix->GetRMS()<< " \\\\" << endl;
    texsummary<<"Pedestal Error Mean : "<< " \\\\" << endl;
    texsummary<<"BPIX : "<<TOT_ErrorPedestalBPix->GetMean()<<" $\\pm$ "<<TOT_ErrorPedestalBPix->GetRMS()<< " \\\\" << endl;
    texsummary<<"FPIX : "<<TOT_ErrorPedestalFPix->GetMean()<<" $\\pm$ "<<TOT_ErrorPedestalFPix->GetRMS()<< " \\\\" << endl<< " \\\\" << endl;
    texsummary.close();
  }
  
  
  cout << setprecision(9) << "Barrel: " << bmin << " to " << bmax << endl;
  cout << "Endcap: " << fmin << " to " << fmax << endl;
  
  
	
  file->Close();
  output->Close(); 
  
}//End of main function


//Get ROC if from a pixel position
int rocId(int col, int row){
 int rocRow = row/80;
 int rocCol = col/52;
 int rocId = rocCol + rocRow*8;
 return rocId;
}

//write TH1F
void write(TH1F* histo , Bool_t logY , Bool_t largePad , TString statOptions){
  TString t=histo->GetName();
  t=t+RunNumber;
  int canv_y = 500 ;
  int canv_x = canv_y;
  if(largePad) canv_x = 2 * canv_y;
  TCanvas* c1 = new TCanvas("c1" , "c1" , canv_x , canv_y);
  c1->cd();
  
  if(logY) gPad->SetLogy(1);
  gStyle->SetOptStat(statOptions);
  
  histo->Draw();
  histo->Write();
  if(PrintImages){
    c1->Print(t+".png","png");
    if(PrintAllFormats){
      c1->Print(t+".eps","eps");
      c1->Print(t+".root","root");
    }
  }
  
  delete c1;
}

//Write TH2F
void write(TH2F* histo , Bool_t logZ){
  TString t=histo->GetName();
  t=t+RunNumber;
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.15);
  if(logZ) gPad->SetLogz(1);
  
  histo->Draw("colz");
  histo->Write();
  if(PrintImages){
    c1->Print(t+".png","png");
    if(PrintAllFormats){
      c1->Print(t+".eps","eps");
      c1->Print(t+".root","root");
    }
  }
  delete c1;
}

//Write pixel, column and ROC level on same canvas
void writeAllLevels(TH1F* pixel,TH1F* column,TH1F* roc,TString name){
  if(pixel->Integral()!=0)  pixel-> Scale(1./pixel->Integral());
  if(column->Integral()!=0) column->Scale(1./column->Integral());
  if(roc->Integral()!=0)    roc->   Scale(1./roc->Integral());
  
  column->SetLineColor(2);
  roc->SetLineColor(4);
  
  double Max=1;
  (pixel->GetMaximum() <= column->GetMaximum() )? Max=column->GetMaximum():Max=pixel->GetMaximum();
  if(Max<roc->GetMaximum()) Max=roc->GetMaximum();
  pixel->SetMaximum(1.1*Max);
  
  pixel->GetYaxis()->SetTitle("a.u.");
  
  TLegend* leg = new TLegend(0.65,0.7,0.85,0.9);
  leg->AddEntry(pixel,"Pixel Level","l");
  leg->AddEntry(column,"Column Level","l");
  leg->AddEntry(roc,"Roc Level","l");
  leg->SetFillColor(10);
  
  gStyle->SetOptStat(0);
  TCanvas* canv = new TCanvas(name,name);
  canv->cd();
  
  gStyle->SetOptStat(0);
  pixel->Draw();
  column->Draw("same");
  roc->Draw("same");
  leg->Draw("same");
  
  canv->Write();
  if(PrintImages){
    canv->Print(name+RunNumber+".png","png");
    if(PrintAllFormats){
      canv->Print(name+RunNumber+".eps","eps");
      canv->Print(name+RunNumber+".root","root");
    }
  }
}

//Get all modules from file
vector<TString> makedirlist(TFile* file,TString comparestring)
{
// make a loop over all plots
  file->cd();
  TList *list = file->GetListOfKeys();
  //Int_t nkeys = file->GetNkeys();
  
  std::vector<TString> dirlist;
  TDirectory *dir = file->GetDirectory("DQMData");
  if(dir==0) return dirlist;

  cout<<"Starting to find Directories"<<endl;
  std::vector<TString> keylist;
  std::vector<TString> notdonelist;
  std::vector<int> nsubdirs;
  //TDirectory *dirsav = dir;
  list = dir->GetListOfKeys();
  int ikey=0;
  //int localkey=0;
  //int ntimes=0;

  //*******************************
  //Get First Directory
  cout<<"Get First Directory"<<endl;
  for(ikey=0;ikey<list->GetEntries();  ikey++)
  {
    TKey *thekey = (TKey*)list->At(ikey);
    if(thekey==0) continue;
      
    TString keyname=thekey->GetName();
    TString keytype=thekey->GetClassName();
    keytype.ReplaceAll(" ","");
    if(keyname=="EventInfo") continue;
    if(keytype=="TDirectoryFile")
    {
      TString dirname=dir->GetPath();
      dirname+="/";
      dirname+=keyname;
      dirname.Remove(0,dirname.Last(':')+1);   
      dir=file->GetDirectory(dirname);
      list=dir->GetListOfKeys();
      if(dirname.Contains(comparestring)) dirlist.push_back(dirname);
      else
      {
	notdonelist.push_back(dirname);
	nsubdirs.push_back(-1);
      }
    }
  }
  
  //*******************************
  //Get List of Directory
  cout<<"Get List of Directory"<<endl;
  unsigned int nempty=0;
  while(nempty!=notdonelist.size())
  {
    for(unsigned int idir=0; idir<notdonelist.size(); ++idir)
    {
      if(nsubdirs[idir]==0) continue;
      dir = file->GetDirectory(notdonelist[idir]); 
      list= dir->GetListOfKeys();
      int ndirectories=0;
      for(ikey=0;ikey<list->GetEntries();  ikey++)
      {
	      TKey *thekey = (TKey*)list->At(ikey);
	      if(thekey==0) continue;
	      TString keyname=thekey->GetName();
	      TString keytype=thekey->GetClassName();
	      keytype.ReplaceAll(" ","");
	      if(keytype=="TDirectoryFile")
	      {
	        TString dirname=dir->GetPath();
	        dirname+="/";
	        dirname+=keyname;
                dirname.Remove(0,dirname.Last(':')+1);
	        ndirectories++;
	        if(dirname.Contains(comparestring)) dirlist.push_back(dirname);
	        else
	        {
	          notdonelist.push_back(dirname);
	          nsubdirs.push_back(-1);
	        }
	      }
      }
      nsubdirs[idir]=ndirectories;
    }
    // count number of done dirs;
    nempty=0;
    for(unsigned int i=0; i<nsubdirs.size(); i++) if(nsubdirs[i]!=-1)	nempty++;
  }

  return dirlist;

}


int getDetId(TString hname){
  return hname.Remove(0,hname.Last('_')+1).Atoi();
}


void getBarrelPos(TString& dirname , int& layer , int& ladder , int& module){
  if(! dirname.Contains("Barrel")){
    cout << "module not in barrel" << endl;
    return ;
  }
  
  TString tmp = dirname , tmp2 = dirname ;
  //getting layer
  module=tmp.Remove(0,tmp.Last('_')+1).Atoi();
  
  tmp2.Remove(tmp2.Last('_') , tmp2.Length());
  tmp = tmp2;
  
  ladder=tmp.Remove(0,tmp.Last('_')+1).Remove(2,tmp.Length()).Atoi();
  
  tmp2.Remove(tmp2.Last('_') , tmp2.Length());
  tmp = tmp2;
  
  layer=tmp.Remove(0,tmp.Last('_')+1).Remove(1,tmp.Length()).Atoi();
  
  if(dirname.Contains("Shell_m"))
    module = -1 * module;
  if(dirname.Contains("Shell_mI") || dirname.Contains("Shell_pI"))
    ladder = -1 * ladder;
}
