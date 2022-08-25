#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "FairTrackParam.h"
#include "CbmStsKFTrackFitter.h"
#include "CbmStsTrack.h"
#include "BmnEventHeader.h"
#include "BmnTrigDigit.h"
#include "BmnZDCDigit.h"
#include "BmnTrack.h"
#include "BmnGemTof700IdentifiableTrack.h"
#include "CbmVertex.h"
#include "BmnGemDchTrack.h"
#include "BmnTofHit.h"
#include "BmnGemTrack.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TSystem.h"

#include <TStopwatch.h>
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include <vector>
#include <fstream>
#include <map>
#include </eos/nica/bmn/users/kovachev/efficiency/cbmroot/bmnroot-root6-csconly/bmndst/DstEventHeader.h>
#include "CbmKF.h"
#include "BmnKalmanFilter.h"
#include "CbmKFTrack.h"

struct
{
    Float_t e;
    Float_t X_digZDC;
    Float_t Y_digZDC;
    Double_t dist;

} Zdc_cell[105];

using namespace std;
using namespace ROOT::Math;



void resid(){

Int_t StartNRun = 3756;
Int_t FinishNRun = 3757;
Int_t ch;
Float_t e;

/*Bool_t compare_tracks(pair<Double_t,>a, pair<Double_t,FairHit>b){
  return(a.first < b.first);
 }*/


TObjArray *Hlist = new TObjArray();
//dist,resx,resy,
  TH1F *hdist_l = new TH1F("hdist_l","dist linear",100,0,160);
 // TH1F *hdist_k = new TH1F("hdist_k","Kalman linear",100,0,160);
   // TH2F *henZDCTtr[105];
    TH1D *hxresid[105];
    TH1D *hyresid[105];
    TH2F *hxyresid[105];

    TH1D *hxresid_line[105];
    TH1D *hyresid_line[105];
    TH2F *hxyresid_line[105];
   // TH1D *dist[105];

    for (Int_t i = 1; i < 105; i++){
     
      //TString en_ZDCTtr("en_ZDCTtr");
      TString hx_resid("hx_resid");
      TString hy_resid("hy_resid");
      TString hxy_resid("hxy_resid");

     // TString dist_line("dist_line");
      TString hx_resid_line_extr("hx_resid_line_extr");
      TString hy_resid_line_extr("hy_resid_line_extr");
      TString hxy_resid_line_extr("hxy_resid_line_extr");
    
     // en_ZDCTtr += to_string(i);
      hx_resid += to_string(i);
      hy_resid += to_string(i);
      hxy_resid += to_string(i);
      
      hx_resid_line_extr += to_string(i);
      hy_resid_line_extr += to_string(i);
      hxy_resid_line_extr += to_string(i);
    // dist_line += to_string(i);

     // henZDCTtr[i] = new TH2F(en_ZDCTtr, "Energy release in the ZDC channel", 100, 0, 20,100,0,20);
      hxresid[i] = new TH1D(hx_resid,"Kalman  extrapolate Residual Xtr-Xmod",100,-200,200);
      hyresid[i] = new TH1D(hy_resid,"Kalman extrapolate Residual Ytr-Ymod",100,-200,200);
      hxyresid[i] = new TH2F(hxy_resid, "Kalman extrapolate Residual XY",100,-200,200,100,-200,200);

      hxresid_line[i] = new TH1D(hx_resid_line_extr,"Linear extrapolate Residual Xtr-Xmod",100,-200,200);
      hyresid_line[i] = new TH1D(hy_resid_line_extr,"Linear extrapolate Residual Ytr-Ymod",100,-200,200);
      hxyresid_line[i] = new TH2F(hxy_resid_line_extr, "Linear extrapolate Residual XY",100,-200,200,100,-200,200);
     // dist[i] = new TH1D(dist_line,"Distance", 100, 0, 150);
    // henZDCTtr[i]->SetFillColor(5);
    //henZDCTtr[i]->SetOption("BOX");
      //Hlist->Add(henZDCTtr[i]);
     /* Hlist->Add(hxresid[i]);
      Hlist->Add(hyresid[i]);
      Hlist->Add(hxyresid[i]);*/
      Hlist->Add(hxresid_line[i]);
      Hlist->Add(hyresid_line[i]);
      Hlist->Add(hxyresid_line[i]);
     // Hlist->Add(dist[i]);
      }

      for (Int_t filenum = StartNRun; filenum <= FinishNRun; filenum = filenum + 1)
      {

          TString fname = Form("/nica/mpd6/ks_alish/pid_ZDC_r7/BmnGemTof700IdentifiableTracks%d_ZDC.root", filenum);
        // TString fname = Form("/nica/mpd6/ks_alish/pid_ZDC_r7/BmnGemTof700IdentifiableTracks%d_ZDC.root", filenum);
          if (gSystem->AccessPathName(fname))
          {
              cout << " file not exist !!! " << fname << endl;
              continue;
          }

          TChain *rec = new TChain("bmndata");
          rec->Add(fname);

          DstEventHeader *EventHeaderIdentified = NULL;
          TClonesArray *tracks = NULL;
          CbmVertex *primaryVertex = NULL;
          TClonesArray *ZdcDigit = NULL;

          rec->SetBranchAddress("EventHeader.", &EventHeaderIdentified);
          rec->SetBranchAddress("tracks", &tracks);
          rec->SetBranchAddress("PrimaryVertex", &primaryVertex);
          rec->SetBranchAddress("ZdcDigit", &ZdcDigit);

        //  FairTrackParam par = FairTrackParam();
          FairTrackParam DchTrackLastParam = FairTrackParam();

          Int_t events = rec->GetEntries();
          Double_t gemX, gemY, gemZ, gemTx, gemTy;
          Double_t gemX_k, gemY_k, gemZ_k, gemTx_k, gemTy_k;
          Double_t gemX_ZDCPlane;
          Double_t gemY_ZDCPlane;

            for (Int_t ev = 0; ev < events; ev++)
          {    
              
              if ((ev % 10000) == 0)
                
              cout << " Event= " << ev << "/" << events << endl;
              rec->GetEntry(ev);

              for (Int_t iTrack0 = 0; iTrack0 < tracks->GetEntriesFast(); iTrack0++)
              {

                  BmnGemTof700IdentifiableTrack *identifiableTracks = (BmnGemTof700IdentifiableTrack *)tracks->At(iTrack0);
                  //if (identifiableTracks->GetTrack()->GetDchTrack()->GetUniqueID() != 100) continue;
                  DchTrackLastParam = *(identifiableTracks->GetTrack()->GetDchTrack()->GetParamLast());

                 // if (identifiableTracks->GetTrack()->GetDchTrack()->GetUniqueID() != 100) continue;
                  gemX = DchTrackLastParam.GetX();
                  gemY = DchTrackLastParam.GetY();
                  gemZ = DchTrackLastParam.GetZ();
                  gemTx = DchTrackLastParam.GetTx();
                  gemTy = DchTrackLastParam.GetTy();
                  Double_t Z0zdc = 1000.0;
                 // Double_t Z0zdc = 955.603;

                 //cout << " X, Y track before extrapolate  " <<  gemX << ", " << gemY << " gemZ" << gemZ << endl;
                
                  // Linearly extrapolate GEM - DCH track to the ZDC cell
                  gemX_ZDCPlane = gemX + (Z0zdc - gemZ) * gemTx;
                  gemY_ZDCPlane = gemY + (Z0zdc - gemZ) * gemTy;

               //  if(abs(gemX_ZDCPlane) > 75 || abs(gemY_ZDCPlane) > 45) continue;
            //   if(abs(gemX_ZDCPlane) < 22.5 || abs(gemY_ZDCPlane) < 22.5 && abs(gemX_ZDCPlane) > 75.5 || abs(gemY_ZDCPlane) > 45.3) continue;

                  //cout<<"Linear extrapolate X, Y  "<<  gemX_ZDCPlane << ", " << gemY_ZDCPlane << endl;
                  // проверка внешних границ(check out)
                  if (abs(gemX_ZDCPlane) > 75.5 || abs(gemY_ZDCPlane) > 45)
                    continue;
                  //проверка внутренних границ(check in)
                  if (abs(gemX_ZDCPlane) < 22.5 && abs(gemY_ZDCPlane) < 22.5)
                    continue;
                  // multimapa<Dist,ii> dist_mod;
                  multimap<Double_t, Int_t> dist_mod;
                //multimap<номер модуля, multimap<дистанция,номер трека> trmodule;
                  map<Int_t, multimap<Double_t, Int_t>> trmodule;
                  for (Int_t ii = 0; ii < ZdcDigit->GetEntries(); ii++)
                  {
                    BmnZDCDigit *p = (BmnZDCDigit *)ZdcDigit->At(ii);
                    ch = p->GetChannel();
                    if (ch < 1 || ch > 104)
                      continue;

                    // coordinate centr module of the ZDC
                   Zdc_cell[ch].X_digZDC = p->GetX() / 10;
                   Zdc_cell[ch].Y_digZDC = p->GetY() / 10;

                   // cout<< "Get zdc centr module"<< Zdc_cell[ch].X_digZDC <<" ,  "<< Zdc_cell[ch].Y_digZDC << endl;
                   Double_t dx_line = gemX_ZDCPlane - Zdc_cell[ch].X_digZDC;
                   Double_t dy_line = gemY_ZDCPlane - Zdc_cell[ch].Y_digZDC;

                   Double_t dist = TMath::Sqrt(dx_line * dx_line + dy_line * dy_line);

                   cout << "dx ,dy residual" << dx_line << "  ,  " << dy_line << " dist =  " << dist << endl;
                   trmodule[ch].insert(make_pair(dist,ii));
                   

                   hxresid_line[ch]->Fill(dx_line);
                   hyresid_line[ch]->Fill(dy_line);
                   hxyresid_line[ch]->Fill(dx_line, dy_line);

                   // dist[ch]->Fill(Zdc_cell[ch].dist);
                  }
                  cout<< " iTrack  " << iTrack0 << "  track module size = " << trmodule.size() << endl;
                 Double_t min_dist = 99999.9;
                 Int_t min_chan = 0;
                 Int_t min_ii = 0;
                  // хоббит ищет минимальный дистанс, удачки ему
                  for(auto it = trmodule.begin(); it != trmodule.end(); it++ ){
                    Int_t chan = it->first;
                    Double_t dist0 = it->second.begin()->first;
                    Int_t ii0 = it->second.begin()->second;
                    if(dist0 < min_dist){

                      min_dist = dist0;
                      min_chan = chan;
                      min_ii = ii0; 
                    }


                  }
                  cout << "  min_dist = " << min_dist << "  min_chan  "<<min_chan << "  min_ii  " << min_ii << endl; 

              } // end loop by good tracks pid tof700
          }     // ev
      }
      TFile *outfile = new TFile("ZDC_extrapolate_resid_line_beamzone_Z.root", "UPDATE");
      Hlist->Write();
      outfile->Write();
      outfile->Close();
      cout << "file write" << endl;
}