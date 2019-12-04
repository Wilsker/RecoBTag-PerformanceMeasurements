#include "TTbarEventAnalysis.h"
#include "TLorentzVector.h"

using namespace std;

//
void TTbarEventAnalysis::prepareOutput(TString outFile)
{

  //prepare output file
  outF_=TFile::Open(outFile,"RECREATE");

  //init KIN tree
  kinTree_=new TTree("kin","kinematics analysis");
  kinTree_->SetDirectory(outF_);
  kinTree_->Branch("EventInfo",    eventInfo_,         "EventInfo[3]/I");
  kinTree_->Branch("ttbar_chan",    &ttbar_chan_,      "ttbar_chan/I");
  kinTree_->Branch("npvn",    &npv_,      "npv/I");
  kinTree_->Branch("flavour",        jetFlavour_,      "flavour/I");
  kinTree_->Branch("jetmult",       &jetmult_,         "jetmult/I");
  kinTree_->Branch("jetpt",          jetPt_,           "jetpt/F");
  kinTree_->Branch("jeteta",         jetEta_,          "jeteta/F");
  kinTree_->Branch("jetrank",       &jetrank_,          "jetrank/I");
  kinTree_->Branch("close_mlj",      close_mlj_,       "close_mlj[5]/F");
  kinTree_->Branch("close_deta",    &close_deta_,      "close_deta/F");
  kinTree_->Branch("close_dphi",    &close_dphi_,      "close_dphi/F");
  kinTree_->Branch("close_ptrel",   &close_ptrel_,     "close_ptrel/F");
  kinTree_->Branch("close_lj2ll_deta",    &close_lj2ll_deta_,      "close_lj2ll_deta/F");
  kinTree_->Branch("close_lj2ll_dphi",    &close_lj2ll_dphi_,      "close_lj2ll_dphi/F");
  kinTree_->Branch("far_mlj",       &far_mlj_,         "far_mlj/F");
  kinTree_->Branch("far_deta",      &far_deta_,        "far_deta/F");
  kinTree_->Branch("far_dphi",      &far_dphi_,        "far_dphi/F");
  kinTree_->Branch("far_ptrel",     &far_ptrel_,       "far_ptrel/F");
  kinTree_->Branch("far_lj2ll_deta",    &far_lj2ll_deta_,      "far_lj2ll_deta/F");
  kinTree_->Branch("far_lj2ll_dphi",    &far_lj2ll_dphi_,      "far_lj2ll_dphi/F");
  kinTree_->Branch("j2ll_deta",     &j2ll_deta_,       "j2ll_deta/F");
  kinTree_->Branch("j2ll_dphi",     &j2ll_dphi_,       "j2ll_dphi/F");
  kinTree_->Branch("kindisc",        kinDisc_,         "kindisc[5]/F");
  kinTree_->Branch("jp",             &(jp_[0]),              "jp/F");
  kinTree_->Branch("svhe",           &(svhe_[0]),            "svhe/F");
  kinTree_->Branch("csv",            &(csv_[0]),             "csv/F");
  kinTree_->Branch("DeepCSVb",               &(DeepCSVb_[0]),        "DeepCSVb/F");
  kinTree_->Branch("DeepCSVc",               &(DeepCSVc_[0]),        "DeepCSVc/F");
  kinTree_->Branch("DeepCSVl",               &(DeepCSVl_[0]),        "DeepCSVl/F");
  kinTree_->Branch("DeepCSVbb",               &(DeepCSVbb_[0]),        "DeepCSVbb/F");
  kinTree_->Branch("DeepCSVcc",               &(DeepCSVcc_[0]),        "DeepCSVcc/F");
  kinTree_->Branch("DeepCSVbN",               &(DeepCSVbN_[0]),        "DeepCSVbN/F");
  kinTree_->Branch("DeepCSVcN",               &(DeepCSVcN_[0]),        "DeepCSVcN/F");
  kinTree_->Branch("DeepCSVlN",               &(DeepCSVlN_[0]),        "DeepCSVlN/F");
  kinTree_->Branch("DeepCSVbbN",               &(DeepCSVbbN_[0]),        "DeepCSVbbN/F");
  kinTree_->Branch("DeepCSVccN",               &(DeepCSVccN_[0]),        "DeepCSVccN/F");
  kinTree_->Branch("DeepCSVbP",               &(DeepCSVbP_[0]),        "DeepCSVbP/F");
  kinTree_->Branch("DeepCSVcP",               &(DeepCSVcP_[0]),        "DeepCSVcP/F");
  kinTree_->Branch("DeepCSVlP",               &(DeepCSVlP_[0]),        "DeepCSVlP/F");
  kinTree_->Branch("DeepCSVbbP",               &(DeepCSVbbP_[0]),        "DeepCSVbbP/F");
  kinTree_->Branch("DeepCSVccP",               &(DeepCSVccP_[0]),        "DeepCSVccP/F");
  kinTree_->Branch("DeepCSVBDisc",               &(DeepCSVBDisc_[0]),        "DeepCSVBDisc/F");
  kinTree_->Branch("DeepCSVBDiscN",               &(DeepCSVBDiscN_[0]),        "DeepCSVBDiscN/F");
  kinTree_->Branch("DeepCSVBDiscP",               &(DeepCSVb_[0]),        "DeepCSVb/F");
  kinTree_->Branch("DeepCSVCvsLDisc",               &(DeepCSVCvsLDisc_[0]),        "DeepCSVCvsLDisc/F");
  kinTree_->Branch("DeepCSVCvsLDiscN",               &(DeepCSVCvsLDiscN_[0]),        "DeepCSVCvsLDiscN/F");
  kinTree_->Branch("DeepCSVCvsLDiscP",               &(DeepCSVCvsLDiscP_[0]),        "DeepCSVCvsLDiscP/F");
  kinTree_->Branch("DeepCSVCvsBDisc",               &(DeepCSVCvsBDisc_[0]),        "DeepCSVCvsBDisc/F");
  kinTree_->Branch("DeepCSVCvsBDiscN",               &(DeepCSVCvsBDiscN_[0]),        "DeepCSVCvsBDiscN/F");
  kinTree_->Branch("DeepCSVCvsBDiscP",               &(DeepCSVCvsBDiscP_[0]),        "DeepCSVCvsBDiscP/F");

  kinTree_->Branch("DeepFlavourBDisc",            &(DeepFlavourBDisc_[0]),           "DeepFlavourBDisc/F");
  kinTree_->Branch("DeepFlavourCvsLDisc",            &(DeepFlavourCvsLDisc_[0]),           "DeepFlavourCvsLDisc/F");
  kinTree_->Branch("DeepFlavourCvsBDisc",            &(DeepFlavourCvsBDisc_[0]),           "DeepFlavourCvsBDisc/F");
  kinTree_->Branch("DeepFlavourB",            &(DeepFlavourB_[0]),           "DeepFlavourB/F");
  kinTree_->Branch("DeepFlavourBB",            &(DeepFlavourBB_[0]),           "DeepFlavourBB/F");
  kinTree_->Branch("DeepFlavourLEPB",            &(DeepFlavourLEPB_[0]),           "DeepFlavourLEPB/F");

  //kinTree_->Branch("weight",         weight_,          "weight[15]/F");
  kinTree_->Branch("weight",         weight_,          "weight[26]/F");

  ftmTree_=new TTree("ftm","flavour tag matching");
  ftmTree_->SetDirectory(outF_);
  ftmTree_->Branch("EventInfo",      eventInfo_,  "EventInfo[3]/I");
  ftmTree_->Branch("ttbar_chan",    &ttbar_chan_, "ttbar_chan/I");
  ftmTree_->Branch("jetmult",       &jetmult_,    "jetmult/I");
  ftmTree_->Branch("flavour",        jetFlavour_, "flavour[2]/I");
  ftmTree_->Branch("jetpt",          jetPt_,      "jetpt[2]/F");
  ftmTree_->Branch("jeteta",         jetEta_,     "jeteta[2]/F");
  ftmTree_->Branch("jp",             jp_,         "jp[2]/F");
  ftmTree_->Branch("svhe",           svhe_,       "svhe[2]/F");
  ftmTree_->Branch("csv",            csv_,        "csv[2]/F");

  ftmTree_->Branch("DeepCSVb",               DeepCSVb_,        "DeepCSVb[2]/F");
  ftmTree_->Branch("DeepCSVc",               DeepCSVc_,        "DeepCSVc[2]/F");
  ftmTree_->Branch("DeepCSVl",               DeepCSVl_,        "DeepCSVl[2]/F");
  ftmTree_->Branch("DeepCSVbb",               DeepCSVbb_,        "DeepCSVbb[2]/F");
  ftmTree_->Branch("DeepCSVcc",               DeepCSVcc_,        "DeepCSVcc[2]/F");
  ftmTree_->Branch("DeepCSVbN",               DeepCSVbN_,        "DeepCSVbN[2]/F");
  ftmTree_->Branch("DeepCSVcN",               DeepCSVcN_,        "DeepCSVcN[2]/F");
  ftmTree_->Branch("DeepCSVlN",               DeepCSVlN_,        "DeepCSVlN[2]/F");
  ftmTree_->Branch("DeepCSVbbN",               DeepCSVbbN_,        "DeepCSVbbN[2]/F");
  ftmTree_->Branch("DeepCSVccN",               DeepCSVccN_,        "DeepCSVccN[2]/F");
  ftmTree_->Branch("DeepCSVbP",               DeepCSVbP_,        "DeepCSVbP[2]/F");
  ftmTree_->Branch("DeepCSVcP",               DeepCSVcP_,        "DeepCSVcP[2]/F");
  ftmTree_->Branch("DeepCSVlP",               DeepCSVlP_,        "DeepCSVlP[2]/F");
  ftmTree_->Branch("DeepCSVbbP",               DeepCSVbbP_,        "DeepCSVbbP[2]/F");
  ftmTree_->Branch("DeepCSVccP",               DeepCSVccP_,        "DeepCSVccP[2]/F");

  ftmTree_->Branch("DeepCSVBDisc",               DeepCSVBDisc_,        "DeepCSVBDisc[2]/F");
  ftmTree_->Branch("DeepCSVBDiscN",               DeepCSVBDiscN_,        "DeepCSVBDiscN[2]/F");
  ftmTree_->Branch("DeepCSVBDiscP",               DeepCSVBDiscP_,        "DeepCSVBDiscP[2]/F");
  ftmTree_->Branch("DeepCSVCvsLDisc",               DeepCSVCvsLDisc_,        "DeepCSVCvsLDisc[2]/F");
  ftmTree_->Branch("DeepCSVCvsLDiscN",               DeepCSVCvsLDiscN_,        "DeepCSVCvsLDiscN[2]/F");
  ftmTree_->Branch("DeepCSVCvsLDiscP",               DeepCSVCvsLDiscP_,        "DeepCSVCvsLDiscP[2]/F");
  ftmTree_->Branch("DeepCSVCvsBDisc",               DeepCSVCvsBDisc_,        "DeepCSVCvsBDisc[2]/F");
  ftmTree_->Branch("DeepCSVCvsBDiscN",               DeepCSVCvsBDiscN_,        "DeepCSVCvsBDiscN[2]/F");
  ftmTree_->Branch("DeepCSVCvsBDiscP",               DeepCSVCvsBDiscP_,        "DeepCSVCvsBDiscP[2]/F");

  ftmTree_->Branch("kindisc",        kinDisc_,    "kindisc[2]/F");
  //ftmTree_->Branch("weight",         weight_,     "weight[15]/F");
  ftmTree_->Branch("weight",         weight_,     "weight[26]/F");

  //prepare histograms
  std::map<TString,TH1F *> baseHistos;
  baseHistos["npvinc" ]  = new TH1F("npvinc", ";N_{PV,good}-N_{HS};Events",              50, 0, 50);
  baseHistos["npv"    ]  = new TH1F("npv",    ";N_{PV,good}-N_{HS};Events",              50, 0, 50);
  baseHistos["rho"    ]  = new TH1F("rho",    ";#rho [GeV];Events",                      50, 0, 30);
  baseHistos["mll"    ]  = new TH1F("mll",    ";Dilepton invariant mass [GeV];Events",   20, 0, 300);
  baseHistos["mllinc" ]  = new TH1F("mllinc", ";Dilepton invariant mass [GeV];Events",   20, 0, 300);
  baseHistos["met"    ]  = new TH1F("met",    ";Missing transverse energy [GeV];Events", 15, 0, 300);
  baseHistos["njets"  ]  = new TH1F("njets",  ";Jet multiplicity;Events;",               6,  2, 8);
  baseHistos["leadjpt"]  = new TH1F("leadjpt",";Leading jet p_{T} [GeV];Events;",        14,30,300);
  baseHistos["leadlpt"]  = new TH1F("leadlpt",";Leading lepton p_{T} [GeV];Events;",     9,20,200);
  baseHistos["trailjpt"] = new TH1F("trailjpt",";Trailing jet p_{T} [GeV];Events;",      14,30,300);
  baseHistos["traillpt"] = new TH1F("traillpt",";Trailing lepton p_{T} [GeV];Events;",   9,20,200);
  baseHistos["leadjeta"]    = new TH1F("leadjeta",    ";Pseudo-rapidity; Jets",              25, 0, 2.5);
  baseHistos["trailjeta"]   = new TH1F("trailjeta",    ";Pseudo-rapidity; Jets",              25, 0, 2.5);
  baseHistos["evsel"]    = new TH1F("evsel",   ";Event selection;Events;",               4,0,4);
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(1,"#geq 2j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(2,"=2j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(3,"=3j");
  baseHistos["evsel"]->GetXaxis()->SetBinLabel(4,"=4j");
  baseHistos["jp"]=new TH1F("jp",";Jet probability;Jets",50,0,3);
  baseHistos["svhe"]=new TH1F("svhe",";Simple secondary vertex (HE);Jets",50,0,6);
  baseHistos["csv"]=new TH1F("csv",";Combined secondary vertex (IVF);Jets",50,0,1.1);
  baseHistos["tche"]=new TH1F("tche",";Track Counting High Efficiency;Jets",50,-20,50);
  baseHistos["jetseltrk"]=new TH1F("jetseltrk",";Selected track multiplicity;Jets",15,0,15);
  baseHistos["jp_leadkin"]=new TH1F("jp_leadkin",";Jet probability;Jets",50,0,3);
  baseHistos["svhe_leadkin"]=new TH1F("svhe_leadkin",";Simple secondary vertex (HE);Jets",50,0,6);
  baseHistos["csv_leadkin"]=new TH1F("csv_leadkin",";Combined secondary vertex (IVF);Jets",50,0,1.1);
  baseHistos["DeepCSVb"]=new TH1F("DeepCSVb",";DeepCSV b ;Jets",50,0,1.10);
  baseHistos["DeepCSVc"]=new TH1F("DeepCSVc",";DeepCSV c;Jets",50,0,1.10);
  baseHistos["DeepCSVl"]=new TH1F("DeepCSVl",";DeepCSV l;Jets",50,0,1.10);
  baseHistos["DeepCSVbb"]=new TH1F("DeepCSVbb",";DeepCSV bb;Jets",50,0,1.10);
  baseHistos["DeepCSVcc"]=new TH1F("DeepCSVcc",";DeepCSV cc;Jets",50,0,1.10);
  baseHistos["DeepCSVbN"]=new TH1F("DeepCSVbN",";DeepCSV bN;Jets",50,0,1.10);
  baseHistos["DeepCSVcN"]=new TH1F("DeepCSVcN",";DeepCSV cN;Jets",50,0,1.10);
  baseHistos["DeepCSVlN"]=new TH1F("DeepCSVlN",";DeepCSV lN;Jets",50,0,1.10);
  baseHistos["DeepCSVbbN"]=new TH1F("DeepCSVbbN",";DeepCSV bbN;Jets",50,0,1.10);
  baseHistos["DeepCSVccN"]=new TH1F("DeepCSVccN",";DeepCSV ccN;Jets",50,0,1.10);
  baseHistos["DeepCSVbP"]=new TH1F("DeepCSVbP",";DeepCSV bP;Jets",50,0,1.10);
  baseHistos["DeepCSVcP"]=new TH1F("DeepCSVcP",";DeepCSV cP;Jets",50,0,1.10);
  baseHistos["DeepCSVlP"]=new TH1F("DeepCSVlP",";DeepCSV lP;Jets",50,0,1.10);
  baseHistos["DeepCSVbbP"]=new TH1F("DeepCSVbbP",";DeepCSV bbP;Jets",50,0,1.10);
  baseHistos["DeepCSVccP"]=new TH1F("DeepCSVccP",";DeepCSV ccP;Jets",50,0,1.10);
  baseHistos["DeepCSVBDisc"]=new TH1F("DeepCSVBDisc",";DeepCSV B Disc. ;Jets",50,0,1.10);
  baseHistos["DeepCSVBDiscN"]=new TH1F("DeepCSVBDiscN",";DeepCSV B Disc. N ;Jets",50,0,1.10);
  baseHistos["DeepCSVBDiscP"]=new TH1F("DeepCSVBDiscP",";DeepCSV B Disc. P ;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsLDisc"]=new TH1F("DeepCSVCvsLDisc",";DeepCSV CvsL Disc. ;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsLDiscN"]=new TH1F("DeepCSVCvsLDiscN",";DeepCSV CvsL Disc. N ;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsLDiscP"]=new TH1F("DeepCSVCvsLDiscP",";DeepCSV CvsL Disc. P ;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsBDisc"]=new TH1F("DeepCSVCvsBDisc",";DeepCSV CvsB Disc. ;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsBDiscN"]=new TH1F("DeepCSVCvsBDiscN",";DeepCSV CvsB Disc. N ;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsBDiscP"]=new TH1F("DeepCSVCvsBDiscP",";DeepCSV CvsB Disc. P ;Jets",50,0,1.10);
  baseHistos["DeepCSVb_leadkin"]=new TH1F("DeepCSVb_leadkin",";DeepCSV b leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVc_leadkin"]=new TH1F("DeepCSVc_leadkin",";DeepCSV c leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVl_leadkin"]=new TH1F("DeepCSVl_leadkin",";DeepCSV l leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVbb_leadkin"]=new TH1F("DeepCSVbb_leadkin",";DeepCSV bb leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVcc_leadkin"]=new TH1F("DeepCSVcc_leadkin",";DeepCSV cc leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVbN_leadkin"]=new TH1F("DeepCSVbN_leadkin",";DeepCSV b N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVcN_leadkin"]=new TH1F("DeepCSVcN_leadkin",";DeepCSV c N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVlN_leadkin"]=new TH1F("DeepCSVlN_leadkin",";DeepCSV l N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVbbN_leadkin"]=new TH1F("DeepCSVbbN_leadkin",";DeepCSV bb N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVccN_leadkin"]=new TH1F("DeepCSVccN_leadkin",";DeepCSV cc N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVbP_leadkin"]=new TH1F("DeepCSVbP_leadkin",";DeepCSV b P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVcP_leadkin"]=new TH1F("DeepCSVcP_leadkin",";DeepCSV c P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVlP_leadkin"]=new TH1F("DeepCSVlP_leadkin",";DeepCSV l P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVbbP_leadkin"]=new TH1F("DeepCSVbbP_leadkin",";DeepCSV bb P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVccP_leadkin"]=new TH1F("DeepCSVccP_leadkin",";DeepCSV cc P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVBDisc_leadkin"]=new TH1F("DeepCSVBDisc_leadkin",";DeepCSV B Disc. leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVBDiscN_leadkin"]=new TH1F("DeepCSVBDiscN_leadkin",";DeepCSV B Disc. N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVBDiscP_leadkin"]=new TH1F("DeepCSVBDiscP_leadkin",";DeepCSV B Disc. P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsLDisc_leadkin"]=new TH1F("DeepCSVCvsLDisc_leadkin",";DeepCSV CvsL Disc. leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsLDiscN_leadkin"]=new TH1F("DeepCSVCvsLDiscN_leadkin",";DeepCSV CvsL Disc. N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsLDiscP_leadkin"]=new TH1F("DeepCSVCvsLDiscP_leadkin",";DeepCSV CvsL Disc. P leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsBDisc_leadkin"]=new TH1F("DeepCSVCvsBDisc_leadkin",";DeepCSV CvsB Disc. leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsBDiscN_leadkin"]=new TH1F("DeepCSVCvsBDiscN_leadkin",";DeepCSV CvsB Disc. N leadkin;Jets",50,0,1.10);
  baseHistos["DeepCSVCvsBDiscP_leadkin"]=new TH1F("DeepCSVCvsBDiscP_leadkin",";DeepCSV CvsB Disc. P leadkin;Jets",50,0,1.10);

  baseHistos["DeepFlavourBDisc"]=new TH1F("DeepFlavourBDisc",";DeepFlavour B Disc. ;Jets",50,0,1.10);
  baseHistos["DeepFlavourCvsLDisc"]=new TH1F("DeepFlavourCvsLDisc",";DeepFlavour CvsL Disc. ;Jets",50,0,1.10);
  baseHistos["DeepFlavourCvsBDisc"]=new TH1F("DeepFlavourCvsBDisc",";DeepFlavour CvsB Disc. ;Jets",50,0,1.10);
  baseHistos["DeepFlavourB"]=new TH1F("DeepFlavourB",";DeepFlavourB ;Jets",50,0,1.10);
  baseHistos["DeepFlavourBB"]=new TH1F("DeepFlavourBB",";DeepFlavourBB ;Jets",50,0,1.10);
  baseHistos["DeepFlavourLEPB"]=new TH1F("DeepFlavourLEPB",";DeepFlavourLEPB ;Jets",50,0,1.10);

  baseHistos["DeepFlavourBDisc_leadkin"]=new TH1F("DeepFlavourBDisc_leadkin",";DeepFlavour B Disc. leadkin;Jets",50,0,1.10);
  baseHistos["DeepFlavourCvsLDisc_leadkin"]=new TH1F("DeepFlavourCvsLDisc_leadkin",";DeepFlavour CvsL Disc. leadkin;Jets",50,0,1.10);
  baseHistos["DeepFlavourCvsBDisc_leadkin"]=new TH1F("DeepFlavourCvsBDisc_leadkin",";DeepFlavour CvsB Disc. leadkin;Jets",50,0,1.10);
  baseHistos["DeepFlavourB_leadkin"]=new TH1F("DeepFlavourB_leadkin",";DeepFlavourB leadkin;Jets",50,0,1.10);
  baseHistos["DeepFlavourBB_leadkin"]=new TH1F("DeepFlavourBB_leadkin",";DeepFlavourBB leadkin;Jets",50,0,1.10);
  baseHistos["DeepFlavourLEPB_leadkin"]=new TH1F("DeepFlavourLEPB_leadkin",";DeepFlavourLEPB leadkin;Jets",50,0,1.10);

  baseHistos["tche_leadkin"]=new TH1F("tche_leadkin",";Track Counting High Efficiency;Jets",50,-20,50);
  baseHistos["jetseltrk_leadkin"]=new TH1F("jetseltrk_leadkin",";Selected track multiplicity;Jets",15,0,15);
  baseHistos["flavour"] = new TH1F("flavour",     ";Jet flavour;Jets",                   4,  0, 4 );
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(1,"unmatched");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(2,"udsg");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(3,"c");
  baseHistos["flavour"]->GetXaxis()->SetBinLabel(4,"b");
  baseHistos["close_mlj"]   = new TH1F("close_mlj",   ";M(lepton,jet) [GeV]; Jets",          50, 0, 250);
  baseHistos["close_deta"]  = new TH1F("close_deta",  ";#Delta#eta(lepton,jet); Jets",       50, 0, 4);
  baseHistos["close_dphi"]  = new TH1F("close_dphi",  ";#Delta#phi(lepton,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["close_ptrel"] = new TH1F("close_ptrel", ";p_{T}^{rel}(lepton,jet) [GeV];Jets", 50, 0, 1);
  baseHistos["close_lj2ll_deta"] = new TH1F("close_lj2ll_deta",  ";#Delta#eta(lj,ll); Jets",       50, 0, 4);
  baseHistos["close_lj2ll_dphi"]  = new TH1F("close_l2jll_dphi",  ";#Delta#phi(lj,ll) [rad]; Jets", 50, 0, 3.15);
  baseHistos["far_mlj"]     = new TH1F("far_mlj",     ";M(lepton,jet) [GeV]; Jets",          50, 0, 250);
  baseHistos["far_deta"]    = new TH1F("far_deta",    ";#Delta#eta(lepton,jet); Jets",       50, 0, 4);
  baseHistos["far_dphi"]    = new TH1F("far_dphi",    ";#Delta#phi(lepton,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["far_ptrel"]   = new TH1F("far_ptrel",   ";p_{T}^{rel}(lepton,jet) [GeV];Jets", 50, 0, 1);
  baseHistos["far_lj2ll_deta"] = new TH1F("far_lj2ll_deta",  ";#Delta#eta(lepton+jet,ll); Jets",       50, 0, 4);
  baseHistos["far_lj2ll_dphi"]  = new TH1F("far_l2jll_dphi",  ";#Delta#phi(lj,ll) [rad]; Jets", 50, 0, 3.15);
  baseHistos["j2ll_deta"]    = new TH1F("j2ll_deta",    ";#Delta#eta(ll,jet); Jets",       50, 0, 4);
  baseHistos["j2ll_dphi"]    = new TH1F("j2ll_dphi",    ";#Delta#phi(lll,jet) [rad]; Jets", 50, 0, 3.15);
  baseHistos["kindisc"]     = new TH1F("kindisc",     ";Kinematics discriminator;Jets",      100, -1, 1);

  //replicate histos per channel
  TString ch[]={"emu","ll","zll"};
  for(size_t i=0; i<sizeof(ch)/sizeof(TString); i++)
  {
    for(std::map<TString,TH1F *>::iterator it=baseHistos.begin(); it!=baseHistos.end(); it++)
    {
      TString tag=ch[i]+"_"+it->first;
      histos_[tag]=(TH1F *)it->second->Clone(tag);
      histos_[tag]->Sumw2();
      histos_[tag]->SetDirectory(outF_);
    }
  }

  histos_["puwgtnorm"] = new TH1F("puwgtnorm", ";puwgtnorm;Events",              4, 0, 4);
  histos_["puwgtnorm"]->Sumw2();
  histos_["puwgtnorm"]->SetDirectory(outF_);
}

//
Int_t TTbarEventAnalysis::processFile(TString inFile,TH1F *xsecWgt, Bool_t isData)
{
  cout << "xsecWgt = " << xsecWgt << endl;
  //loop over events
  TFile *inF=TFile::Open(inFile);
  TTree *tree=(TTree *)inF->Get("btagana/ttree");
  Int_t nentries=tree->GetEntriesFast();
  std::cout << "...opening " << inFile << " -> analysing " << nentries << " events -> " << outF_->GetName();
  std::cout << std::endl;

  if (nentries == 0){
    inF->Close();
    return 0;
  }

  //prepare reader
  std::vector<Float_t> tmvaVars( tmvaVarNames_.size(), 0. );
  if(weightsDir_!=""){
    tmvaReader_=new TMVA::Reader( "!Color:!Silent" );
    for(size_t ivar=0; ivar<tmvaVarNames_.size(); ivar++)
    tmvaReader_->AddVariable( tmvaVarNames_[ivar], &tmvaVars[ivar] );

    TString jranks[]={"leading",  "others",  "subleading" };
    for(size_t i=0; i<sizeof(jranks)/sizeof(TString); i++){
      tmvaReader_->BookMVA("BDT_"+jranks[i], weightsDir_+"/"+jranks[i]+"/TMVAClassification_BDT.weights.xml");
    }
  }

  //prepare to read the tree (for jets only interested in a couple of variables)
  struct MyEventInfoBranches_t
  {
    Int_t Run,Evt,LumiBlock,nPV;
    Int_t   ttbar_chan, ttbar_trigWord, ttbar_metfilterWord;
    Int_t   ttbar_nl, ttbar_lid[10], ttbar_lgid[10], ttbar_lch[10];
    Float_t ttbar_lpt[10], ttbar_leta[10], ttbar_lphi[10], ttbar_lm[10];
    Float_t ttbar_metpt,ttbar_metphi;
    Int_t   ttbar_nw;
    Int_t nPU;
    Float_t nPUtrue;
    Float_t ttbar_w[1095];
    Int_t nJet;
    Float_t Jet_pt[100],Jet_genpt[100],Jet_area[100],Jet_jes[100],Jet_eta[100],Jet_phi[100],Jet_mass[100];
    Float_t Jet_Svx[100],Jet_CombIVF[100],Jet_Proba[100],Jet_Ip2P[100];
    Float_t Jet_DeepCSVb[100], Jet_DeepCSVc[100], Jet_DeepCSVl[100], Jet_DeepCSVbN[100], Jet_DeepCSVcN[100], Jet_DeepCSVlN[100];
    Float_t Jet_DeepCSVBDisc[100],Jet_DeepCSVBDiscN[100],Jet_DeepCSVCvsLDisc[100],Jet_DeepCSVCvsLDiscN[100],Jet_DeepCSVCvsBDisc[100],Jet_DeepCSVCvsBDiscN[100];
    Float_t Jet_DeepFlavourBDisc[100], Jet_DeepFlavourCvsLDisc[100], Jet_DeepFlavourCvsBDisc[100];
    Float_t Jet_DeepFlavourB[100];

    Int_t Jet_nseltracks[100];
    Int_t Jet_flavour[100];
  };
  MyEventInfoBranches_t ev;
  tree->SetBranchAddress("Run"        , &ev.Run        );
  tree->SetBranchAddress("Evt"        , &ev.Evt        );
  tree->SetBranchAddress("LumiBlock"  , &ev.LumiBlock  );
  tree->SetBranchAddress("nPV"        , &ev.nPV        );
  tree->SetBranchAddress("nPU"        , &ev.nPU        );
  tree->SetBranchAddress("nPUtrue",     &ev.nPUtrue );
  tree->SetBranchAddress("ttbar_chan" , &ev.ttbar_chan);
  tree->SetBranchAddress("ttbar_metfilterWord", &ev.ttbar_metfilterWord);
  tree->SetBranchAddress("ttbar_trigWord", &ev.ttbar_trigWord);
  tree->SetBranchAddress("ttbar_nl"   ,  &ev.ttbar_nl);
  tree->SetBranchAddress("ttbar_lpt"  ,   ev.ttbar_lpt);
  tree->SetBranchAddress("ttbar_leta" ,   ev.ttbar_leta);
  tree->SetBranchAddress("ttbar_lphi" ,   ev.ttbar_lphi);
  tree->SetBranchAddress("ttbar_lm"   ,   ev.ttbar_lm);
  tree->SetBranchAddress("ttbar_lid"  ,   ev.ttbar_lid);
  tree->SetBranchAddress("ttbar_lgid" ,   ev.ttbar_lgid);
  tree->SetBranchAddress("ttbar_lch"  ,   ev.ttbar_lch);
  tree->SetBranchAddress("ttbar_metpt",  &ev.ttbar_metpt);
  tree->SetBranchAddress("ttbar_metphi", &ev.ttbar_metphi);
  tree->SetBranchAddress("ttbar_nw",     &ev.ttbar_nw);
  tree->SetBranchAddress("ttbar_w",      ev.ttbar_w);
  tree->SetBranchAddress("nJet",            &ev.nJet);
  tree->SetBranchAddress("Jet_pt",          ev.Jet_pt);
  if(!isData){
    tree->SetBranchAddress("Jet_genpt",       ev.Jet_genpt);
  }
  tree->SetBranchAddress("Jet_area",        ev.Jet_area);
  tree->SetBranchAddress("Jet_jes",         ev.Jet_jes);
  tree->SetBranchAddress("Jet_eta",         ev.Jet_eta);
  tree->SetBranchAddress("Jet_phi",         ev.Jet_phi);
  tree->SetBranchAddress("Jet_mass",        ev.Jet_mass);
  tree->SetBranchAddress("Jet_Svx",         ev.Jet_Svx);
  tree->SetBranchAddress("Jet_CombIVF",     ev.Jet_CombIVF);
  tree->SetBranchAddress("Jet_Proba",       ev.Jet_Proba);
  tree->SetBranchAddress("Jet_Ip2P",        ev.Jet_Ip2P);
  tree->SetBranchAddress("Jet_nseltracks",  ev.Jet_nseltracks);
  tree->SetBranchAddress("Jet_flavour",     ev.Jet_flavour);
  tree->SetBranchAddress("Jet_DeepCSVb",  ev.Jet_DeepCSVb);
  tree->SetBranchAddress("Jet_DeepCSVc",  ev.Jet_DeepCSVc);
  tree->SetBranchAddress("Jet_DeepCSVl",  ev.Jet_DeepCSVl);
  tree->SetBranchAddress("Jet_DeepCSVbN",  ev.Jet_DeepCSVbN);
  tree->SetBranchAddress("Jet_DeepCSVcN",  ev.Jet_DeepCSVcN);
  tree->SetBranchAddress("Jet_DeepCSVlN",  ev.Jet_DeepCSVlN);
  tree->SetBranchAddress("Jet_DeepCSVBDisc", ev.Jet_DeepCSVBDisc);
  tree->SetBranchAddress("Jet_DeepCSVBDiscN", ev.Jet_DeepCSVBDiscN);
  tree->SetBranchAddress("Jet_DeepCSVCvsLDisc", ev.Jet_DeepCSVCvsLDisc);
  tree->SetBranchAddress("Jet_DeepCSVCvsLDiscN", ev.Jet_DeepCSVCvsLDiscN);
  tree->SetBranchAddress("Jet_DeepCSVCvsBDisc", ev.Jet_DeepCSVCvsBDisc);
  tree->SetBranchAddress("Jet_DeepCSVCvsBDiscN", ev.Jet_DeepCSVCvsBDiscN);

  tree->SetBranchAddress("Jet_DeepFlavourBDisc", ev.Jet_DeepFlavourBDisc);
  tree->SetBranchAddress("Jet_DeepFlavourCvsLDisc", ev.Jet_DeepFlavourCvsLDisc);
  tree->SetBranchAddress("Jet_DeepFlavourCvsBDisc", ev.Jet_DeepFlavourCvsBDisc);

  tree->SetBranchAddress("Jet_DeepFlavourB", ev.Jet_DeepFlavourB);

  int Event_i = 0;

  cout << "nentries: " << nentries << endl;

  for(Int_t i=Event_i; i<nentries; i++){

    tree->GetEntry(i);

    //progress bar
    //if(i%100==0) std::cout << "\r[ " << int(100.*i/nentries) << "/100 ] to completion" << std::flush;

    //generator level weights
    Float_t genWgt=ev.ttbar_nw==0 ? 1.0 : ev.ttbar_w[0];
    Float_t qcdScaleLo(1.0),qcdScaleHi(1.0),hdampLo(1.0),hdampHi(1.0);
    double isrRedHi=1;
    double fsrRedHi=1;
    double isrRedLo=1;
    double fsrRedLo=1;
    double isrDefHi=1;
    double fsrDefHi=1;
    double isrDefLo=1;
    double fsrDefLo=1;
    double isrConHi=1;
    double fsrConHi=1;
    double isrConLo=1;
    double fsrConLo=1;
    if(readTTJetsGenWeights_ && ev.ttbar_nw>17){
      // Weight * [sum(weights for given systematic) / sum(weights nominal)
      // i.e. renormalise to same number of events as nominal.
      // N.B. for GetBinContent(X) 'X' must be weight index + 1 due to way histogram is filled (0th bin always underflow).
      qcdScaleLo=ev.ttbar_w[9]*(xsecWgt->GetBinContent(10)/xsecWgt->GetBinContent(1));
      qcdScaleHi=ev.ttbar_w[5]*(xsecWgt->GetBinContent(6)/xsecWgt->GetBinContent(1));
      //hdampLo=ev.ttbar_w[ev.ttbar_nw-17]*(xsecWgt->GetBinContent(ev.ttbar_nw-17+1)/xsecWgt->GetBinContent(1));
      //hdampHi=ev.ttbar_w[ev.ttbar_nw-9]*(xsecWgt->GetBinContent(ev.ttbar_nw-9+1)/xsecWgt->GetBinContent(1));
      hdampLo=ev.ttbar_w[ev.ttbar_nw-17]*(xsecWgt->GetBinContent(ev.ttbar_nw-17+1)/xsecWgt->GetBinContent(1));
      hdampHi=ev.ttbar_w[ev.ttbar_nw-9]*(xsecWgt->GetBinContent(ev.ttbar_nw-9+1)/xsecWgt->GetBinContent(1));

      // >>> PSWeights <<<
      // Vector of weight to be used instead of old ISR/FSR varied alternative samples.
      // First weight (weightID= 1081) corresponds to central ME weight value.
      // The remaining 12 values (weightIDs = 1082 to 1093) correspond to the PS weights in the following order (ISR up, FSR up, ISR down, FSR down) x 3 sets, i.e.
      // 1082 = isrRedHi isr:muRfac=0.707, 1083 = fsrRedHi fsr:muRfac=0.707, 1084 = isrRedLo isr:muRfac=1.414, 1085 = fsrRedLo fsr:muRfac=1.414,
      // 1086 = isrDefHi isr:muRfac=0.5,   1087 = fsrDefHi fsr:muRfac=0.5,   1088 = isrDefLo isr:muRfac=2.0,   1089 = fsrDefLo fsr:muRfac=2.0,
      // 1090 = isrConHi isr:muRfac=0.25,  1091 = fsrConHi fsr:muRfac=0.25,  1092 = isrConLo isr:muRfac=4.0,   1093 = fsrConLo fsr:muRfac=4.0


      /*std::cout << "--- ttbar_nw " << ev.ttbar_nw << " ---" << std::endl;
      cout << "--- number of bins in xsecWgt: " << xsecWgt.GetNbinsX() << "---" << endl;
      std::cout << "Nominal sum of weights (xsecWgt) = " << xsecWgt->GetBinContent(1) << std::endl;
      std::cout << "qcdScaleLo : " << ev.ttbar_w[9]*(xsecWgt->GetBinContent(10)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "qcdScaleHi : " << ev.ttbar_w[5]*(xsecWgt->GetBinContent(6)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "hdampLo : " << ev.ttbar_w[ev.ttbar_nw-31]*(xsecWgt->GetBinContent(ev.ttbar_nw-31+1)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "hdampHi : " << ev.ttbar_w[ev.ttbar_nw-23]*(xsecWgt->GetBinContent(ev.ttbar_nw-23+1)/xsecWgt->GetBinContent(1)) << std::endl;
      // Store generator weights
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "isrRedHi event weight = " << ev.ttbar_w[1082] << std::endl;
      std::cout << "isrRedHi sum of event weights = " << xsecWgt->GetBinContent(1083) << std::endl;
      std::cout << "isrRedHi renormalisation factor w.r.t. nominal = " << xsecWgt->GetBinContent(1083)/xsecWgt->GetBinContent(1) << std::endl;
      std::cout << "Stored isrRedHi weight = " << ev.ttbar_w[1082]*(xsecWgt->GetBinContent(1083)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "isrRedLo event weight = " << ev.ttbar_w[1084] << std::endl;
      std::cout << "isrRedLo sum of event weights = " << xsecWgt->GetBinContent(1083) << std::endl;
      std::cout << "isrRedLo renormalisation factor w.r.t. nominal = " << xsecWgt->GetBinContent(1085)/xsecWgt->GetBinContent(1) << std::endl;
      std::cout << "Stored isrRedLo weight = " << ev.ttbar_w[1084]*xsecWgt->GetBinContent(1085)/xsecWgt->GetBinContent(1) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "fsrRedHi event weight = " << ev.ttbar_w[1083] << std::endl;
      std::cout << "fsrRedHi sum of event weights = " << xsecWgt->GetBinContent(1084) << std::endl;
      std::cout << "fsrRedHi renormalisation factor w.r.t. nominal = " << xsecWgt->GetBinContent(1084)/xsecWgt->GetBinContent(1) << std::endl;
      std::cout << "Stored fsrRedHi weight =  " << ev.ttbar_w[1083]*(xsecWgt->GetBinContent(1084)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "fsrRedLo event weight = " << ev.ttbar_w[1085] << std::endl;
      std::cout << "fsrRedLo sum of event weights = " << xsecWgt->GetBinContent(1086) << std::endl;
      std::cout << "fsrRedLo renormaliation factor w.r.t. nominal = " << xsecWgt->GetBinContent(1086)/xsecWgt->GetBinContent(1) << std::endl;
      std::cout << "Stored fsrRedLo weight = " << ev.ttbar_w[1085]*(xsecWgt->GetBinContent(1086)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "isrDefHi: " << ev.ttbar_w[1086]*(xsecWgt->GetBinContent(1087)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "fsrDefHi: " << ev.ttbar_w[1087]*(xsecWgt->GetBinContent(1088)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "isrDefLo: " << ev.ttbar_w[1088]*(xsecWgt->GetBinContent(1089)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "fsrDefLo: " << ev.ttbar_w[1089]*(xsecWgt->GetBinContent(1090)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "isrConHi: " << ev.ttbar_w[1090]*(xsecWgt->GetBinContent(1091)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "fsrConHi: " << ev.ttbar_w[1091]*(xsecWgt->GetBinContent(1092)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "isrConLo: " << ev.ttbar_w[1092]*(xsecWgt->GetBinContent(1093)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "fsrConLo: " << ev.ttbar_w[1093]*(xsecWgt->GetBinContent(1094)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
      std::cout << "Anything else????: " << ev.ttbar_w[1094]*(xsecWgt->GetBinContent(1095)/xsecWgt->GetBinContent(1)) << std::endl;
      std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;*/

      isrRedHi = ev.ttbar_w[1082]*(xsecWgt->GetBinContent(1083)/xsecWgt->GetBinContent(1));
      fsrRedHi = ev.ttbar_w[1083]*(xsecWgt->GetBinContent(1084)/xsecWgt->GetBinContent(1));
      isrRedLo = ev.ttbar_w[1084]*(xsecWgt->GetBinContent(1085)/xsecWgt->GetBinContent(1));
      fsrRedLo = ev.ttbar_w[1085]*(xsecWgt->GetBinContent(1086)/xsecWgt->GetBinContent(1));
      isrDefHi = ev.ttbar_w[1086]*(xsecWgt->GetBinContent(1087)/xsecWgt->GetBinContent(1));
      fsrDefHi = ev.ttbar_w[1087]*(xsecWgt->GetBinContent(1088)/xsecWgt->GetBinContent(1));
      isrDefLo = ev.ttbar_w[1088]*(xsecWgt->GetBinContent(1089)/xsecWgt->GetBinContent(1));
      fsrDefLo = ev.ttbar_w[1089]*(xsecWgt->GetBinContent(1090)/xsecWgt->GetBinContent(1));
      isrConHi = ev.ttbar_w[1090]*(xsecWgt->GetBinContent(1091)/xsecWgt->GetBinContent(1));
      fsrConHi = ev.ttbar_w[1091]*(xsecWgt->GetBinContent(1092)/xsecWgt->GetBinContent(1));
      isrConLo = ev.ttbar_w[1092]*(xsecWgt->GetBinContent(1093)/xsecWgt->GetBinContent(1));
      fsrConLo = ev.ttbar_w[1093]*(xsecWgt->GetBinContent(1094)/xsecWgt->GetBinContent(1));
    }
    //pileup weights
    Float_t puWgtLo(1.0), puWgtNom(1.0), puWgtHi(1.0);
    if(!isData){
      if(puWgtGr_)     puWgtNom = puWgtGr_->Eval(ev.nPU);
      if(puWgtDownGr_) puWgtLo  = puWgtDownGr_->Eval(ev.nPU);
      if(puWgtUpGr_)   puWgtHi  = puWgtUpGr_->Eval(ev.nPU);
      if(puWgtNom<0)   {/*cout << "puWgtNom =: " << puWgtNom << " setting to 0" << endl; */puWgtNom = 0;}
      if(puWgtLo <0)   puWgtLo  = 0;
      if(puWgtHi <0)   puWgtHi  = 0;
    }
    histos_["puwgtnorm" ]->Fill(0.,1.0);
    histos_["puwgtnorm" ]->Fill(1.,puWgtNom);
    histos_["puwgtnorm" ]->Fill(2.,puWgtLo);
    histos_["puwgtnorm" ]->Fill(3.,puWgtHi);

    //
    //CHANNEL ASSIGNMENT
    //
    if(ev.ttbar_nl<2 || ev.nJet<2) continue;
    ev.ttbar_chan=ev.ttbar_lid[0]*ev.ttbar_lch[0]*ev.ttbar_lid[1]*ev.ttbar_lch[1];

    TString ch("");
    if(ev.ttbar_chan==-11*13) ch="emu";
    if(ev.ttbar_chan==-11*11 || ev.ttbar_chan==-13*13) ch="ll";
    if(ch=="") continue;

    //
    //TRIGGER
    //
    bool hasTrigger( triggerBits_.size()==0  ? true : false);
    for(size_t ibit=0; ibit<triggerBits_.size(); ibit++){
      if(triggerBits_[ibit].second!=ev.ttbar_chan) continue;
      hasTrigger |= ((ev.ttbar_trigWord>>triggerBits_[ibit].first) & 1);
    }
    if(!hasTrigger) continue;

    //trigger efficiency weight
    Float_t trigWgtLo(1.0), trigWgtNom(1.0), trigWgtHi(1.0);
    if(!isData){
      std::pair<float,float> eff=getTriggerEfficiency(ev.ttbar_lid[0],ev.ttbar_lpt[0],ev.ttbar_leta[0],
        ev.ttbar_lid[1],ev.ttbar_lpt[1],ev.ttbar_leta[1],
        ev.ttbar_chan);
        trigWgtLo=eff.first-eff.second;
        trigWgtNom=eff.first;
        trigWgtHi=eff.first+eff.second;
      }

      //lepton selection efficiency
      Float_t lepSelEffLo(1.0), lepSelEffNom(1.0), lepSelEffHi(1.0);
      if(!isData){
        for(size_t il=0; il<2; il++)
        {
          std::pair<float,float> lepSF = getLeptonSelectionEfficiencyScaleFactor(ev.ttbar_lid[il],ev.ttbar_lpt[il],ev.ttbar_leta[il]);
          lepSelEffLo  *= (lepSF.first-lepSF.second);
          lepSelEffNom *= lepSF.first;
          lepSelEffHi  *= (lepSF.first+lepSF.second);
        }
      }

      //dilepton invariant mass
      std::vector<TLorentzVector> lp4;
      for(Int_t il=0; il<ev.ttbar_nl; il++){
        lp4.push_back( TLorentzVector(0,0,0,0) );
        lp4[il].SetPtEtaPhiM(ev.ttbar_lpt[il],ev.ttbar_leta[il],ev.ttbar_lphi[il],0.);
      }

      TLorentzVector dilepton(lp4[0]+lp4[1]);
      Float_t mll=dilepton.M();
      if(mll<12) continue;
      if(mll<90) continue;
      if(lp4[0].Pt()<25 || lp4[1].Pt()<25) continue;

      //nominal event weight
      Float_t evWgt(1.0);
      if(!isData){
        evWgt *= puWgtNom*trigWgtNom*lepSelEffNom*genWgt;
        if(xsecWgt) {
          evWgt *= xsecWgt->GetBinContent(1);
        }
      }
      histos_[ch+"_npvinc"]->Fill(ev.nPV-1,evWgt);
      npv_=ev.nPV;

      //
      //JET/MET SELECTION
      //
      Int_t jetCount[5]={0,0,0,0,0};
      std::vector<Int_t> selJets;
      std::vector<std::vector<Float_t> > selJetsKINDisc;
      std::vector< std::vector<TLorentzVector> > selJetsP4;
      std::vector< std::vector< std::vector<LJKinematics_t> > > selJetsLJKinematics;
      for(Int_t ij=0; ij<ev.nJet; ij++){
        //convert to P4
        TLorentzVector jp4(0,0,0,0);
        jp4.SetPtEtaPhiM(ev.Jet_pt[ij],ev.Jet_eta[ij],ev.Jet_phi[ij],ev.Jet_mass[ij]);

        //cross clean jets wrt to leptons
        Float_t minDRlj(9999.);
        for(size_t il=0; il<2; il++) minDRlj = TMath::Min( (Float_t)minDRlj, (Float_t)lp4[il].DeltaR(jp4) );
        if(minDRlj<0.4) continue;

        //update jet energy scale/resolution
        Float_t jrawsf=1./ev.Jet_jes[ij];
        Float_t jarea=ev.Jet_area[ij];


        // update JES+JER for this jet
        std::vector<float> jesSF(3,1.0);
        jecUnc_->setJetEta(fabs(jp4.Eta()));
        jecUnc_->setJetPt(jp4.Pt());
        float unc = jecUnc_->getUncertainty(true);
        jesSF[1]=(1.+fabs(unc));
        jesSF[2]=(1.-fabs(unc));

        std::vector<float> jerSF(3);

        TLorentzVector oldjp4(jp4);
        if (!isData){
          Float_t genjpt=ev.Jet_genpt[ij];
          jerSF = getJetResolutionScales(jesSF[0]*jp4.Pt(), jp4.Eta(), genjpt);
          jp4 = jp4*jesSF[0]*jerSF[0];
        }
        else{
          jerSF[0] = 1.;
          jerSF[1] = 1.;
          jerSF[2] = 1.;
          jp4 = jp4*jesSF[0];
        }

        // apply energy shifts according to systematic variation
        Bool_t canBeSelected(false);
        std::vector<TLorentzVector> varjp4;
        std::vector<Float_t> varkindisc;
        std::vector< std::vector<LJKinematics_t> > varLJKinematics;
        for(Int_t iSystVar=0; iSystVar<5; iSystVar++){
          varjp4.push_back( jp4 );

          if(iSystVar==1) varjp4[iSystVar] *= jesSF[1]/jesSF[0];
          if(iSystVar==2) varjp4[iSystVar] *= jesSF[2]/jesSF[0];
          if(iSystVar==3) varjp4[iSystVar] *= jerSF[1]/jerSF[0];
          if(iSystVar==4) varjp4[iSystVar] *= jerSF[2]/jerSF[0];

          //prepare variables for MVA
          std::vector< LJKinematics_t > ljkinematics;
          for(Int_t il=0; il<2; il++)
          {
            LJKinematics_t iljkin;
            iljkin.dr         = lp4[il].DeltaR(varjp4[iSystVar]);
            iljkin.dphi       = fabs(lp4[il].DeltaPhi(varjp4[iSystVar]));
            iljkin.deta       = fabs(lp4[il].Eta()-varjp4[iSystVar].Eta());
            iljkin.ptrel      = ROOT::Math::VectorUtil::Perp(lp4[il].Vect(),varjp4[iSystVar].Vect().Unit())/lp4[il].P();
            TLorentzVector ljP4(lp4[il]+varjp4[iSystVar]);
            iljkin.mlj        = ljP4.M();
            iljkin.lj2ll_deta = fabs(ljP4.Eta()-dilepton.Eta());
            iljkin.lj2ll_dphi = fabs(ljP4.DeltaPhi(dilepton));
            ljkinematics.push_back(iljkin);
          }
          sort(ljkinematics.begin(),ljkinematics.end(),sortLJKinematicsByDR);
          varLJKinematics.push_back(ljkinematics);
          
          //evaluate the MVA
          Float_t kindisc(0.0);
          if(tmvaReader_)
          {
            for(size_t ivar=0; ivar<tmvaVarNames_.size(); ivar++){
              if( tmvaVarNames_[ivar].Contains("j2ll_") ){
                if(tmvaVarNames_[ivar]=="j2ll_deta") tmvaVars[ivar]=fabs(varjp4[iSystVar].Eta()-dilepton.Eta());
                if(tmvaVarNames_[ivar]=="j2ll_phi")  tmvaVars[ivar]=fabs(varjp4[iSystVar].DeltaPhi(dilepton));
              }
              else{
                int ljidx( tmvaVarNames_[ivar].Contains("close") ? 0 : 1);
                if(tmvaVarNames_[ivar].Contains("_dr"))    tmvaVars[ivar]=ljkinematics[ljidx].dr;
                if(tmvaVarNames_[ivar].Contains("_dphi")){
                  if(tmvaVarNames_[ivar].Contains("lj2ll_")) tmvaVars[ivar]=ljkinematics[ljidx].lj2ll_dphi;
                  else                                       tmvaVars[ivar]=ljkinematics[ljidx].dphi;
                }
                if(tmvaVarNames_[ivar].Contains("_deta")){
                  if(tmvaVarNames_[ivar].Contains("lj2ll_")) tmvaVars[ivar]=ljkinematics[ljidx].lj2ll_deta;
                  else                                       tmvaVars[ivar]=ljkinematics[ljidx].deta;
                }
                if(tmvaVarNames_[ivar].Contains("_ptrel")) tmvaVars[ivar]=ljkinematics[ljidx].ptrel;
                if(tmvaVarNames_[ivar].Contains("_mlj"))   tmvaVars[ivar]=ljkinematics[ljidx].mlj;
              }
            }

            TString methodPFix("_others");
            if(selJets.size()==0) methodPFix="_leading";
            else if(selJets.size()==1) methodPFix="_subleading";
            varkindisc.push_back( tmvaReader_->EvaluateMVA("BDT"+methodPFix) );
          }

          //check if can be selected for this variation
          if(varjp4[iSystVar].Pt()<30 || TMath::Abs(varjp4[iSystVar].Eta())>2.5) continue;
          canBeSelected=true;
          jetCount[iSystVar]++;
        }

        //add jet if it is selectable
        if(!canBeSelected) continue;
        selJets.push_back(ij);
        selJetsP4.push_back(varjp4);
        selJetsLJKinematics.push_back( varLJKinematics );
        if(tmvaReader_) selJetsKINDisc.push_back(varkindisc);
      }

      //
      // base selection and n-1 plots
      //
      bool zCand( (ch.Contains("ll") && TMath::Abs(mll-91)<15) ? true : false );
      bool passMet( (ch.Contains("emu") || ev.ttbar_metpt>40) ?  true : false);
      bool passJets(selJets.size()>=2 ? true : false);

      if(!passJets) continue;
      if(zCand)  ch="z"+ch;
      histos_[ch+"_mllinc"]->Fill(mll,evWgt);
      histos_[ch+"_met"]->Fill(ev.ttbar_metpt,evWgt);

      if(!passMet) continue;
      histos_[ch+"_evsel"]->Fill(0.,evWgt);
      if(selJets.size()<5)   histos_[ch+"_evsel"]->Fill(selJets.size()-1,evWgt);
      histos_[ch+"_npv"]->Fill(ev.nPV-1,evWgt);
      histos_[ch+"_mll"]->Fill(mll,evWgt);
      histos_[ch+"_njets"]->Fill(selJets.size(),evWgt);
      histos_[ch+"_leadjpt"]->Fill(selJetsP4[0][0].Pt(),evWgt);
      histos_[ch+"_leadjeta"]->Fill((selJetsP4[0][0].Eta()),evWgt);
      histos_[ch+"_leadlpt"]->Fill(lp4[0].Pt(),evWgt);
      histos_[ch+"_trailjpt"]->Fill(selJetsP4[1][0].Pt(),evWgt);
      histos_[ch+"_trailjeta"]->Fill(fabs(selJetsP4[1][0].Eta()),evWgt);
      histos_[ch+"_traillpt"]->Fill(lp4[1].Pt(),evWgt);

      std::vector<float> leadingkindisc(2,-9999);
      std::vector<int> leadingkindiscIdx(2,-1);
      for(size_t ij=0; ij<selJets.size(); ij++){
        Int_t jetIdx(selJets[ij]);
        histos_[ch+"_close_mlj"]->Fill(selJetsLJKinematics[ij][0][0].mlj,evWgt);
        histos_[ch+"_close_deta"]->Fill(selJetsLJKinematics[ij][0][0].deta,evWgt);
        histos_[ch+"_close_dphi"]->Fill(selJetsLJKinematics[ij][0][0].dphi,evWgt);
        histos_[ch+"_close_ptrel"]->Fill(selJetsLJKinematics[ij][0][0].ptrel,evWgt);
        histos_[ch+"_close_lj2ll_deta"]->Fill(selJetsLJKinematics[ij][0][0].lj2ll_deta,evWgt);
        histos_[ch+"_close_lj2ll_dphi"]->Fill(selJetsLJKinematics[ij][0][0].lj2ll_dphi,evWgt);
        histos_[ch+"_far_mlj"]->Fill(selJetsLJKinematics[ij][0][1].mlj,evWgt);
        histos_[ch+"_far_deta"]->Fill(selJetsLJKinematics[ij][0][1].deta,evWgt);
        histos_[ch+"_far_dphi"]->Fill(selJetsLJKinematics[ij][0][1].dphi,evWgt);
        histos_[ch+"_far_ptrel"]->Fill(selJetsLJKinematics[ij][0][1].ptrel,evWgt);
        histos_[ch+"_far_lj2ll_deta"]->Fill(selJetsLJKinematics[ij][0][1].lj2ll_deta,evWgt);
        histos_[ch+"_far_lj2ll_dphi"]->Fill(selJetsLJKinematics[ij][0][1].lj2ll_dphi,evWgt);
        histos_[ch+"_j2ll_deta"]->Fill(fabs(selJetsP4[ij][0].Eta()-dilepton.Eta()),evWgt);
        histos_[ch+"_j2ll_dphi"]->Fill(fabs(selJetsP4[ij][0].DeltaPhi(dilepton)),evWgt);
        if(tmvaReader_) histos_[ch+"_kindisc"]->Fill(selJetsKINDisc[ij][0],evWgt);
        histos_[ch+"_jp"]->Fill(ev.Jet_Proba[jetIdx],evWgt);
        histos_[ch+"_svhe"]->Fill(ev.Jet_Svx[jetIdx],evWgt);
        histos_[ch+"_csv"]->Fill(ev.Jet_CombIVF[jetIdx],evWgt);
        histos_[ch+"_tche"]->Fill(ev.Jet_Ip2P[jetIdx],evWgt);
        histos_[ch+"_jetseltrk"]->Fill(ev.Jet_nseltracks[jetIdx],evWgt);
        histos_[ch+"_DeepCSVb"]->Fill(ev.Jet_DeepCSVb[jetIdx],evWgt);
        histos_[ch+"_DeepCSVc"]->Fill(ev.Jet_DeepCSVc[jetIdx],evWgt);
        histos_[ch+"_DeepCSVl"]->Fill(ev.Jet_DeepCSVl[jetIdx],evWgt);
        histos_[ch+"_DeepCSVbN"]->Fill(ev.Jet_DeepCSVbN[jetIdx],evWgt);
        histos_[ch+"_DeepCSVcN"]->Fill(ev.Jet_DeepCSVcN[jetIdx],evWgt);
        histos_[ch+"_DeepCSVlN"]->Fill(ev.Jet_DeepCSVlN[jetIdx],evWgt);
        histos_[ch+"_DeepCSVBDisc"]->Fill(ev.Jet_DeepCSVBDisc[jetIdx],evWgt);
        histos_[ch+"_DeepCSVBDiscN"]->Fill(ev.Jet_DeepCSVBDiscN[jetIdx],evWgt);
        histos_[ch+"_DeepCSVCvsLDisc"]->Fill(ev.Jet_DeepCSVCvsLDisc[jetIdx],evWgt);
        histos_[ch+"_DeepCSVCvsLDiscN"]->Fill(ev.Jet_DeepCSVCvsLDiscN[jetIdx],evWgt);
        histos_[ch+"_DeepCSVCvsBDisc"]->Fill(ev.Jet_DeepCSVCvsBDisc[jetIdx],evWgt);
        histos_[ch+"_DeepCSVCvsBDiscN"]->Fill(ev.Jet_DeepCSVCvsBDiscN[jetIdx],evWgt);

        histos_[ch+"_DeepFlavourBDisc"]->Fill(ev.Jet_DeepFlavourBDisc[jetIdx],evWgt);
        histos_[ch+"_DeepFlavourCvsLDisc"]->Fill(ev.Jet_DeepFlavourCvsLDisc[jetIdx],evWgt);
        histos_[ch+"_DeepFlavourCvsBDisc"]->Fill(ev.Jet_DeepFlavourCvsBDisc[jetIdx],evWgt);
        histos_[ch+"_DeepFlavourB"]->Fill(ev.Jet_DeepFlavourB[jetIdx],evWgt);


        Int_t flavBin(0),partonFlav(abs(ev.Jet_flavour[jetIdx]));
        if(partonFlav==21 || (partonFlav>0 && partonFlav<4)) flavBin=1;
        if(partonFlav==4) flavBin=2;
        if(partonFlav==5) flavBin=3;
        histos_[ch+"_flavour"]->Fill(flavBin,evWgt);

        //rank jets by kinematics discriminator
        if(tmvaReader_)
        {
          if(selJetsKINDisc[ij][0]>leadingkindisc[0])
          {
            leadingkindisc[1]=leadingkindisc[0];     leadingkindiscIdx[1]=leadingkindiscIdx[0];
            leadingkindisc[0]=selJetsKINDisc[ij][0]; leadingkindiscIdx[0]=jetIdx;
          }
          else if(selJetsKINDisc[ij][0]>leadingkindisc[1])
          {
            leadingkindisc[1]=selJetsKINDisc[ij][0]; leadingkindiscIdx[1]=jetIdx;
          }
        }
      }

      //control b-tagging quantities for the most promising jets in the KIN discriminator
      if(tmvaReader_)
      {
        for(size_t ij=0; ij<2; ij++)
        {
          size_t jetIdx=leadingkindiscIdx[ij];
          histos_[ch+"_jp_leadkin"]->Fill(ev.Jet_Proba[jetIdx],evWgt);
          histos_[ch+"_svhe_leadkin"]->Fill(ev.Jet_Svx[jetIdx],evWgt);
          histos_[ch+"_csv_leadkin"]->Fill(ev.Jet_CombIVF[jetIdx],evWgt);
          histos_[ch+"_tche_leadkin"]->Fill(ev.Jet_Ip2P[jetIdx],evWgt);
          histos_[ch+"_jetseltrk_leadkin"]->Fill(ev.Jet_nseltracks[jetIdx],evWgt);
          histos_[ch+"_DeepCSVb_leadkin"]->Fill(ev.Jet_DeepCSVb[jetIdx],evWgt);
          histos_[ch+"_DeepCSVc_leadkin"]->Fill(ev.Jet_DeepCSVc[jetIdx],evWgt);
          histos_[ch+"_DeepCSVl_leadkin"]->Fill(ev.Jet_DeepCSVl[jetIdx],evWgt);
          histos_[ch+"_DeepCSVbN_leadkin"]->Fill(ev.Jet_DeepCSVbN[jetIdx],evWgt);
          histos_[ch+"_DeepCSVcN_leadkin"]->Fill(ev.Jet_DeepCSVcN[jetIdx],evWgt);
          histos_[ch+"_DeepCSVlN_leadkin"]->Fill(ev.Jet_DeepCSVlN[jetIdx],evWgt);
          histos_[ch+"_DeepCSVBDisc_leadkin"]->Fill(ev.Jet_DeepCSVBDisc[jetIdx],evWgt);
          histos_[ch+"_DeepCSVBDiscN_leadkin"]->Fill(ev.Jet_DeepCSVBDiscN[jetIdx],evWgt);
          histos_[ch+"_DeepCSVCvsLDisc_leadkin"]->Fill(ev.Jet_DeepCSVCvsLDisc[jetIdx],evWgt);
          histos_[ch+"_DeepCSVCvsLDiscN_leadkin"]->Fill(ev.Jet_DeepCSVCvsLDiscN[jetIdx],evWgt);
          histos_[ch+"_DeepCSVCvsBDisc_leadkin"]->Fill(ev.Jet_DeepCSVCvsBDisc[jetIdx],evWgt);
          histos_[ch+"_DeepCSVCvsBDiscN_leadkin"]->Fill(ev.Jet_DeepCSVCvsBDiscN[jetIdx],evWgt);
          histos_[ch+"_DeepFlavourBDisc_leadkin"]->Fill(ev.Jet_DeepFlavourBDisc[jetIdx],evWgt);
          histos_[ch+"_DeepFlavourCvsLDisc_leadkin"]->Fill(ev.Jet_DeepFlavourCvsLDisc[jetIdx],evWgt);
          histos_[ch+"_DeepFlavourCvsBDisc_leadkin"]->Fill(ev.Jet_DeepFlavourCvsBDisc[jetIdx],evWgt);
          histos_[ch+"_DeepFlavourB_leadkin"]->Fill(ev.Jet_DeepFlavourB[jetIdx],evWgt);
        }
      }

      //
      //prepare to store trees
      //
      eventInfo_[0]=ev.Run;
      eventInfo_[1]=ev.Evt;
      eventInfo_[2]=ev.LumiBlock;

      jetmult_=selJets.size();
      ttbar_chan_=ev.ttbar_chan;
      if(zCand) ttbar_chan_+= 230000;

      //weights for systematic uncertainties
      for(Int_t iSystVar=0; iSystVar<5; iSystVar++)
      {
        Float_t selWeight(jetCount[iSystVar]>=2 ? 1.0 : 0.0);
        weight_[iSystVar]=evWgt*selWeight;
      }
      weight_[5] = puWgtNom>0 ? evWgt*puWgtLo/puWgtNom : evWgt;
      weight_[6] = puWgtLo>0  ? evWgt*puWgtHi/puWgtNom : evWgt;
      weight_[7] = evWgt*trigWgtLo/trigWgtNom;
      weight_[8] = evWgt*trigWgtHi/trigWgtNom;
      weight_[9] = evWgt*lepSelEffLo/lepSelEffNom;
      weight_[10]= evWgt*lepSelEffHi/lepSelEffNom;
      weight_[11]= evWgt*qcdScaleLo/genWgt;
      weight_[12]= evWgt*qcdScaleHi/genWgt;
      weight_[13]= evWgt*hdampLo/genWgt;
      weight_[14]= evWgt*hdampHi/genWgt;
      weight_[15]= evWgt*isrRedLo/genWgt;
      weight_[16]= evWgt*isrRedHi/genWgt;
      weight_[17]= evWgt*fsrRedLo/genWgt;
      weight_[18]= evWgt*fsrRedHi/genWgt;
      weight_[19]= evWgt*isrDefLo/genWgt;
      weight_[20]= evWgt*isrDefHi/genWgt;
      weight_[21]= evWgt*fsrDefLo/genWgt;
      weight_[22]= evWgt*fsrDefHi/genWgt;
      weight_[23]= evWgt*isrConLo/genWgt;
      weight_[24]= evWgt*isrConHi/genWgt;
      weight_[25]= evWgt*fsrConLo/genWgt;
      weight_[26]= evWgt*fsrConHi/genWgt;

      //fill trees
      for(size_t ij=0; ij<selJets.size(); ij++)
      {
        Int_t jetIdx(selJets[ij]);

        jetrank_ = ij;
        jetFlavour_[0] = ev.Jet_flavour[jetIdx];
        jetPt_[0]      = selJetsP4[ij][0].Pt();
        jetEta_[0]     = selJetsP4[ij][0].Eta();
        for(size_t iSystVar=0; iSystVar<5; iSystVar++)
        {
          close_mlj_[iSystVar] = selJetsLJKinematics[ij][iSystVar][0].mlj;
          if(tmvaReader_) kinDisc_[iSystVar] = selJetsKINDisc[ij][iSystVar];
          else            kinDisc_[iSystVar] = -999;
        }
        close_deta_ =selJetsLJKinematics[ij][0][0].deta;
        close_dphi_ =selJetsLJKinematics[ij][0][0].dphi;
        close_ptrel_=selJetsLJKinematics[ij][0][0].ptrel;
        close_lj2ll_deta_ = selJetsLJKinematics[ij][0][0].lj2ll_deta;
        close_lj2ll_dphi_ = selJetsLJKinematics[ij][0][0].lj2ll_dphi;

        far_mlj_    =selJetsLJKinematics[ij][0][1].mlj;
        far_deta_   =selJetsLJKinematics[ij][0][1].deta;
        far_dphi_   =selJetsLJKinematics[ij][0][1].dphi;
        far_ptrel_  =selJetsLJKinematics[ij][0][1].ptrel;
        far_lj2ll_deta_ = selJetsLJKinematics[ij][0][1].lj2ll_deta;
        far_lj2ll_dphi_ = selJetsLJKinematics[ij][0][1].lj2ll_dphi;

        j2ll_deta_  =fabs(selJetsP4[ij][0].Eta()-dilepton.Eta());
        j2ll_dphi_  =fabs(selJetsP4[ij][0].DeltaPhi( dilepton ));

        jp_[0]=ev.Jet_Proba[jetIdx];
        svhe_[0]=ev.Jet_Svx[jetIdx];
        csv_[0]=ev.Jet_CombIVF[jetIdx];
        DeepCSVb_[0]=ev.Jet_DeepCSVb[jetIdx];
        DeepCSVc_[0]=ev.Jet_DeepCSVc[jetIdx];
        DeepCSVl_[0]=ev.Jet_DeepCSVl[jetIdx];
        DeepCSVbN_[0]=ev.Jet_DeepCSVbN[jetIdx];
        DeepCSVcN_[0]=ev.Jet_DeepCSVcN[jetIdx];
        DeepCSVlN_[0]=ev.Jet_DeepCSVlN[jetIdx];

        DeepCSVBDisc_[0]=ev.Jet_DeepCSVBDisc[jetIdx];
        DeepCSVBDiscN_[0]=ev.Jet_DeepCSVBDiscN[jetIdx];
        DeepCSVCvsLDisc_[0]=ev.Jet_DeepCSVCvsLDisc[jetIdx];
        DeepCSVCvsLDiscN_[0]=ev.Jet_DeepCSVCvsLDiscN[jetIdx];
        DeepCSVCvsBDisc_[0]=ev.Jet_DeepCSVCvsBDisc[jetIdx];
        DeepCSVCvsBDiscN_[0]=ev.Jet_DeepCSVCvsBDiscN[jetIdx];

        DeepFlavourBDisc_[0]=ev.Jet_DeepFlavourBDisc[jetIdx];
        DeepFlavourCvsLDisc_[0]=ev.Jet_DeepFlavourCvsLDisc[jetIdx];
        DeepFlavourCvsBDisc_[0]=ev.Jet_DeepFlavourCvsBDisc[jetIdx];

        DeepFlavourB_[0]=ev.Jet_DeepFlavourB[jetIdx];

        kinTree_->Fill();
      }

      //FtM tree is filled with the two leading jets in the KIN discriminator
      if(tmvaReader_)
      {
        for(size_t ij=0; ij<2; ij++)
        {
          size_t jetIdx=leadingkindiscIdx[ij];
          jetFlavour_[ij] = ev.Jet_flavour[jetIdx];
          jetPt_[ij]      = ev.Jet_pt[jetIdx];
          jetEta_[ij]     = ev.Jet_eta[jetIdx];
          jp_[ij]         = ev.Jet_Proba[jetIdx];
          svhe_[ij]       = ev.Jet_Svx[jetIdx];
          csv_[ij]        = ev.Jet_CombIVF[jetIdx];
          kinDisc_[ij]    = leadingkindisc[ij];
          DeepCSVb_[ij] = ev.Jet_DeepCSVb[jetIdx];
          DeepCSVc_[ij] = ev.Jet_DeepCSVc[jetIdx];
          DeepCSVl_[ij] = ev.Jet_DeepCSVl[jetIdx];
          DeepCSVbN_[ij] = ev.Jet_DeepCSVbN[jetIdx];
          DeepCSVcN_[ij] = ev.Jet_DeepCSVcN[jetIdx];
          DeepCSVlN_[ij] = ev.Jet_DeepCSVlN[jetIdx];

          DeepCSVBDisc_[ij] = ev.Jet_DeepCSVBDisc[jetIdx];
          DeepCSVBDiscN_[ij] = ev.Jet_DeepCSVBDiscN[jetIdx];
          DeepCSVCvsLDisc_[ij] = ev.Jet_DeepCSVCvsLDisc[jetIdx];
          DeepCSVCvsLDiscN_[ij] = ev.Jet_DeepCSVCvsLDiscN[jetIdx];
          DeepCSVCvsBDisc_[ij] = ev.Jet_DeepCSVCvsBDisc[jetIdx];
          DeepCSVCvsBDiscN_[ij] = ev.Jet_DeepCSVCvsBDiscN[jetIdx];

          DeepFlavourBDisc_[ij] = ev.Jet_DeepFlavourBDisc[jetIdx];
          DeepFlavourCvsLDisc_[ij] = ev.Jet_DeepFlavourCvsLDisc[jetIdx];
          DeepFlavourCvsBDisc_[ij] = ev.Jet_DeepFlavourCvsBDisc[jetIdx];
          DeepFlavourB_[ij] = ev.Jet_DeepFlavourB[jetIdx];
          //std::cout << ij << " " <<  jetFlavour_[ij] << " "
          //<< jetPt_[ij] << " " << csv_[ij]  << " " << kinDisc_[ij]
          //		<< std::endl;
        }
        ftmTree_->Fill();
      }
    }
    //all done with this file
    inF->Close();
    return 1;
  }


std::pair<float,float> TTbarEventAnalysis::getTriggerEfficiency(int id1,float pt1,float eta1,int id2,float pt2,float eta2,int ch)
{
  std::pair<float,float>res(1.0,0.0);
  if(ch==-11*13) { res.first=1.0; res.second=0.05; }
  if(ch==-11*11) { res.first=1.0; res.second=0.05; }
  if(ch==-13*13) { res.first=1.0; res.second=0.05; } // this is a single muon trigger
 return res;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#Electron_Reconstruction_Scale_Fa
//https://soffi.web.cern.ch/soffi/EGM-ID/SF-17Nov2017-MCv2-IDv1-020618/Electrons/egammaEffi.txt_egammaPlots_runBCDEF_passingMedium94X.pdf
std::pair<float,float> TTbarEventAnalysis::getLeptonSelectionEfficiencyScaleFactor(int id,float pt,float eta)
{
  std::pair<float,float>res(1.0,0.0);

  //electrons
  if(abs(id)==11)
  {
    if (2.0 < eta < 2.5)
    {
      if (10.0<pt<20.0)      { res.first=1.034; res.second=0.028; }
      else if (20.0<pt<35.0) { res.first=0.926; res.second=0.018; }
      else if (35.0<pt<50.0) { res.first=0.942; res.second=0.005; }
      else if (50.0<pt<100.0) { res.first=0.954; res.second=0.010; }
      else if (100.0<pt<200.0) { res.first=1.003; res.second=0.0015; }
      else if (200.0<pt<500.0) { res.first=1.007; res.second=0.059; }
    }
    else if (1.566 < eta < 2.0)
    {
      if (10.0<pt<20.0)      { res.first=0.996; res.second=0.019; }
      else if (20.0<pt<35.0) { res.first=0.935; res.second=0.024; }
      else if (35.0<pt<50.0) { res.first=0.957; res.second=0.003; }
      else if (50.0<pt<100.0) { res.first=0.971; res.second=0.013; }
      else if (100.0<pt<200.0) { res.first=1.009; res.second=0.021; }
      else if (200.0<pt<500.0) { res.first=0.998; res.second=0.022; }
    }
    else if (0.8 < eta < 1.444)
    {
      if (10.0<pt<20.0)      { res.first=1.042; res.second=0.0143; }
      else if (20.0<pt<35.0) { res.first=0.962; res.second=0.002; }
      else if (35.0<pt<50.0) { res.first=0.964; res.second=0.002; }
      else if (50.0<pt<100.0) { res.first=0.970; res.second=0.002; }
      else if (100.0<pt<200.0) { res.first=0.993; res.second=0.011; }
      else if (200.0<pt<500.0) { res.first=1.018; res.second=0.022; }
    }
    else if (0.0 < eta < 0.8)
    {
      if (10.0<pt<20.0)      { res.first=1.053; res.second=0.029; }
      else if (20.0<pt<35.0) { res.first=0.963; res.second=0.017; }
      else if (35.0<pt<50.0) { res.first=0.966; res.second=0.0025; }
      else if (50.0<pt<100.0) { res.first=0.970; res.second=0.0074; }
      else if (100.0<pt<200.0) { res.first=0.998; res.second=0.0065; }
      else if (200.0<pt<500.0) { res.first=0.990; res.second=0.0030; }
    }
    else if (-0.8 < eta < 0.0)
    {
      if (10.0<pt<20.0)      { res.first=1.053; res.second=0.029; }
      else if (20.0<pt<35.0) { res.first=0.958; res.second=0.017; }
      else if (35.0<pt<50.0) { res.first=0.963; res.second=0.0027; }
      else if (50.0<pt<100.0) { res.first=0.966; res.second=0.0075; }
      else if (100.0<pt<200.0) { res.first=0.993; res.second=0.0065; }
      else if (200.0<pt<500.0) { res.first=0.991; res.second=0.031; }
    }
    else if (-1.444 < eta < -0.8)
    {
      if (10.0<pt<20.0)      { res.first= 1.034; res.second= 0.015; }
      else if (20.0<pt<35.0) { res.first= 0.960; res.second= 0.024; }
      else if (35.0<pt<50.0) { res.first= 0.959; res.second= 0.0022; }
      else if (50.0<pt<100.0) { res.first= 0.963; res.second= 0.0099; }
      else if (100.0<pt<200.0) { res.first= 0.998; res.second= 0.015; }
      else if (200.0<pt<500.0) { res.first= 0.999; res.second= 0.024; }
    }
    else if (-2.0 < eta < -1.566)
    {
      if (10.0<pt<20.0)      { res.first= 0.966; res.second= 0.019; }
      else if (20.0<pt<35.0) { res.first= 0.936; res.second= 0.024; }
      else if (35.0<pt<50.0) { res.first= 0.959; res.second= 0.0027; }
      else if (50.0<pt<100.0) { res.first= 0.955; res.second= 0.013; }
      else if (100.0<pt<200.0) { res.first= 0.994; res.second= 0.020; }
      else if (200.0<pt<500.0) { res.first= 0.988; res.second= 0.022; }
    }
    else if (-2.5 < eta < -2.0)
    {
      if (10.0<pt<20.0)      { res.first= 1.005; res.second= 0.029; }
      else if (20.0<pt<35.0) { res.first= 0.926; res.second= 0.014; }
      else if (35.0<pt<50.0) { res.first= 0.941; res.second= 0.0047; }
      else if (50.0<pt<100.0) { res.first= 0.955; res.second= 0.010; }
      else if (100.0<pt<200.0) { res.first= 0.986; res.second= 0.015; }
      else if (200.0<pt<500.0) { res.first= 0.980; res.second= 0.063; }
    }
    else if (1.444 < fabs(eta) < 2.566)
    {
        cout << "Electron should not be in calorimeter crack!" << endl;
    }
  }

  //muons
  // SFs and Uncs for NUM_TightID_DEN_genTracks
  // taken from https://twiki.cern.ch/twiki/pub/CMS/MuonReferenceEffs2017/RunBCDEF_SF_ID.json
  if (abs(id)==13)
  {
    if (fabs(eta)<0.9)
    {
      if (pt<25)      { res.first=0.9910777627756951; res.second=0.0034967203087024274; }
      else if (pt<30) { res.first=0.987410468262084; res.second=0.0018975082960536634; }
      else if (pt<40) { res.first=0.9907753279135898; res.second=0.0003977503197109705; }
      else if (pt<50) { res.first=0.9892483588952047; res.second=0.00032329941312374114; }
      else if (pt<60) { res.first=0.9855545160334763; res.second=0.0008602730194340781; }
      else            { res.first=0.9898057377093389; res.second=0.0016204327419630327; }
    }
    else if(fabs(eta)<1.2)
    {
      if (pt<25)      { res.first=0.9927389275515244; res.second=0.004743737066128366; }
      else if (pt<30) { res.first=0.985063939762512; res.second=0.02187876084241925; }
      else if (pt<40) { res.first=0.9865359464182247; res.second=0.0006903526352667042; }
      else if (pt<50) { res.first=0.984913093101493;  res.second=0.02013503091494561; }
      else if (pt<60) { res.first=0.9839056384760008; res.second=0.00159172336926836; }
      else            { res.first=0.984060403143468; res.second=0.012127878049119219; }
    }
    else if(fabs(eta)<2.1)
    {
      if (pt<25)      { res.first=0.9924252719877384; res.second=0.007785952740242002; }
      else if (pt<30) { res.first=0.9890884461284933; res.second=0.014905356047753367; }
      else if (pt<40) { res.first=0.9946469069883841; res.second=0.012424223689919973; }
      else if (pt<50) { res.first=0.9926528825155183; res.second=0.009931678115496967; }
      else if (pt<60) { res.first=0.9906364222943529; res.second=0.000971321379850273; }
      else            { res.first=0.9920464322143979; res.second=0.0021353964567237746; }
    }
    else
    {
      if (pt<25)      { res.first=0.9758095839531763; res.second=0.004399315121784104; }
      else if (pt<30) { res.first=0.9745153594179884;    res.second=0.0027111009825340473; }
      else if (pt<40) { res.first=0.9787410500158746; res.second=0.001003557787201416; }
      else if (pt<50) { res.first=0.978189122919501; res.second=0.0011230605941385397; }
      else if (pt<60) { res.first=0.9673568416097894; res.second=0.0037006525169638958; }
      else            { res.first=0.9766311856731202; res.second=0.00862666466882855; }
    }
  }
  return res;
}

//Sources
//  Assuming nominal JER but uncertainties from Run I
//  https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
std::vector<float> TTbarEventAnalysis::getJetResolutionScales(float pt, float eta, float genjpt)
{
  std::vector<float> res(3,1.0);

  float ptSF(1.0), ptSF_err(0.06);
  if(TMath::Abs(eta)<0.522){
    //ptSF=1.1432;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2));
    //ptSF_err =  0.0222;
  }
  else if(TMath::Abs(eta)<0.783){
    //ptSF=1.1815;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2));
    //ptSF_err = 0.0484;
  }
  else if(TMath::Abs(eta)<1.131){
    //ptSF=1.0989;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2));
    //ptSF_err = 0.0456;
  }
  else if(TMath::Abs(eta)<1.305){
    //ptSF=1.1137;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.1397;
  }
  else if(TMath::Abs(eta)<1.740){
    //ptSF=1.1307;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.1470;
  }
  else if(TMath::Abs(eta)<1.930){
    //ptSF=1.1600;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.0976;
  }
  else if(TMath::Abs(eta)<2.043){
    //ptSF=1.2393;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.1909;
  }
  else if(TMath::Abs(eta)<2.322){
    //ptSF=1.2604;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.1501;
  }
  else if(TMath::Abs(eta)<2.500){
    //ptSF=1.4085;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.2020;
  }
  else if(TMath::Abs(eta)<2.853){
    //ptSF=1.9909;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.5684;
  }
  else if(TMath::Abs(eta)<2.964){
    //ptSF=2.2923;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.3743;
  }
  else if(TMath::Abs(eta)<3.139){
    //ptSF=1.2696;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2));
    //ptSF_err = 0.1089;
  }
  else{
    //ptSF=1.1542;
    ptSF=1.0;
    ptSF_err = TMath::Sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2));
    //ptSF_err = 0.1524;
  }

  res[0] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF)*(pt-genjpt)))/pt;
  res[1] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF-ptSF_err)*(pt-genjpt)))/pt;
  res[2] = TMath::Max((Float_t)0.,(Float_t)(genjpt+(ptSF+ptSF_err)*(pt-genjpt)))/pt;

  return res;
}


//
void TTbarEventAnalysis::finalizeOutput()
{
  //dump results to file
  outF_->cd();

  //pileup weighting screws up a bit normalization - fix it a posteriori
  float puwgtSF(histos_["puwgtnorm" ]->GetBinContent(1)/histos_["puwgtnorm" ]->GetBinContent(2));

  for(std::map<TString,TH1F *>::iterator it = histos_.begin(); it != histos_.end(); it++){
    if(it->first!="puwgtnorm")
    it->second->Scale(puwgtSF);
    it->second->Write();
  }
  kinTree_->Write();
  if(ftmTree_) ftmTree_->Write();
  outF_->Close();
}
