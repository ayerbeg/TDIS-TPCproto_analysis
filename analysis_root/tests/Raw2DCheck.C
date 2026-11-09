#include "../TPC.C"
#include <TKey.h>
#include <chrono>
// End-to-end fixture: mapping, hit cuts, exact upper edges, both-axis flow
// corners, saved canvases and separate cleaned output. No real run is edited.
void Raw2DCheck() {
    using namespace PersonalTPC;
    auto folder=fs::temp_directory_path()/("tpc-raw2d-check-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder);
    struct Cleanup {fs::path path;~Cleanup(){fs::remove_all(path);}} cleanup{folder};
    const int run=999001;
    {
        TFile file((folder/"ntuple999001.root").c_str(),"RECREATE");
        TTree tree("ntuple","Synthetic pulses");
        int event,channel,hit,time,charge,pedestal=0;
        tree.Branch("Event",&event);tree.Branch("Channel",&channel);tree.Branch("Hit",&hit);tree.Branch("Time",&time);tree.Branch("Charge",&charge);
        tree.Branch("Pedestal",&pedestal);
        for(const auto& row:std::vector<Hit>{{0,152,1,0,-1},{1,150,2,1,0},{2,0,1,2,1},{3,575,3,3,2},{1,152,2,4,3},{1,15,1,5,4}}) {
            event=row.event;channel=row.channel;hit=row.hit;time=row.time;charge=row.charge;tree.Fill();
        }
        tree.Write();file.Close();
    }
    auto cuts=folder/"cuts.ini",veto=folder/"veto.csv",pairs=folder/"pairs.csv";
    {std::ofstream file(cuts);file<<"[global]\nreject = Channel == 152\n";}
    {std::ofstream file(veto);file<<"run,event\n";}
    {
        std::ofstream file(pairs);
        file<<"Name,X,Y,XBinWidth,XBins,XMin,XMax,YBinWidth,YBins,YMin,YMax,LogZ\n";
        for(const auto& row:std::vector<std::vector<std::string>>{
            {"Event_vs_Pad","Pad","Event","","","","","","2","0","2","false"},
            {"Time_vs_Charge","Charge","Time","","4","0","4","","4","0","4","true"}}) {
            for(size_t i=0;i<row.size();++i)file<<(i?",":"")<<row[i];file<<'\n';
        }
    }
    Config config(folder.c_str(),(folder/"results").c_str(),veto.c_str(),cuts.c_str(),(folder/"missing-study.ini").c_str());
    config.plots2DPath=pairs;
    Raw2D(config,run,false);Raw2D(config,run,true);
    for(bool cleaned:{false,true}) {
        auto out=config.output/"raw/run999001";if(cleaned)out/="cleaned";out/="2d";
        TFile file((out/"raw_2d_histograms.root").c_str(),"READ");
        auto pad=dynamic_cast<TH2D*>(file.Get("raw2d_Event_vs_Pad"));
        auto charge=dynamic_cast<TH2D*>(file.Get("raw2d_Time_vs_Charge"));
        Check(pad && charge,"Missing saved TH2s");
        Check(pad->GetEntries()==(cleaned?4:6),"Mapping or cuts lost/added hits");
        Check(pad->Integral()==(cleaned?2:4),"Wrong visible event/pad count");
        Check(pad->GetBinContent(0,3)==1,"Unmapped/event-overflow corner lost");
        Check(charge->GetBinContent(5,5)==1,"Both-axis overflow corner lost");
        Check(charge->GetBinContent(4,5)==(cleaned?0:1),"Exact upper Time edge was included");
        Check(charge->Integral()==3,"Wrong visible charge/time count");
        if(cleaned)Check(pad->Integral(1,1,0,3)==0,"Electronics exclusion did not reach mapped Pad0");
        std::string stem="TPC_raw2d_999001"+std::string(cleaned?"_clean":"");
        Check(file.Get(stem.c_str()) && file.Get((stem+"_Event_vs_Pad").c_str()) && file.Get((stem+"_Time_vs_Charge").c_str()),"Missing summary/individual canvases");
        file.Close();
    }
    // Exercise the public unified entry and both dimension-selection flags.
    for(const auto& mode:std::vector<std::string>{"both","1d","2d"}) {
        auto output=folder/("unified_"+mode);
        std::string command=mode=="both"?"raw":"raw --"+mode;
        TPC(command.c_str(),run,-1,folder.c_str(),output.c_str(),veto.c_str(),cuts.c_str(),(folder/"missing-study.ini").c_str());
        auto one=output/"raw/run999001/raw_histograms.root";
        auto two=output/"raw/run999001/2d/raw_2d_histograms.root";
        Check(fs::exists(one)==(mode!="2d"),"Wrong 1D group selection");
        Check(fs::exists(two)==(mode!="1d"),"Wrong 2D group selection");
        for(const auto& item:std::vector<std::pair<fs::path,int>>{{one,8},{two,5}})if(fs::exists(item.first)) {
            TFile file(item.first.c_str(),"READ");int count=0;
            for(auto object:*file.GetListOfKeys())if(std::string(static_cast<TKey*>(object)->GetClassName())=="TCanvas")++count;
            Check(count==item.second,"Wrong number of summary/individual canvases");file.Close();
        }
    }
    std::cout<<"PASS: native unified groups/flags, mapping/cuts, exact edges, flow corners and saved canvases."<<std::endl;
}
