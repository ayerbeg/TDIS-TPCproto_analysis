#include "../TPCDetector.h"
#include <chrono>
// Run: root -l -b -q 'analysis_root/tests/RunCutsCheck.C()'
// Temporary fixtures do not change the working manual-cut database.
void RunCutsCheck() {
    using namespace PersonalTPC;
    auto path=fs::temp_directory_path()/("tpc-native-cuts-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".ini");
    auto studyPath=path;studyPath+=".study.ini";
    struct Cleanup {fs::path p;~Cleanup(){fs::remove(p);}} cleanup{path};
    Cleanup cleanupStudy{studyPath};
    {std::ofstream f(studyPath);f<<"[global]\nreject = Channel >= 285\n";}
    {std::ofstream f(path);f<<"[global]\nreject = Time >= 228\n[run 01616]\nreject = Channel in 40..44 && Time >= 90 && Time < 95\nreject = Event in 1000..1002\nreject = Charge < -1.5\n";}
    Config config("","","",path.c_str(),studyPath.c_str());Geometry geometry(config);
    Check(!PassHit({1,40,1,90,100},config,false,1616),"Correlated lower endpoint");
    Check(!PassHit({1,44,1,94,100},config,false,1616),"Correlated inclusive channel endpoint");
    Check(PassHit({1,44,1,95,100},config,false,1616),"Time upper endpoint must pass");
    Check(PassHit({1,40,1,90,100},config,false,1615),"Rules leaked to other run");
    Check(!PassHit({1002,10,1,30,100},config,false,1616),"Event inclusive range endpoint");
    Check(PassHit({1,335,1,30,100},config,false,1615),"Event profile excludes large pads");
    Check(!PassHit({1,335,1,30,100},config,true,1615),"Voltage profile omitted");
    Check(!PassHit({1,10,1,30,-2},config,false,1616),"Negative threshold");
    {std::ofstream f(path);f<<"[global]\nreject = FADC == 153\nreject = Pad == 24\n";}
    Config aliases("","","",path.c_str());
    Check(!PassHit({1,152,1,30,100},aliases),"FADC conversion");
    Check(!PassHit({1,15,1,30,100},aliases),"Pad conversion");
    {std::ofstream f(path);f<<"[run 1616]\nreject = Channel == 42 && Time == 20\n";}
    Config pulses("","","",path.c_str());RunData data;data.run=1616;
    data.hits={{1,42,1,20,100},{1,42,2,40,200}};data.events[1]={0,1};
    auto hits=DisplayHits(data,1,pulses,geometry,{});
    Check(hits.size()==1 && hits[0].time==40,"First hit was reduced before DB cuts");
    {std::ofstream f(path);f<<"[global]\n";}
    Config empty("","","",path.c_str(),path.c_str());
    Check(PassHit({1,335,1,500,100},empty,true,1616),"Hidden duplicate global or study cuts");
    auto missingStudyPath=path;missingStudyPath+=".missing";
    Config independent("","","",path.c_str(),missingStudyPath.c_str());
    Check(PassHit({1,335,1,30,100},independent,false,1616),"Detector cleaning read an independent study DB");
    bool missingStudyFailed=false;
    try{PassHit({1,335,1,30,100},independent,true,1616);}catch(...){missingStudyFailed=true;}
    Check(missingStudyFailed,"Voltage study ignored a missing study DB");
    {std::ofstream f(path);f<<"[global]\nreject = Pedestal > 100\n";}
    Config missing("","","",path.c_str());bool failed=false;
    try{missing.Cuts().Validate(1616,{"Event","Channel","Hit","Time","Charge"});}catch(...){failed=true;}
    Check(failed,"Missing branch ignored");
    {std::ofstream f(path);f<<"[run1616]\nreject = Channnel == 1\n";}
    failed=false;try{CutDatabase bad(path);}catch(...){failed=true;}Check(failed,"Invalid syntax accepted");
    {std::ofstream f(path);f<<"[voltage]\nreject = Channel >= 285\n";}
    failed=false;try{CutDatabase bad(path);}catch(...){failed=true;}Check(failed,"Study block accepted in detector DB");
    std::cout<<"PASS: native global/run cuts, independent study DB, ranges, correlations, aliases, pulse order and validation."<<std::endl;
}
