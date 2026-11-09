#ifndef PERSONAL_TPC_DETECTOR_H
#define PERSONAL_TPC_DETECTOR_H
// C++17 + ROOT, with no Python dependency.
// This file handles raw input, selection, pad geometry, and event noise.
// TPC.C handles menus, plots, and exports. Manual rules are in the shared
// config/run_cuts.ini; settings.csv controls automatic-noise thresholds.
// Only the independent voltage command adds config/voltage_study_cuts.ini.
// All paths are relative to this folder unless explicitly supplied to TPC().
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <cctype>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PersonalTPC {
namespace fs = std::filesystem;
inline fs::path SourceDirectory() {
    // Normalize ../ before taking parent paths (also when included by tests).
    static const fs::path directory = fs::absolute(__FILE__).lexically_normal().parent_path();
    return directory;
}
inline std::string Trim(std::string text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    return first == std::string::npos ? "" : text.substr(first, text.find_last_not_of(" \t\r\n")-first+1);
}
// RFC-style CSV fields: quoted commas and doubled quotes are supported.
// Spaces/tabs around headers and values are ignored; commas remain required.
inline std::vector<std::string> CSV(const std::string& line) {
    std::vector<std::string> result;
    std::string field; bool quoted = false;
    for (size_t i=0;i<line.size();++i) {
        const char ch=line[i];
        if (ch=='"') {
            if (quoted && i+1<line.size() && line[i+1]=='"') {field+='"';++i;}
            else quoted=!quoted;
        } else if (ch==',' && !quoted) {result.push_back(Trim(field));field.clear();}
        else field+=ch;
    }
    result.push_back(Trim(field));return result;
}
inline std::string QuoteCSV(const std::string& text) {
    std::string out="\"";for(char c:text) {if(c=='"')out+='"';out+=c;}return out+'"';
}
inline std::set<int> IntegerList(const std::string& text) {
    std::set<int> result;std::stringstream ss(text);std::string item;
    while(std::getline(ss,item,';'))if(!Trim(item).empty())result.insert(std::stoi(item));
    return result;
}

// ---- Shared manual rejection DB (config/run_cuts.ini) ----
// Lines describe values to REMOVE. Lines are OR-ed, conditions joined by
// && on one line are AND-ed. No arbitrary C++/ROOT expression is executed.
inline std::string Lower(std::string text) {
    std::transform(text.begin(),text.end(),text.begin(),[](unsigned char c){return std::tolower(c);});return text;
}
struct CutCondition {
    std::string field,op;std::vector<double> values;bool range=false;
    bool Matches(double x)const {
        if(op=="in" || op=="not in") {
            bool found=range?(x>=values[0] && x<=values[1]):std::find(values.begin(),values.end(),x)!=values.end();
            return op=="not in"?!found:found;
        }
        double y=values[0];if(op=="==")return x==y;if(op=="!=")return x!=y;
        if(op=="<")return x<y;if(op=="<=")return x<=y;if(op==">")return x>y;return x>=y;
    }
};
struct CutRule {std::string section,expression,note;int line;std::vector<CutCondition> conditions;};
class CutDatabase {
public:
    fs::path path;std::string text;std::vector<CutRule> rules;
    explicit CutDatabase(const fs::path& filename):path(fs::absolute(filename)) {
        std::ifstream file(path);if(!file)throw std::runtime_error("Cannot read cut database: "+path.string());
        text=std::string(std::istreambuf_iterator<char>(file),{});
        std::stringstream lines(text);std::string original,section;std::set<std::string> sections;int line=0;
        while(std::getline(lines,original)) {
            ++line;auto hash=original.find('#');std::string code=Trim(original.substr(0,hash));
            std::string note=hash==std::string::npos?"":Trim(original.substr(hash+1));if(code.empty())continue;
            try {
                if(code[0]=='[') {
                    std::smatch m;
                    if(!std::regex_match(code,m,std::regex(R"(\[\s*(global|run\s+[0-9]+)\s*\])",std::regex::icase)))
                        throw std::runtime_error("Expected [global] or [run NNNN]; study cuts belong in a separate file");
                    section=Lower(Trim(m[1]));
                    if(section.substr(0,3)=="run")section="run "+std::to_string(std::stoi(section.substr(3)));
                    if(!sections.insert(section).second)throw std::runtime_error("Duplicate block ["+section+"]");
                    continue;
                }
                if(section.empty())throw std::runtime_error("Put rules inside a named block");
                auto equal=code.find('=');
                if(equal==std::string::npos || Lower(Trim(code.substr(0,equal)))!="reject")throw std::runtime_error("Expected reject = CONDITION");
                CutRule rule{section,Trim(code.substr(equal+1)),note,line,{}};
                size_t begin=0;
                while(true){auto end=rule.expression.find("&&",begin);rule.conditions.push_back(ParseCondition(Trim(rule.expression.substr(begin,end==std::string::npos?end:end-begin))));
                    if(end==std::string::npos)break;begin=end+2;}
                rules.push_back(std::move(rule));
            } catch(const std::exception& e) {throw std::runtime_error(path.string()+":"+std::to_string(line)+": "+e.what());}
        }
    }
    static double Number(const std::string& value) {
        auto text=Trim(value);size_t consumed=0;
        if(!std::regex_match(text,std::regex(R"([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)")))
            throw std::runtime_error("Expected finite numeric cut value: "+text);
        double result=std::stod(text,&consumed);
        if(consumed!=text.size() || !std::isfinite(result))throw std::runtime_error("Invalid numeric cut value");return result;
    }
    static CutCondition ParseCondition(const std::string& term) {
        std::smatch m;
        if(!std::regex_match(term,m,std::regex(R"(([A-Za-z_][A-Za-z_0-9]*)\s*(not\s+in|in|==|!=|<=|>=|<|>)\s*(.+))")))
            throw std::runtime_error("Invalid condition: "+term);
        const std::map<std::string,std::string> names={{"event","Event"},{"channel","Channel"},{"hit","Hit"},{"time","Time"},{"charge","Charge"},{"pedestal","Pedestal"},{"blockhead","Blockhead"},{"fadc","FADC"},{"pad","Pad"}};
        auto found=names.find(Lower(m[1]));if(found==names.end())throw std::runtime_error("Unknown field: "+m[1].str());
        std::string op=m[2],rhs=Trim(m[3]);if(op.substr(0,3)=="not")op="not in";
        CutCondition c{found->second,op,{},false};
        if(op=="in" || op=="not in") {
            auto dots=rhs.find("..");
            if(dots!=std::string::npos){c.values={Number(rhs.substr(0,dots)),Number(rhs.substr(dots+2))};c.range=true;if(c.values[0]>c.values[1])throw std::runtime_error("Range lower endpoint exceeds upper endpoint");}
            else {if(rhs.empty() || rhs.back()==',')throw std::runtime_error("Empty list value");std::stringstream ss(rhs);std::string value;while(std::getline(ss,value,','))c.values.push_back(Number(value));}
        } else c.values={Number(rhs)};
        return c;
    }
    bool Active(const CutRule& rule,int run)const {
        return rule.section=="global" || (run>=0 && rule.section=="run "+std::to_string(run));
    }
    template<class Getter> bool Accept(int run,Getter get)const {
        for(const auto& rule:rules)if(Active(rule,run)) {
            bool bad=true;for(const auto& c:rule.conditions)if(!c.Matches(get(c.field))){bad=false;break;}
            if(bad)return false;
        }
        return true;
    }
    void Validate(int run,const std::set<std::string>& branches)const {
        for(const auto& rule:rules)if(Active(rule,run))for(const auto& c:rule.conditions) {
            std::string field=(c.field=="FADC" || c.field=="Pad")?"Channel":c.field;
            if(!branches.count(field))throw std::runtime_error(path.string()+":"+std::to_string(rule.line)+": missing input branch "+field);
        }
    }
    std::string Describe(int run=-1)const {
        std::string result;for(const auto& rule:rules)if(Active(rule,run)){if(!result.empty())result+=" && ";result+="NOT("+rule.expression+")";}
        return result.empty()?"all rows":result;
    }
    void Snapshot(const fs::path& destination)const {std::ofstream out(destination);out<<text;}
};

struct Config {
    fs::path source=SourceDirectory(), input=source.parent_path()/"data"/"root", output;
    fs::path vetoOverride,cutDBPath,binningPath,studyCutDBPath,plots2DPath;
    mutable std::shared_ptr<CutDatabase> database,studyDatabase;
    std::map<int,int> padNumbers;
    std::map<std::string,std::string> values;
    int timeMax, noiseChannels, moduleChannels;
    int moduleSynchronous, moduleHitsPerChannel, channelsPerModule, readoutModules;
    std::set<int> excludedRuns;
    Config(const char* data="",const char* out="",const char* veto="",const char* cuts="",const char* studyCuts="") {
        // Organized workspace: data/root/ for input, results/root/ for
        // output. Retain support for sharing beside a flat ROOT dataset.
        if(!fs::is_directory(input))input=source.parent_path();
        if(data && *data)input=fs::absolute(data);
        auto project=(input.filename()=="root" && input.parent_path().filename()=="data")
                     ?input.parent_path().parent_path():input;
        output=(out && *out)?fs::absolute(out):project/"results"/"root";
        if(veto && *veto)vetoOverride=fs::absolute(veto);
        cutDBPath=(cuts && *cuts)?fs::absolute(cuts):source.parent_path()/"config"/"run_cuts.ini";
        studyCutDBPath=(studyCuts && *studyCuts)?fs::absolute(studyCuts):source.parent_path()/"config"/"voltage_study_cuts.ini";
        binningPath=source.parent_path()/"config"/"raw_binning.csv";
        plots2DPath=source.parent_path()/"config"/"raw_2d_plots.csv";
        std::ifstream mapping(source/"auxiliary"/"pad_map.csv");std::string mappingLine;std::getline(mapping,mappingLine);
        while(std::getline(mapping,mappingLine)){auto row=CSV(mappingLine);if(row.size()>=3)padNumbers[std::stoi(row[2])]=std::stoi(row[0]);}
        std::ifstream file(source/"settings.csv");
        if(!file)throw std::runtime_error("Cannot read analysis_root/settings.csv");
        std::string line;std::getline(file,line);
        while(std::getline(file,line)){auto row=CSV(line);if(row.size()>=2)values[row[0]]=row[1];}
        timeMax=std::stoi(values.at("time_max"));
        noiseChannels=std::stoi(values.at("noise_channels"));
        moduleChannels=std::stoi(values.at("module_channels"));
        moduleSynchronous=std::stoi(values.at("module_synchronous"));
        moduleHitsPerChannel=std::stoi(values.at("module_hits_per_channel"));
        channelsPerModule=std::stoi(values.at("channels_per_module"));
        readoutModules=std::stoi(values.at("readout_modules"));
        excludedRuns=IntegerList(values.at("excluded_runs"));
        if(timeMax<=0 || channelsPerModule<=0 || readoutModules<=0)throw std::runtime_error("Invalid settings.csv values");
    }
    fs::path Auxiliary(const char* name)const{return source/"auxiliary"/name;}
    fs::path VetoFile()const {
        if(!vetoOverride.empty())return vetoOverride;
        auto generated=output/"voltage"/"rejected_events.csv";
        return fs::exists(generated)?generated:Auxiliary("rejected_events.csv");
    }
    const CutDatabase& Cuts()const {
        if(!database)database=std::make_shared<CutDatabase>(cutDBPath);return *database;
    }
    // Only the independent voltage workflow requests this extra database.
    // Event displays/raw cleaning can run even if this study file is absent.
    const CutDatabase& StudyCuts()const {
        if(!studyDatabase)studyDatabase=std::make_shared<CutDatabase>(studyCutDBPath);
        return *studyDatabase;
    }
    std::string CutText(bool voltage=false,int run=-1)const {
        std::string text=Cuts().Describe(run);
        if(voltage)text+=" && "+StudyCuts().Describe(run);
        return text;
    }

};

struct Pad {
    int number,fadc,channel,row,col;std::string section;
    double x,y,width,height;
    double CenterX()const{return x+width/2;}double CenterY()const{return y+height/2;}
};
class Geometry {
public:
    std::vector<Pad> pads;
    std::map<int,int> byChannel;
    double width=0,height=0;
    explicit Geometry(const Config& config) {
        std::ifstream file(config.Auxiliary("pad_map.csv"));
        if(!file)throw std::runtime_error("Cannot read auxiliary/pad_map.csv");
        std::string line;std::getline(file,line);std::set<int> ids;
        while(std::getline(file,line)) {
            auto r=CSV(line);if(r.size()<10)continue;
            Pad p{std::stoi(r[0]),std::stoi(r[1]),std::stoi(r[2]),std::stoi(r[4]),std::stoi(r[5]),r[3],
                  std::stod(r[6]),std::stod(r[7]),std::stod(r[8]),std::stod(r[9])};
            if(p.number<0 || p.number>=410 || p.fadc!=p.channel+1 || p.width<=0 || p.height<=0 ||
               !ids.insert(p.number).second || byChannel.count(p.channel))throw std::runtime_error("Invalid/duplicate pad mapping");
            byChannel[p.channel]=p.number;pads.push_back(p);
            width=std::max(width,p.x+p.width);height=std::max(height,p.y+p.height);
        }
        std::sort(pads.begin(),pads.end(),[](const Pad&a,const Pad&b){return a.number<b.number;});
        if(pads.size()!=410)throw std::runtime_error("Expected 410 pads, numbered 0-409");
        for(size_t i=0;i<pads.size();++i)if(pads[i].number!=int(i))throw std::runtime_error("Pad numbering has holes");
    }
    int Number(int channel)const {auto it=byChannel.find(channel);return it==byChannel.end()?-1:it->second;}
};

// A Hit is one RAW pulse row. Keep these rows, including Hit index, for any
// later track analysis. No aggregation, time conversion, or z inference here.
struct Hit {int event,channel,hit,time,charge;int pedestal=0,blockhead=0;};
struct RunData {
    int run;std::vector<Hit> hits;std::map<int,std::vector<size_t>> events;std::set<std::string> branches;
};
inline RunData ReadRun(int run,const Config& config) {
    auto path=config.input/("ntuple"+std::to_string(run)+".root");
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
    if(!file || file->IsZombie())throw std::runtime_error("Cannot open "+path.string());
    auto tree=dynamic_cast<TTree*>(file->Get("ntuple"));
    if(!tree)throw std::runtime_error("No ntuple tree in "+path.string());
    for(const char* branch:{"Event","Channel","Hit","Time","Charge"})
        if(!tree->GetBranch(branch))throw std::runtime_error(std::string("Missing branch ")+branch);
    RunData data;data.run=run;data.hits.reserve(tree->GetEntries());
    for(const char* name:{"Event","Channel","Hit","Time","Charge","Pedestal","Blockhead"})
        if(tree->GetBranch(name))data.branches.insert(name);
    TTreeReader reader(tree);
    TTreeReaderValue<Int_t> event(reader,"Event"),channel(reader,"Channel"),hit(reader,"Hit"),time(reader,"Time"),charge(reader,"Charge");
    std::unique_ptr<TTreeReaderValue<Int_t>> pedestal,blockhead;
    if(data.branches.count("Pedestal"))pedestal=std::make_unique<TTreeReaderValue<Int_t>>(reader,"Pedestal");
    if(data.branches.count("Blockhead"))blockhead=std::make_unique<TTreeReaderValue<Int_t>>(reader,"Blockhead");
    while(reader.Next()) {
        size_t index=data.hits.size();data.hits.push_back({*event,*channel,*hit,*time,*charge,pedestal?**pedestal:0,blockhead?**blockhead:0});
        data.events[*event].push_back(index);
    }
    return data;
}
inline double HitValue(const Hit& h,const std::string& field,const Config& config) {
    if(field=="Event")return h.event;if(field=="Channel")return h.channel;if(field=="Hit")return h.hit;
    if(field=="Time")return h.time;if(field=="Charge")return h.charge;
    if(field=="Pedestal")return h.pedestal;if(field=="Blockhead")return h.blockhead;
    if(field=="FADC")return h.channel+1;
    if(field=="Pad"){auto it=config.padNumbers.find(h.channel);return it==config.padNumbers.end()?-1:it->second;}
    throw std::runtime_error("Unknown branch: "+field);
}
inline bool PassHit(const Hit& hit,const Config& config,bool voltage=false,int run=-1) {
    auto get=[&](const std::string& field){return HitValue(hit,field,config);};
    if(!config.Cuts().Accept(run,get))return false;
    return !voltage || config.StudyCuts().Accept(run,get);
}
inline void ValidateCuts(const RunData& data,const Config& config,bool voltage=false) {
    config.Cuts().Validate(data.run,data.branches);
    if(voltage)config.StudyCuts().Validate(data.run,data.branches);
}
inline std::set<int> ReadVetoes(int run,const Config& config) {
    auto path=config.VetoFile();std::ifstream file(path);
    if(!file)throw std::runtime_error("Missing veto ledger: "+path.string());
    std::set<int> result;std::string line;std::getline(file,line);
    while(std::getline(file,line)){auto r=CSV(line);if(r.size()>=2 && std::stoi(r[0])==run)result.insert(std::stoi(r[1]));}
    return result;
}
inline std::vector<Hit> DisplayHits(const RunData& data,int event,const Config& config,
                                    const Geometry& geometry,const std::set<int>& vetoes) {
    // FIRST passing pulse/channel is a display convention only. Voltage
    // spectra count all passing rows; future fitting should use raw hits.
    std::vector<Hit> result;std::set<int> seen;
    auto it=data.events.find(event);if(it==data.events.end() || vetoes.count(event))return result;
    for(size_t index:it->second) {
        const auto& h=data.hits[index];
        if(PassHit(h,config,false,data.run) && geometry.Number(h.channel)>=0 && seen.insert(h.channel).second)result.push_back(h);
    }
    return result;
}
struct ModuleEvidence {int module,channels,hits,synchronous;};
struct NoiseInfo {
    int event=0,distinct=0,totalHits=0,baseHits=0;
    bool detectorWide=false;std::vector<ModuleEvidence> modules;
    bool Vetoed()const{return detectorWide || !modules.empty();}
};
inline NoiseInfo ClassifyEvent(const RunData& data,int event,const std::vector<size_t>& indices,const Config& config) {
    // BEFORE all hit cuts: repeated hits on one channel cannot inflate
    // detector occupancy. Module groups are electronics blocks, not rows.
    struct Module {std::set<int> channels;int hits=0;};
    std::vector<Module> modules(config.readoutModules);std::set<int> all;
    NoiseInfo result;result.event=event;result.totalHits=indices.size();
    for(size_t index:indices) {
        const auto& h=data.hits[index];all.insert(h.channel);
        if(PassHit(h,config,true,data.run))++result.baseHits;
        if(h.channel<0 || h.channel>=config.channelsPerModule*config.readoutModules)
            throw std::runtime_error("Channel outside configured readout range; edit settings.csv");
        auto& module=modules[h.channel/config.channelsPerModule];module.channels.insert(h.channel);++module.hits;
    }
    result.distinct=all.size();result.detectorWide=result.distinct>=config.noiseChannels;
    for(int m=0;m<config.readoutModules;++m) {
        const auto& module=modules[m];if(int(module.channels.size())<config.moduleChannels)continue;
        std::map<int,std::set<int>> channelsAtTime;
        for(size_t index:indices){const auto& h=data.hits[index];if(h.channel/config.channelsPerModule==m)channelsAtTime[h.time].insert(h.channel);}
        int simultaneous=0;for(const auto& item:channelsAtTime)simultaneous=std::max(simultaneous,int(item.second.size()));
        if(simultaneous>=config.moduleSynchronous || module.hits>=config.moduleHitsPerChannel*int(module.channels.size()))
            result.modules.push_back({m,int(module.channels.size()),module.hits,simultaneous});
    }
    return result;
}
inline std::string ModuleJSON(const NoiseInfo& info) {
    std::ostringstream out;out<<"[";
    for(size_t i=0;i<info.modules.size();++i) {
        const auto& m=info.modules[i];if(i)out<<", ";
        out<<"{\"module\": "<<m.module<<", \"channels\": "<<m.channels<<", \"hits\": "<<m.hits<<", \"synchronous_channels\": "<<m.synchronous<<"}";
    }
    return out.str()+"]";
}
inline void Check(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
}
#endif
