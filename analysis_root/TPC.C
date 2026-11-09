// Native C++/ROOT interface. This does not call or modify the Python tools.
//
// From the workspace root (input files are under data/root/):
//   root -l 'analysis_root/TPC.C()'                       // menu
//   root -l 'analysis_root/TPC.C("event",1616,5584)'      // browse
//   root -l -b -q 'analysis_root/TPC.C("save",1616,5584)' // save 2D/3D
//   root -l -b -q 'analysis_root/TPC.C("examples")'
//   root -l -b -q 'analysis_root/TPC.C("voltage")'
//   root -l -b -q 'analysis_root/TPC.C("check")'
// For compilation, use TPC.C+(...) instead of TPC.C(...).
// Optional arguments after event: inputDirectory, outputDirectory, vetoCSV,
// cutsDB, studyCutsDB. General cleaning lives in config/run_cuts.ini;
// the independent voltage study adds config/voltage_study_cuts.ini.
// Edit drawing functions here for plot styles;
// TPCDetector.h isolates raw hits, mapping, selection and noise algorithms.
//
// READING GUIDE / MAIN EDIT POINTS
//   * ../config/run_cuts.ini: shared global and per-run manual rejection rules.
//   * ../config/voltage_study_cuts.ini: independent voltage-study restrictions.
//   * ../config/raw_binning.csv: raw histogram axes and log scales.
//   * ../config/raw_2d_plots.csv: 2D pairs, axis overrides and colour scales.
//   * settings.csv: automatic-noise thresholds and voltage histogram limits.
//   * TPCDetector.h: raw input, mapping, PassHit(), DisplayHits(), ClassifyEvent().
//   * Section 2 below: pad colors, labels, margins and 2D/3D drawing.
//   * Sections 6-8: voltage histograms, normalization, run choice and exports.
//   * Section 10: the public TPC(...) command dispatcher (start here to see
//     which function each menu/command invokes).
//
// DATA CONVENTIONS TO KEEP IN MIND
//   ROOT Channel starts at 0; printed fADC = Channel + 1; Pad runs 0..409.
//   Medium pads 0..249 are below large pads 250..409.
//   Event displays keep the first passing hit/channel and include both sizes.
//   Voltage histograms count ALL accepted pulses and add the voltage-only
//   Channel limit. Their normalization is by accepted HITS, not events.
//   Time and Charge stay in native branch units. No track fit or drift-z
//   conversion is performed by this macro.

#include "TPCDetector.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2Poly.h>
#include <TBox.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TColor.h>
#include <TGraph.h>
#include <TLine.h>
#include <TArrow.h>
#include <TNamed.h>
#include <TTree.h>
#include <sys/select.h>
#include <limits>
#include <unistd.h>


void set_color_env(){   
  //******************************************************************
  //code to improve the color palette
  //from the root web page and color codes from:
  //http://ultrahigh.org/2007/08/20/making-pretty-root-color-palettes/
  //******************************************************************
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
  //******************************************************************
  

}





namespace PersonalTPC {


    // =========================================================================
    // SECTION 1: CANVAS AND OUTPUT UTILITIES
    // =========================================================================

    // ROOT uses object names as identifiers. A fresh suffix prevents a new
    // histogram from silently replacing an earlier one in the same session.
    inline std::string Unique(const std::string& stem) {
        static int n = 0;
        return stem + "_" + std::to_string(++n);
    }

    // Reuse the named window instead of opening a new one for every event.
    // Clear() removes its previous drawing before the next event is shown.
    inline TCanvas* Canvas(const std::string& name, const std::string& title, int w = 1100, int h = 900) {
        auto c = dynamic_cast<TCanvas*>(gROOT->GetListOfCanvases()->FindObject(name.c_str()));
        if (c) c->Clear();
        else c = new TCanvas(name.c_str(), title.c_str(), w, h);
        c->SetTitle(title.c_str());
        return c;
    }

    // stem is a full path WITHOUT an extension. Update() forces ROOT to
    // finish painting axes/palettes before saving the same view as PNG/PDF.
    inline void SaveCanvas(TCanvas* canvas, const fs::path& stem) {
        fs::create_directories(stem.parent_path());
        canvas->Update();
        canvas->SaveAs((stem.string() + ".png").c_str());
        canvas->SaveAs((stem.string() + ".pdf").c_str());
    }

    // Convert selected hits into one value per PAD for drawing:
    // charge=true uses Charge; false uses Time. The selection/first-hit
    // decision has already been made by DisplayHits() in TPCDetector.h.
    inline std::map<int, double> Values(const std::vector<Hit>& hits, const Geometry& geometry, bool charge) {
        std::map<int, double> values;
        for (const auto& h : hits) values[geometry.Number(h.channel)] = charge ? h.charge : h.time;
        return values;
    }


    // =========================================================================
    // SECTION 2: DRAW THE PHYSICAL PAD PLANE
    // =========================================================================

    // Draw on the CURRENT pad (gPad), so this works in a standalone canvas
    // or in one panel of a comparison. values is keyed by pad number.
    // threeD=false: true-size rectangles and labels; true: ROOT LEGO bars.
    inline void DrawPhysical(const Geometry& geometry, const std::map<int, double>& values, const std::string& title,
        bool threeD = false) {
        double maximum = 1;
        for (const auto& item : values) maximum = std::max(maximum, item.second);

        // Margins are fractions of the drawing area. The right margin makes
        // room for the color scale; the bottom margin makes room for x labels.
        gPad->SetRightMargin(.15);
        gPad->SetBottomMargin(.12);
        if (!threeD) {
            // TH2Poly preserves each pad's actual 4/5-mm footprint. Missing
            // values draw as 0; that is not a declaration that the pad is dead.
            auto h = new TH2Poly(Unique("pad_plane").c_str(), (title + ";x (mm);y (mm)").c_str(), 0., geometry.width,
                0., geometry.height);

            // Keep this display histogram independent of an open ROOT file.
            // kCanDelete lets the canvas clean it up when its drawing is cleared.
            h->SetDirectory(nullptr);
            h->SetBit(kCanDelete);
            for (const auto& pad : geometry.pads) {
                int bin = h->AddBin(pad.x, pad.y, pad.x + pad.width, pad.y + pad.height);
                auto it = values.find(pad.number);
                h->SetBinContent(bin, it == values.end() ? 0 : it->second);
            }
            h->SetMinimum(0);
            h->SetMaximum(maximum);
            h->Draw("COLZ");
            for (const auto& pad : geometry.pads) {
                TBox border(pad.x, pad.y, pad.x + pad.width, pad.y + pad.height);
                border.SetFillStyle(0);
                border.SetLineColor(kGray + 1);
                border.SetLineWidth(1);
                border.DrawClone();
                auto it = values.find(pad.number);
                double val = it == values.end() ? 0 : it->second;
                TLatex label;
                label.SetTextAlign(22);

                // Change this normalized text size to make pad numbers larger/smaller.
                label.SetTextSize(.011);
                label.SetTextColor(val > maximum / 2 ? kBlack : kWhite);
                label.DrawLatex(pad.CenterX(), pad.CenterY(), std::to_string(pad.number).c_str());
            }
        } else {
            // ROOT LEGO requires uniform bins. Rasterize each pad onto a 1-mm
            // display grid; this repeats HEIGHT, not hits. Never integrate this
            // rendering histogram for physics or treat Time bars as drift z.
            int nx = std::lround(geometry.width), ny = std::lround(geometry.height);
            auto h = new TH2D(Unique("pad_lego").c_str(), (title + ";x (mm);y (mm);Branch value").c_str(), nx, 0,
                geometry.width, ny, 0, geometry.height);
            h->SetDirectory(nullptr);
            h->SetBit(kCanDelete);
            for (const auto& pad : geometry.pads) {
                auto it = values.find(pad.number);
                double val = it == values.end() ? 0 : it->second;

                // The tiny inward offsets avoid assigning a shared edge to the
                // neighboring pad when FindBin() converts mm coordinates into bins.
                int x1 = h->GetXaxis()->FindBin(pad.x + .001), x2 = h->GetXaxis()->FindBin(pad.x + pad.width - .001);
                int y1 = h->GetYaxis()->FindBin(pad.y + .001), y2 = h->GetYaxis()->FindBin(pad.y + pad.height - .001);
                for (int ix = x1; ix <= x2; ++ix)
                    for (int iy = y1; iy <= y2; ++iy)
                        h->SetBinContent(ix, iy, val);
            }
            h->SetMinimum(0);
            h->SetMaximum(maximum);
            h->Draw("LEGO2Z");
        }
        gPad->Modified();
        gPad->Update();
    }

    // Open/update a Charge canvas and a Time canvas for the SAME selected
    // event. Explicit mode names keep 2D and 3D windows/exports separate.
    inline std::vector<TCanvas*> EventCanvases(int run, int event, const std::vector<Hit>& hits, const Geometry& geometry,
        bool threeD) {
        std::vector<TCanvas*> result;
        for (bool charge : { true, false }) {
            std::string kind = charge ? "charge" : "time", mode = threeD ? "3d" : "2d";
            auto c = Canvas("TPC_" + kind + "_run" + std::to_string(run) + "_" + mode, "TPC event");
            c->cd();
            DrawPhysical(geometry, Values(hits, geometry, charge), kind + " - Run " + std::to_string(run) + " - Event "
                + std::to_string(event), threeD);
            result.push_back(c);
        }
        return result;
    }


    // =========================================================================
    // SECTION 3: SAVE ONE EVENT: FIGURES, DISPLAYED-HIT CSV AND ROOT TREE
    // =========================================================================

    // Each CSV row corresponds to one displayed hit. x/y are PAD CENTERS;
    // width/height describe its footprint, not a fitted track uncertainty.
    inline void WriteHitCSV(std::ostream& out, int run, const std::vector<Hit>& hits, const Geometry& geometry) {
        for (const auto& h : hits) {
            const auto& p = geometry.pads.at(geometry.Number(h.channel));
            out << run << ',' << h.event << ',' << p.number << ',' << p.fadc << ',' << h.channel << ',' << h.charge
                << ',' << h.time << ',' << h.hit << ',' << p.CenterX() << ',' << p.CenterY() << ',' << p.width
                << ',' << p.height << '\n';
        }
    }

    // Workflow: read a run -> load its vetoes -> select one event -> export.
    // event<0 means choose the first event with accepted mapped hits.
    // This function writes only under the configured OUTPUT directory.
    inline std::vector<Hit> SaveEvent(const Config& config, const Geometry& geometry, int run, int event) {

        // Raw data is read-only. Cuts/vetoes and first-hit-per-channel reduction
        // are applied by DisplayHits(), not by changing the input tree.
        auto data = ReadRun(run, config);
        ValidateCuts(data, config);
        auto vetoes = ReadVetoes(run, config);
        if (event < 0)
            for (const auto& e : data.events)
                if (!DisplayHits(data, e.first, config, geometry, vetoes).empty()) {
                    event = e.first;
                    break;
                }
        auto hits = DisplayHits(data, event, config, geometry, vetoes);
        if (hits.empty()) throw std::runtime_error("Requested event is vetoed, absent, or has no selected mapped hits");
        auto out = config.output / "events" / ("run" + std::to_string(run));
        fs::create_directories(out);
        std::string tag = "event_" + std::to_string(run) + "_" + std::to_string(event);

        // RECREATE replaces an earlier export of THIS event in the output
        // folder. The original ntupleNNNN.root is never opened for writing.
        auto rootPath = out / (tag + ".root");
        TFile file(rootPath.c_str(), "RECREATE");
        TNamed("selection", config.CutText(false, run).c_str()).Write();
        TNamed("cuts_db", config.cutDBPath.c_str()).Write();
        config.Cuts().Snapshot(out / (tag + "_cuts.ini"));
        TNamed("hit_policy", "First passing hit per Channel; display only. Time has native units, not calibrated z.").Write();
        // ROOT TTree and CSV export exactly the displayed hits. Their transverse
        // pad centers/widths are useful debugging data, not reconstructed tracks.
        TTree tree("display_hits", "First passing mapped hit per channel");

        // The local tree owns its memory separately from the output TFile.
        // Write() below persists it without relying on the file to own the tree.
        tree.SetDirectory(nullptr);
        int runValue = run, eventValue = event, pad, fadc, channel, hit, time, charge;
        double x, y, width, height;
        tree.Branch("run", &runValue);
        tree.Branch("event", &eventValue);
        tree.Branch("pad", &pad);
        tree.Branch("fadc", &fadc);
        tree.Branch("channel", &channel);
        tree.Branch("hit", &hit);
        tree.Branch("time", &time);
        tree.Branch("charge", &charge);
        tree.Branch("x_mm", &x);
        tree.Branch("y_mm", &y);
        tree.Branch("width_mm", &width);
        tree.Branch("height_mm", &height);
        for (const auto& h : hits) {
            const auto& p = geometry.pads.at(geometry.Number(h.channel));
            pad = p.number;
            fadc = p.fadc;
            channel = h.channel;
            hit = h.hit;
            time = h.time;
            charge = h.charge;
            x = p.CenterX();
            y = p.CenterY();
            width = p.width;
            height = p.height;
            tree.Fill();
        }
        tree.Write();

        // Always save both 2D and 3D views, each for Charge and Time. The same
        // selected hits feed all four views; drawing mode does not change cuts.
        for (bool threeD : { false, true }) {
            auto canvases = EventCanvases(run, event, hits, geometry, threeD);
            for (size_t i = 0; i < canvases.size(); ++i) {
                std::string kind = i == 0 ? "charge" : "time";
                auto directory = threeD ? out / "3d" : out;
                SaveCanvas(canvases[i], directory / (tag + "_" + kind));
                file.cd();
                canvases[i]->Write((kind + (threeD ? "_3d" : "_2d")).c_str());
            }
        }
        std::ofstream csv(out / (tag + "_hits.csv"));
        csv << "run,event,pad,fadc,root_channel,charge,time,hit,x_mm,y_mm,width_mm,height_mm\n";
        WriteHitCSV(csv, run, hits, geometry);
        file.Close();
        int large = 0;
        for (const auto& h : hits) large += geometry.Number(h.channel) >= 250;
        std::cout << "Run " << run << ", event " << event << ": " << hits.size() - large << " medium and "
            << large << " large pads. Saved in " << out << std::endl;
        return hits;
    }


    // =========================================================================
    // SECTION 4: INTERACTIVE EVENT NAVIGATION
    // =========================================================================

    // A blocking getline alone would freeze ROOT window updates. Poll stdin
    // while processing GUI events so windows remain usable (including 3D
    // rotation). End-of-input is treated as q, just like an explicit quit.
    inline std::string ReadLineWithGUI() {
        while (true) {
            fd_set fds;
            FD_ZERO(&fds);
            FD_SET(STDIN_FILENO, &fds);
            timeval timeout { 0, 50000 };
            int ready = select(STDIN_FILENO + 1, &fds, nullptr, nullptr, &timeout);
            if (ready < 0) throw std::runtime_error("Cannot read terminal input");
            if (ready > 0) break;
            gSystem->ProcessEvents();
        }
        std::string line;
        if (!std::getline(std::cin, line)) return "q";
        return Trim(line);
    }

    // Build a list of eligible event IDs once, then move within that list.
    // Enter/p/jump cannot land on a noise-vetoed or empty event. Pressing s
    // invokes SaveEvent(), so interactive and batch exports use one policy.
    inline void Browse(const Config& config, const Geometry& geometry, int run, int event, bool threeD = false) {
        if (gROOT->IsBatch()) throw std::runtime_error("Use save in batch mode; event browsing needs interactive ROOT");
        auto data = ReadRun(run, config);
        ValidateCuts(data, config);
        auto vetoes = ReadVetoes(run, config);
        std::vector<int> accepted;
        for (const auto& item : data.events)
            if (!DisplayHits(data, item.first, config, geometry, vetoes).empty())
                accepted.push_back(item.first);
        if (accepted.empty()) throw std::runtime_error("No accepted mapped events");
        size_t position = 0;
        if (event >= 0) {
            auto found = std::find(accepted.begin(), accepted.end(), event);
            if (found == accepted.end()) throw std::runtime_error("Requested event is absent or vetoed");
            position = found - accepted.begin();
        }
        while (true) {
            auto hits = DisplayHits(data, accepted[position], config, geometry, vetoes);
            EventCanvases(run, accepted[position], hits, geometry, threeD);
            std::cout << "Event " << accepted[position] << " (" << position + 1 << '/' << accepted.size() << "): "
                << hits.size() << " pads. Enter=next, p=previous, event number=jump, s=save, q=quit: " << std::flush;
            auto line = ReadLineWithGUI();
            if (line == "q") break;
            if (line.empty()) {
                if (position + 1 < accepted.size()) ++position;
            } else if (line == "p") {
                if (position) --position;
            } else if (line == "s") SaveEvent(config, geometry, run, accepted[position]);
            else {
                try {
                    int target = std::stoi(line);
                    auto it = std::find(accepted.begin(), accepted.end(), target);
                    if (it != accepted.end()) position = it - accepted.begin();
                    else std::cout << "Event not accepted\n";
                } catch (...) {
                    std::cout << "Use Enter, p, s, q, or an event number\n";
                }
            }
        }
    }


    // =========================================================================
    // SECTION 5: REPRODUCE THE TWO REFERENCE EXAMPLES
    // =========================================================================

    // These are deliberately chosen illustrations, not an automatic track
    // selection: run1615/event42535 and run1616/event5584. Change the two
    // SaveEvent calls and the corresponding titles/filenames together if
    // you want a different reference pair.
    inline void Examples(const Config& config, const Geometry& geometry) {
        auto first = SaveEvent(config, geometry, 1615, 42535), second = SaveEvent(config, geometry, 1616, 5584);
        auto c = Canvas("TPC_examples", "Examples: 1615 and 1616", 1700, 1400);
        c->Divide(2, 2);
        c->cd(1);
        DrawPhysical(geometry, Values(first, geometry, true), "Run 1615, event 42535: Charge");
        c->cd(2);
        DrawPhysical(geometry, Values(first, geometry, false), "Run 1615, event 42535: Time");
        c->cd(3);
        DrawPhysical(geometry, Values(second, geometry, true), "Run 1616, event 5584: Charge");
        c->cd(4);
        DrawPhysical(geometry, Values(second, geometry, false), "Run 1616, event 5584: Time");
        auto directory = config.output / "comparisons";
        SaveCanvas(c, directory / "examples_1615_1616");
        std::ofstream csv(directory / "example_hits_1615_1616.csv");
        csv << "run,event,pad,fadc,root_channel,charge,time,hit,x_mm,y_mm,width_mm,height_mm\n";
        WriteHitCSV(csv, 1615, first, geometry);
        WriteHitCSV(csv, 1616, second, geometry);
        config.Cuts().Snapshot(directory / "examples_1615_1616_cuts.ini");
    }

    // =========================================================================
    // SECTION 6: PER-RUN NUMERICAL ANALYSIS FOR THE VOLTAGE SCAN
    // =========================================================================

    // One point/event for the run1610 diagnostic. "hits" means raw pulse
    // rows that pass voltage hit cuts; "channels" is BEFORE those cuts.
    struct AuditPoint {
        int event, channels, hits;
        bool detectorWide, vetoed;
    };
    // Collected numerical results for ONE run. Keep this separate from
    // drawing: plots may use copies/clones, while these histograms remain
    // the authoritative accepted-hit counts written to the ROOT output.
    // time: all accepted hits (including unmapped electronics).
    // pad/timePad: only mapped hits. occupancy: all events before hit cuts.
    // previous: detector-wide veto only, retained for the historical audit.
    struct Summary {
        int run = 0, voltage = -1;
        long long totalEvents = 0, acceptedEvents = 0, totalHits = 0, baseHits = 0, acceptedHits = 0, mappedHits = 0;
        bool selected = false;
        std::vector<NoiseInfo> vetoes;
        std::vector<AuditPoint> audit;
        std::shared_ptr<TH1D> time, pad, occupancy, previous;
        std::shared_ptr<TH2D> timePad;
    };

    // Analyze each event in this order:
    //   1) ClassifyEvent() on ALL raw hits (no Time/channel cuts yet).
    //   2) PassHit(..., true) applies the voltage-specific hit cuts.
    //   3) Reject the whole event if either noise criterion flagged it.
    //   4) Fill spectra once per accepted RAW pulse, without deduplication.
    inline Summary Analyze(int run, int voltage, const Config& config, const Geometry& geometry) {
        auto data = ReadRun(run, config);
        ValidateCuts(data, config, true);
        Summary s;
        s.run = run;
        s.voltage = voltage;
        s.totalEvents = data.events.size();
        s.totalHits = data.hits.size();
        auto tag = "run" + std::to_string(run);

        // Integer-centered Time bins have edges -0.5, 0.5, ..., timeMax-0.5.
        // The strict Time < timeMax cut excludes the upper boundary itself.
        s.time = std::make_shared<TH1D>(("time_" + tag).c_str(), ";Time (native branch units);Accepted hits", config.timeMax,
            -.5, config.timeMax - .5);
        s.previous = std::make_shared<TH1D>(("previous_" + tag).c_str(), "", config.timeMax, -.5, config.timeMax - .5);
        s.pad = std::make_shared<TH1D>(("pad_" + tag).c_str(), ";Pad number (0-based);Accepted hits", 410, -.5,
            409.5);
        int channels = config.readoutModules * config.channelsPerModule;
        s.occupancy = std::make_shared<TH1D>(("occupancy_" + tag).c_str(), ";Distinct channels/event;Events", channels
            + 1, -.5, channels + .5);
        s.timePad = std::make_shared<TH2D>(("time_pad_" + tag).c_str(), ";Pad number (0-based);Time (native units)",
            410, -.5, 409.5, config.timeMax, -.5, config.timeMax - .5);

        // Detach all numerical histograms from the current file/directory.
        // Their shared_ptr owners keep them alive until the scan is exported.
        s.time->SetDirectory(nullptr);
        s.previous->SetDirectory(nullptr);
        s.pad->SetDirectory(nullptr);
        s.occupancy->SetDirectory(nullptr);
        s.timePad->SetDirectory(nullptr);
        for (const auto& item : data.events) {
            auto noise = ClassifyEvent(data, item.first, item.second, config);
            s.occupancy->Fill(noise.distinct);
            s.baseHits += noise.baseHits;
            if (noise.Vetoed()) s.vetoes.push_back(noise);
            if (run == 1610) s.audit.push_back({ item.first, noise.distinct, noise.baseHits, noise.detectorWide,
                noise.Vetoed() });

            // Track whether this event contributes at least one hit. This counts
            // accepted EVENTS separately from the accepted HIT total.
            bool hasAccepted = false;
            for (size_t index : item.second) {
                const auto& h = data.hits[index];

                // true selects the VOLTAGE profile, including its extra Channel limit.
                // Event displays use the other profile through DisplayHits().
                if (!PassHit(h, config, true, run)) continue;

                // The audit reference excludes only detector-wide noise. The final
                // spectrum below also excludes module-only noise.
                if (!noise.detectorWide) s.previous->Fill(h.time);
                if (noise.Vetoed()) continue;
                // One accepted pulse = one histogram entry, including repeated
                // hits on the same pad. Unmapped channels stay in time spectra.
                s.time->Fill(h.time);
                ++s.acceptedHits;
                hasAccepted = true;
                int pad = geometry.Number(h.channel);
                if (pad >= 0) {
                    s.pad->Fill(pad);
                    s.timePad->Fill(pad, h.time);
                    ++s.mappedHits;
                }
            }
            if (hasAccepted) ++s.acceptedEvents;
        }
        Check(s.time->GetBinContent(0) == 0 && s.time->GetBinContent(config.timeMax + 1) == 0, "Accepted Time outside histogram range");
        return s;
    }


    // =========================================================================
    // SECTION 7: VOLTAGE-SCAN PLOTS AND NOISE AUDIT
    // =========================================================================

    // Legend event totals are distinct Event values before cuts, not the
    // number of histogram entries. The spectra themselves count raw hits.
    inline std::string Label(const Summary& s) {
        return std::to_string(s.voltage) + " V | run " + std::to_string(s.run) + " | " + std::to_string(s.totalEvents)
            + " events";
    }

    // normalized=false: compare raw hit counts/exposure.
    // normalized=true: divide each spectrum by its accepted-hit integral.
    // pads=false uses Time; true uses continuous pad-number occupancy.
    inline TCanvas* Overlay(const std::vector<Summary>& summaries, const std::string& name, bool normalized, bool pads,
        const Config& config) {
        auto c = Canvas("TPC_" + name, name, 1400, 800);
        c->cd();
        c->SetRightMargin(.34);
        c->SetBottomMargin(.15);
        std::vector<TH1D*> drawn;
        double maximum = 0;

        // Edit these colors for the run curves. LineStyle(1) below deliberately
        // keeps all curves solid; one distinct color identifies each voltage.
        const char* colors[] = { "#1559D4", "#D00000", "#008525", "#9A00B5", "#111111", "#E07000", "#009C9C", "#8B4513" };
        for (const auto& s : summaries) if (s.selected) {

            // Draw clones so normalization and style edits cannot modify the
            // authoritative histograms stored in Summary.
            auto h = static_cast<TH1D*>((pads ? s.pad : s.time)->Clone(Unique("overlay").c_str()));
            h->SetDirectory(nullptr);
            h->SetBit(kCanDelete);
            if (normalized && h->Integral() > 0) h->Scale(1. / h->Integral());
            h->SetLineColor(TColor::GetColor(colors[drawn.size() % 8]));
            h->SetLineWidth(2);
            h->SetLineStyle(1);
            maximum = std::max(maximum, h->GetMaximum());
            drawn.push_back(h);
        }

        // Legend coordinates are normalized to the canvas (0..1), not Time
        // or pad coordinates. The large right margin reserves this space.
        auto legend = new TLegend(.68, .60, .99, .94);
        legend->SetTextSize(.023);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetBit(kCanDelete);
        size_t index = 0;
        for (const auto& s : summaries) if (s.selected) {
            auto h = drawn[index++];
            h->SetTitle((name + ";" + (pads ? "Pad number (0-based)" : "Time (native units)") + ";" + (normalized ? "Fraction of accepted hits" : "Accepted hits")).c_str());
            h->SetMinimum(0);
            h->SetMaximum(maximum * 1.08);
            h->Draw(index == 1 ? "HIST" : "HIST SAME");
            legend->AddEntry(h, Label(s).c_str(), "l");
        }
        legend->Draw();
        TLatex text;
        text.SetNDC();
        text.SetTextSize(.022);
        text.DrawLatex(.1, .035, config.CutText(true).c_str());
        c->Modified();
        return c;
    }

    // This diagnostic shows ALL events before hit cuts. A log y-axis makes
    // rare flooding visible. The dashed vertical line is only the global
    // detector-wide threshold; veto totals also include module flooding.
    inline TCanvas* NoiseOccupancy(const std::vector<Summary>& summaries, const Config& config) {
        int rows = (summaries.size() + 3) / 4;
        auto c = Canvas("TPC_noise_occupancy", "Event occupancy before cuts", 1600, rows * 350);
        c->Divide(4, rows);
        for (size_t i = 0; i < summaries.size(); ++i) {
            c->cd(i + 1);
            gPad->SetLogy();
            const auto& s = summaries[i];
            s.occupancy->SetMinimum(.5);
            s.occupancy->SetMaximum(std::max(2., s.occupancy->GetMaximum() * 2));
            s.occupancy->SetTitle(("Run " + std::to_string(s.run) + " | vetoed " + std::to_string(s.vetoes.size())
                + ";Distinct channels/event;Events").c_str());
            s.occupancy->DrawCopy("HIST");
            TLine line(config.noiseChannels - .5, .5, config.noiseChannels - .5, s.occupancy->GetMaximum());
            line.SetLineColor(kRed + 1);
            line.SetLineStyle(2);
            line.DrawClone();
        }
        return c;
    }

    // Run-specific historical audit, not a generic classifier or track fit.
    // Left: mutually exclusive kept/global/module/special-event groups.
    // Right: identical hit cuts, comparing global-only versus all vetoes.
    // Event10743 is annotated because its synchronous module burst has no
    // hits passing the voltage channel selection. Its ID/text are specific
    // to this dataset; update them if repurposing the audit for another run.
    inline TCanvas* Audit1610(const Summary& s) {
        auto c = Canvas("TPC_run1610_audit", "Run 1610 noise audit", 1600, 650);
        c->Divide(2, 1);
        c->cd(1);
        double ymax = 1;
        for (const auto& p : s.audit) ymax = std::max(ymax, double(p.hits));
        auto frame = gPad->DrawFrame(0, -ymax * .03, 520, ymax * 1.15, "Final event classification;Distinct channels/event (before cuts);Hits passing time/channel cuts");
        (void) frame;
        std::vector<TGraph*> groups;
        for (int i = 0; i < 4; ++i) {
            auto graph = new TGraph();
            graph->SetName(Unique("events").c_str());
            graph->SetBit(kCanDelete);
            groups.push_back(graph);
        }
        for (const auto& p : s.audit) {
            int group = !p.vetoed ? 0 : p.detectorWide ? 1 : p.event == 10743 ? 3 : 2;
            groups[group]->SetPoint(groups[group]->GetN(), p.channels, p.hits);
        }

        // Graph index -> class: 0 kept, 1 detector-wide veto, 2 module-only
        // veto, 3 the separately annotated event10743. Markers distinguish them.
        const int colors[] = { kBlue + 1, kRed + 1, kOrange + 7, kViolet + 1 }, markers[] = { 20, 5, 22, 29 };
        const char* names[] = { "Kept (passes noise checks)", "Vetoed: detector-wide flooding", "Vetoed: module flooding",
            "Vetoed: synchronous event10743" };
        auto legend = new TLegend(.12, .70, .79, .94);
        legend->SetTextSize(.025);
        legend->SetBorderSize(0);
        legend->SetBit(kCanDelete);
        for (int i = 0; i < 4; ++i) {
            auto g = groups[i];
            g->SetMarkerColor(colors[i]);
            g->SetMarkerStyle(markers[i]);
            g->SetMarkerSize(i == 0 ? .3 : 1.3);
            g->Draw("P SAME");
            legend->AddEntry(g, (std::string(names[i]) + " (" + std::to_string(g->GetN()) + ")").c_str(), "p");
        }
        legend->Draw();
        if (groups[3]->GetN()) {
            TLatex note;
            note.SetTextSize(.025);
            note.SetTextColor(kViolet + 1);
            note.DrawLatex(140, ymax * .16, "Event10743: 61 synchronous channels in module5");
            note.DrawLatex(140, ymax * .11, "0 hits pass the voltage time/channel cuts");
        }
        c->cd(2);
        auto previous = static_cast<TH1D*>(s.previous->Clone(Unique("previous_audit").c_str()));
        previous->SetDirectory(nullptr);
        previous->SetBit(kCanDelete);
        previous->SetLineColor(kRed + 1);
        previous->SetLineWidth(2);
        previous->SetTitle("Effect of module-flooding vetoes;Time (native units);Accepted hits");
        previous->Draw("HIST");
        auto final = static_cast<TH1D*>(s.time->Clone(Unique("final_audit").c_str()));
        final->SetDirectory(nullptr);
        final->SetBit(kCanDelete);
        final->SetLineColor(kBlue + 1);
        final->SetLineWidth(2);
        final->Draw("HIST SAME");
        auto l = new TLegend(.45, .78, .94, .94);
        l->SetTextSize(.027);
        l->SetBorderSize(0);
        l->SetBit(kCanDelete);
        l->AddEntry(previous, "Detector-wide veto only", "l");
        l->AddEntry(final, "All noise vetoes (final)", "l");
        l->Draw();
        return c;
    }


    // =========================================================================
    // SECTION 8: FULL VOLTAGE WORKFLOW: READ, SELECT RUNS, EXPORT, PLOT
    // =========================================================================

    // This command recomputes noise evidence from raw data; it does NOT
    // merely reuse the bundled veto list. Its exported rejected_events.csv
    // becomes the preferred ledger for later event browsing/saving.
    inline void Voltage(const Config& config, const Geometry& geometry) {
        // Edit YOUR copy of the run/voltage table in analysis_root/auxiliary/.
        // Conflicting voltage assignments for one run are excluded, never guessed.

        // 8.1 Read run/voltage assignments from your auxiliary table. One run
        // may not have two different voltages; such conflicts get voltage=-1.
        std::ifstream table(config.Auxiliary("datafiles-voltage.txt"));
        if (!table) throw std::runtime_error("Missing voltage table");
        std::map<int, std::set<int>> assignments;
        std::string line;
        while (std::getline(table, line)) {
            line = line.substr(0, line.find('#'));
            std::stringstream ss(line);
            int run, v;
            if (ss >> run >> v) assignments[run].insert(v);
        }
        std::vector<Summary> summaries;
        for (const auto& item : assignments) {
            int run = item.first;
            if (!fs::exists(config.input / ("ntuple" + std::to_string(run) + ".root"))) {
                std::cout << "Missing run " << run << '\n';
                continue;
            }
            int voltage = item.second.size() == 1 ? *item.second.begin() : -1;
            if (voltage < 0) std::cout << "Run " << run << " has conflicting voltages; it will not enter overlays\n";
            auto summary = Analyze(run, voltage, config, geometry);
            std::cout << "Run " << run << ": " << summary.totalEvents << " events, " << summary.vetoes.size()
                << " vetoes, " << summary.acceptedHits << " accepted hits\n";
            summaries.push_back(std::move(summary));
        }
        if (summaries.empty()) throw std::runtime_error("No listed ROOT files found");

        // 8.2 Choose one eligible file per voltage. Tie-break order: most raw
        // events, then most clean events, then the lower run number.
        std::map<int, size_t> winners;
        // Most EVENTS means distinct Event values before cuts, not hit entries.
        for (size_t i = 0; i < summaries.size(); ++i) {
            auto& s = summaries[i];
            if (s.voltage < 0 || config.excludedRuns.count(s.run)) continue;
            auto it = winners.find(s.voltage);
            if (it == winners.end()) winners[s.voltage] = i;
            else {
                const auto& old = summaries[it->second];
                auto clean = s.totalEvents - s.vetoes.size(), oldClean = old.totalEvents - old.vetoes.size();
                if (s.totalEvents > old.totalEvents || (s.totalEvents == old.totalEvents && (clean > oldClean
                    || (clean == oldClean && s.run < old.run)))) it->second = i;
            }
        }
        if (winners.empty()) throw std::runtime_error("No eligible unambiguous voltage runs");
        for (auto& item : winners) summaries[item.second].selected = true;
        // Plotting and export order follow voltage, matching the Python overlay.

        // 8.3 Use ascending voltage for overlays and selected histogram tables.
        // The summary/veto CSVs are separately ordered by run below.
        std::stable_sort(summaries.begin(), summaries.end(), [](const Summary& a, const Summary& b) {
            if (a.selected != b.selected) return a.selected > b.selected;
            return a.selected ? a.voltage < b.voltage : a.run < b.run;
        });
        auto out = config.output / "voltage";
        fs::create_directories(out);

        // 8.4 Keep human/machine-readable CSVs beside the figures. The long
        // header strings define the exported columns; their order is intentional.
        std::ofstream summaryCSV(out / "run_summary.csv"), vetoCSV(out / "rejected_events.csv"), histCSV(out / "histograms.csv"),
            padCSV(out / "pad_histograms.csv");
        summaryCSV << "run,voltage,total_events,clean_events,accepted_events,total_hits,base_cut_hits,noise_events,additional_module_noise_events,noise_hits_after_base_cuts,accepted_hits,selected,mapped_accepted_hits,unmapped_accepted_hits\n";
        vetoCSV << "run,event,distinct_channels,total_hits,hits_passing_base_cuts,detector_wide,flooded_modules\n";
        histCSV << "run,voltage,time,count\n";
        padCSV << "run,voltage,pad,count\n";
        // CSV rows remain ordered by run; histogram rows follow selected voltage.
        std::vector<const Summary*> runOrder;
        for (const auto& s : summaries) runOrder.push_back(&s);
        std::sort(runOrder.begin(), runOrder.end(), [](const Summary* a, const Summary* b) {
            return a->run < b->run;
        });
        for (const auto* ptr : runOrder) {
            const auto& s = *ptr;
            long long noiseHits = 0, extra = 0;
            for (const auto& n : s.vetoes) {
                noiseHits += n.baseHits;
                extra += !n.detectorWide;
                vetoCSV << s.run << ',' << n.event << ',' << n.distinct << ',' << n.totalHits << ',' << n.baseHits
                    << ',' << (n.detectorWide ? "True" : "False") << ',' << QuoteCSV(ModuleJSON(n)) << '\n';
            }
            summaryCSV << s.run << ',';
            if (s.voltage >= 0) summaryCSV << s.voltage;
            summaryCSV << ',' << s.totalEvents << ',' << s.totalEvents - s.vetoes.size() << ',' << s.acceptedEvents
                << ',' << s.totalHits << ',' << s.baseHits << ',' << s.vetoes.size() << ',' << extra << ',' << noiseHits
                << ',' << s.acceptedHits << ',' << (s.selected ? "True" : "False") << ',' << s.mappedHits << ','
                << s.acceptedHits - s.mappedHits << '\n';
        }

        // 8.5 Native ROOT output: numerical TH1D/TH2D objects plus canvases.
        // Write only selected runs into this file; all run summaries and veto
        // evidence remain in the CSVs, including unselected/conflicting runs.
        auto rootPath = out / "time_histograms.root";
        TFile file(rootPath.c_str(), "RECREATE");
        TNamed("selection", config.CutText(true).c_str()).Write();
        TNamed("cuts_db", config.cutDBPath.c_str()).Write();
        config.Cuts().Snapshot(out / "cuts_used.ini");
        TNamed("voltage_study_cuts_db", config.studyCutDBPath.c_str()).Write();
        config.StudyCuts().Snapshot(out / "voltage_study_cuts_used.ini");
        TNamed("settings_file", (config.source / "settings.csv").c_str()).Write();
        for (const auto& s : summaries) if (s.selected) {
            file.cd();
            s.time->Write();
            s.pad->Write();
            s.timePad->Write();
            s.occupancy->Write();
            for (int b = 1; b <= config.timeMax; ++b) histCSV << s.run << ',' << s.voltage << ',' << b - 1
                << ',' << std::llround(s.time->GetBinContent(b)) << '\n';
            for (int b = 1; b <= 410; ++b) padCSV << s.run << ',' << s.voltage << ',' << b - 1 << ',' << std::llround(s.pad->GetBinContent(b))
                << '\n';
        }

        // 8.6 Overlay options are (normalized, pads). These three views share
        // the same selected runs and cuts; only representation/normalization changes.
        for (auto spec : std::vector<std::pair<std::string, std::pair<bool, bool>>> {
                 { "time_overlay_counts",     { false, false } },
                 { "time_overlay_normalized", { true,  false } },
                 { "hits_by_pad",             { false, true  } }
             }) {
            auto c = Overlay(summaries, spec.first, spec.second.first, spec.second.second, config);
            SaveCanvas(c, out / spec.first);
            file.cd();
            c->Write();
        }

        // 8.7 Save the full pre-cut occupancy diagnostic, then one Time-vs-Pad
        // panel for each selected voltage and the run1610 audit.
        auto occupancy = NoiseOccupancy(summaries, config);
        SaveCanvas(occupancy, out / "noise_occupancy");
        file.cd();
        occupancy->Write();
        int selected = winners.size();
        auto timePad = Canvas("TPC_time_vs_pad", "Time versus continuous pad number", 1300, selected * 320);
        timePad->Divide(1, selected);
        int position = 0;
        for (const auto& s : summaries) if (s.selected) {
            timePad->cd(++position);
            gPad->SetRightMargin(.12);
            s.timePad->SetTitle(Label(s).c_str());
            s.timePad->DrawCopy("COLZ");
        }
        SaveCanvas(timePad, out / "time_vs_pad");
        file.cd();
        timePad->Write();
        for (const auto& s : summaries) if (s.run == 1610) {
            auto audit = Audit1610(s);
            SaveCanvas(audit, out / "run1610_noise_audit");
            file.cd();
            audit->Write();
        }
        file.Close();
        std::cout << "Voltage outputs: " << out << "\nSelected:";
        for (const auto& item : winners) std::cout << ' ' << item.first << " V";
        std::cout << std::endl;
    }



    // =========================================================================
    // SECTION 8B: RAW BRANCH INSPECTION (NO CUTS BY DEFAULT)
    // =========================================================================
    // Axis choices live in config/raw_binning.csv. The same per-branch
    // definition can later be reused for either axis of a 2D histogram.
    struct RawAxis {
        int bins;double low,high;bool logY;
    };

    // Thousands separators make large annotation counts easier to read.
    inline std::string CountText(long long count) {
        std::string text=std::to_string(count);
        for(int position=int(text.size())-3;position>0;position-=3)
            text.insert(position,",");
        return text;
    }

    inline RawAxis RawBinning(const std::vector<double>& values,const std::string& branch,const Config& config,
                             const std::map<std::string,std::string>& overrides={}) {
        std::ifstream file(config.binningPath);
        if(!file)throw std::runtime_error("Cannot read raw binning file: "+config.binningPath.string());
        std::string line;std::getline(file,line);std::vector<std::string> settings(6);
        while(std::getline(file,line)){auto row=CSV(line);if(row.size()>=6 && row[0]==branch){settings=row;break;}}
        auto option=[&](const std::string& key){auto it=overrides.find(key);return it==overrides.end()?std::string():it->second;};
        if(!option("BinWidth").empty() && option("Bins").empty())settings[2].clear();
        for(const auto& entry:std::map<std::string,int>{{"BinWidth",1},{"Bins",2},{"Min",3},{"Max",4}})
            if(!option(entry.first).empty())settings[entry.second]=option(entry.first);
        double low=0,high=0;
        if(!values.empty()){auto bounds=std::minmax_element(values.begin(),values.end());low=*bounds.first;high=*bounds.second;}
        if(branch=="Event" || branch=="Channel" || branch=="Hit" || branch=="Time" || branch=="Blockhead")low=std::min(0.,low);
        if(branch=="Channel")high=std::max(high,double(config.channelsPerModule*config.readoutModules-1));
        low=std::floor(low)-.5;high=std::ceil(high)+.5;
        if(!settings[3].empty())low=CutDatabase::Number(settings[3]);
        if(!settings[4].empty())high=CutDatabase::Number(settings[4]);
        if(high<=low)throw std::runtime_error("Histogram edges must increase: "+branch);
        int bins=0;
        if(!settings[2].empty()) {
            double count=CutDatabase::Number(settings[2]);
            if(count!=std::floor(count) || count<1 || count>1000000)throw std::runtime_error("Bins must be 1..1000000");
            bins=int(count);
        } else {
            double width=settings[1].empty()?1.:CutDatabase::Number(settings[1]);
            if(width<=0)throw std::runtime_error("BinWidth must be positive");
            double count=std::ceil((high-low)/width);
            if(count<1 || count>1000000)throw std::runtime_error("Too many bins; edit config/raw_binning.csv for "+branch);
            bins=int(count);high=low+bins*width;
        }
        auto scale=Lower(settings[5]);
        if(!scale.empty() && scale!="true" && scale!="false" && scale!="1" && scale!="0")throw std::runtime_error("LogY must be true/false");
        return {bins,low,high,scale=="true" || scale=="1"};
    }

    inline void ShowCuts(const Config& config,int run) {
        const auto& db=config.Cuts();std::cout<<"Manual rejection DB: "<<db.path<<'\n';
        for(const auto& rule:db.rules)if(db.Active(rule,run))
            std::cout<<'['<<rule.section<<"] line "<<rule.line<<": REJECT "<<rule.expression<<(rule.note.empty()?"":" # "+rule.note)<<'\n';
        std::cout<<"\nRules listed only; this command does not open plots.\n"
                 <<"To apply detector cuts and show 1D + 2D: TPC(\"raw-clean\","<<run<<");\n"
                 <<"Keep any custom input/DB overrides when plotting.\n";
    }

    inline void Raw(const Config& config,int run,bool cleaned=false,const std::string& onlyBranch="",const RunData* sharedData=nullptr) {
        // Mapping always converts electronics Channel to physical Pad.
        // It does not reject rows: unconnected channels become Pad=-1
        // (underflow). Physics cuts/vetoes run only when cleaned=true.
        RunData ownedData;if(!sharedData)ownedData=ReadRun(run,config);
        const RunData& data=sharedData?*sharedData:ownedData;
        std::set<int> vetoes;
        if(cleaned){ValidateCuts(data,config);vetoes=ReadVetoes(run,config);}
        std::vector<std::string> branches;
        // Six mapped analysis quantities form the overview. Channel is an
        // extra, standalone electronics diagnostic; Blockhead is omitted.
        const std::string requested=onlyBranch;
        for(const char* name:{"Event","Pad","Hit","Time","Charge","Pedestal","Channel"}) {
            const std::string field=name;
            bool available=data.branches.count(field=="Pad"?"Channel":field);
            if(available && (requested.empty() || requested==field))branches.emplace_back(field);
        }
        if(branches.empty())throw std::runtime_error("No matching raw branch: "+onlyBranch);
        auto out=config.output/"raw"/("run"+std::to_string(run));if(cleaned)out/="cleaned";
        // A focused inspection must not replace the full run's ROOT file.
        if(!requested.empty())out/="branches_"+requested;
        fs::create_directories(out);
        std::vector<const Hit*> selected;
        for(const auto& h:data.hits)if(!cleaned || (!vetoes.count(h.event) && PassHit(h,config,false,run)))selected.push_back(&h);
        auto unmapped=std::count_if(data.hits.begin(),data.hits.end(),[&](const Hit& h){return HitValue(h,"Pad",config)<0;});
        auto selectedUnmapped=std::count_if(selected.begin(),selected.end(),[&](const Hit* h){return HitValue(*h,"Pad",config)<0;});
        int summaryCount=requested.empty()?int(std::count_if(branches.begin(),branches.end(),[](const std::string& b){return b!="Channel";})):0;
        const std::string canvasStem="TPC_raw_"+std::to_string(run)+(cleaned?"_clean":"");
        TCanvas* canvas=nullptr;
        if(summaryCount) {
            int columns=std::min(3,summaryCount),rows=(summaryCount+columns-1)/columns;
            canvas=Canvas(canvasStem,cleaned?"After DB cuts and event vetoes":"UNCUT raw data: all hit rows",columns*500,rows*360);
            canvas->Divide(columns,rows);
        }
        TFile output((out/"raw_histograms.root").c_str(),"RECREATE");
        TNamed("selection",cleaned?config.CutText(false,run).c_str():"UNCUT: all raw hit rows; no vetoes").Write();
        auto mappingPath=config.source/"auxiliary/pad_map.csv";
        TNamed("mapping",(mappingPath.string()+": Pad from ROOT Channel; FADC=Channel+1; Pad=-1 is unmapped (underflow)").c_str()).Write();
        fs::copy_file(mappingPath,out/"pad_map_used.csv",fs::copy_options::overwrite_existing);
        std::ofstream summary(out/"raw_summary.csv");
        summary<<"branch,total_rows,selected_rows,bins,low,high,underflow,overflow,log_y,selection,unmapped_rows,selected_unmapped_rows\n";
        for(size_t i=0;i<branches.size();++i) {
            const auto& branch=branches[i];std::vector<double> values;values.reserve(selected.size());
            for(const auto* h:selected)values.push_back(HitValue(*h,branch,config));
            // Use the uncut extent for both views: removed tail bins remain
            // visible as zeros instead of shrinking the cleaned axis.
            std::vector<double> rawAxisValues;rawAxisValues.reserve(data.hits.size());
            for(const auto& h:data.hits)rawAxisValues.push_back(HitValue(h,branch,config));
            auto axis=RawBinning(rawAxisValues,branch,config);
            // Put run and selection on the exported figure as well as in
            // the filename, so an isolated image remains identifiable.
            std::string title = "Run " + std::to_string(run) + " | "
                + (cleaned ? "CLEANED" : "UNCUT") + " | " + branch;
            std::string axisLabel=branch=="Pad"?"Pad number (0-based)":branch+" (native values)";
            if(branch=="Channel")axisLabel="ROOT Channel (0-based, fADC = Channel + 1)";
            TH1D histogram(("raw_"+branch).c_str(),(title+";"+axisLabel+";Raw Hits").c_str(),axis.bins,axis.low,axis.high);
            histogram.SetDirectory(nullptr);
            for(double value:values)histogram.Fill(value);
            auto under=std::llround(histogram.GetBinContent(0)),over=std::llround(histogram.GetBinContent(axis.bins+1));
            std::vector<std::string> annotations={"Underflow: "+CountText(under)+"   Overflow: "+CountText(over)};
            if(branch=="Pad" || branch=="Channel")
                annotations.push_back("Hits on unconnected channels: "+CountText(selectedUnmapped));
            // DrawCopy gives each canvas its own histogram, so closing the
            // output file or zooming one view does not change another view.
            auto draw = [&](TCanvas* target, int pad) {
                if(pad)target->cd(pad);else target->cd();
                gPad->SetBottomMargin(.14);
                gPad->SetLogy(axis.logY && histogram.GetMaximum()>0);
                if(axis.logY && histogram.GetMaximum()>0)histogram.SetMinimum(.7);
                histogram.SetLineColor(kBlue+1);
                histogram.SetLineWidth(1);
                histogram.DrawCopy("HIST");
                TLatex note;note.SetNDC();note.SetTextSize(.028);
                for(size_t line=0;line<annotations.size();++line)
                    note.DrawLatex(.13,.84-.045*line,annotations[line].c_str());
            };
            if(canvas && branch!="Channel")draw(canvas,int(i+1));
            auto individual=Canvas(canvasStem+"_"+branch,title,1000,700);
            draw(individual,0);
            SaveCanvas(individual,out/("raw_"+branch));
            output.cd();individual->Write();
            // All ROOT bins, including the last, are [low,high). A user
            // plotting range never silently discards under/overflow counts.
            std::ofstream bins(out/(branch+"_bins.csv"));bins<<"low_inclusive,high_exclusive,count\n";bins<<std::setprecision(17);
            for(int b=1;b<=axis.bins;++b)bins<<histogram.GetXaxis()->GetBinLowEdge(b)<<','<<histogram.GetXaxis()->GetBinUpEdge(b)<<','<<std::llround(histogram.GetBinContent(b))<<'\n';
            summary<<branch<<','<<data.hits.size()<<','<<selected.size()<<','<<axis.bins<<','<<axis.low<<','<<axis.high<<','<<under<<','<<over<<','<<(axis.logY?"True":"False")<<','<<(cleaned?"clean":"raw")<<','<<unmapped<<','<<selectedUnmapped<<'\n';
            output.cd();histogram.Write();
        }
        if(canvas){SaveCanvas(canvas,out/"raw_overview");output.cd();canvas->Write();}
        output.Close();
        if(cleaned)config.Cuts().Snapshot(out/"cuts_used.ini");
        std::cout<<"Run "<<run<<": "<<data.hits.size()<<" raw rows, "<<selected.size()<<" plotted rows, selection="<<(cleaned?"clean":"raw")<<". Saved in "<<out<<std::endl;
    }

    // =========================================================================
    // SECTION 8C: TWO-DIMENSIONAL RAW/CLEANED INSPECTION
    // =========================================================================
    struct Raw2DPlot {
        std::string name,x,y;
        bool logZ;
        std::map<std::string,std::string> xOptions,yOptions;
    };

    // The CSV is the edit point for axis orientation, binning and colour
    // scale. Empty overrides inherit the established 1D axis settings.
    inline std::vector<Raw2DPlot> Load2DPlots(const Config& config) {
        std::ifstream input(config.plots2DPath);
        if(!input)throw std::runtime_error("Cannot read 2D plot definitions: "+config.plots2DPath.string());
        std::string line;std::getline(input,line);auto header=CSV(line);
        std::vector<Raw2DPlot> plots;std::set<std::string> names;int number=1;
        while(std::getline(input,line)) {
            ++number;if(Trim(line).empty())continue;
            try {
                auto cells=CSV(line);std::map<std::string,std::string> row;
                for(size_t i=0;i<header.size();++i)row[header[i]]=i<cells.size()?cells[i]:"";
                Raw2DPlot plot;plot.name=row["Name"];
                if(!std::regex_match(plot.name,std::regex("[A-Za-z_][A-Za-z_0-9]*")) || !names.insert(plot.name).second)
                    throw std::runtime_error("Plot names must be unique filename-safe identifiers");
                plot.x=CutDatabase::ParseCondition(row["X"]+" == 0").field;
                plot.y=CutDatabase::ParseCondition(row["Y"]+" == 0").field;
                auto scale=Lower(row["LogZ"]);
                if(scale!="true" && scale!="false" && scale!="1" && scale!="0")throw std::runtime_error("LogZ must be true/false");
                plot.logZ=scale=="true" || scale=="1";
                for(const char* key:{"BinWidth","Bins","Min","Max"}) {
                    plot.xOptions[key]=row[std::string("X")+key];
                    plot.yOptions[key]=row[std::string("Y")+key];
                }
                plots.push_back(plot);
            } catch(const std::exception& error) {
                throw std::runtime_error(config.plots2DPath.string()+":"+std::to_string(number)+": "+error.what());
            }
        }
        if(plots.empty())throw std::runtime_error("No 2D plot definitions");
        return plots;
    }

    inline std::string Raw2DLabel(const std::string& field,const RawAxis& axis) {
        std::string label=field=="Pad"?"Pad number (0-based)":field=="Event"?"Event ID":field=="Hit"?"Hit (pulse index)":field+" (native values)";
        std::ostringstream text;text<<label<<" [bin width "<<(axis.high-axis.low)/axis.bins<<']';return text.str();
    }

    inline void Raw2D(const Config& config,int run,bool cleaned=false,const std::string& onlyPlot="",const RunData* sharedData=nullptr) {
        RunData ownedData;if(!sharedData)ownedData=ReadRun(run,config);
        const RunData& data=sharedData?*sharedData:ownedData;
        std::set<int> vetoes;
        if(cleaned){ValidateCuts(data,config);vetoes=ReadVetoes(run,config);}
        std::vector<const Hit*> selected;
        for(const auto& hit:data.hits)
            if(!cleaned || (!vetoes.count(hit.event) && PassHit(hit,config,false,run)))selected.push_back(&hit);
        auto unmapped=std::count_if(selected.begin(),selected.end(),[&](const Hit* h){return HitValue(*h,"Pad",config)<0;});
        auto definitions=Load2DPlots(config);
        if(!onlyPlot.empty()) {
            definitions.erase(std::remove_if(definitions.begin(),definitions.end(),[&](const Raw2DPlot& p){return p.name!=onlyPlot;}),definitions.end());
            if(definitions.empty())throw std::runtime_error("Unknown 2D plot: "+onlyPlot);
        }
        auto out=config.output/"raw"/("run"+std::to_string(run));
        if(cleaned)out/="cleaned";
        out/="2d";
        if(!onlyPlot.empty())out/="plots_"+onlyPlot;
        fs::create_directories(out);
        fs::copy_file(config.Auxiliary("pad_map.csv"),out/"pad_map_used.csv",fs::copy_options::overwrite_existing);
        fs::copy_file(config.plots2DPath,out/"raw_2d_plots_used.csv",fs::copy_options::overwrite_existing);
        fs::copy_file(config.binningPath,out/"raw_binning_used.csv",fs::copy_options::overwrite_existing);
        const std::string profile=cleaned?"CLEANED":"UNCUT";
        const std::string stem="TPC_raw2d_"+std::to_string(run)+(cleaned?"_clean":"");
        TCanvas* summary=nullptr;
        if(onlyPlot.empty()) {
            int columns=std::min(2,int(definitions.size())),rows=(definitions.size()+columns-1)/columns;
            summary=Canvas(stem,"Run "+std::to_string(run)+" | "+profile+" | 2D summary",650*columns,500*rows);
            summary->Divide(columns,rows);
        }
        TFile file((out/"raw_2d_histograms.root").c_str(),"RECREATE");
        TNamed("selection",cleaned?config.CutText(false,run).c_str():"UNCUT: all raw hits; no vetoes").Write();
        TNamed("mapping","Pad from established mapping; Pad=-1 for unconnected channels, retained in flow cells").Write();
        if(cleaned){TNamed("cuts_db",config.cutDBPath.c_str()).Write();config.Cuts().Snapshot(out/"cuts_used.ini");}
        std::ofstream report(out/"raw_2d_summary.csv");
        report<<"name,x,y,total_rows,selected_rows,x_bins,y_bins,x_low,x_high,y_low,y_high,entries,visible_hits,outside_hits,nonfinite_hits,x_underflow,x_overflow,y_underflow,y_overflow,selected_unmapped_hits,log_z\n";
        report<<std::setprecision(17);
        for(size_t index=0;index<definitions.size();++index) {
            const auto& plot=definitions[index];
            for(const auto& field:{plot.x,plot.y}) {
                auto branch=(field=="Pad" || field=="FADC")?"Channel":field;
                if(!data.branches.count(branch))throw std::runtime_error(plot.name+": missing branch "+branch);
            }
            // Axis extents come from uncut data, even in cleaned mode.
            std::vector<double> rawX,rawY;rawX.reserve(data.hits.size());rawY.reserve(data.hits.size());
            for(const auto& hit:data.hits){rawX.push_back(HitValue(hit,plot.x,config));rawY.push_back(HitValue(hit,plot.y,config));}
            auto x=RawBinning(rawX,plot.x,config,plot.xOptions),y=RawBinning(rawY,plot.y,config,plot.yOptions);
            if(size_t(x.bins+2)*size_t(y.bins+2)>4000000)
                throw std::runtime_error(plot.name+": more than 4 million 2D cells; increase BinWidth or narrow limits in raw_2d_plots.csv");
            std::string title="Run "+std::to_string(run)+" | "+profile+" | "+plot.y+" vs "+plot.x;
            std::string labels=title+";"+Raw2DLabel(plot.x,x)+";"+Raw2DLabel(plot.y,y)+";Raw Hits";
            TH2D histogram(("raw2d_"+plot.name).c_str(),labels.c_str(),x.bins,x.low,x.high,y.bins,y.low,y.high);
            histogram.SetDirectory(nullptr);histogram.SetStats(false);histogram.SetOption("COLZ");
            long long nonfinite=0;
            for(const auto* hit:selected) {
                double xv=HitValue(*hit,plot.x,config),yv=HitValue(*hit,plot.y,config);
                if(std::isfinite(xv) && std::isfinite(yv))histogram.Fill(xv,yv);else ++nonfinite;
            }
            auto visible=std::llround(histogram.Integral(1,x.bins,1,y.bins));
            auto outside=std::llround(histogram.GetEntries())-visible;
            std::vector<std::string> notes={"Outside plotted ranges: "+CountText(outside)};
            if(plot.x=="Pad" || plot.y=="Pad")notes.push_back("Hits on unconnected channels: "+CountText(unmapped));
            double maximum=histogram.GetMaximum();bool logZ=plot.logZ && maximum>0;
            if(logZ){histogram.SetMinimum(1);histogram.SetMaximum(std::max(2.,maximum));}
            auto draw=[&](TCanvas* canvas,int pad) {
                if(pad)canvas->cd(pad);else canvas->cd();
                gPad->SetLeftMargin(.13);gPad->SetRightMargin(.17);gPad->SetBottomMargin(.15);gPad->SetTopMargin(.19);gPad->SetLogz(logZ);
                histogram.DrawCopy("COLZ");
                TLatex text;text.SetNDC();text.SetTextSize(.028);
                for(size_t line=0;line<notes.size();++line)text.DrawLatex(.13,.91-.045*line,notes[line].c_str());
            };
            if(summary)draw(summary,int(index+1));
            auto individual=Canvas(stem+"_"+plot.name,title,1000,800);draw(individual,0);
            SaveCanvas(individual,out/("raw2d_"+plot.name));file.cd();histogram.Write();individual->Write();

            // Sparse export includes all NONZERO cells, including flows.
            // Low/high are infinite on the open ends of under/overflow bins.
            auto bounds=[](const RawAxis& axis,int bin) {
                double low=bin==0?-std::numeric_limits<double>::infinity():axis.low+(bin-1)*(axis.high-axis.low)/axis.bins;
                double high=bin==axis.bins+1?std::numeric_limits<double>::infinity():axis.low+bin*(axis.high-axis.low)/axis.bins;
                return std::make_pair(low,high);
            };
            std::ofstream cells(out/(plot.name+"_bins.csv"));cells<<std::setprecision(17);
            cells<<"x_bin,y_bin,x_low_inclusive,x_high_exclusive,y_low_inclusive,y_high_exclusive,count\n";
            for(int iy=0;iy<=y.bins+1;++iy)for(int ix=0;ix<=x.bins+1;++ix) {
                auto count=std::llround(histogram.GetBinContent(ix,iy));if(!count)continue;
                auto xb=bounds(x,ix),yb=bounds(y,iy);
                cells<<ix<<','<<iy<<','<<xb.first<<','<<xb.second<<','<<yb.first<<','<<yb.second<<','<<count<<'\n';
            }
            report<<plot.name<<','<<plot.x<<','<<plot.y<<','<<data.hits.size()<<','<<selected.size()<<','<<x.bins<<','<<y.bins<<','<<x.low<<','<<x.high<<','<<y.low<<','<<y.high<<','<<std::llround(histogram.GetEntries())<<','<<visible<<','<<outside<<','<<nonfinite<<','
                <<std::llround(histogram.Integral(0,0,0,y.bins+1))<<','<<std::llround(histogram.Integral(x.bins+1,x.bins+1,0,y.bins+1))<<','<<std::llround(histogram.Integral(0,x.bins+1,0,0))<<','<<std::llround(histogram.Integral(0,x.bins+1,y.bins+1,y.bins+1))<<','<<unmapped<<','<<(plot.logZ?"True":"False")<<'\n';
        }
        if(summary){SaveCanvas(summary,out/"raw2d_overview");file.cd();summary->Write();}
        file.Close();
        std::cout<<"Run "<<run<<": "<<profile<<", "<<selected.size()<<" hits; "<<definitions.size()<<" 2D views saved to "<<out<<std::endl;
    }

    // =========================================================================
    // SECTION 9: USER INFORMATION, QUICK CHECKS AND THE TEXT MENU
    // =========================================================================

    // Print the actual paths and active profiles for this invocation.
    // Useful after editing settings or when running from another directory.
    inline void Info(const Config& config, const Geometry& geometry) {
        std::cout << "C++/ROOT analysis\nSettings: " << config.source / "settings.csv" << "\nInput: "
            << config.input << "\nOutput: " << config.output << "\nMapping: " << geometry.pads.size() << " pads, medium 0-249 below large 250-409\nEvent cuts: "
            << config.CutText() << "\nCut DB: " << config.cutDBPath << "\nIndependent voltage-study DB: " << config.studyCutDBPath << "\nVeto ledger: " << config.VetoFile()
            << std::endl;
    }

    // Regression fixtures for the present detector/setup. This is a quick
    // consistency check, not a proof that a physics selection is optimal.
    // Deliberate geometry/cut changes may require updating the fixtures.
    inline void SelfCheck(const Config& config, const Geometry& geometry) {
        for (auto anchor : std::map<int, int> { { 153, 0 }, { 151, 1 }, { 149, 2 }, { 147, 3 }, { 16, 24 }, { 163,
            25 }, { 505, 250 }, { 576, 395 }, { 568, 399 }, { 431, 409 } }) Check(geometry.Number(anchor.first - 1) == anchor.second,
            "Wiring anchor mismatch");
        Check(geometry.Number(0) == -1 && geometry.Number(1) == -1 && geometry.Number(2) == -1, "Unconnected channels mapped unexpectedly");
        // The working DB is user-editable. Validate its schema, rather than
        // requiring historical Time/Channel exclusions to stay enabled.
        // Fixed cut-boundary fixtures live in tests/RunCutsCheck.C.
        config.Cuts().Validate(0,{"Event","Channel","Hit","Time","Charge","Pedestal","Blockhead"});

        // Create a synthetic synchronous module burst without needing an input
        // ROOT file. Noise must be recognized before ordinary hit selection.
        RunData synthetic;
        synthetic.run = 0;
        for (int c = 0; c < config.moduleChannels; ++c) {
            synthetic.hits.push_back({ 1, c, 1, 40, 100 });
            synthetic.events[1].push_back(c);
        }
        auto noise = ClassifyEvent(synthetic, 1, synthetic.events.at(1), config);
        Check(noise.Vetoed(), "Synchronous module burst not vetoed");

        // For the current 100x80-mm plane, every 1-mm cell must belong to
        // exactly one pad. This detects gaps, overlap, and incorrect stacking.
        std::vector<int> coverage(100 * 80);
        for (const auto& p : geometry.pads)
            for (int y = std::lround(p.y); y < std::lround(p.y + p.height); ++y)
                for (int x = std::lround(p.x); x < std::lround(p.x + p.width); ++x) {
                    Check(x >= 0 && x < 100 && y >= 0 && y < 80, "Pad outside board");
                    ++coverage[y * 100 + x];
                }
        Check(std::all_of(coverage.begin(), coverage.end(), [](int n) {
            return n == 1;
        }), "Pad footprints have gaps or overlap");
        std::cout << "PASS: wiring anchors, disconnected channels, detector cuts, separate study restrictions, module noise and complete 410-pad board coverage."
            << std::endl;
    }

    // Translate simple keyboard choices into the SAME command names used
    // by TPC("...", ...). Suggested event IDs are conveniences only; they
    // are not event cuts and do not affect the voltage scan.
    inline std::string Menu(int& run, int& event) {
        std::cout << "\nTPC ROOT analysis\n 1 Browse events (2D)\n 2 Save one event (2D + 3D)\n 3 Recreate 1615/1616 examples\n 4 Independent voltage study\n 5 Show settings/paths\n 6 Check mapping/cuts\n 7 Inspect raw data (1D + 2D)\n 8 Show manual rules for a run\n 0 Quit\nChoose [0]: "
            << std::flush;
        std::string choice;
        if (!std::getline(std::cin, choice)) return "quit";
        choice = Trim(choice);
        if (choice.empty() || choice == "0") return "quit";
        if (choice == "1" || choice == "2" || choice == "7" || choice == "8") {
            std::cout << "Run [1616]: " << std::flush;
            std::string text;
            std::getline(std::cin, text);
            run = Trim(text).empty() ? 1616 : std::stoi(text);
            if (choice == "7") return "raw";
            if (choice == "8") return "cuts";
            int suggested = run == 1615 ? 42535 : run == 1616 ? 5584 : run == 1610 ? 687 : -1;
            std::cout << "Event [" << (suggested >= 0 ? std::to_string(suggested) : "first accepted") << "]: "
                << std::flush;
            std::getline(std::cin, text);
            event = Trim(text).empty() ? suggested : std::stoi(text);
            return choice == "1" ? "event" : "save";
        }
        if (choice == "3") return "examples";
        if (choice == "4") return "voltage";
        if (choice == "5") return "info";
        if (choice == "6") return "check";
        throw std::runtime_error("Choose 0-8");
    }
}


// =========================================================================
// SECTION 10: PUBLIC ROOT ENTRY POINT / COMMAND DISPATCH
// =========================================================================

// Arguments, in order:
//   command: menu, event, event3d, save, examples, voltage, raw, raw:Branch,
//            raw-clean, raw2d, raw2d-clean, cuts, info, check.
//   raw/raw-clean show BOTH dimensions. Add --1d or --2d inside command
//   to select a group, e.g. TPC("raw --1d",1616). --both is explicit default.
//   run/event: event commands use these; event=-1 chooses the first eligible.
//   inputDirectory: folder containing ntupleNNNN.root (empty = default).
//   outputDirectory: generated output root (empty = results/root/).
//   vetoCSV: explicit automatic-event ledger (empty = generated or snapshot).
//   cutsDB: manual DB override (empty = ../config/run_cuts.ini).
//   studyCutsDB: independent voltage-study DB override (voltage command only).
//
// Config is reconstructed on each call, so settings.csv edits take effect
// immediately. Manual cuts/binning are shared through config/; ROOT noise
// settings and source remain independent from the Python implementation.
void TPC(const char* command = "menu", int run = 1616, int event = -1, const char* inputDirectory = "", const char* outputDirectory = "",
    const char* vetoCSV = "", const char* cutsDB = "", const char* studyCutsDB = "") {
    using namespace PersonalTPC;

    // Load configuration and the established pad map once per command.
    Config config(inputDirectory, outputDirectory, vetoCSV, cutsDB, studyCutsDB);
    Geometry geometry(config);

    // Apply your palette before any histogram/canvas is drawn. Edit the
    // stops/RGB arrays in set_color_env() above to personalize it. Calling
    // SetPalette(kViridis) afterwards would replace your custom palette.
    // Explicit voltage curve colours are defined separately in Overlay().
    gStyle->SetOptStat(0);
    set_color_env();
    std::string action = command ? command : "menu";
    // A single entry point and small flags keep routine inspection simple.
    // Existing raw2d and raw:Branch shortcuts remain available.
    std::istringstream arguments(action);
    arguments>>action;
    std::string flag,dimensions="both";
    bool dimensionalFlag=false;
    while(arguments>>flag) {
        if(dimensionalFlag || (flag!="--1d" && flag!="--2d" && flag!="--both"))
            throw std::runtime_error("Choose one inspection flag: --1d, --2d, or --both");
        dimensions=flag=="--1d"?"1d":flag=="--2d"?"2d":"both";
        dimensionalFlag=true;
    }
    if(dimensionalFlag && action!="raw" && action!="raw-clean")
        throw std::runtime_error("Dimension flags apply to raw or raw-clean commands");
    if (action == "menu") action = Menu(run, event);
    if (action == "quit") return;
    if (action == "info") {
        Info(config, geometry);
        return;
    }
    if (action == "check") {
        SelfCheck(config, geometry);
        return;
    }
    if (action == "event" || action == "event3d") {
        Browse(config, geometry, run, event, action == "event3d");
        return;
    }
    if (action == "cuts") { ShowCuts(config, run); return; }
    if(action=="raw" || action=="raw-clean") {
        // Read the ntuple once. Both drawing groups share input and the
        // same Config, including the cached detector-cut database.
        auto data=ReadRun(run,config);bool cleaned=action=="raw-clean";
        if(dimensions!="2d")Raw(config,run,cleaned,"",&data);
        if(dimensions!="1d")Raw2D(config,run,cleaned,"",&data);
        return;
    }
    if (action == "raw2d" || action.rfind("raw2d:", 0) == 0) {
        Raw2D(config, run, false, action == "raw2d" ? "" : action.substr(6)); return;
    }
    if (action == "raw2d-clean" || action.rfind("raw2d-clean:", 0) == 0) {
        Raw2D(config, run, true, action == "raw2d-clean" ? "" : action.substr(12)); return;
    }
    if (action.rfind("raw:", 0) == 0) {
        Raw(config, run, false, action.substr(4)); return;
    }
    if (action.rfind("raw-clean:", 0) == 0) {
        Raw(config, run, true, action.substr(10)); return;
    }
    // Batch exports do not require windows. Restore the caller's ROOT mode
    // afterwards so a later interactive command can still open canvases.

    // RAII guard: even if export throws, the destructor restores the
    // previous batch setting instead of leaving the user session headless.
    struct BatchGuard {
        bool previous;
        BatchGuard() : previous(gROOT->IsBatch()) {
            gROOT->SetBatch(true);
        }
        ~BatchGuard() {
            gROOT->SetBatch(previous);
        }
    } guard;
    if (action == "save") SaveEvent(config, geometry, run, event);
    else if (action == "examples") Examples(config, geometry);
    else if (action == "voltage") Voltage(config, geometry);
    else throw std::runtime_error("Commands: menu, event, event3d, save, examples, voltage, raw, raw:Branch, raw-clean, raw2d, raw2d-clean, cuts, info, check");


   
}
