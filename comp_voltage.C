// if is_peaks_found = true, take the peak positions saved in the
// txt file and fit for gain and ped, then calibrate spectrum
// if is_peaks_found = false, fit peaks on the fly and fit for 
// gain and ped, then calibrate spectrum
bool is_peaks_found = true; 
bool overwrite_peaks_file = false; // default to be false
bool is_peaks_found_height = true; // default to be false
bool overwrite_peaks_file_height = false; // default to be false

//is_peaks_found = false;
//is_peaks_found_height = false;

void comp_voltage(const char* filelist=""){
    
    cout << " ======== is_peaks_found set to " << is_peaks_found << endl;
    cout << " ======== overwrite_peaks_file set to " << overwrite_peaks_file << endl;

    // load peak finding script
    gROOT->ProcessLine(".L peak_finding.C");


    gStyle->SetOptStat(0);

    TH1F* h_int[20]; // integral of the pulse
    TH1F* h_int_calib[20]; // calibrated
    TH1F* h_int_calib_fr[20]; // full range calibrated
    TH1F* h_height[20]; // height of the pulse
    TH1F* h_height_calib[20]; // calibrated
    TH1F* h_height_calib_fr[20]; // full range calibrated
    TH1F* h_ped_rms[20]; // pedestal rms 
    TH1F* h_ped_mean[20]; // pedestal rms 
    TH1F* h_sig_rms[20]; // signal rms 
    TH1F* h_sig_mean[20]; // signal mean 
    TH1F* h_peak_pos[20]; // position of peak
    TH1F* h_pulse_wid[20]; // width of pulse

    double gain_int[20];
    double ped_int[20];
    double gain_height[20];
    double ped_height[20];

    TLegend* lg = new TLegend(0.25, 0.50, 0.80, 0.90);
    lg->SetFillColor(0);
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);

    const char* listname="list.txt";
    if (filelist!="") listname = filelist;
    ifstream infile(listname); // file list
    TString name;
    int npt = 0;    

    TCanvas* c1 = new TCanvas();
    TCanvas* c2 = new TCanvas();
    TCanvas* c3 = new TCanvas();
    c3->Divide(3,2);
    //TCanvas* c8 = new TCanvas();
    TCanvas* c4 = new TCanvas();
    TCanvas* c5 = new TCanvas();
    TCanvas* c6 = new TCanvas();
    TCanvas* c7 = new TCanvas();

    FILE* pfile_int = NULL;
    FILE* pfile_height = NULL;

    TString last_datedir = "";
    while(true){
        string line;
        if(!infile.good()) break; // one more empty line will be loaded
        getline(infile, line);
        if(line.empty()) continue;
        
        istringstream sline(line);
        if(sline.peek()=='#'){
            //cout << "skip line " << line.data() << endl;
            continue;
        }
        sline >> name;
        cout << npt << " " << name.Data() << endl;

        int pos = name.Last('/');
        TString datedir = name(0, pos);

        TString BV = decode_bias_from_name(name);

        if(last_datedir!="" && datedir!=last_datedir){
            fclose(pfile_int);
            fclose(pfile_height);
            pfile_int = NULL;
            pfile_height = NULL;
        }

        if(pfile_int==NULL){
            pfile_int = fopen(Form("root/%s/gain_int.txt", datedir.Data()), "w");
        }
        if(pfile_height==NULL){
            pfile_height = fopen(Form("root/%s/gain_height.txt", datedir.Data()), "w");
        }

        last_datedir = datedir;

        // plot for each file
        TString filename = Form("root/%s/charge.root", 
                                name.Data());
        TFile* file = new TFile(filename, "read");
        //h[npt] = (TH1F*)file->Get("h");
    
        TTree* wf_info = (TTree*)file->Get("wf_info");

        // temp canvas to form hist
        const char* histname;
        //TCut cut = "sig_rms<1.8";
        //TCut cut = "peak_pos>238";
        //TCut cut = "peak_pos>190 && peak_pos<250";
        //TCut cut = "(min_pos-peak_pos)>10 && (min_pos-peak_pos)<35&&Entry$<100000";
        TCut cut = "";
        cout << "******* CUT used: " << cut.GetTitle() << endl;
        // plot figures
        //h_int[npt] = create_hist(wf_info,   "charge_near_peak_local", "h_int",      npt, 4200, -1, 25, "Integral (A.U.)", cut, "", -0.1, 5);
        h_int[npt] = create_hist(wf_info,   "charge", "h_int",      npt, 4200, -1, 25, "Integral (A.U.)", cut);
        
        //h_height[npt] = create_hist(wf_info, "height_local", "h_height",   npt, 4200, -0.1, 0.5, "Pulse height [V]", cut, "", 0.0, 1);
        //h_height[npt] = create_hist(wf_info, "height_local", "h_height",   npt, 4200, -0.1, 0.5, "Pulse height [V]", cut, "");
        h_height[npt] = create_hist(wf_info, "height", "h_height",   npt, 4200, -0.1, 0.5, "Pulse height [V]", cut, "");
        
        h_ped_rms[npt] = create_hist(wf_info,"ped_rms", "h_ped_rms", npt, 1000, 0, 0.02, "Baseline RMS [V]");

        h_ped_mean[npt] = create_hist(wf_info,"ped", "h_ped_mean", npt, 1000, -0.01, 0.05, "Baseline mean [V]");

        h_sig_rms[npt] = create_hist(wf_info,"sig_rms", "h_sig_rms", npt, 1000, 0, 0.02, "Signal RMS [V]");

        h_sig_mean[npt] = create_hist(wf_info,"sig_mean", "h_sig_mean", npt, 1000, -0.01, 0.05, "Signal mean [V]");

        h_peak_pos[npt] = create_hist(wf_info,"peak_pos", "h_peak_pos", npt, 1024, 0, 1024, "Peak position");

        h_pulse_wid[npt] = create_hist(wf_info,"(min_pos-peak_pos)", "h_pulse_wid", npt, 2048, -1024, 1024, "pulse width");

        double res[2];
        c1->cd();
        if(npt==0) h_int[npt]->Draw("hist");
        else h_int[npt]->Draw("histsame");
        lg->AddEntry(h_int[npt], name.Data(), "l");
        // calibration using fit with single gaus
        //fit_single_gaus(h_int[npt], res);
        //fprintf(pfile_int, "%f %f\n", res[0], res[1]);
        //gain_int[npt] = res[0];

        // calibration with peak finding using multiple peaks
        const char* peak_int_filename = Form("root/%s/peak_pos_int.txt", name.Data());
        if(!is_peaks_found){
            if(!gSystem->AccessPathName(peak_int_filename)){
                // file already existed
                if(!overwrite_peaks_file){
                    cout << peak_int_filename << " already existed, exit ..." << endl;
                    cout << "SET overwrite_peaks_file = true to ignore." << endl;
                    cout << "If you do want to find peaks again, SET overwrite_peaks_file = true to ignore." << endl;
                    cout << "Otherwise SET is_found_peaks = true" << endl;
                    return;
                }
            }
            get_peaks(h_int[npt], peak_int_filename);
        }
        // read peak pos from peak_int_filename and fit gain
        get_gain_calib(peak_int_filename, gain_int[npt], ped_int[npt]);
        fprintf(pfile_int, "%s %f %f\n", BV.Data(), gain_int[npt], ped_int[npt]);
        
        c2->cd();
        if(npt==0) h_height[npt]->Draw("hist");
        else h_height[npt]->Draw("histsame");
        cout << "Mean " << h_height[npt]->GetMean() << endl;
        //fit_single_gaus(h_height[npt], res);
        //fprintf(pfile_height, "%f %f\n", res[0], res[1]);
        //gain_height[npt] = res[0];

        // calibration with peak finding using multiple peaks
        const char* peak_height_filename = Form("root/%s/peak_pos_height.txt", name.Data());
        if(!is_peaks_found_height){
            if(!gSystem->AccessPathName(peak_height_filename)){
                // file already existed
                if(!overwrite_peaks_file_height){
                    cout << peak_height_filename << " already existed, exit ..." << endl;
                    cout << "If you do want to find peaks again, SET overwrite_peaks_file = true to ignore." << endl;
                    cout << "Otherwise SET is_found_peaks = true" << endl;
                    return;
                }
            }
            get_peaks(h_height[npt], peak_height_filename);
        }
        get_gain_calib(peak_height_filename, gain_height[npt], ped_height[npt]);
        fprintf(pfile_height, "%s %f %f\n", BV.Data(), gain_height[npt], ped_height[npt]);

        c3->cd();
        c3->cd(1);
        if(npt==0) h_ped_rms[npt]->Draw("hist");
        else h_ped_rms[npt]->Draw("histsame");

        c3->cd(2);
        if(npt==0) h_ped_mean[npt]->Draw("hist");
        else h_ped_mean[npt]->Draw("histsame");

        c3->cd(3);
        if(npt==0) h_sig_rms[npt]->Draw("hist");
        else h_sig_rms[npt]->Draw("histsame");

        c3->cd(4);
        if(npt==0) h_sig_mean[npt]->Draw("hist");
        else h_sig_mean[npt]->Draw("histsame");

        c3->cd(5);
        if(npt==0) h_peak_pos[npt]->Draw("hist");
        else h_peak_pos[npt]->Draw("histsame");

        c3->cd(6);
        if(npt==0) h_pulse_wid[npt]->Draw("hist");
        else h_pulse_wid[npt]->Draw("histsame");

        // draw calibrated spectrum
        h_int_calib[npt] = create_hist(wf_info, Form("(charge-%f)/%f", ped_int[npt], gain_int[npt]), "h_int_calib", npt, 1500, 0, 15, "Integral Charge (PE)", cut, BV.Data());
        c4->cd();
        if(npt==0) h_int_calib[npt]->Draw("hist");
        else h_int_calib[npt]->Draw("histsame");

        h_height_calib[npt] = create_hist(wf_info, Form("(height-%f)/%f", ped_height[npt], gain_height[npt]), "h_height_calib", npt, 300, 0, 15, "Height Charge (PE)", cut, BV.Data());
        //h_height_calib[npt] = create_hist(wf_info, Form("(height)/0.02", ped_height[npt], gain_height[npt]), "h_height_calib", npt, 300, 0, 15, "Height Charge (PE)", cut, BV.Data());
        c5->cd();
        if(npt==0) h_height_calib[npt]->Draw("hist");
        else h_height_calib[npt]->Draw("histsame");

        h_int_calib_fr[npt] = create_hist(wf_info, Form("(charge-%f)/%f", ped_int[npt], gain_int[npt]), "h_int_calib_fr", npt, 640, 0, 160, "Integral Charge (PE)", cut, BV.Data());
        c6->cd();
        if(npt==0) h_int_calib_fr[npt]->Draw("hist");
        else h_int_calib_fr[npt]->Draw("histsame");

        h_height_calib_fr[npt] = create_hist(wf_info, Form("(height-%f)/%f", ped_height[npt], gain_height[npt]), "h_height_calib_fr", npt, 320, 0, 160, "Height Charge (PE)", cut, BV.Data());
        //h_height_calib_fr[npt] = create_hist(wf_info, Form("(height)/0.02", ped_height[npt], gain_height[npt]), "h_height_calib_fr", npt, 320, 0, 160, "Height Charge (PE)", cut, BV.Data());
        c7->cd();
        if(npt==0) h_height_calib_fr[npt]->Draw("hist");
        else h_height_calib_fr[npt]->Draw("histsame");


        npt++; // total number of points

    }


    c1->cd();
    lg->Draw();
    gPad->SetLogy();

    c2->cd();
    lg->Draw();
    gPad->SetLogy();

    for(int i=0;i<4;i++){
        c3->cd(i+1);
        gPad->SetLogy();
        lg->Draw();
    }

    //c8->cd();
    //gPad->SetLogy();
    //lg->Draw();

    c4->cd();
    gPad->SetLogy();
    gPad->SetGridx();
    lg->Draw();
    
    c5->cd();
    gPad->SetLogy();
    gPad->SetGridx();
    lg->Draw();

    c6->cd();
    gPad->SetLogy();
    lg->Draw();

    c7->cd();
    gPad->SetLogy();
    lg->Draw();

    c1->Print(Form("root/%s/integral.pdf", datedir.Data()));
    c2->Print(Form("root/%s/height.pdf", datedir.Data()));
    c3->Print(Form("root/%s/ped_sig_mean_rms.pdf", datedir.Data()));
    //c8->Print(Form("root/%s/ped_mean.pdf", datedir.Data()));
    c4->Print(Form("root/%s/calib_integral.pdf", datedir.Data()));
    c5->Print(Form("root/%s/calib_height.pdf", datedir.Data()));
    c6->Print(Form("root/%s/calib_integral_full.pdf", datedir.Data()));
    c7->Print(Form("root/%s/calib_height_full.pdf", datedir.Data()));

    fclose(pfile_int);
    fclose(pfile_height);

    // save output hist
    TFile* outfile = new TFile(Form("root/%s/spectrum.root", datedir.Data()), "recreate");
    for(int i=0;i<npt;i++){
        //h_int[i]->Write();
        //h_height[i]->Write();
        h_int_calib[i]->Write();
        h_int_calib_fr[i]->Write();
        h_height_calib[i]->Write();
        h_height_calib_fr[i]->Write();
    }
    outfile->Close();

}

TString decode_bias_from_name(TString filename){
    
    TObjArray* objArray = filename.Tokenize("_");
    int n = objArray->GetEntries();
    for(int i=0;i<n;i++){
        TObjString* objstr = (TObjString*)objArray->At(i);
        TString str = objstr->GetString();
        if(str.EndsWith("V")){
            int pos = str.Last('V');
            TString vol = str(0, pos);
            cout << "Bias voltage: " << vol << endl;
            return vol;
        }
    }


}

TH1F* create_hist(TTree* tree,
                  const char* branchname,
                  const char* histname,
                  int npt,
                  int nbins,
                  double x_low,
                  double x_high,
                  const char* xtitle,
                  const char* cut="",
                  const char* title="",
                  double norm_range=-1e9,
                  double norm_range_end=-1e9){

    gROOT->SetBatch(true);
    TCanvas* c = new TCanvas();
    const char* histname_full = Form("%s_%d", histname, npt);
    TH1F* h = new TH1F(histname_full, cut, nbins, x_low, x_high);
    tree->Draw(Form("%s>>%s", branchname, histname_full), cut, "");
    h->SetTitle(title);
    format_hist(h, npt, xtitle);
    if(norm_range>-1e9){
        int start_bin = h->FindBin(norm_range);
        int end_bin = nbins+1;
        if(norm_range_end>-1e9) end_bin = h->FindBin(norm_range_end);
        cout << start_bin << " " << end_bin << " " << h->Integral(start_bin, end_bin) << endl;
        //cout << "Scale " << 1./h->Integral(start_bin, end_bin) << endl;
        h->Scale(1./h->Integral(start_bin, end_bin));
    }
    c->Close();
    gROOT->SetBatch(false);
    return h;

}

void format_hist(TH1F* h, int npt, const char* xtitle){

    h->Sumw2();
    //h->Scale(1./h->GetEntries());
    int color = npt+1;
    if (color>=5) color+=1; // skip yellow
    h->SetLineColor(color);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitle(xtitle);

}

void fit_single_gaus(TH1F* h, 
                     double* res = NULL /* fit results to be returned */
                    ){

    TString name = h->GetName();
    double xmin = 0;
    if(name.Contains("h_int")) xmin=0.2;
    if(name.Contains("h_height")) xmin=0.02;
    

    h->GetXaxis()->SetRangeUser(xmin, h->GetXaxis()->GetXmax());
    int bin = h->GetMaximumBin();
    double x = h->GetBinCenter(bin);

    double entries = h->GetBinContent(bin);
    int low_bin=1; 
    int high_bin = h->GetNbinsX();
    // get fit range lower limit
    for(int i=bin-1;i>low_bin;i--){
        double v = h->GetBinContent(i);
        if(v<entries/4){
            low_bin = i;
            break;
        }
    }
    // get fit range upper limit
    for(int i=bin+1;i<high_bin;i++){
        double v = h->GetBinContent(i);
        if(v<entries/4){
            high_bin = i;
            break;
        }
    }

    double fit_low= h->GetBinLowEdge(low_bin);
    double fit_up = h->GetBinLowEdge(high_bin+1);

    cout << "Fit window: [" << fit_low << ", " << fit_up << "]" << endl;
    
    h->Fit("gaus", "", "", fit_low, fit_up);
    TF1* gaus = h->GetListOfFunctions()->FindObject("gaus");
    gaus->Draw("same");

    if(res!=NULL){
        res[0] = gaus->GetParameter(1);
        res[1] = gaus->GetParError(1);
        //res[1] = gaus->GetParameter(2); // save sigma
    }

    return;
}

//int get_peak_loc(TH1F* h){
//
//}

void fit_double_gaus(TH1F* h){
    
    double fit_low = 2.5;
    double fit_mid = 5.8;
    double fit_up = 6.5;

    TF1* func = new TF1("func", "[0]*exp(-pow((x-[1])/[2],2)/2)+[3]*exp(-pow((x-[4])/[5],2)/2)", fit_low, fit_up);
    //TF1* func = new TF1("func", "[0]*exp(-x/[1])+[2]*exp(-pow((x-[3])/[4],2)/2)", fit_low, fit_up);
    TF1* gaus1 = new TF1("gaus1", "[0]*exp(-pow((x-[1])/[2],2)/2)", fit_low, fit_mid);
    TF1* gaus2 = new TF1("gaus2", "[0]*exp(-pow((x-[1])/[2],2)/2)", fit_mid, fit_up);
    double best_guess[6];
    gaus1->SetParameter(0, 1);
    gaus1->SetParameter(1, 2.5);
    gaus1->SetParameter(2, 0.2);

    gaus2->SetParameter(0, 0.5);
    gaus2->SetParameter(1, 6.0);
    gaus2->SetParameter(2, 0.3);

    h->Fit(gaus1, "R");
    best_guess[0] = gaus1->GetParameter(0);
    best_guess[1] = gaus1->GetParameter(1);
    best_guess[2] = gaus1->GetParameter(2);
    h->Fit(gaus2, "R");
    best_guess[3] = gaus2->GetParameter(0);
    best_guess[4] = gaus2->GetParameter(1);
    best_guess[5] = gaus2->GetParameter(2);
    
    func->SetParameters(best_guess);
    h->Fit(func);
    func->Draw("same");
    gaus1->Draw("same");
    gaus2->Draw("same");
    
    double* p = func->GetParameters();
    cout << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " " << p[4] << " " << p[5] << endl;

}
