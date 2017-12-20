void comp_voltage(bool calib=false){

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

    double gain_int[20];
    double gain_height[20];

    TLegend* lg = new TLegend(0.25, 0.50, 0.80, 0.90);
    lg->SetFillColor(0);
    lg->SetFillStyle(0);
    lg->SetBorderSize(0);

    ifstream infile("list.txt");
    TString name;
    int npt = 0;    

    TCanvas* c1 = new TCanvas();
    TCanvas* c2 = new TCanvas();
    TCanvas* c3 = new TCanvas();
    c3->Divide(2,2);
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
        //TCut cut = "ped_mean<0.01";
        TCut cut = "";
        // plot figures
        h_int[npt] = create_hist(wf_info,   "charge", "h_int",      npt, 440, -1, 10, "Integral (A.U.)", cut);
        
        h_height[npt] = create_hist(wf_info,"height", "h_height",   npt, 4440, -0.1, 10, "Pulse height [V]", cut);
        
        h_ped_rms[npt] = create_hist(wf_info,"ped_rms", "h_ped_rms", npt, 1000, 0, 0.02, "Baseline RMS [V]");

        h_ped_mean[npt] = create_hist(wf_info,"ped_mean", "h_ped_mean", npt, 1000, -0.01, 0.05, "Baseline mean [V]");

        h_sig_rms[npt] = create_hist(wf_info,"sig_rms", "h_sig_rms", npt, 1000, 0, 0.02, "Signal RMS [V]");

        h_sig_mean[npt] = create_hist(wf_info,"sig_mean", "h_sig_mean", npt, 1000, -0.01, 0.05, "Signal mean [V]");

        double res[2];
        c1->cd();
        if(npt==0) h_int[npt]->Draw("hist");
        else h_int[npt]->Draw("histsame");
        lg->AddEntry(h_int[npt], name.Data(), "l");
        fit_single_gaus(h_int[npt], res);
        fprintf(pfile_int, "%f %f\n", res[0], res[1]);
        gain_int[npt] = res[0];
        
        c2->cd();
        if(npt==0) h_height[npt]->Draw("hist");
        else h_height[npt]->Draw("histsame");
        fit_single_gaus(h_height[npt], res);
        //fprintf(pfile_height, "%f %f\n", res[0], res[1]);
        //gain_height[npt] = res[0];

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

        // draw calibrated spectrum
        h_int_calib[npt] = create_hist(wf_info, Form("charge/%f", gain_int[npt]), "h_int_calib", npt, 1500, 0, 15, "Integral Charge (PE)", cut);
        c4->cd();
        if(npt==0) h_int_calib[npt]->Draw("hist");
        else h_int_calib[npt]->Draw("histsame");

        h_height_calib[npt] = create_hist(wf_info, Form("height/%f", gain_height[npt]), "h_height_calib", npt, 300, 0, 15, "Height Charge (PE)", cut);
        c5->cd();
        if(npt==0) h_height_calib[npt]->Draw("hist");
        else h_height_calib[npt]->Draw("histsame");

        h_int_calib_fr[npt] = create_hist(wf_info, Form("charge/%f", gain_int[npt]), "h_int_calib_fr", npt, 640, 0, 160, "Integral Charge (PE)", cut);
        c6->cd();
        if(npt==0) h_int_calib_fr[npt]->Draw("hist");
        else h_int_calib_fr[npt]->Draw("histsame");

        h_height_calib_fr[npt] = create_hist(wf_info, Form("height/%f", gain_height[npt]), "h_height_calib_fr", npt, 320, 0, 160, "Height Charge (PE)", cut);
        c7->cd();
        if(npt==0) h_height_calib_fr[npt]->Draw("hist");
        else h_height_calib_fr[npt]->Draw("histsame");

        npt++; // total number of points

    }


    c1->cd();
    lg->Draw();

    c2->cd();
    lg->Draw();

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
        h_int_calib[i]->Write();
        h_int_calib_fr[i]->Write();
    }
    outfile->Close();

}

TH1F* create_hist(TTree* tree,
                  const char* branchname,
                  const char* histname,
                  int npt,
                  int nbins,
                  double x_low,
                  double x_high,
                  const char* xtitle,
                  const char* cut=""){

    gROOT->SetBatch(true);
    TCanvas* c = new TCanvas();
    const char* histname_full = Form("%s_%d", histname, npt);
    TH1F* h = new TH1F(histname_full, cut, nbins, x_low, x_high);
    tree->Draw(Form("%s>>%s", branchname, histname_full), cut);
    format_hist(h, npt, xtitle);
    c->Close();
    gROOT->SetBatch(false);
    return h;

}

void format_hist(TH1F* h, int npt, const char* xtitle){

    h->Sumw2();
    //h->Scale(1./h->GetEntries());
    h->SetLineColor(npt+1);
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
