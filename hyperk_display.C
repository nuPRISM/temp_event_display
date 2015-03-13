{
gROOT->ProcessLine(".L $WCSIMDIR/libWCSimRoot.so");
//gSystem->AddIncludePath(" -I$FITQUN_ROOT");
gSystem->AddIncludePath(" -I$WCSIMDIR/include/");
gROOT->ProcessLine(".L hyperk_esd_html_summary.C+");
gROOT->ProcessLine(".L THKGamma.C+");
gROOT->ProcessLine(".L fitQunDisplay.C+");
gROOT->ProcessLine(".L Picker.C+");
gROOT->ProcessLine(".x hyperk_esd.C");
}

