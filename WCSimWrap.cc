#include "WCSimWrap.h"
#include <TMath.h>
// initialize static values
WCSimWrap *     WCSimWrap::staticthis = NULL;
TFile *         WCSimWrap::fIn = NULL;
TTree *         WCSimWrap::fTIn = NULL;
TTree *         WCSimWrap::fTGeo = NULL;

int                WCSimWrap::fNEv = 0;
int                WCSimWrap::fCurEntry = -1;
WCSimRootEvent*    WCSimWrap::fEv = NULL;
int                WCSimWrap::fNSubEv = 0;
WCSimRootTrigger * WCSimWrap::fTrig[MAXNSUBEV] = {
  NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
  NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
  NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL };
double WCSimWrap::fTOffset[MAXNSUBEV]={
  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
  0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

WCSimRootGeom * WCSimWrap::fGeo     = NULL;

Int_t           WCSimWrap::fNPMT    = 0;
Float_t **      WCSimWrap::fPMTpos  = NULL;
Float_t **      WCSimWrap::fPMTdir  = NULL;
Float_t *       WCSimWrap::fPMTQE   = NULL;
Int_t   *       WCSimWrap::fPMTtype = NULL;
Float_t *       WCSimWrap::fPMTR    = NULL;
Float_t *       WCSimWrap::fPMTphi  = NULL;
Int_t   *       WCSimWrap::fPMTnum  = NULL;
Int_t   *       WCSimWrap::fPMTloc  = NULL;

//=============================================================================
// WCSimWrap() Constructor does nothing.
// Everything happens when calling Get() method
WCSimWrap::WCSimWrap(){

  return;
}

//=============================================================================
/// ~WCSimWrap() Destructor closes file, unsets instance and deletes 
///  saved goemetry info 
WCSimWrap::~WCSimWrap(){
  staticthis = 0;
  if (fIn!=0){
    if (fIn->IsOpen()) fIn->Close();
  }
  fNEv = 0;
  fNSubEv = 0;
  if ( fPMTpos != NULL ){
    for (int ipmt=0; ipmt<fNPMT; ipmt++){
      if ( fPMTpos[ipmt] != NULL) delete [] fPMTpos[ipmt];
      if ( fPMTdir[ipmt] != NULL) delete [] fPMTdir[ipmt];
    }
    delete [] fPMTpos;
    delete [] fPMTdir;
  }
  if (fPMTQE!=NULL)   delete [] fPMTQE;
  if (fPMTtype!=NULL) delete [] fPMTtype;
  if (fPMTR!=NULL)    delete [] fPMTR;
  if (fPMTphi!=NULL)  delete [] fPMTphi;
  if (fPMTnum!=NULL)  delete [] fPMTnum;
  if (fPMTloc!=NULL)  delete [] fPMTloc;
  fNPMT = 0;
  return;
}

//=============================================================================
/// Method to open and prepare reading from a
/// WCSim input file by the name in afilename
/// and load the first Tree entry        
/// Return values: 0 = all okay
///               -1 = missing file
///               -2 = could not open file
///               -3 = missing event tree
///               -4 = missing geometry tree
///               -5 = error reading geometry
///               -6 = error reading first event
///               -7 = requested event number not found
int WCSimWrap::SetInput( char* afilename ){

  std::cout<<"<WCSimWrap::SetInput> "<<afilename<<std::endl;
  fTOffset[0] = 0.0;
  
  // setup wcsim input file
  fIn = new TFile(afilename,"read");
  if (!fIn) {
    std::cout<<"<WCSimWrap::SetInput> Unable to open wcsim input file: "<<afilename<<std::endl;
    return -1;
  }
  if (fIn->IsOpen()==kFALSE){
    std::cout<<"<WCSimWrap::SetInput> wcsim input file not open: "<<afilename<<std::endl;
    return -2;
  }
  fTIn = (TTree*) fIn->Get("wcsimT");
  if (!fTIn){
    std::cout<<"<WCSimWrap::SetInput> Unable to find wcsimT TTree in "<<afilename<<std::endl;
    return -3;
  }
  fTGeo = (TTree*) fIn->Get("wcsimGeoT");
  if (!fTGeo){
    std::cout<<"<WCSimWrap::SetInput> Unable to find wcsimGeoT TTree in "<<afilename<<std::endl;
    return -4;
  }

  // set up event storage
  fEv = new WCSimRootEvent();
  fTIn->SetBranchAddress("wcsimrootevent",&fEv);
  fTIn->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  fNEv = fTIn->GetEntries();
  fNSubEv = 0;

  // set up geometry storage
  int istat = ReadGeom();
  if (istat<0) return istat;

  // Load first entry
  return LoadEntry();
} 

//=============================================================================
/// ReadGeom() is a helper method to read in the geometry
int WCSimWrap::ReadGeom(){

  if (fGeo!=NULL){
    std::cout<<"<WCSimWrap::ReadGeom> WARNING fGeo!=NULL, Reloading it anyway, but maybe memory leak"<<std::endl;
  } else{
    fGeo = new WCSimRootGeom();
  }
  fTGeo->SetBranchAddress("wcsimrootgeom",&fGeo);
  // first "event" in branch has the geometry
  fTGeo->GetEntry(0);
  
  fNPMT = fGeo->GetWCNumPMT();

  // now set up storage 
  fPMTpos  = new Float_t * [ fNPMT ];
  fPMTdir  = new Float_t * [ fNPMT ]; 
  fPMTQE   = new Float_t   [ fNPMT ];
  fPMTtype = new Int_t     [ fNPMT ];
  fPMTR    = new Float_t   [ fNPMT ];
  fPMTphi  = new Float_t   [ fNPMT ];
  fPMTnum  = new Int_t     [ fNPMT ];
  fPMTloc  = new Int_t     [ fNPMT ];
  for (int ipmt=0; ipmt<fNPMT; ipmt++) fPMTpos[ipmt] = new Float_t [3];
  for (int ipmt=0; ipmt<fNPMT; ipmt++) fPMTdir[ipmt] = new Float_t [3];
  
  // fill the arrays of geometry info
  WCSimRootPMT pmt;
  for (int ipmt=0; ipmt<fNPMT; ipmt++){    
    pmt = fGeo->GetPMT(ipmt);
    //std::cout<<"<WCSimWrap::ReadGeom> ipmt="<<ipmt<<" tubeno="<<pmt.GetTubeNo()<<std::endl;

    fPMTpos[ipmt][0] = pmt.GetPosition(0);
    fPMTpos[ipmt][1] = pmt.GetPosition(1);
    fPMTpos[ipmt][2] = pmt.GetPosition(2);
    fPMTdir[ipmt][0] = pmt.GetOrientation(0);
    fPMTdir[ipmt][1] = pmt.GetOrientation(1);
    fPMTdir[ipmt][2] = pmt.GetOrientation(2);
    fPMTQE[ipmt] = 1.0; // currently not in wcsim output file!
    fPMTtype[ipmt] = 1; // currently not in wcsim output file!
    fPMTR[ipmt] = TMath::Sqrt( fPMTpos[ipmt][0]*fPMTpos[ipmt][0] + fPMTpos[ipmt][1]*fPMTpos[ipmt][1] );
    fPMTphi[ipmt] = TMath::ATan2(  fPMTpos[ipmt][1], fPMTpos[ipmt][0] );
    fPMTnum[ipmt] = pmt.GetTubeNo(); 
    fPMTloc[ipmt] = pmt.GetCylLoc();
  }
  return 0;
}
  
//=============================================================================
/// Method to load next event
/// Return values: as per Set Input
int WCSimWrap::LoadNextEvent(){
  fCurEntry++;
  LoadEntry( fCurEntry );
  return 0;
}

//=============================================================================
/// Method to load file to a particular TTree entry number
/// Return values: as per Set Input
int WCSimWrap::LoadEntry(int aentry){
  fCurEntry = aentry;
  if ( fTIn == NULL ) return -3;
  int anbytes = fTIn->GetEntry(fCurEntry);
  if (anbytes<=0) return -7;
  fNSubEv = std::min(fEv->GetNumberOfSubEvents(),MAXNSUBEV);
  if (fNSubEv==0){
    fNSubEv = 1; // ? not sure if this is right
    fTrig[0] = fEv->GetTrigger(0);
  } else {
    // figure out time offsets
    int aTime = 0.0;
    int aTimeZero = 0.0;
    for (int isub=0; isub<std::min(MAXNSUBEV,fNSubEv); isub++){
      fTrig[isub] = fEv->GetTrigger(isub);
      aTime = fTrig[isub]->GetHeader()->GetDate();
      std::cout<<"Subevent "<<isub<<" Date="<<aTime<<" dT="<<float(aTime-aTimeZero)<<std::endl;
      if ( isub>0 ){
	fTOffset[isub] = float(aTime-aTimeZero);
      } else {
	aTimeZero = aTime;
      }
 	
    }
  }

  std::cout<<"<WCSimWrap::LoadEntry> "<<aentry<<" with "<<fNSubEv<<" sub events"<<std::endl;
  return 0;
}

//=============================================================================
/// Method to load a particular event number
/// Return values: as per Set Input
int WCSimWrap::LoadEvent(int aevent){
  // first check if by chance the loaded entry is the one looking for
  if (aevent == fTrig[0]->GetHeader()->GetEvtNum()) return 0;

  for (fCurEntry = 0; fCurEntry<fNEv; fCurEntry++){
    fTIn->GetEntry(fCurEntry);
    fTrig[0] = fEv->GetTrigger(0);
    if ( aevent == fTrig[0]->GetHeader()->GetEvtNum() ){
      return LoadEntry( fCurEntry );
    }
  }
  return -7;
}

