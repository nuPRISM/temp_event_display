#ifndef WCSimWrap_h_
#define WCSimWrap_h_

#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"

const int MAXNSUBEV = 30; //< Maximum number of sub events to allow

/// Singleton class to provide global access
/// to the data and geometry info from WCSim input
/// author: Blair Jamieson, Jan 2013
class WCSimWrap {
 private:
  WCSimWrap();                       //< Constructor method for this class
  int ReadGeom();                    //< Helper method to read geometry

  static WCSimWrap* staticthis;  //< The instance of this class

  static TFile * fIn;            //< Input WCSim file pointer
  static TTree * fTIn;           //< Input WCSim data and truth TTree pointer
  static TTree * fTGeo;          //< Input WCSim pmt Geometry TTree pointer
 
  static int fNEv;               //< Number of events in current file
  static int fCurEntry;          //< Current tree entry number loaded
  static WCSimRootEvent *fEv;    //< Pointer to current super event

  static int fNSubEv;            //< Number of subevents in current super event
  /// Pointers to up to MAXNSUBEV subevents
  static WCSimRootTrigger *fTrig[MAXNSUBEV]; 
  static double fTOffset[MAXNSUBEV];

  static WCSimRootGeom * fGeo;   //< Pointer to detector geometry
  // save copies of geometry here for faster access
  static Int_t fNPMT;            //< Number of PMTs
  static Float_t ** fPMTpos;     //< PMT position Index by [fNPMT][x,y,z]
  static Float_t ** fPMTdir;     //< PMT normal to front face Index by [fNPMT][dx,dy,dz]
  static Float_t * fPMTQE;       //< Relative PMT QE (~1)
  static Int_t   * fPMTtype;     //< PMT Type (0 old, 1 new, -1 missing)
  static Float_t * fPMTR;        //< PMT R (calculated)
  static Float_t * fPMTphi;      //< PMT phi angle (calculated)
  static Int_t   * fPMTnum;      //< PMT tube number
  static Int_t   * fPMTloc;      //< PMT location (endcap1, wall, endcap2)

 public:
  /// Method to get the instance of this class
  /// Optionally can set the filename to try to read in
  static WCSimWrap* Get( char * afilename){
    if (!staticthis) {
      std::cout<<"<WCSimWrap::Get> Initializing (w/file)..."<<std::endl;
      staticthis = new WCSimWrap();
    }
    staticthis->SetInput( afilename );
    return staticthis;
  }
  static WCSimWrap* Get(){
    if (!staticthis) {
      std::cout<<"<WCSimWrap::Get> Initializing (no file)..."<<std::endl;
      staticthis = new WCSimWrap();
    }
    return staticthis;
  }
  
  /// Destructor closes file, unsets instance and deletes saved goemetry info 
  ~WCSimWrap(); 

  /// Method to open and prepare reading from a
  /// WCSim input file by the name in afilename
  /// and load the first Tree entry        
  /// Return values: 0 = all okay
  ///               -1 = missing file
  ///               -2 = missing event tree
  ///               -3 = missing event tree
  ///               -4 = missing geometry tree
  ///               -5 = error reading geometry
  ///               -6 = error reading first event
  ///               -7 = requested event number not found
  int SetInput( char* afilename); 
  
  /// Method to load next event
  /// Return values: as per Set Input
  int LoadNextEvent();

  /// Method to load file to a particular TTree entry number
  /// Return values: as per Set Input
  int LoadEntry(int aentry=0);

  /// Method to load a particular event number
  /// Return values: as per Set Input
  int LoadEvent(int aevent=0);

  /// Get Number of Events in File
  int NEvt(){ return fNEv; }

  /// Get Number Sub Events of current event
  int NSubEvt(){ return fNSubEv; }
  
  /// Get pointer to Current Event Object
  WCSimRootEvent * Evt(){ return fEv; }

  /// Get pointer to Sub Event
  WCSimRootTrigger * SubEvt(int asubindex){ 
    // if there are no sub events, sub event 0 is always the "event"
    if ( fNSubEv==0 && asubindex==0 ){
      //std::cout<<"<WCSimWrap::SubEvt("<<asubindex<<") ="<<fTrig[asubindex]<<std::endl;
      return fTrig[0];
    } else if (asubindex>=0 && asubindex<fNSubEv){
      //std::cout<<"<WCSimWrap::SubEvt("<<asubindex<<") ="<<fTrig[asubindex]<<std::endl;
      return fTrig[asubindex];
    } else {
      //std::cout<<"<WCSimWrap::SubEvt> out of range 0<="<<asubindex<<"<"<<fNSubEv<<std::endl;
      return NULL;
    }
  }

  double GetTOffset( int asubindex){
    if (asubindex>=0 && asubindex<fNSubEv){
      //std::cout<<"<WCSimWrap::SubEvt("<<asubindex<<") ="<<fTrig[asubindex]<<std::endl;
      return fTOffset[asubindex];
    } else {
      std::cout<<"<WCSimWrap::GetTOffset> subindex  out of range 0<="<<asubindex<<"<"<<fNSubEv<<std::endl;
      return 0.0;
    }
  }

  /// Get pointer to pointers to Sub Events
  WCSimRootTrigger ** SubEvts(){ return fTrig; }

  /// Get pointer to the Geometry
  WCSimRootGeom * Geo(){ return fGeo; }

  // Helper methods to access the Geometry
  /// Get number of PMTs in Geometry
  Int_t NPMT(){ return fNPMT; }

  /// Get the cylindrical tank size (cylinder radius)
  Float_t CylRad(){ return fGeo->GetWCCylRadius(); }
  
  /// Get cylindrical tank size (length)
  Float_t CylLen(){ return fGeo->GetWCCylLength(); }

  /// Get PMT radius
  Float_t PMTRad(){ return fGeo->GetWCPMTRadius(); }

  /// Get PMT position for the ith PMT
  Float_t PMTpos(int ipmt,int ixyz){ return fPMTpos[ipmt][ixyz]; }

  /// Get PMT positions array
  Float_t ** PMTposs(){ return fPMTpos; }

  /// Get PMT normal to front face for ith PMT
  Float_t PMTdir(int ipmt, int ixyz){ return fPMTdir[ipmt][ixyz]; }

  /// Get PMT normal to front face array
  Float_t ** PMTdirs(){ return fPMTdir; }

  /// Get PMT relative Quantum Efficiency for ith PMT
  Float_t PMTQE(int ipmt){ return fPMTQE[ipmt]; }

  /// Get PMT QE array
  Float_t * PMTQEs(){ return fPMTQE; }

  /// Get PMT type (0 old, 1 new, -1 missing) for ith PMT
  Int_t PMTtype(int ipmt){ return fPMTtype[ipmt]; }

  /// Get PMT type array
  Int_t * PMTtypes(){ return fPMTtype; }

  /// Get PMT R for ith PMT
  Float_t PMTR(int ipmt){ return fPMTR[ipmt]; }

  /// Get PMT R array
  Float_t * PMTRs(){ return fPMTR; }

  /// Get PMT phi for ith PMT
  Float_t PMTphi(int ipmt){ return fPMTphi[ipmt]; }

  /// Get PMT phi array
  Float_t * PMTphis(){ return fPMTphi; }

  /// Get PMT number for the ith PMT
  Float_t PMTnum(int ipmt){ return fPMTnum[ipmt]; }

  /// Get PMt number array
  Int_t * PMTnums(){ return fPMTnum; }

  /// Get PMT location type (endcap1, wall, endcap2) for ith PMT
  Int_t PMTloc(int ipmt){ return fPMTloc[ipmt]; }

  //// Get PMT location type array
  Int_t * PMTlocs(){ return fPMTloc; }

  

  

};


#endif
