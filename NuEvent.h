////////// ///////// ///////// //////// ///////// ///////// ///////// 72
////////////////////////////////////////////////////////////////////////
// Contact: Jeff Hartnell j.j.hartnell@rl.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef NUEVENT_H
#define NUEVENT_H

#include <cassert>
#include <string>
#include <vector>

#include "TObject.h"
namespace minos{
class NuEvent : public TObject
{
 public:
  NuEvent();
  virtual ~NuEvent() {};

  //non-const methods
  void Reset();

  //////////////////////////////////////////////////////////////////
  /// EVERYTIME A TRUTH VARIABLE IS ADDED TO THIS CLASS IT MUST
  /// ALSO BE ADDED TO NuMCEvent
  ///////////////////////////////////////////////////////////////////

  /// \name Book-keeping quantities
  //@{
  Int_t index; ///< the n'th NuEvent object to be written to the tree
  Int_t entry; ///< the n'th snarl to be processed
  //@}

  /// \name Snarl/run based quantities
  //@{
  Int_t run; ///< fHeader.fRun
  Int_t subRun; ///< fHeader.fSubRun
  Int_t runType; ///< fHeader.fRunType
  Int_t errorCode; ///< fHeader.fErrorCode
  Int_t snarl; ///< actual DAQ snarl number
  Int_t trigSrc; ///< fHeader.fTrigSrc
  Int_t timeFrame; ///< fHeader.fTimeFrame
  Int_t remoteSpillType; ///< fHeader.fRemoteSpillType

  Int_t detector; ///< fHeader.fVldContext.fDetector
  Int_t simFlag; ///< fHeader.fVldContext.fSimFlag
  Int_t timeSec; ///< fHeader.fVldContext.fTimeStamp.fSec
  Int_t timeNanoSec; ///< fHeader.fVldContext.fTimeStamp.fNanoSec
  Double_t timeSeconds; ///< VldTimeStamp.GetSeconds() sec+nanoSec/1e9

  Double_t trigtime; ///< evthdr.trigtime
  Double_t medianTime; ///< the median time in event (with respect to timeFrame start)
  Double_t timeEvtMin; ///< earliest time in event (with respect to
		       ///timeNanoSec). This comes from
		       ///NtpSRStrip::time1 (or an average of time1
		       ///and time0 in the FD), which is a
		       ///charge-weighted average of the times of the
		       ///digits in the earliest strip.
  Double_t timeEvtMax; ///< latest time in event (with respect to timeNanoSec)
  Double_t crateT0; ///< Crate t0
  Int_t sgate_10mhz;
  Int_t sgate_53mhz;
  Int_t rollover_53mhz;
  Int_t rollover_last_53mhz;
  Int_t timingFiducial;

  Int_t nearestSpillSec; ///< time of nearest spill (s)
  Int_t nearestSpillNanosec; ///< time of nearest spill (ns);
  Double_t timeToNearestSpill; ///< time to nearest spill (uses SpillTimeFinder)

  Int_t planeEvtHdrBeg; ///< evthdr.plane.beg (different than evt.plane.beg!)
  Int_t planeEvtHdrEnd; ///< evthdr.plane.end (different than evt.plane.end!)
  Double_t snarlPulseHeight; ///< evthdr.ph.sigcor
  //@}

  /// \name Data quality variables
  //@{
  Bool_t isGoodDataQuality; ///< uses DataUtil::IsGoodData()
  Bool_t isGoodDataQualityRUN; ///< uses DataUtil::IsGoodDataRUN()
  Bool_t isGoodDataQualityCOIL; ///< uses DataUtil::IsGoodDataCOIL()
  Bool_t isGoodDataQualityHV; ///< uses DataUtil::IsGoodDataHV()
  Bool_t isGoodDataQualityGPS; ///< uses DataUtil::IsGoodDataGPS()

  Int_t numActiveCrates; ///<  number of active readout crates ('cratemask')
  Int_t numTimeFrames; ///<  number of timeframes, gives run duration
  Int_t numGoodSnarls; ///<  number of 'good' snarls during run
  Float_t snarlRateMedian; ///<  median snarl rate during run
  Float_t snarlRateMax; ///<  maximum snarl rate during run

  Int_t deltaSecToSpillGPS; ///< time from trigger to nearest spill (s)
  Int_t deltaNanoSecToSpillGPS; ///< time from trigger to nearest spill (ns)
  Int_t gpsError; ///< GPS error (ns)
  Int_t gpsSpillType; ///<  type of spill

  Bool_t coilIsOk; ///< CoilTools::IsOk(vldc);
  Bool_t coilIsReverse; ///< CoilTools::IsReverse(vldc);
  Float_t coilCurrent; ///< detector B field coil current (not always filled)

  Bool_t isLI; ///< LISieve::IsLI()
  Int_t litag; ///< CC analysis definition: reco.FDRCBoundary(nu);
  Double_t litime; ///< evthdr.litime (==-1 if not LI)
  //@}

  /// \name Reconstructed quantities
  /// These variables are here to add a level of indirection, they don't store
  /// information directly from the sntp files, but rather they can be set at
  /// analysis-time to what the user wants them to contain
  /// %e.g. \ref shwEn could be set to the deweighted, energy corrected,
  /// CC shower energy if the analysis selected the event as a CC event.
  /// Alternatively, if the analysis
  /// selected the event as %NC then \ref shwEn could be set to the linear %NC,
  /// no energy correction, %NC shower energy.
  /// Don't rely on the defaults, but they typically assume you
  /// have a CC event and want linear, energy corrected values of the energy

  //@{

  /// \brief Best reconstructed energy.
  ///
  /// Always set to \ref energyCC by \ref NuReco::GetEvtEnergy. Unless
  /// \ref applyEnergyShifts is set in which case \ref NuReco::ApplyReweights
  /// sets it to \ref energyRw.
  Float_t energy;

  /// \brief Reconstructed energy assuming this is a CC event
  ///
  /// Always set to the sum of \ref trkEn and \ref shwEnCC by
  /// \ref NuReco::GetEvtEnergy.
  Float_t energyCC;

  /// \brief Reconstructed energy assuming this is an %NC event
  ///
  /// Always set to \ref shwEnNC by NuReco::GetEvtEnergy
  Float_t energyNC;

  /// \brief Reconstructed energy for RAF events
  ///
  /// Always set to \ref trkEn by NuReco::GetEvtEnergy
  Float_t energyRM;

  /// \brief Best reconstructed track energy
  ///
  /// Set in \ref NuReco::CalcTrkVariables. Zero if not \ref trkExists.
  /// Set to \ref trkEnRange if \ref usedRange, or \ref trkEnCurv if
  /// \ref usedCurv. These should never both be false.
  /// If \ref applyEnergyShifts then \ref NuReco::ApplyReweights sets this to
  /// \ref trkEnRw.
  Float_t trkEn;

  /// \brief Reconstructed track energy calculated from range
  ///
  /// Set by \ref NuReco::GetTrackEnergyFromRange. This applies energy
  /// corrections to \ref trkMomentumRange and then sums in the muon mass.
  Float_t trkEnRange;

  /// \brief Reconstructed track energy calculated from curvature
  ///
  /// Set by \ref NuReco::GetTrackEnergyFromCurv. This applies energy
  /// corrections to 1/\ref qp and then sums in the muon mass.
  Float_t trkEnCurv;

  /// \brief Best reconstructed shower energy
  ///
  /// Always set to \ref shwEnCC by \ref NuReco::CalcShwVariables.
  /// If \ref applyEnergyShifts then \ref NuReco::ApplyReweights sets this to
  /// \ref shwEnRw.
  Float_t shwEn;

  /// \brief Reconstructed shower energy assuming this is a CC event
  ///
  /// Calculated in \ref NuReco::GetShowerEnergyCC. Returns \ref shwEnkNN
  /// unless \ref hornIsReverse in which case it applies energy corrections
  /// to \ref shwEnLinCCNoCor
  Float_t shwEnCC;

  /// \brief Reconstructed shower energy assuming this is an %NC event
  ///
  /// Calculated in \ref NuReco::GetShowerEnergyNC. Applies energy corrections
  /// to \ref shwEnLinNCNoCor.
  Float_t shwEnNC;

  /// \brief Estimated energy resolution of this event in GeV
  ///
  /// Set by \ref NuReco::CalcResolution to the quadrature sum of
  /// \ref trkResolution and \ref shwResolution.
  Float_t resolution;
  //@}

  /// \name Convenience flags
  /// These 3 flags are just for convenience when playing with DSTs
  /// they are ABSOLUTELY NOT TO BE USED FOR ANALYSIS
  //@{
  Bool_t isInFidVolCC; ///< Is the event in the CC 2010 fiducial volume
  Bool_t isGoodTrackReclamation; ///< Pass trkFitPass or ND reclamation
  Bool_t isNCClean; ///< Pass %NC cleaning (for all events, not just NCs)
  //@}

  /// \name Basic kinematics
  //@{
  Float_t y; ///< reco'd y = (Enu-Emu)/Enu
  Float_t q2; ///< reco'd q2 = -q^{2} = 2 Enu Emu (1 - cos(theta_mu))
  Float_t x; ///< reco'd x = Q2 / (2M* Ehad)  {M=nucleon mass}
  Float_t w2; ///< reco'd w2 = M^{2} + q^{2} + 2M(Enu-Emu)
  Float_t dirCosNu; ///< reco'd direction cosine of track
  //@}

  /// \name Event info extracted from the ntuples
  //@{
  Int_t evt; ///< reco'd event number, evt.index
  Int_t slc; ///< reco'd slice number, evt.slc
  Int_t nevt; ///<  Number of reconstructed events in the snarl
  Int_t ndigitEvt; ///< evt.ndigit
  Int_t nstripEvt; ///< evt.nstrip
  Int_t nshw; ///< evt.nshower
  Int_t ntrk; ///< evt.ntrack
  Int_t primshw; ///< evt.primshw
  Int_t primtrk; ///< evt.primtrk
  Float_t rawPhEvt; ///< evt.ph.raw
  Float_t evtphsigcor; ///< evt.ph.sigcor
  Float_t evtphsigmap; ///< evt.ph.sigmap
  Int_t planeEvtN; ///< number of planes in event
  Int_t planeEvtNu; ///< evt.plane.nu
  Int_t planeEvtNv; ///< evt.plane.nv
  //@}

  //Information from Andy Blake's track vertex timing fix. This is
  //event-based information
  Bool_t vtxFitPassFail; //Does this event pass the fit? (This
			 //basically requires there to be a track.)
  Double_t vtxFitTime; //The best fit vertex time.
  Double_t vtxFitTimeError; //The error in the reconstructed vertex time.
  Double_t vtxFitChi2DoF; //The chi2/DoF from the fit.


  /// \name RO PID variables
  /// These RO PID variables store the values obtained using the
  /// \ref NtpSREvent. Rustem's code defines which track to use by using
  /// the number of "active" planes
  //@{
  Float_t roIDEvt; ///< RO's PID variable (got using the evt)
  Float_t knn01TrkActivePlanesEvt; ///< number of active planes in trk
  Float_t knn10TrkMeanPhEvt; ///< average ph per plane in trk
  Float_t knn20LowHighPhEvt; ///< average of low ph strips / average of high ph strips
  Float_t knn40TrkPhFracEvt; ///< fraction of ph in trk
  Float_t roIDNuMuBarEvt; ///< RO's PID NuMuBar selection (0 or 1)
  Float_t relativeAngleEvt; ///< RO's track angle relative to muon dir.
  //@}

  /// \name  Jasmine Ma/Ratchford new PID variables (Event)
  //@{
  Float_t jmIDEvt;
  Float_t jmTrackPlaneEvt;
  Float_t jmEndPhEvt;
  Float_t jmMeanPhEvt;
  Float_t jmScatteringUEvt;
  Float_t jmScatteringVEvt;
  Float_t jmScatteringUVEvt;
  Float_t jmEventknnIDEvt;
  Float_t jmEventknn208Evt; ///< Summed ph in path end
  Float_t jmEventknn207Evt; ///< Standard devation of ph in path end
  Float_t jmEventknn206Evt; ///< Mean ph in path
  Float_t jmEventknn205Evt; ///< "BVD" the value of the Dijkstra algorithm
                            ///< (proportional to the probability that it is a track)
  Float_t jmEventknn204Evt; ///< Size of the path found by the Dijkstra algorithm
  //@}

  /// \name Event vertex and end
  //@{
  Float_t xEvtVtx; ///< evt.vtx.x
  Float_t yEvtVtx; ///< evt.vtx.y
  Float_t zEvtVtx; ///< evt.vtx.z
  Float_t uEvtVtx; ///< evt.vtx.u
  Float_t vEvtVtx; ///< evt.vtx.v
  Float_t tEvtVtx; ///< evt.vtx.t
  Int_t planeEvtVtx; ///< evt.vtx.plane
  Int_t planeEvtBeg; ///< evt.plane.beg
  Int_t planeEvtBegu; ///< evt.plane.begu
  Int_t planeEvtBegv; ///< evt.plane.begv

  Float_t xEvtEnd; ///< evt.end.x
  Float_t yEvtEnd; ///< evt.end.y
  Float_t zEvtEnd; ///< evt.end.z
  Float_t uEvtEnd; ///< evt.end.u
  Float_t vEvtEnd; ///< evt.end.v
  Int_t planeEvtEnd; ///< evt.plane.end
  Int_t planeEvtEndu; ///< evt.plane.endu
  Int_t planeEvtEndv; ///< evt.plane.endv
  //@}

  /// \name Best track
  //@{
  Bool_t trkExists; ///< simply state if track exists
  Int_t trkIndex; ///< trk.index, position in TClonesArray in sntp file
  Int_t ndigitTrk; ///< trk.ndigit
  Int_t nstripTrk; ///< trk.nstrip

  /// \brief Equivalent to \ref trkEnRange
  Float_t trkEnCorRange;

  /// \brief Equivalent to \ref trkEnCurv
  Float_t trkEnCorCurv;

  /// \brief Reconstructed shower energy within 1m of the track vertex
  ///
  /// Calculated by \ref NuReco::GetShowerEnergyNearTrack. This is the sum of
  /// \ref NtpSRShowerPulseHeight::linCCgev for all reconstructed showers with
  /// vertex within 1m of the track vertex (and times within -25ns,+75ns).
  Float_t trkShwEnNear;

  /// \brief Reconstructed deweighted shower energy within 1m of the track vertex
  ///
  /// As for \ref trkShwEnNear but using \ref NtpSRShowerPulseHeight::wtCCgev.
  /// -1 if there are no showers within 1m.
  Float_t trkShwEnNearDW;

  /// \brief Uncorrected track momentum from range
  ///
  /// Copied from the value in \ref NtpSRMomentum::range.
  Float_t trkMomentumRange;

  /// \brief Estimated resolution on the track energy in GeV
  ///
  /// Calculated in \ref EnergyResolution::MuonResolution based on \ref trkEn,
  /// \ref sigqp and \ref usedRange.
  Float_t trkResolution;

  Int_t containedTrk; ///< trk.contained
  Int_t trkfitpass; ///< trk.fit.pass
  Float_t trkvtxdcosz; ///< trk.vtx.dcosz
  Float_t trkvtxdcosy; ///< trk.vtx.dcosy
  Int_t trknplane; ///< trk.nplane
  Int_t charge; ///< +1 or -1: simply derived from qp (from track fit)
  Float_t qp; ///< track Q/P from fit (no EnCor)
  Float_t qp_rangebiased; ///< track Q/P from fit (with range bias)
  Float_t sigqp; ///< track sigma Q/P from fit
  Float_t qp_sigqp; ///< qp / sigqp
  Float_t chi2; ///< track chi2 of fit
  Float_t ndof; ///< track ndof of fit
  Float_t qpFraction; ///< trk.stpfitqp[i], npositive/nstrip
  Int_t trkVtxUVDiffPl; ///< trk.plane.begu-trk.plane.begv;
  Int_t trkLength; ///< abs(trk.plane.end-trk.plane.beg+1);
  Int_t planeTrkNu; ///< trk.plane.nu: number of u planes hit
  Int_t planeTrkNv; ///< trk.plane.nv: number of v planes hit
  Int_t ntrklike; ///< the number of trk-like planes
  Float_t trkphsigcor; ///< trk.ph.sigcor
  Float_t trkphsigmap; ///< trk.ph.sigmap
  Int_t trkIdMC; ///< true identity of the track
  Float_t trkds; ///< total track path length from vertex to end (m) trk.ds

  Int_t trkfitpassSA; ///< variables from the SA track fitter
  Float_t trkvtxdcoszSA; ///< fitsa.fit.dcosz
  Int_t chargeSA; ///< definition same as variable without SA postfix
  Float_t qpSA; ///< definition same as variable without SA postfix
  Float_t sigqpSA; ///< definition same as variable without SA postfix
  Float_t chi2SA; ///< definition same as variable without SA postfix
  Float_t ndofSA; ///< definition same as variable without SA postfix
  Float_t probSA; ///< definition same as variable without SA postfix
  Float_t xTrkVtxSA; ///< calculated from u&v
  Float_t yTrkVtxSA; ///< calculated from u&v
  Float_t zTrkVtxSA; ///< fitsa.fit.z
  Float_t uTrkVtxSA; ///< fitsa.fit.u
  Float_t vTrkVtxSA; ///< fitsa.fit.v

  Float_t jitter; ///< from MajCInfo
  Float_t jPID; ///< from MajCInfo
  Float_t majC; ///< from MajCInfo
  Float_t smoothMajC; ///< from MajCInfo

  //The best RO PID (roID) variables are located with the other PIDs
  //see below

  Float_t xTrkVtx; ///< trk.vtx.x
  Float_t yTrkVtx; ///< trk.vtx.y
  Float_t zTrkVtx; ///< trk.vtx.z
  Float_t uTrkVtx; ///< trk.vtx.u
  Float_t vTrkVtx; ///< trk.vtx.v
  Float_t tTrkVtx;
  Int_t planeTrkVtx; ///< trk.vtx.plane
  Int_t planeTrkBeg; ///< trk.plane.beg
  Int_t planeTrkBegu; ///< trk.plane.begu
  Int_t planeTrkBegv; ///< trk.plane.begv

  //Some track timing information from the NtpSRTrackTime object
  Float_t inverseBetaTrk; ///< trk.time.cdtds, which is 1/beta
  Double_t t0Trk; ///< trk.time.t0: the offset from time fit used to
		  ///determine beta
  Float_t chi2TimeTrk; ///< trk.time.chi2: chi2 of fit done on track direction
		       ///determination
  UShort_t ndigitTimeTrk; ///< trk.time.ndigit: number of digits used
			  ///to determine the track direction from
			  ///timing
  Double_t forwardRMSTrk; ///< trk.time.forwardRMS: RMS for forward
			  ///time fit
  Int_t forwardNDOFTrk; ///< trk.time.forwardNDOF: Digits used for
			///forward time fit

  /// \name Track strip geometry information
  //@{
  Int_t stripTrkBeg; ///< Strip number of the first strip hit
  Int_t stripTrkBegu; ///< Strip number of the first strip hit in a u plane
  Int_t stripTrkBegv; ///< Strip number of the first strip hit in a v plane
  Short_t stripTrkEnd; ///< Strip number of the last strip hit
  Short_t stripTrkEndu; ///< Strip number of the last strip hit in a u plane
  Short_t stripTrkEndv; ///< Strip number of the last strip hit in a u plane
  Bool_t stripTrkBegIsu; ///< True if the first strip hit in this track was in a u plane
  Bool_t stripTrkEndIsu; ///< True if the last strip hit in this track was in a u plane
  //@}

  /// \name RAF (rock and anti-fiducial) variables
  //@{

  /*!
    This is set to -2 if not at the far detector and set to -1 if there 
    is no track.  It is 0 for tracks with their vertex in the fiducial 
    volume.  For tracks with their vertex outside the fiducial volume, 
    it is a number from 1-18 found by the sum of:

    1 if the vertex is near the outside radial edge <br />
    2 if the vertex is near the coil hole

    4 if the vertex is at the south end of supermodule 1 <br />
    8 if the vertex is at the north end of supermodule 1 <br />
    12 if the vertex is at the south end of supermodule 2 <br />
    16 if the vertex is at the north end of supermodule 2 <br />

  */
  Int_t regionTrkVtx; ///< The region of the detector the track vertex is in

  /*!
    This is calculated only for the far detector.  If not in the far 
    detector or there is no track, it's set to -1.

    Octants boundaries are at the corners.  Counting starts with 0 at 
    +x (west vertical edge) and goes counter-clockwise (looking south) 
    around to to 7 (lower west diagonal edge).  While this variable was
    intended for use only in the anti-fiducial volume near the edge, 
    it is defined for all positions.

    Note that the detector is not a regular octagon and so the corners
    are not at exactly (n+0.5)*360/8 degrees.  Because the scintillator
    boundary has a complex shape, the exact locations of the corners
    is arbitrary on the level of a few centimeters.

    In the 2010 CC analysis, only the parity of this value is used.
    i.e. if it is even, this is a horizontal/vertical edge and if it is
    odd, it is a diagonal edge.
  */
  Char_t edgeRegionTrkVtx; ///< The octant of the detector that the track vertex is in

  Char_t edgeRegionTrkEnd; ///< Same as edgeRegionTrkVtx but for the track end

  /*!
    <strong>I don't think this is directly used in any analysis and is 
    trivially calcuable with TMath::ATan2(yTrkVtx, xTrkVtx).  Remove?</strong>
  */
  Float_t phiTrkVtx; ///< Location of the track vertex in Phi = atan2(vtx.x / vtx.y)

  /*!
    <strong>I don't think this is directly used in any analysis and is 
    trivially calcuable with TMath::ATan2(yTrkEnd, xTrkEnd).  Remove?</strong>
  */
  Float_t phiTrkEnd; ///< Location of the track end in Phi = atan2(vtx.x / vtx.y)

  /*!
    If there is no track, this is set to -999.  If this is not the far
    detector, or the track vertex is not in one of the octants of the 
    far detector defined by a diagonal edge, this is set to -10.

    While this is almost certainly only interesting for anti-fiducial
    events, it is set even if the track vertex is far from the edge so
    that it doesn't matter if the fiducial cut is moved.

    This gives the number of strips from the edge, not the raw strip
    number, so on any edge if the outermost strip is the first hit,
    this is set to 0.

    In the 2010 CC analysis, a track is considered rock-like if this
    is 0 or 1.
  */   
  Short_t parallelStripTrkVtx; ///< On a diagonal edge, the strip number of the first track hit

  /*!
    To account for alignment differences between the MC and reality, 
    or to apply an alignment systematic, parallelStripTrkVtx may be 
    modified from what the raw MC would give.

    As of the 2010 analysis, this is not done since it is unclear how 
    to handle it.  Instead an analysis strategy was used that doesn't 
    depend strongly on the precise alignment.  So this should always 
    be equal to parallelStripTrkVtx at the moment. <strong>It seems 
    unlikely that this will ever change and so probably this variable
    should be removed.</strong>
  */
  Short_t parallelStripTrkVtxNoShift; ///< parallelStripTrkVtx without post-MC corrections

  Short_t parallelStripTrkEnd; ///< Same as parallelStripTrkVtx, but for the track end

  /*!
    If there is no track, this is not the far detector, or the track 
    vertex is not in one of the octants of the far detector defined by 
    a diagonal edge, this is set to -10.
   
    While this is almost certainly only interesting for anti-fiducial 
    events, it is set even if the track vertex is far from the edge so 
    that it doesn't matter if the fiducial cut is moved.

    If the first hit is on a strip that is perpendicular to the edge 
    (these overhang the strips parallel to the edge), it is set to 1.  
    It's set to 0 if this hit is on a parallel strip. 

    This could be useful for rock/detector discrimination since the 
    overhanging portion of the perpendicular strips are the farthest 
    extent of scintillator on these edges.  However, as the exact 
    extent of these strips is not very clear, this was not used in the 
    2010 CC analysis. <strong>It seems unlikely that anyone will 
    start using this variable, so should we remove it?</strong>
   */
  Char_t stripTrkBegPerpFlag; ///< Whether the first hit on a diagonal edge is on a perpendicular (overhanging) strip

  Char_t stripTrkEndPerpFlag; ///< Same as stripTrkBegPerpFlag, but for the last hit in the track


  /*!
    If there is no track, this is not the far detector, or the track 
    vertex is not in one of the octants of the far detector defined by 
    a horizontal or vertical edge, this is set to -999.

    This number gives the distance from the edge in terms of strips 
    found from looking at the first U and V hits in a track.  It is 
    defined so that the same number means the same distance on all 
    four horizontal/vertical edges.  It has an arbitrary offset.  (On 
    the bottom edge it is the U strip number plus the V strip number, 
    except a +1 crept in somehow and was retained for consistency.)  
    The definition is:

    On the west edge, 192 - U strip + V strip<br />
    On the top edge,  383 - U strip - V strip<br />
    On the east edge, 192 + U strip - V strip<br />
    On the bottom edge, 1 + U strip + V strip

    The intent of this variable is to give a more reliable and/or 
    easy-to-understand-the-systematics-of measure of distance from the 
    edge that does not rely on the track fit.  

    While this is almost certainly only interesting for anti-fiducial 
    events, it is set even if the track vertex is far from the edge so 
    that it doesn't matter if the fiducial cut is moved.

    In the 2010 CC analysis, a track was considered rock-like if 
    stripHoveNumTrkVtx <= 57.
  */
  Short_t stripHoveNumTrkVtx; ///< HOrizontal/VErtical strip number --- a measure of distance from the edge

  /*!
    stripHoveNumTrkVtx must be modified before analysis since alignment
    variation in the real detector isn't modeled in the daikon 
    geometry.  Essentially a smearing is added to account for the 
    alignment variations; this smeared value is stored in 
    stripHoveNumTrkVtx and the original value is stored here.  
   */
  Short_t stripHoveNumTrkVtxNoShift; ///< stripHoveNumTrkVtx without post-MC corrections

  Short_t stripHoveNumTrkEnd; ///< stripHoveNumTrkVtx but for the track end

  /*!
    <strong>This was considered for use in RAF analysis and then dropped.  Remove?</strong>
  */ 
  Float_t cosPrTrkVtx; ///< Cosine of angle between theta and radial momentum

  //@}

  /// \name More best track variables copied from the SNTP

  Float_t xTrkEnd; ///< trk.end.x
  Float_t yTrkEnd; ///< trk.end.y
  Float_t zTrkEnd; ///< trk.end.z
  Float_t uTrkEnd; ///< trk.end.u
  Float_t vTrkEnd; ///< trk.end.v
  Int_t planeTrkEnd; ///< trk.plane.end
  Int_t planeTrkEndu; ///< trk.plane.endu
  Int_t planeTrkEndv; ///< trk.plane.endv

  Float_t drTrkFidall; ///< trk.fidall.dr
  Float_t dzTrkFidall; ///< trk.fidall.dz
  Float_t drTrkFidvtx; ///< trk.fidbeg.dr
  Float_t dzTrkFidvtx; ///< trk.fidbeg.dz
  Float_t drTrkFidend; ///< trk.fidend.dr
  Float_t dzTrkFidend; ///< trk.fidend.dz
  Float_t traceTrkFidall; ///< trk.fidall.trace
  Float_t traceTrkFidvtx; ///< trk.fidvtx.trace
  Float_t traceTrkFidend; ///< trk.fidend.trace

  //@}

  /// \name Best shower
  //@{
  Bool_t shwExists; ///< simply state if shower exists
  Int_t ndigitShw; ///< shw.ndigit
  Int_t nstripShw; ///< shw.nstrip
  Int_t nplaneShw; ///< shw.plane.n

  Float_t shwEnCor; ///< Equivalent to \ref shwEnLinCCCor
  Float_t shwEnNoCor; ///< Copied from \ref NtpSRShowerPulseHeight::linCCgev
  Float_t shwEnMip; ///< Copied from \ref NtpSRStripPulseHeight::mip
  Float_t shwEnLinCCNoCor; ///< shw.shwph.linCCgev

  /// \brief Calorimetric shower energy assuming a CC event
  ///
  /// Calculated in \ref NuReco::GetShowerEnergyCor.
  /// This is the shower energy to use if you don't want kNN energy
  Float_t shwEnLinCCCor;
  Float_t shwEnWtCCNoCor; ///< Copied from \ref NtpSRShowerPulseHeight::wtCCgev
  Float_t shwEnWtCCCor; ///< Deweighted shower energy with energy corrections
  Float_t shwEnLinNCNoCor; ///< Copied from \ref NtpSRShowerPulseHeight::linNCgev

  /// Shower energy with energy corrections assuming this is an %NC event
  Float_t shwEnLinNCCor;
  Float_t shwEnWtNCNoCor; ///< Copied from \ref NtpSRShowerPulseHeight::wtNCgev
  Float_t shwEnWtNCCor; ///< Deweighted %NC shower energy with corrections

  /// \brief Estimated resolution on the shower energy in GeV
  ///
  /// Calculated by \ref EnergyResolution::ShowerResolution based on
  /// \ref shwEnkNN.
  Float_t shwResolution; ///< Cambridge resolution, shower only

  /// \brief Shower energy as calculated by the kNN algorithm
  ///
  /// Includes corrections. Calculated in \ref NuReco::GetShowerEnergykNN.
  Float_t shwEnkNN;

  /// \brief Uncorrected shower energy as calculated by the kNN algorithm
  ///
  /// Calculated in \ref NuReco::GetShowerEnergykNN.
  Float_t shwEnkNNNoCor;

  /// \brief Shower energy resolution estimate from the kNN method
  ///
  /// Calculated in \ref NuReco::GetShowerEnergykNN as the standard-deviation
  /// of the true energies of the k nearest neighbours.
  Float_t shwEnReskNN;

  Int_t planeShwBeg; ///< shw.plane.beg
  Int_t planeShwEnd; ///< shw.plane.end
  Int_t planeShwMax; ///< Shower plane with maximum pulseheight
  Float_t xShwVtx; ///< shw.vtx.x
  Float_t yShwVtx; ///< shw.vtx.y
  Float_t zShwVtx; ///< shw.vtx.z
  //@}

  /// \brief MINOS+ kNN shower energy variables
  /// shwEnkNNPlus and shwEnkNNPlusNoCor are computed in NtupleUtils/NuReco.cxx
  /// shwRmsT, shwPctWidth, ..., planeAllShwMax are computed in NtupleUtils/NuExtraction.cxx
  /// junting@physics.utexas.edu
  //@{
  Float_t shwEnkNNPlus;       ///< kNN shower energy with a reselection of variables and mahalanobis metric (for MINOS+ only)
  Float_t shwEnkNNPlusNoCor;  ///< non-corrected version of shwEnkNNPlus
  ///// topology variables for the best shower
  Float_t shwRmsT;            ///< RMS of transverse distribution
  Float_t shwPctWidth;        ///< sum of widths with 10%, 20%, ..., 90% containment
  Float_t shwPctLength;       ///< sum of widths with 10%, 20%, ..., 90% containment
  Float_t shwSumHitLength;    ///< sum of hit z position weighted by mip
  Float_t shwSumHitWidth;     ///< sum of hit transverse position weighted by mip
  ///// kinematics variables for the best shower
  Float_t avgShwHitMip;       ///< average hit mip
  Float_t maxShwPlaneMip;     ///< max mip in a plane
  Float_t maxShwHitMip;       ///< max mip of a hit
  ///// variables calculated from hits in all showers
  Float_t shwRmsTAllShws;          ///< all-shower version of shwRmsT
  Float_t shwPctWidthAllShws;      ///< all-shower version of shwPctWidth 
  Float_t shwPctLengthAllShws;     ///< all-shower version of shwPctLength
  Float_t avgShwHitMipAllShws;     ///< all-shower version of avgShwHitMipAllShws
  Float_t shwSumHitLengthAllShws;  ///< all-shower version of shwSumHitLengthAllShws
  Float_t shwSumHitWidthAllShws;   ///< all-shower version of shwSumHitWidthAllShws
  Float_t maxShwPlaneMipAllShws;   ///< all-shower version of maxShwPlaneMipAllShws
  Float_t maxShwHitMipAllShws;     ///< all-shower version of maxShwHitMipAllShws
  Float_t nplaneAllShws;           ///< all-shower version of nplaneShw
  Float_t planeAllShwMax;          ///< all-shower version of planeShwMax
  //@}

  /// \name Primary track
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t trkExists1;
  Int_t trkIndex1;
  Int_t ndigitTrk1;
  Int_t nstripTrk1;
  Float_t trkEnCorRange1;
  Float_t trkEnCorCurv1;
  Float_t trkShwEnNear1;
  Float_t trkShwEnNearDW1;
  Float_t trkMomentumRange1;
  Int_t containedTrk1;
  Int_t trkfitpass1;
  Float_t trkvtxdcosz1;
  Float_t trkvtxdcosy1;
  Int_t trknplane1;
  Int_t charge1;
  Float_t qp1;
  Float_t qp_rangebiased1;
  Float_t sigqp1;
  Float_t qp_sigqp1;
  Float_t chi21;
  Float_t ndof1;
  Float_t qpFraction1;
  Int_t trkVtxUVDiffPl1;
  Int_t trkLength1;
  Int_t planeTrkNu1;
  Int_t planeTrkNv1;
  Int_t ntrklike1;
  Float_t trkphsigcor1;
  Float_t trkphsigmap1;
  Int_t trkIdMC1;
  Float_t trkds1;

  Int_t trkfitpassSA1;
  Float_t trkvtxdcoszSA1;
  Int_t chargeSA1;
  Float_t qpSA1;
  Float_t sigqpSA1;
  Float_t chi2SA1;
  Float_t ndofSA1;
  Float_t probSA1;
  Float_t xTrkVtxSA1;
  Float_t yTrkVtxSA1;
  Float_t zTrkVtxSA1;
  Float_t uTrkVtxSA1;
  Float_t vTrkVtxSA1;

  Float_t jitter1;
  Float_t jPID1;
  Float_t majC1;
  Float_t smoothMajC1;

  Float_t roID1;
  Float_t knn01TrkActivePlanes1;
  Float_t knn10TrkMeanPh1;
  Float_t knn20LowHighPh1;
  Float_t knn40TrkPhFrac1;
  Float_t roIDNuMuBar1;
  Float_t relativeAngle1;

  //----------------------------
  //  Jasmine Ma/Ratchford new PID variables (track1)
  //----------------------------
  Float_t jmID1;
  Float_t jmTrackPlane1;
  Float_t jmEndPh1;
  Float_t jmMeanPh1;
  Float_t jmScatteringU1;
  Float_t jmScatteringV1;
  Float_t jmScatteringUV1;
  Float_t xTrkVtx1;
  Float_t yTrkVtx1;
  Float_t zTrkVtx1;
  Float_t uTrkVtx1;
  Float_t vTrkVtx1;
  Float_t tTrkVtx1;
  Int_t planeTrkVtx1;
  Int_t planeTrkBeg1;
  Int_t planeTrkBegu1;
  Int_t planeTrkBegv1;
  //Track timing information
  Float_t inverseBetaTrk1;
  Double_t t0Trk1;
  Float_t chi2TimeTrk1;
  UShort_t ndigitTimeTrk1;
  Double_t forwardRMSTrk1;
  Int_t forwardNDOFTrk1;
  // Strip geometry information
  Int_t stripTrkBeg1;
  Int_t stripTrkBegu1;
  Int_t stripTrkBegv1;
  Short_t stripTrkEnd1;
  Short_t stripTrkEndu1;
  Short_t stripTrkEndv1;
  Bool_t stripTrkBegIsu1;
  Bool_t stripTrkEndIsu1;
  Int_t regionTrkVtx1;
  Char_t edgeRegionTrkVtx1;
  Char_t edgeRegionTrkEnd1;
  Float_t phiTrkVtx1;
  Float_t phiTrkEnd1;

  // More, entirely derivative, geometry variables for the rock analysis
  Short_t parallelStripTrkVtx1;
  Short_t parallelStripTrkEnd1;
  Char_t stripTrkBegPerpFlag1;
  Char_t stripTrkEndPerpFlag1; 
  Short_t stripHoveNumTrkVtx1;
  Short_t stripHoveNumTrkEnd1;

  Float_t xTrkEnd1;
  Float_t yTrkEnd1;
  Float_t zTrkEnd1;
  Float_t uTrkEnd1;
  Float_t vTrkEnd1;
  Int_t planeTrkEnd1;
  Int_t planeTrkEndu1;
  Int_t planeTrkEndv1;

  Float_t drTrkFidall1;
  Float_t dzTrkFidall1;
  Float_t drTrkFidvtx1;
  Float_t dzTrkFidvtx1;
  Float_t drTrkFidend1;
  Float_t dzTrkFidend1;
  Float_t traceTrkFidall1;
  Float_t traceTrkFidvtx1;
  Float_t traceTrkFidend1;
  Float_t cosPrTrkVtx1;

  //@}

  /// \name Primary shower
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t shwExists1;
  Int_t ndigitShw1;
  Int_t nstripShw1;
  Int_t nplaneShw1;
  Float_t shwEnCor1;
  Float_t shwEnNoCor1;
  Float_t shwEnLinCCNoCor1;
  Float_t shwEnLinCCCor1;
  Float_t shwEnWtCCNoCor1;
  Float_t shwEnWtCCCor1;
  Float_t shwEnLinNCNoCor1;
  Float_t shwEnLinNCCor1;
  Float_t shwEnWtNCNoCor1;
  Float_t shwEnWtNCCor1;
  Float_t shwEnMip1;
  Int_t planeShwBeg1;
  Int_t planeShwEnd1;
  Int_t planeShwMax1;
  Float_t xShwVtx1;
  Float_t yShwVtx1;
  Float_t zShwVtx1;

  //@}

  /// \name Second track
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t trkExists2;
  Int_t trkIndex2;
  Int_t ndigitTrk2;
  Int_t nstripTrk2;
  Float_t trkEnCorRange2;
  Float_t trkEnCorCurv2;
  Float_t trkShwEnNear2;
  Float_t trkShwEnNearDW2;
  Float_t trkMomentumRange2;
  Int_t containedTrk2;
  Int_t trkfitpass2;
  Float_t trkvtxdcosz2;
  Float_t trkvtxdcosy2;
  Int_t trknplane2;
  Int_t charge2;
  Float_t qp2;
  Float_t qp_rangebiased2;
  Float_t sigqp2;
  Float_t qp_sigqp2;
  Float_t chi22;
  Float_t ndof2;
  Float_t qpFraction2;
  Int_t trkVtxUVDiffPl2;
  Int_t trkLength2;
  Int_t planeTrkNu2;
  Int_t planeTrkNv2;
  Int_t ntrklike2;
  Float_t trkphsigcor2;
  Float_t trkphsigmap2;
  Int_t trkIdMC2;
  Float_t trkds2;

  Int_t trkfitpassSA2;
  Float_t trkvtxdcoszSA2;
  Int_t chargeSA2;
  Float_t qpSA2;
  Float_t sigqpSA2;
  Float_t chi2SA2;
  Float_t ndofSA2;
  Float_t probSA2;
  Float_t xTrkVtxSA2;
  Float_t yTrkVtxSA2;
  Float_t zTrkVtxSA2;
  Float_t uTrkVtxSA2;
  Float_t vTrkVtxSA2;

  Float_t jitter2;
  Float_t jPID2;
  Float_t majC2;
  Float_t smoothMajC2;

  Float_t roID2;
  Float_t knn01TrkActivePlanes2;
  Float_t knn10TrkMeanPh2;
  Float_t knn20LowHighPh2;
  Float_t knn40TrkPhFrac2;
  Float_t roIDNuMuBar2;
  Float_t relativeAngle2;

  //----------------------------
  //  Jasmine Ma/Ratchford new PID variables (track2)
  //----------------------------
  Float_t jmID2;
  Float_t jmTrackPlane2;
  Float_t jmEndPh2;
  Float_t jmMeanPh2;
  Float_t jmScatteringU2;
  Float_t jmScatteringV2;
  Float_t jmScatteringUV2;

  Float_t xTrkVtx2;
  Float_t yTrkVtx2;
  Float_t zTrkVtx2;
  Float_t uTrkVtx2;
  Float_t vTrkVtx2;
  Float_t tTrkVtx2;
  Int_t planeTrkVtx2;
  Int_t planeTrkBeg2;
  Int_t planeTrkBegu2;
  Int_t planeTrkBegv2;
  //Track timing information
  Float_t inverseBetaTrk2;
  Double_t t0Trk2;
  Float_t chi2TimeTrk2;
  UShort_t ndigitTimeTrk2;
  Double_t forwardRMSTrk2;
  Int_t forwardNDOFTrk2;
  //Strip geometry information
  Int_t stripTrkBeg2;
  Int_t stripTrkBegu2;
  Int_t stripTrkBegv2;
  Short_t stripTrkEnd2;
  Short_t stripTrkEndu2;
  Short_t stripTrkEndv2;
  Bool_t stripTrkBegIsu2;
  Bool_t stripTrkEndIsu2;
  Int_t regionTrkVtx2;
  Char_t edgeRegionTrkVtx2;
  Char_t edgeRegionTrkEnd2;
  Float_t phiTrkVtx2;
  Float_t phiTrkEnd2;

  // More, entirely derivative, geometry variables for the rock analysis
  Short_t parallelStripTrkVtx2;
  Short_t parallelStripTrkEnd2;
  Char_t stripTrkBegPerpFlag2;
  Char_t stripTrkEndPerpFlag2;
  Short_t stripHoveNumTrkVtx2;
  Short_t stripHoveNumTrkEnd2;

  Float_t xTrkEnd2;
  Float_t yTrkEnd2;
  Float_t zTrkEnd2;
  Float_t uTrkEnd2;
  Float_t vTrkEnd2;
  Int_t planeTrkEnd2;
  Int_t planeTrkEndu2;
  Int_t planeTrkEndv2;

  Float_t drTrkFidall2;
  Float_t dzTrkFidall2;
  Float_t drTrkFidvtx2;
  Float_t dzTrkFidvtx2;
  Float_t drTrkFidend2;
  Float_t dzTrkFidend2;
  Float_t traceTrkFidall2;
  Float_t traceTrkFidvtx2;
  Float_t traceTrkFidend2;
  Float_t cosPrTrkVtx2;

  //@}

  /// \name Second shower
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t shwExists2;
  Int_t ndigitShw2;
  Int_t nstripShw2;
  Int_t nplaneShw2;
  Float_t shwEnCor2;
  Float_t shwEnNoCor2;
  Float_t shwEnLinCCNoCor2;
  Float_t shwEnLinCCCor2;
  Float_t shwEnWtCCNoCor2;
  Float_t shwEnWtCCCor2;
  Float_t shwEnLinNCNoCor2;
  Float_t shwEnLinNCCor2;
  Float_t shwEnWtNCNoCor2;
  Float_t shwEnWtNCCor2;
  Float_t shwEnMip2;
  Int_t planeShwBeg2;
  Int_t planeShwEnd2;
  Int_t planeShwMax2;
  Float_t xShwVtx2;
  Float_t yShwVtx2;
  Float_t zShwVtx2;

  //@}

  /// \name Third track
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t trkExists3;
  Int_t trkIndex3;
  Int_t ndigitTrk3;
  Int_t nstripTrk3;
  Float_t trkEnCorRange3;
  Float_t trkEnCorCurv3;
  Float_t trkShwEnNear3;
  Float_t trkShwEnNearDW3;
  Float_t trkMomentumRange3;
  Int_t containedTrk3;
  Int_t trkfitpass3;
  Float_t trkvtxdcosz3;
  Float_t trkvtxdcosy3;
  Int_t trknplane3;
  Int_t charge3;
  Float_t qp3;
  Float_t qp_rangebiased3;
  Float_t sigqp3;
  Float_t qp_sigqp3;
  Float_t chi23;
  Float_t ndof3;
  Float_t qpFraction3;
  Int_t trkVtxUVDiffPl3;
  Int_t trkLength3;
  Int_t planeTrkNu3;
  Int_t planeTrkNv3;
  Int_t ntrklike3;
  Float_t trkphsigcor3;
  Float_t trkphsigmap3;
  Int_t trkIdMC3;
  Float_t trkds3;

  Int_t trkfitpassSA3;
  Float_t trkvtxdcoszSA3;
  Int_t chargeSA3;
  Float_t qpSA3;
  Float_t sigqpSA3;
  Float_t chi2SA3;
  Float_t ndofSA3;
  Float_t probSA3;
  Float_t xTrkVtxSA3;
  Float_t yTrkVtxSA3;
  Float_t zTrkVtxSA3;
  Float_t uTrkVtxSA3;
  Float_t vTrkVtxSA3;

  Float_t jitter3;
  Float_t jPID3;
  Float_t majC3;
  Float_t smoothMajC3;

  Float_t roID3;
  Float_t knn01TrkActivePlanes3;
  Float_t knn10TrkMeanPh3;
  Float_t knn20LowHighPh3;
  Float_t knn40TrkPhFrac3;
  Float_t roIDNuMuBar3;
  Float_t relativeAngle3;

  //----------------------------
  //  Jasmine Ma/Ratchford new PID variables (track3)
  //----------------------------
  Float_t jmID3;
  Float_t jmTrackPlane3;
  Float_t jmEndPh3;
  Float_t jmMeanPh3;
  Float_t jmScatteringU3;
  Float_t jmScatteringV3;
  Float_t jmScatteringUV3;

  Float_t xTrkVtx3;
  Float_t yTrkVtx3;
  Float_t zTrkVtx3;
  Float_t uTrkVtx3;
  Float_t vTrkVtx3;
  Float_t tTrkVtx3;
  Int_t planeTrkVtx3;
  Int_t planeTrkBeg3;
  Int_t planeTrkBegu3;
  Int_t planeTrkBegv3;
  //Track timing information
  Float_t inverseBetaTrk3;
  Double_t t0Trk3;
  Float_t chi2TimeTrk3;
  UShort_t ndigitTimeTrk3;
  Double_t forwardRMSTrk3;
  Int_t forwardNDOFTrk3;
  // Strip geometry
  Int_t stripTrkBeg3;
  Int_t stripTrkBegu3;
  Int_t stripTrkBegv3;
  Short_t stripTrkEnd3;
  Short_t stripTrkEndu3;
  Short_t stripTrkEndv3;
  Bool_t stripTrkBegIsu3;
  Bool_t stripTrkEndIsu3;
  Int_t regionTrkVtx3;
  Char_t edgeRegionTrkVtx3;
  Char_t edgeRegionTrkEnd3;
  Float_t phiTrkVtx3;
  Float_t phiTrkEnd3;

  // More, entirely derivative, geometry variables for the rock analysis
  Short_t parallelStripTrkVtx3;
  Short_t parallelStripTrkEnd3;
  Char_t stripTrkBegPerpFlag3;
  Char_t stripTrkEndPerpFlag3;
  Short_t stripHoveNumTrkVtx3;
  Short_t stripHoveNumTrkEnd3;

  Float_t xTrkEnd3;
  Float_t yTrkEnd3;
  Float_t zTrkEnd3;
  Float_t uTrkEnd3;
  Float_t vTrkEnd3;
  Int_t planeTrkEnd3;
  Int_t planeTrkEndu3;
  Int_t planeTrkEndv3;

  Float_t drTrkFidall3;
  Float_t dzTrkFidall3;
  Float_t drTrkFidvtx3;
  Float_t dzTrkFidvtx3;
  Float_t drTrkFidend3;
  Float_t dzTrkFidend3;
  Float_t traceTrkFidall3;
  Float_t traceTrkFidvtx3;
  Float_t traceTrkFidend3;
  Float_t cosPrTrkVtx3;

  //@

  /// \name Third shower
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t shwExists3;
  Int_t ndigitShw3;
  Int_t nstripShw3;
  Int_t nplaneShw3;
  Float_t shwEnCor3;
  Float_t shwEnNoCor3;
  Float_t shwEnLinCCNoCor3;
  Float_t shwEnLinCCCor3;
  Float_t shwEnWtCCNoCor3;
  Float_t shwEnWtCCCor3;
  Float_t shwEnLinNCNoCor3;
  Float_t shwEnLinNCCor3;
  Float_t shwEnWtNCNoCor3;
  Float_t shwEnWtNCCor3;
  Float_t shwEnMip3;
  Int_t planeShwBeg3;
  Int_t planeShwEnd3;
  Int_t planeShwMax3;
  Float_t xShwVtx3;
  Float_t yShwVtx3;
  Float_t zShwVtx3;

  //@}

  /// \name Fourth shower
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t shwExists4;
  Int_t ndigitShw4;
  Int_t nstripShw4;
  Int_t nplaneShw4;
  Float_t shwEnCor4;
  Float_t shwEnNoCor4;
  Float_t shwEnLinCCNoCor4;
  Float_t shwEnLinCCCor4;
  Float_t shwEnWtCCNoCor4;
  Float_t shwEnWtCCCor4;
  Float_t shwEnLinNCNoCor4;
  Float_t shwEnLinNCCor4;
  Float_t shwEnWtNCNoCor4;
  Float_t shwEnWtNCCor4;
  Float_t shwEnMip4;
  Int_t planeShwBeg4;
  Int_t planeShwEnd4;
  Int_t planeShwMax4;
  Float_t xShwVtx4;
  Float_t yShwVtx4;
  Float_t zShwVtx4;

  //@}

  /// \name Fifth shower
  /// For documentation see the same variables without the suffix
  //@{

  Bool_t shwExists5;
  Int_t ndigitShw5;
  Int_t nstripShw5;
  Int_t nplaneShw5;
  Float_t shwEnCor5;
  Float_t shwEnNoCor5;
  Float_t shwEnLinCCNoCor5;
  Float_t shwEnLinCCCor5;
  Float_t shwEnWtCCNoCor5;
  Float_t shwEnWtCCCor5;
  Float_t shwEnLinNCNoCor5;
  Float_t shwEnLinNCCor5;
  Float_t shwEnWtNCNoCor5;
  Float_t shwEnWtNCCor5;
  Float_t shwEnMip5;
  Int_t planeShwBeg5;
  Int_t planeShwEnd5;
  Int_t planeShwMax5;
  Float_t xShwVtx5;
  Float_t yShwVtx5;
  Float_t zShwVtx5;

  //@}

  /// \name Other info calculated
  //@{
  Float_t rEvtVtx; ///< Distance of event vertex from (x,y) = (0,0)
  Float_t rEvtEnd; ///< Distance of event end from (x,y) = (0,0)
  Float_t distToEdgeEvtVtx; ///< event vertex distance to edge of scintillator
  Int_t evtVtxUVDiffPl; ///< Number of planes between first U plane of event and first V plane of event

  Float_t rTrkVtx; ///< Distance of track vertex from (x,y) = (0,0)
  Float_t rTrkEnd; ///< Distance of track end from (x,y) = (0,0)
  Float_t sigqp_qp; ///< best track sigma Q/P / Q/P from fit
  Float_t chi2PerNdof; ///< track's chi2/ndof of fit
  Float_t prob; ///< track's probability of fit

  Int_t containmentFlag; ///< containment flag to use
  Int_t containmentFlagCC0093Std; ///< CC 0.9e20 pot containment flag
  Int_t containmentFlagCC0250Std; ///< CC 2.5e20 pot containment flag
  Int_t containmentFlagPitt; ///< pitt specific flag: 1-4, up/down stop/exit
  Int_t usedRange; ///< 1=used range, 0=not
  Int_t usedCurv; ///< 1=used curvature, 0=not
  //@}


  /// \name Weights
  //@{
 
  /*!
     This is the SKZP beam weight (beamWeight) plus a few other small
     corrections.  By itself, it does <strong>not</strong> give a
     proper prediction of the far detector spectrum.    
  */
  Float_t rw; ///< The final weight applied to the event. 
  Float_t fluxErr; ///< The error on the flux from SKZP.
  Float_t rwActual; ///< This is the weight as it leaves NuReco::ApplyReweights.  Any later changes (i.e. systematics) only affect rw.
  Float_t generatorWeight; ///< weight factor from generator
  Float_t detectorWeight; ///< SKZP detector weight to use as default

  /*!
    Daikon mismodels the coil by using "coil material" for neutrino
    interactions, but steel and scintillator for particle propagation.

    As a good approximation, this can be fixed by weighting events
    in which the true neutrino interaction is in the coil by the raio
    of densities.   This is set to 1 if the true neutrino interaction
    is not in the coil and that ratio if it is.

    This contributes to the total weight rw.
  */
  Float_t coilCorrWeight; ///< weight to account for daikon coil mismodeling

  /*!
    Daikon uses an old cross section model.  Newer models predict
    somewhat smaller cross sections for light elements relative to
    heavier elements. This weight applies a correction that roughly
    converts our monte carlo to the newer model.  The correction
    grows as the element moves farther from iron.  

    Since cross sections for detector events cancels near-to-far, this
    correction is only done for rock events and is set to 1 for 
    detector events.

    This contributes to the total weight rw.
  */
  Float_t rockxsecCorrWeight; ///< weight updating cross section model for rock

  /// \brief Weight applied to this event by the 2010 CC analysis
  ///
  /// This weight is intended for plot-making purposes only. 
  /// It is not included in \ref rw
  ///
  /// The analysis weights can't be known at DST-making time, so this weight
  /// defaults to one. The files should be rewritten so this weight represents
  /// the weight applied to this event by the beam-matrix extrapolation.
  /// Including this weight along with \ref rw should give fully correct
  /// nominal Far Detector predictions. Best-fit systematics are more 
  /// complicated and not included here.
  Float_t anaWeightCC2010;

  /*! 
    Daikon04 and earlier: average POT weighted weights, e.g. (1.27*RunI + 1.23*RunII)/2.50
    Daikon07 and later: Just the correct weights for this MC event
 */
  Float_t trkEnWeight; ///< SKZP weight to apply to \ref trkEn
  Float_t shwEnWeight; ///< SKZP weight to apply to \ref shwEn
  Float_t beamWeight; ///< SKZP weight for beam (%e.g. hadron production)
  Float_t fluxErrHadProdAfterTune; ///< flux error, kHadProdAfterTune
  Float_t fluxErrTotalErrorPreTune; ///< flux error, kTotalErrorPreTune
  Float_t fluxErrTotalErrorAfterTune; ///< flux error kTotalErrorAfterTune
  Float_t detectorWeightNMB; ///< SKZP detector weight (%e.g. xsec)
  Float_t detectorWeightNM; ///< NuMu rather than NuMuBar above

  /// \brief Reconstructed energy with SKZP detector weights applied
  ///
  /// Set to the sum of \ref trkEnRw and \ref shwEnRw by
  /// \ref NuReco::ApplyReweights
  Float_t energyRw;

  /// \brief Reconstructed energy without detector weights
  ///
  /// Equivalent to \ref energy when \ref applyEnergyShifts is not set
  Float_t energyNoRw;

  /// \brief Reconstructed track energy with SKZP detector weights applied
  ///
  /// Set to \ref trkEnNoRw times \ref trkEnWeight by
  /// \ref NuReco::ApplyReweights
  Float_t trkEnRw;

  /// \brief Reconstructed track energy without detector weights
  ///
  /// Equivalent to \ref trkEn when \ref applyEnergyShifts is not set
  Float_t trkEnNoRw;

  /// \brief Reconstructed shower energy with SKZP detector weights applied
  ///
  /// Set to \ref shwEnNoRw times \ref shwEnWeight by
  /// \ref NuReco::ApplyReweights
  Float_t shwEnRw;

  /// \brief Reconstructed shower energy without detector weights
  ///
  /// Equivalent to \ref shwEn when \ref applyEnergyShifts is not set
  Float_t shwEnNoRw;
  //@}

  /// \name PIDs
  //@{
  Float_t dpID; ///< DP's PID variable
  Float_t abID; ///< AB's PID variable
  Float_t roID; ///< RO's PID variable
  Float_t knn01TrkActivePlanes; ///< number of active planes in trk
  Float_t knn10TrkMeanPh; ///< average ph per plane in trk
  Float_t knn20LowHighPh; ///< av of low ph strips/av of high ph strips
  Float_t knn40TrkPhFrac; ///< fraction of ph in trk
  Float_t roIDNuMuBar; ///< RO's PID variable for NuMuBar selection (0 or 1)
  Float_t relativeAngle; ///< RO's track angle relative to muon dir.
  Float_t poID; ///< PO's PID variable
  Float_t poIDKin; ///< PO's PID variable (kinematics only)
  //@}

  /// \name MINOS+ PIDs
  /// junting@physics.utexas.edu
  //@{
  Float_t roIDPlus; ///< RO's PID made with mahalanobis metric (for MINOS+ only)
  //@}

  /// \name Jasmine Ma/Ratchford variables
  //@{
  Float_t jmID;
  Float_t jmTrackPlane;
  Float_t jmEndPh;
  Float_t jmMeanPh;
  Float_t jmScatteringU;
  Float_t jmScatteringV;
  Float_t jmScatteringUV;
  Float_t jmEventknnID;
  Float_t jmEventknn208;
  Float_t jmEventknn207;
  Float_t jmEventknn206;
  Float_t jmEventknn205;
  Float_t jmEventknn204;
  //@}

  /// \name NC variables
  /// Comments from %NC code
  //@{

  // event
  Float_t closeTimeDeltaZ; ///< distance in z to event closest in time
  Int_t edgeActivityStrips; ///< number of strips in partially instrumented region in 40ns time window
  Float_t edgeActivityPH; ///< summed PH in partially instrumented region in 40ns time window
  Int_t oppEdgeStrips; ///< number of strips in opposite edge region in 40ns time window
  Float_t oppEdgePH; ///< summed PH in opposite edge region in 40ns time window
  Float_t vtxMetersToCoilEvt; ///< distance to the coil hole
  Float_t vtxMetersToCloseEdgeEvt; ///< distance to nearest edge (in XY)of partial plane outline in ND
  Double_t minTimeSeparation; ///< time difference to closest event in time

  Float_t slicePHFraction;
  Int_t maxConsecutivePlanes;

  // shower
  Float_t transverseRMSU; ///< rms of transverse strip positions in U
  Float_t transverseRMSV; ///< rms of transverse strip positions in V

  // track
  Float_t dtdz; ///< gradient of t(z) fit (with c_light=1)
  Float_t endMetersToCloseEdge; //distance to nearest edge (in XY) of full plane outline in ND
  Float_t vtxMetersToCloseEdgeTrk; ///< distance to nearest edge (in XY) of partial plane outline in ND
  Float_t vtxMetersToCoilTrk; ///< distance to the center of the coil hole
  Float_t traceEndZ; ///< delta in z from end to projected exit location

  //@}

  /// \name Beam variables
  //@{

  Float_t pot; ///< pot in current spill
  Float_t potDB; ///< pot in current spill from database
  Float_t potSinceLastEvt; ///< includes pot count for spills with no events
  Float_t potSinceLastEvtGood; ///< as above but also good beam+det spills
  Float_t potSinceLastEvtBad; ///< as above but for bad beam+det spills
  Float_t potSinceLastEvtDB; ///< as above but with database version of pots
  Float_t potSinceLastEvtGoodDB; ///< as above but with database version of pots
  Float_t potSinceLastEvtBadDB; ///< as above but with database version of pots
  Int_t nBatch;
  std::vector<Double_t> potListSinceLastEvtGood;//A list of the
						//spill-by-spill PoT
						//since the last
						//event. Usually only
						//one entry (the
						//current spill), but
						//will contain more
						//than one if there
						//have been spills
						//with no events

  Int_t runPeriod; ///< Major analysis run period, e.g. "Run IV" == 4
  Bool_t hornIsReverse;
  //Conventions/BeamType.h definition
  Int_t beamTypeDB; ///<  From BeamMonSpill
  Int_t beamType; ///< From NuConfig
  Float_t intensity; ///< Only filled for MC events
  Float_t hornCur; ///< BeamMonSpill::fHornCur
  Bool_t goodBeam; ///< BMSpillAna.SelectSpill()
  Bool_t goodBeamSntp; ///< uses the beam data in the sntp file

  //@}


  /// \name Truth variables
  ///
  /// EVERYTIME %A TRUTH VARIABLE IS ADDED TO THIS CLASS IT MUST
  /// ALSO BE ADDED TO \ref NuMCEvent
  //@{

  Bool_t isInFidVolCCMC; ///< Is the true vertex in the CC fiducial volume?

  /// \brief True energy in principle visible in the detector
  ///
  /// For CC events this is equivalent to \ref neuEnMC.
  /// For NC events a factor \ref yMC is multiplied to account for the
  /// unobservable energy carried away by the final-state neutrino.
  Float_t energyMC;

  /// \brief True neutrino energy
  ///
  /// The energy component of the neutrino 4-vector: p4neu[3]
  Float_t neuEnMC;
  Float_t neuPxMC; ///< p4neu[0];
  Float_t neuPyMC; ///< p4neu[1];
  Float_t neuPzMC; ///< p4neu[2];

  Float_t mu1EnMC; ///< p4mu1[3];
  Float_t mu1PxMC; ///< p4mu1[0];
  Float_t mu1PyMC; ///< p4mu1[1];
  Float_t mu1PzMC; ///< p4mu1[2];

  Float_t tgtEnMC; ///< p4tgt[3]
  Float_t tgtPxMC; ///< p4tgt[0];
  Float_t tgtPyMC; ///< p4tgt[1];
  Float_t tgtPzMC; ///< p4tgt[2];

  Int_t zMC; ///< z;
  Int_t aMC; ///< a;
  Int_t nucleusMC; ///< encoding from Mad of the nucleus type
  Int_t initialStateMC; ///< encoding from Mad of the initial state
  Int_t hadronicFinalStateMC; ///< encoding from Mad of the hadronic final state
  Int_t numPreInukeFSprotMC; ///< Number of final state protons (before inuke)
  Int_t numPreInukeFSneutMC; ///< Number of final state neutrons (before inuke)
  Float_t maxMomPreInukeFSprotMC; ///< The highest momentum final state proton (before inuke)
  Float_t maxMomPreInukeFSneutMC; ///< The highest momentum final state neutron (before inuke)

  Float_t yMC; ///< kinematic y value from sntp file (not y position)
  Float_t y2MC; ///< p4shw[3]/(fabs(p4mu1[3])+p4shw[3])
  Float_t xMC; ///< kinematic x (not x position)
  Float_t q2MC; ///< q2
  Float_t w2MC; ///< w2

  /// \brief True energy of the muon. Signed according to true muon charge.
  ///
  /// The energy component of the muon 4-vector: p4mu1[3]
  Float_t trkEnMC;

  /// True leptonic energy: (1-y)*p4neu[3]
  Float_t trkEn2MC;

  /// \brief True energy of the shower
  ///
  /// The energy component of the shower 4-vector: p4shw[3]
  Float_t shwEnMC;

  /// True hadronic energy: y*p4neu[3];
  Float_t shwEn2MC;

  Float_t trkEndEnMC; ///< NtpMCStdHep.dethit[1].pE, particle energy at last scint hit
  Float_t trkStartEnMC; ///< muon energy at first scintillator hit
  Bool_t  trkContainmentMC; ///<  Experimental: MC Containment of primary muon

  Float_t sigma; ///< mc.sigma=cross-section
  Int_t iaction; ///< CC=1, NC=0
  Int_t iresonance; ///< QE=1001, RES=1002, DIS=1003, CPP=1004
  Int_t inu; ///< PDG code of neutrino (positive = particles, negative = anti-particles)
  Int_t inunoosc; ///< PDG code id of neutrino at birth
  Int_t itg; ///< neutrino interaction target

  Float_t vtxxMC; ///< vtxx: x vertex of neutrino interaction
  Float_t vtxyMC; ///< vtxy: y vertex of neutrino interaction
  Float_t vtxzMC; ///< vtxz: z vertex of neutrino interaction
  Float_t vtxuMC; ///< vertex u, calculated from x&y
  Float_t vtxvMC; ///< vertex v, calculated from x&y
  Int_t planeTrkVtxMC; ///< calculated from z
  Float_t rTrkVtxMC; ///< calculated from x&y

  Int_t mc; ///< the index of the object to be used
  Int_t mcTrk; ///< track mc index
  Int_t mcShw; ///< shower mc index
  Int_t mcEvt; ///< event mc index

  Int_t mcTrk1; ///< 1st track mc index
  Int_t mcTrk2; ///< 2nd track mc index
  Int_t mcTrk3; ///< 3rd track mc index

  Int_t mcShw1; ///< 1st shower mc index
  Int_t mcShw2; ///< 2nd shower mc index
  Int_t mcShw3; ///< 3rd shower mc index
  Int_t mcShw4; ///< 4th shower mc index
  Int_t mcShw5; ///< 5th shower mc index
  
  /// \name  Truth flux variables, copied from underlying flux files
  //@{
  
  Float_t Npz; ///< Neutrino momentum (GeV/c) along the \f$z\f$-axis (beam axis)
  Float_t NdxdzNea; ///< Direction slopes for a neutrino forced towards the center of the near detecto
  Float_t NdydzNea; ///< Direction slopes for a neutrino forced towards the center of the near detecto
  Float_t NenergyN; ///<  Energy for a neutrino forced towards the center of the near detector
  Float_t NWtNear; ///<  Weight for a neutrino forced towards the center of the near detector
  Float_t NdxdzFar; ///< Direction slopes for a neutrino forced towards the center of the far detector
  Float_t NdydzFar; ///< Direction slopes for a neutrino forced towards the center of the far detector
  Float_t NenergyF; ///<  Neutrino energy (GeV) for a decay forced to the center of the far detector
  Float_t NWtFar; ///<  Neutrino weight for a decay forced to the center of the far detector
  Int_t Ndecay; ///<  Decay process that produced the neutrino, see \tab{decays}
  Float_t Vx; ///< Neutrino production vertex (cm) 
  Float_t Vy; ///< Neutrino production vertex (cm) 
  Float_t Vz; ///< Neutrino production vertex (cm) 
  Float_t pdPx; ///< Momentum (GeV/c) of the neutrino parent at the neutrino production vertex (parent decay point)
  Float_t pdPy; ///< Momentum (GeV/c) of the neutrino parent at the neutrino production vertex (parent decay point) 
  Float_t pdPz; ///< Momentum (GeV/c) of the neutrino parent at the neutrino production vertex (parent decay point)
  Float_t ppdxdz; ///< Direction of the neutrino parent at its production point (which may be in the target)
  Float_t ppdydz; ///< Direction of the neutrino parent at its production point (which may be in the target)
  Float_t pppz; ///<  z momentum (GeV/c) of the neutrino parent at its production point
  Float_t ppenergy; ///<  Energy (GeV) of the neutrino parent at its production point
  Float_t ppmedium; ///< Code for the material the neutrino parent was produced in (see doc-6316) 
  Float_t ppvx; ///< Production vertex (cm) of the neutrino parent 
  Float_t ppvy; ///< Production vertex (cm) of the neutrino parent 
  Float_t ppvz; ///< Production vertex (cm) of the neutrino parent (used to determine decay pipe events)
  Int_t ptype; ///<  Neutrino parent species (GEANT codes)
  Float_t Necm; ///<  Neutrino energy (GeV) in the center-of-mass frame
  Float_t Nimpwt; ///<  Importance weight of the neutrino
  Float_t tvx; ///< Position (cm) of the neutrino ancestor as it exits target (possibly, but not necessarily, the direct neutrino parent)
  Float_t tvy; ///< Position (cm) of the neutrino ancestor as it exits target (possibly, but not necessarily, the direct neutrino parent)
  Float_t tvz; ///< Position (cm) of the neutrino ancestor as it exits target (possibly, but not necessarily, the direct neutrino parent)
  Float_t tpx; ///< Momentum (GeV/c) of the ancestor as it exits target
  Float_t tpy; ///< Momentum (GeV/c) of the ancestor as it exits target
  Float_t tpz; ///< Momentum (GeV/c) of the ancestor as it exits target
  Int_t tptype; ///<  Species of the ancestor exiting the target  (GEANT codes)
  Int_t tgen; ///<  Neutrino parent generation in cascade. 1 = primary proton, 2 = particles produced by proton interaction, 3 = particles from 
                                                               
  //@}
  
  /// \name  Truth  Intranuke weighting (Jasmine ma)
  //@{

  Int_t   InukeNwts;
  Float_t InukePiCExchgP;   ///< 0
  Float_t InukePiCExchgN;   ///< 1
  Float_t InukePiEScatP;    ///< 2
  Float_t InukePiEScatN;    ///< 3
  Float_t InukePiInEScatP;  ///< 4
  Float_t InukePiInEScatN;  ///< 5
  Float_t InukePiAbsorbP;   ///< 6
  Float_t InukePiAbsorbN;   ///< 7
  Float_t InukePi2PiP;      ///< 8
  Float_t InukePi2PiN;      ///< 9
  Float_t InukeNknockP;     ///< 10
  Float_t InukeNknockN;     ///< 11
  Float_t InukeNNPiP;       ///< 12
  Float_t InukeNNPiN;       ///< 13
  Float_t InukeFormTP;      ///< 14
  Float_t InukeFormTN;      ///< 15
  Float_t InukePiXsecP;     ///< 16
  Float_t InukePiXsecN;     ///< 17
  Float_t InukeNXsecP;      ///< 18
  Float_t InukeNXsecN;      ///< 19
  Float_t InukeNucrad;
  Float_t InukeWrad;

  //@}

  //////////////////////////////////////////////////////////////////
  //EVERYTIME A TRUTH VARIABLE IS ADDED TO THIS CLASS IT MUST
  //ALSO BE ADDED TO NuMCEvent
  ///////////////////////////////////////////////////////////////////

  /// \name Program control variables
  //@{

  Int_t anaVersion; ///< different cuts etc
  Int_t releaseType; ///< Conventions/ReleaseType.h definition
  Int_t recoVersion; ///< Birch/Cedar
  Int_t mcVersion; ///< Carrot/Daikon, used to select pdfs in data case
  Int_t reweightVersion; ///< the beam reweighting to use

  Bool_t useGeneratorReweight; ///< switch to turn off generator reweighting
  std::string sGeneratorConfigName; ///< name of generator rw configuration
  Int_t generatorConfigNo; ///< number of generator rw configuration

  Bool_t useDBForDataQuality; ///< flag to use DB for data quality
  Bool_t useDBForSpillTiming; ///< flag to use DB for spill timing
  Bool_t useDBForBeamInfo; ///< flag to use DB for beam info

  Bool_t cutOnDataQuality; ///< flag to enable cut on data quality
  Bool_t cutOnSpillTiming; ///< flag to enable cut on spill timing
  Bool_t cutOnBeamInfo; ///< flag to enable cut on beam info

  Bool_t applyEnergyShifts; ///< Flag to use energy shifts
  Bool_t applyBeamWeight; ///< Flag to use beam weight
  Bool_t apply1SigmaWeight; ///< Flag to use +1-sigma SKZP error shift
  Bool_t applyDetectorWeight; ///< Flag to use detector weight, %e.g. xsec
  Bool_t applyGeneratorWeight; ///< Flag to use generator weight

  Bool_t calcMajCurv; ///< flag to run majorityCurv or not
  Bool_t calcRoID; ///< flag to run RoID or not
  Bool_t calcJmID; ///<  flag to run JmID or not

  //@}

  ClassDef(NuEvent, 94);
};


//// Accessor functions for the numbered track and shower variables:
//
//// This anonymous namespace prevents us from getting "multiple definition"
//// linker errors.
//namespace{
//
//// For working out the return types. It seems like we should be able to use
//// typeof(NuEvent::var) instead of typeof(dummy.var), but apparently not...
//NuEvent dummy;
//
//#define TrkVarFunc(var)                           \
//typeof(dummy.var##1)& get_##var(NuEvent& nu,      \
//				int idx)          \
//{                                                 \
//  switch(idx){                                    \
//  case 1: return nu.var##1;                       \
//  case 2: return nu.var##2;                       \
//  case 3: return nu.var##3;                       \
//  default:                                        \
//    assert(0 && "Track index out of range");      \
//  }                                               \
//}                                                 \
//                                                  \
//typeof(dummy.var##1) get_##var(const NuEvent& nu, \
//			       int idx)           \
//{                                                 \
//  switch(idx){                                    \
//  case 1: return nu.var##1;                       \
//  case 2: return nu.var##2;                       \
//  case 3: return nu.var##3;                       \
//  default:                                        \
//    assert(0 && "Track index out of range");      \
//  }                                               \
//}                                                 \
//                                                  \
//typeof(dummy.var##1)& get_##var(NuEvent* nu,      \
//				int idx)          \
//{                                                 \
//  switch(idx){                                    \
//  case 1: return nu->var##1;                      \
//  case 2: return nu->var##2;                      \
//  case 3: return nu->var##3;                      \
//  default:                                        \
//    assert(0 && "Track index out of range");      \
//  }                                               \
//}
//
//// eg this creates a function: Bool_t& get_trkExists(NuEvent& nu, int idx);
//TrkVarFunc(trkExists)
//TrkVarFunc(trkIndex)
//TrkVarFunc(ndigitTrk)
//TrkVarFunc(nstripTrk)
//TrkVarFunc(trkEnCorRange)
//TrkVarFunc(trkEnCorCurv)
//TrkVarFunc(trkShwEnNear)
//TrkVarFunc(trkShwEnNearDW)
//TrkVarFunc(trkMomentumRange)
//TrkVarFunc(containedTrk)
//TrkVarFunc(trkfitpass)
//TrkVarFunc(trkvtxdcosz)
//TrkVarFunc(trkvtxdcosy)
//TrkVarFunc(trknplane)
//TrkVarFunc(charge)
//TrkVarFunc(qp)
//TrkVarFunc(qp_rangebiased)
//TrkVarFunc(sigqp)
//TrkVarFunc(qp_sigqp)
//TrkVarFunc(chi2)
//TrkVarFunc(ndof)
//TrkVarFunc(qpFraction)
//TrkVarFunc(trkVtxUVDiffPl)
//TrkVarFunc(trkLength)
//TrkVarFunc(planeTrkNu)
//TrkVarFunc(planeTrkNv)
//TrkVarFunc(ntrklike)
//TrkVarFunc(trkphsigcor)
//TrkVarFunc(trkphsigmap)
//TrkVarFunc(trkIdMC)
//TrkVarFunc(trkds)
//TrkVarFunc(trkfitpassSA)
//TrkVarFunc(trkvtxdcoszSA)
//TrkVarFunc(chargeSA)
//TrkVarFunc(qpSA)
//TrkVarFunc(sigqpSA)
//TrkVarFunc(chi2SA)
//TrkVarFunc(ndofSA)
//TrkVarFunc(probSA)
//TrkVarFunc(xTrkVtxSA)
//TrkVarFunc(yTrkVtxSA)
//TrkVarFunc(zTrkVtxSA)
//TrkVarFunc(uTrkVtxSA)
//TrkVarFunc(vTrkVtxSA)
//TrkVarFunc(jitter)
//TrkVarFunc(jPID)
//TrkVarFunc(majC)
//TrkVarFunc(smoothMajC)
//TrkVarFunc(roID)
//TrkVarFunc(knn01TrkActivePlanes)
//TrkVarFunc(knn10TrkMeanPh)
//TrkVarFunc(knn20LowHighPh)
//TrkVarFunc(knn40TrkPhFrac)
//TrkVarFunc(roIDNuMuBar)
//TrkVarFunc(relativeAngle)
//TrkVarFunc(jmID)
//TrkVarFunc(jmTrackPlane)
//TrkVarFunc(jmEndPh)
//TrkVarFunc(jmMeanPh)
//TrkVarFunc(jmScatteringU)
//TrkVarFunc(jmScatteringV)
//TrkVarFunc(jmScatteringUV)
//TrkVarFunc(xTrkVtx)
//TrkVarFunc(yTrkVtx)
//TrkVarFunc(zTrkVtx)
//TrkVarFunc(uTrkVtx)
//TrkVarFunc(vTrkVtx)
//TrkVarFunc(tTrkVtx)
//TrkVarFunc(planeTrkVtx)
//TrkVarFunc(planeTrkBeg)
//TrkVarFunc(planeTrkBegu)
//TrkVarFunc(planeTrkBegv)
//TrkVarFunc(inverseBetaTrk)
//TrkVarFunc(t0Trk)
//TrkVarFunc(chi2TimeTrk)
//TrkVarFunc(ndigitTimeTrk)
//TrkVarFunc(forwardRMSTrk)
//TrkVarFunc(forwardNDOFTrk)
//TrkVarFunc(stripTrkBeg)
//TrkVarFunc(stripTrkBegu)
//TrkVarFunc(stripTrkBegv)
//TrkVarFunc(stripTrkEnd)
//TrkVarFunc(stripTrkEndu)
//TrkVarFunc(stripTrkEndv)
//TrkVarFunc(stripTrkBegIsu)
//TrkVarFunc(stripTrkEndIsu)
//TrkVarFunc(regionTrkVtx)
//TrkVarFunc(edgeRegionTrkVtx)
//TrkVarFunc(edgeRegionTrkEnd)
//TrkVarFunc(phiTrkVtx)
//TrkVarFunc(phiTrkEnd)
//TrkVarFunc(parallelStripTrkVtx)
//TrkVarFunc(parallelStripTrkEnd)
//TrkVarFunc(stripTrkBegPerpFlag)
//TrkVarFunc(stripTrkEndPerpFlag)
//TrkVarFunc(stripHoveNumTrkVtx)
//TrkVarFunc(stripHoveNumTrkEnd)
//TrkVarFunc(xTrkEnd)
//TrkVarFunc(yTrkEnd)
//TrkVarFunc(zTrkEnd)
//TrkVarFunc(uTrkEnd)
//TrkVarFunc(vTrkEnd)
//TrkVarFunc(planeTrkEnd)
//TrkVarFunc(planeTrkEndu)
//TrkVarFunc(planeTrkEndv)
//TrkVarFunc(drTrkFidall)
//TrkVarFunc(dzTrkFidall)
//TrkVarFunc(drTrkFidvtx)
//TrkVarFunc(dzTrkFidvtx)
//TrkVarFunc(drTrkFidend)
//TrkVarFunc(dzTrkFidend)
//TrkVarFunc(traceTrkFidall)
//TrkVarFunc(traceTrkFidvtx)
//TrkVarFunc(traceTrkFidend)
//TrkVarFunc(cosPrTrkVtx)
//TrkVarFunc(mcTrk)
//
//#undef TrkVarFunc
//
//#define ShwVarFunc(var)                           \
//typeof(dummy.var##1)& get_##var(NuEvent& nu,      \
//				int idx)          \
//{                                                 \
//  switch(idx){                                    \
//  case 1: return nu.var##1;                       \
//  case 2: return nu.var##2;                       \
//  case 3: return nu.var##3;                       \
//  case 4: return nu.var##4;                       \
//  case 5: return nu.var##5;                       \
//  default:                                        \
//    assert(0 && "Shower index out of range");     \
//  }                                               \
//}                                                 \
//                                                  \
//typeof(dummy.var##1) get_##var(const NuEvent& nu, \
//			       int idx)           \
//{                                                 \
//  switch(idx){                                    \
//  case 1: return nu.var##1;                       \
//  case 2: return nu.var##2;                       \
//  case 3: return nu.var##3;                       \
//  case 4: return nu.var##4;                       \
//  case 5: return nu.var##5;                       \
//  default:                                        \
//    assert(0 && "Shower index out of range");     \
//  }                                               \
//}                                                 \
//                                                  \
//typeof(dummy.var##1)& get_##var(NuEvent* nu,      \
//				int idx)          \
//{                                                 \
//  switch(idx){                                    \
//  case 1: return nu->var##1;                      \
//  case 2: return nu->var##2;                      \
//  case 3: return nu->var##3;                      \
//  case 4: return nu->var##4;                      \
//  case 5: return nu->var##5;                      \
//  default:                                        \
//    assert(0 && "Shower index out of range");     \
//  }                                               \
//}
//
//ShwVarFunc(shwExists)
//ShwVarFunc(ndigitShw)
//ShwVarFunc(nstripShw)
//ShwVarFunc(nplaneShw)
//ShwVarFunc(shwEnCor)
//ShwVarFunc(shwEnNoCor)
//ShwVarFunc(shwEnLinCCNoCor)
//ShwVarFunc(shwEnLinCCCor)
//ShwVarFunc(shwEnWtCCNoCor)
//ShwVarFunc(shwEnWtCCCor)
//ShwVarFunc(shwEnLinNCNoCor)
//ShwVarFunc(shwEnLinNCCor)
//ShwVarFunc(shwEnWtNCNoCor)
//ShwVarFunc(shwEnWtNCCor)
//ShwVarFunc(shwEnMip)
//ShwVarFunc(planeShwBeg)
//ShwVarFunc(planeShwEnd)
//ShwVarFunc(planeShwMax)
//ShwVarFunc(xShwVtx)
//ShwVarFunc(yShwVtx)
//ShwVarFunc(zShwVtx)
//ShwVarFunc(mcShw)
//
//#undef ShwVarFunc
//
//} // end namespace
} // end minos namespace
#endif //NUEVENT_H

