////////// ///////// ///////// //////// ///////// ///////// ///////// 72
////////////////////////////////////////////////////////////////////////
// Contact: Jeff Hartnell j.j.hartnell@rl.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef NUMCEVENT_H
#define NUMCEVENT_H

#include "TObject.h"
namespace minos {
    class NuMCEvent : public TObject {
    public:
        NuMCEvent();

        //non-const methods
        void Reset();

        //////////////////////////////////////////////////////////////////
        //EVERYTIME A VARIABLE IS ADDED TO THIS CLASS IT MUST BE ADDED TO
        //THE CODE IN NtupleUtils/NuExtraction::NuMCEventFromNuEvent
        //                        (const NuEvent& nu,NuMCEvent& numc) const
        ///////////////////////////////////////////////////////////////////

        ///////////////////////////
        //book keeping  quantities
        ///////////////////////////

        Int_t index;//the number in the ouput tree of NuMCEvent objects
        Int_t entry;//the n'th snarl to be processed

        /////////////////////////////
        //snarl/run based quantities
        /////////////////////////////

        Int_t run;//fHeader.fRun
        Int_t subRun;//fHeader.fSubRun
        Int_t snarl;//actual DAQ snarl number

        // Flux Info
        Int_t runPeriod;
        Float_t intensity;
        Bool_t hornIsReverse;
        Int_t beamType;//Conventions/BeamType.h definition

        Int_t detector;//fHeader.fVldContext.fDetector
        Int_t simFlag;//fHeader.fVldContext.fSimFlag
        Int_t timeSec;//fHeader.fVldContext.fTimeStamp.fSec
        Int_t timeNanoSec;//fHeader.fVldContext.fTimeStamp.fNanoSec
        Double_t timeSeconds;//VldTimeStamp.GetSeconds() sec+nanoSec/1e9
        Double_t trigtime;//evthdr.trigtime

        Int_t releaseType;//Conventions/ReleaseType.h definition
        Int_t anaVersion;//different cuts etc


        //////////////////
        //truth variables
        //////////////////
        Bool_t isInFidVolCCMC;//Is the true vertex in the CC fiducial volume?

        Float_t energyMC;//this is what could be truely reco'd: nuEn*y for NC

        Float_t neuEnMC;//p4neu[3];
        Float_t neuPxMC;//p4neu[0];
        Float_t neuPyMC;//p4neu[1];
        Float_t neuPzMC;//p4neu[2];

        Float_t mu1EnMC;//p4mu1[3];
        Float_t mu1PxMC;//p4mu1[0];
        Float_t mu1PyMC;//p4mu1[1];
        Float_t mu1PzMC;//p4mu1[2];

        Float_t tgtEnMC;//p4tgt[3]
        Float_t tgtPxMC;//p4tgt[0];
        Float_t tgtPyMC;//p4tgt[1];
        Float_t tgtPzMC;//p4tgt[2];

        Int_t zMC;//z;
        Int_t aMC;//a;
        Int_t nucleusMC;//encoding from Mad of the nucleus type
        Int_t initialStateMC;//encoding from Mad of the initial state
        Int_t hadronicFinalStateMC;//encoding from Mad of the hfs
        Int_t numPreInukeFSprotMC;//Number of final state protons (before inuke)
        Int_t numPreInukeFSneutMC;//Number of final state neutrons (before inuke)
        Float_t maxMomPreInukeFSprotMC;//The highest momentum final state proton
        //(before inuke)
        Float_t maxMomPreInukeFSneutMC;//The highest momentum final state neutron
        //(before inuke)

        Float_t yMC;//y value
        Float_t y2MC;//p4shw[3]/(fabs(p4mu1[3])+p4shw[3]);
        Float_t xMC;//x;
        Float_t q2MC;//q2;
        Float_t w2MC;//w2;

        Float_t trkEnMC;//p4mu1[3] - not proper p4: muon energy (+/- !!!)
        Float_t trkEn2MC;//(1-y)*p4neu[3];
        Float_t shwEnMC;//p4shw[3]
        Float_t shwEn2MC;//y*p4neu[3];

        Float_t trkEndEnMC;//NtpMCStdHep.dethit[1].pE, particle En at last scint hit
        Float_t trkStartEnMC; // muon energy at first scintillator hit
        Bool_t trkContainmentMC;// Experimental: MC Containment of primary muon

        Float_t sigma;//mc.sigma=cross-section
        Int_t iaction;//CC=1, NC=0
        Int_t iresonance;//QE=1001, RES=1002, DIS=1003, CPP=1004
        Int_t inu;//>0 particles, <0 anti-particles
        Int_t inunoosc;//id of neutrino at birth
        Int_t itg;//neutrino interaction target

        Float_t vtxxMC;//vtxx: x vtx of neutrino interaction
        Float_t vtxyMC;//vtxy: y vtx of neutrino interaction
        Float_t vtxzMC;//vtxz: z vtx of neutrino interaction
        Float_t vtxuMC;//vtx u, calculated from x&y
        Float_t vtxvMC;//vtx v, calculated from x&y
        Int_t planeTrkVtxMC;//calculated from z
        Float_t rTrkVtxMC;//calculated from x&y

        Int_t mc;//the index of the object to be used

        //http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/v19/output_gnumi.html
        Float_t Npz;
        Float_t NdxdzNea;
        Float_t NdydzNea;
        Float_t NenergyN;
        Float_t NWtNear;
        Float_t NdxdzFar;
        Float_t NdydzFar;
        Float_t NenergyF;
        Float_t NWtFar;
        Int_t Ndecay;
        Float_t Vx;
        Float_t Vy;
        Float_t Vz;
        Float_t pdPx;
        Float_t pdPy;
        Float_t pdPz;
        Float_t ppdxdz;
        Float_t ppdydz;
        Float_t pppz;
        Float_t ppenergy;
        Float_t ppmedium;
        Float_t ppvx;//mc.flux.ppvx = parent production vtx x, e.g. decay pipe
        Float_t ppvy;//mc.flux.ppvy
        Float_t ppvz;//mc.flux.ppvz
        Int_t ptype;//mc.flux.ptype = pion, kaon, muon, etc
        Float_t Necm;
        Float_t Nimpwt;
        Float_t tvx;
        Float_t tvy;
        Float_t tvz;
        Float_t tpx;//mc.flux.tpx = grand^n parent target momentum x (Fluka)
        Float_t tpy;//mc.flux.tpy
        Float_t tpz;//mc.flux.tpz
        Int_t tptype;//mc.flux.tptype = pion, kaon, muon, etc
        Int_t tgen;//mc.flux.tgen


        /////////
        //weights
        /////////
        Int_t reweightVersion; ///< the beam reweighting to use
        Bool_t applyBeamWeight; ///< flag to use beam weight
        Bool_t apply1SigmaWeight; ///< flag to use +1-sigma SKZP error shift
        Bool_t applyDetectorWeight; ///< flag to use detector weight, e.g. xsec
        Bool_t applyGeneratorWeight; ///< flag to use generator weight

        Float_t rw; ///< The final weight applied to the event.  By default this is just beamWeight (SKZP beam weights).
        Float_t fluxErr; ///< The error on the flux from SKZP.
        Float_t rwActual; ///< This is the weight as it leaves NuReco::ApplyReweights.  Any later changes (i.e. systeamtics) only affect rw.
        Float_t generatorWeight; ///< weight factor from generator
        Float_t detectorWeight; ///< SKZP detector weight to use as default
        Float_t anaWeightCC2010; ///< CC 2010 analysis weight

        Float_t trkEnWeight; ///< SKZP weight applied to trkEn (number close to 1).  By default, this is not used.
        Float_t shwEnWeight; ///< SKZP weight applied to shwEn (number close to 1).  By default, this is not used.
        Float_t beamWeight; ///< SKZP weight for beam (e.g. hadron prod.).  By default, this is used.
        Float_t detectorWeightNMB; ///< SKZP detector weight with an antineutrino selection.  By default, this is not used.
        Float_t detectorWeightNM; ///< SKZP detector weight with a neutrino selection.  By default, this is not used.

        //----------------------------

        Float_t InukeNwts;
        Float_t InukePiCExchgP;
        Float_t InukePiCExchgN;
        Float_t InukePiEScatP;
        Float_t InukePiEScatN;
        Float_t InukePiInEScatP;
        Float_t InukePiInEScatN;
        Float_t InukePiAbsorbP;
        Float_t InukePiAbsorbN;
        Float_t InukePi2PiP;
        Float_t InukePi2PiN;
        Float_t InukeNknockP;
        Float_t InukeNknockN;
        Float_t InukeNNPiP;
        Float_t InukeNNPiN;
        Float_t InukeFormTP;
        Float_t InukeFormTN;
        Float_t InukePiXsecP;
        Float_t InukePiXsecN;
        Float_t InukeNXsecP;
        Float_t InukeNXsecN;
        Float_t InukeNucrad;
        Float_t InukeWrad;
        //-----------------------------

        //////////////////////////////////////////////////////////////////
        //EVERYTIME A VARIABLE IS ADDED TO THIS CLASS IT MUST BE ADDED TO
        //THE CODE IN NtupleUtils/NuExtraction::NuMCEventFromNuEvent
        //                        (const NuEvent& nu,NuMCEvent& numc) const
        //AND         NtupleUtils/NuExtraction::NuEventFromNuMCEvent
        //                        (const NuMCEvent& mc,NuEvent& nu) const
        ///////////////////////////////////////////////////////////////////

        ClassDef(NuMCEvent,
        20);
    };
} // end namespace
#endif //NUMCEVENT_H

