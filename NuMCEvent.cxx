////////// ///////// ///////// //////// ///////// ///////// ///////// 72
////////////////////////////////////////////////////////////////////////
// Contact: Jeff Hartnell (j.j.hartnell@rl.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "CovarianceMatrixFit/modules/minos/NuMCEvent.h"

//......................................................................
namespace minos {
    NuMCEvent::NuMCEvent() {
        this->Reset();
    }

//......................................................................

    void NuMCEvent::Reset() {
        //////////////////////////////////////////////////////////////////
        //EVERYTIME A VARIABLE IS ADDED TO THIS CLASS IT MUST BE ADDED TO
        //THE CODE IN NtupleUtils/NuExtraction::NuMCEventFromNuEvent
        //                        (const NuEvent& nu,NuMCEvent& numc) const
        ///////////////////////////////////////////////////////////////////

        ///////////////////////////
        //book keeping  quantities
        ///////////////////////////
        index = -1;
        entry = -1;

        /////////////////////////////
        //snarl/run based quantities
        /////////////////////////////
        run = -1;
        subRun = -1;
        snarl = -1;

        // Flux Type
        runPeriod = -1;
        intensity = -1;
        hornIsReverse = false;
        beamType = 0;

        detector = -1;
        simFlag = -1;
        timeSec = -1;
        timeNanoSec = -1;
        timeSeconds = -1;
        trigtime = -1;

        releaseType = -1;
        anaVersion = 0;

        //////////////////
        //truth variables
        //////////////////
        isInFidVolCCMC = false;

        energyMC = -1;

        neuEnMC = -1;
        neuPxMC = -1;
        neuPyMC = -1;
        neuPzMC = -1;

        mu1EnMC = -1;
        mu1PxMC = -1;
        mu1PyMC = -1;
        mu1PzMC = -1;

        tgtEnMC = -1;
        tgtPxMC = -1;
        tgtPyMC = -1;
        tgtPzMC = -1;

        zMC = -1;
        aMC = -1;
        nucleusMC = -1;
        initialStateMC = -1;
        hadronicFinalStateMC = -1;
        numPreInukeFSprotMC = -1;
        numPreInukeFSneutMC = -1;
        maxMomPreInukeFSprotMC = -1;
        maxMomPreInukeFSneutMC = -1;

        yMC = -1;
        y2MC = -1;
        xMC = -1;
        q2MC = -1;
        w2MC = -1;

        trkEnMC = -1;
        trkEn2MC = -1;
        shwEnMC = -1;
        shwEn2MC = -1;

        trkEndEnMC = -1;
        trkStartEnMC = -1;
        trkContainmentMC = false;

        sigma = 999999;
        iaction = -1;
        iresonance = -1;
        inu = 0;
        inunoosc = 0;
        itg = 0;

        vtxxMC = -999;
        vtxyMC = -999;
        vtxzMC = -999;
        vtxuMC = -999;
        vtxvMC = -999;
        planeTrkVtxMC = -999;
        rTrkVtxMC = -999;

        mc = -1;

        Npz = -1;
        NdxdzNea = -1;
        NdydzNea = -1;
        NenergyN = -1;
        NWtNear = -1;
        NdxdzFar = -1;
        NdydzFar = -1;
        NenergyF = -1;
        NWtFar = -1;
        Ndecay = -1;
        Vx = -1;
        Vy = -1;
        Vz = -1;
        pdPx = -1;
        pdPy = -1;
        pdPz = -1;
        ppdxdz = -1;
        ppdydz = -1;
        pppz = -1;
        ppenergy = -1;
        ppmedium = -1;
        ppvx = -1;
        ppvy = -1;
        ppvz = -1;
        ptype = -1;
        Necm = -1;
        Nimpwt = -1;
        tvx = -1;
        tvy = -1;
        tvz = -1;
        tpx = -1;
        tpy = -1;
        tpz = -1;
        tptype = -1;
        tgen = -1;

        reweightVersion = 0;
        applyBeamWeight = true;
        apply1SigmaWeight = false;
        applyDetectorWeight = false;
        applyGeneratorWeight = false;

        rw = 1;//default of no reweight
        fluxErr = -1;
        rwActual = 1;//default of no reweight
        generatorWeight = 1;//no weight
        detectorWeight = 1;//no weight
        anaWeightCC2010 = 1; // unknown as yet

        trkEnWeight = 1;//no weight
        shwEnWeight = 1;//no weight
        beamWeight = 1;//no weight
        detectorWeightNMB = 1;//no weight
        detectorWeightNM = 1;//no weight

        InukeNwts = -999;
        InukePiCExchgP = -999;
        InukePiCExchgN = -999;
        InukePiEScatP = -999;
        InukePiEScatN = -999;
        InukePiInEScatP = -999;
        InukePiInEScatN = -999;
        InukePiAbsorbP = -999;
        InukePiAbsorbN = -999;
        InukePi2PiP = -999;
        InukePi2PiN = -999;
        InukeNknockP = -999;
        InukeNknockN = -999;
        InukeNNPiP = -999;
        InukeNNPiN = -999;
        InukeFormTP = -999;
        InukeFormTN = -999;
        InukePiXsecP = -999;
        InukePiXsecN = -999;
        InukeNXsecP = -999;
        InukeNXsecN = -999;
        InukeNucrad = -999;
        InukeWrad = -999;

        //////////////////////////////////////////////////////////////////
        //EVERYTIME A VARIABLE IS ADDED TO THIS CLASS IT MUST BE ADDED TO
        //THE CODE IN NtupleUtils/NuExtraction::NuMCEventFromNuEvent
        //                        (const NuEvent& nu,NuMCEvent& numc) const
        ///////////////////////////////////////////////////////////////////
    }

//......................................................................
} // end namespace