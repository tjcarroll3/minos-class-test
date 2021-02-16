////////// ///////// ///////// //////// ///////// ///////// ///////// 72
////////////////////////////////////////////////////////////////////////
// Contact: Jeff Hartnell (j.j.hartnell@rl.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "CovarianceMatrixFit/modules/minos/NuEvent.h"

namespace minos {

//......................................................................

    NuEvent::NuEvent() {
        this->Reset();
    }

//......................................................................

    void NuEvent::Reset() {
        //////////////////////////////////////////////////////////////////
        /// EVERYTIME A TRUTH VARIABLE IS ADDED TO THIS CLASS IT MUST
        /// ALSO BE ADDED TO NuMCEvent
        ///////////////////////////////////////////////////////////////////

        ///See header file for a description of the variables

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
        runType = -1;
        errorCode = -1;
        snarl = -1;
        trigSrc = -1;
        timeFrame = -1;
        remoteSpillType = -1;

        detector = -1;
        simFlag = -1;
        timeSec = -1;
        timeNanoSec = -1;
        timeSeconds = -1;

        trigtime = -1;
        medianTime = -1;
        timeEvtMin = -999999;
        timeEvtMax = -999999;
        crateT0 = -999999.;
        sgate_10mhz = -999999;
        sgate_53mhz = -999999;
        rollover_53mhz = -999999;
        rollover_last_53mhz = -999999;
        timingFiducial = -999999;

        nearestSpillSec = -1;
        nearestSpillNanosec = 0;
        timeToNearestSpill = -999999;

        planeEvtHdrBeg = -1;
        planeEvtHdrEnd = -1;
        snarlPulseHeight = -1;

        //////////////////////////
        // data quality variables
        //////////////////////////
        isGoodDataQuality = false;
        isGoodDataQualityRUN = false;
        isGoodDataQualityCOIL = false;
        isGoodDataQualityHV = false;
        isGoodDataQualityGPS = false;

        numActiveCrates = 0;
        numTimeFrames = 0;
        numGoodSnarls = 0;
        snarlRateMedian = 0.0;
        snarlRateMax = 0.0;

        coilIsOk = true;//false is safer, but old ntuples don't have it
        coilIsReverse = false;//assume forward (unlikely to happen upon reverse)
        coilCurrent = 0;//will never want 0 most likely

        deltaSecToSpillGPS = -999999;
        deltaNanoSecToSpillGPS = 0;
        gpsError = -1;
        gpsSpillType = -1;

        isLI = false;
        litag = 0;
        litime = -1;

        ///////////////////////////
        //reconstructed quantities
        ///////////////////////////
        energy = -1;
        energyCC = -1;
        energyNC = -1;
        energyRM = -1;
        trkEn = -1;
        trkEnRange = -1;
        trkEnCurv = -1;
        shwEn = -1;
        shwEnCC = -1;
        shwEnNC = -1;
        resolution = -1;

        isInFidVolCC = false;
        isGoodTrackReclamation = false;
        isNCClean = false;

        y = -1;
        q2 = -1;
        x = -1;
        w2 = -1;
        dirCosNu = -999;

        //fvnmb=false;//DEPRECATED
        //fvpitt=false;//DEPRECATED
        //fvcc=false;//DEPRECATED

        ////////////////////////////////////////
        //event info extracted from the ntuples
        evt = -1;
        slc = -1;
        nevt = -1;
        ndigitEvt = -1;
        nstripEvt = -1;
        nshw = -1;
        ntrk = -1;
        primshw = -1;
        primtrk = -1;
        rawPhEvt = -1;
        evtphsigcor = -1;
        evtphsigmap = -1;
        planeEvtN = -1;
        planeEvtNu = -1;
        planeEvtNv = -1;

        vtxFitPassFail = false;
        vtxFitTime = -999.9;
        vtxFitTimeError = -999.9;
        vtxFitChi2DoF = -999.9;

        roIDEvt = -999;
        knn01TrkActivePlanesEvt = -999;
        knn10TrkMeanPhEvt = -999;
        knn20LowHighPhEvt = -999;
        knn40TrkPhFracEvt = -999;
        roIDNuMuBarEvt = -999;
        relativeAngleEvt = -999;
        // jasmine's test variables
        jmIDEvt = -1;
        jmTrackPlaneEvt = -999;
        jmMeanPhEvt = -999;
        jmEndPhEvt = -999;
        jmScatteringUEvt = -999;
        jmScatteringVEvt = -999;
        jmScatteringUVEvt = -999;
        jmEventknnIDEvt = -999;
        jmEventknn208Evt = -999;
        jmEventknn207Evt = -999;
        jmEventknn206Evt = -999;
        jmEventknn205Evt = -999;
        jmEventknn204Evt = -999;


        xEvtVtx = -999;
        yEvtVtx = -999;
        zEvtVtx = -999;
        uEvtVtx = -999;
        vEvtVtx = -999;
        tEvtVtx = -999;
        planeEvtVtx = -1;
        planeEvtBeg = -1;
        planeEvtBegu = -1;
        planeEvtBegv = -1;

        xEvtEnd = -999;
        yEvtEnd = -999;
        zEvtEnd = -999;
        uEvtEnd = -999;
        vEvtEnd = -999;
        planeEvtEnd = -1;
        planeEvtEndu = -1;
        planeEvtEndv = -1;

        /////////////////////////////////////////////////////////
        //these are the variables of the "best" track and shower
        trkExists = false;
        trkIndex = -1;
        ndigitTrk = -1;
        nstripTrk = -1;
        trkEnCorRange = -1;
        trkEnCorCurv = -1;
        trkShwEnNear = -1;
        trkShwEnNearDW = -1;
        trkMomentumRange = -1;
        trkResolution = -1;
        containedTrk = 0;
        trkfitpass = -1;
        trkvtxdcosz = -999;
        trkvtxdcosy = -999;
        trknplane = -999;
        charge = 0;
        qp = -999;
        qp_rangebiased = -999;
        sigqp = -1;
        qp_sigqp = -999;
        chi2 = -1;
        ndof = 0;
        qpFraction = -1;
        trkVtxUVDiffPl = -999;
        trkLength = -1;
        planeTrkNu = -1;
        planeTrkNv = -1;
        ntrklike = -1;
        trkphsigcor = -1;
        trkphsigmap = -1;
        trkIdMC = 0;
        trkds = 0;


        trkfitpassSA = -1;
        trkvtxdcoszSA = -999;
        chargeSA = 0;
        qpSA = -999;
        sigqpSA = -1;
        chi2SA = -1;
        ndofSA = 0;
        probSA = -1;
        xTrkVtxSA = -999;
        yTrkVtxSA = -999;
        zTrkVtxSA = -999;
        uTrkVtxSA = -999;
        vTrkVtxSA = -999;

        jitter = -1;
        jPID = -999;
        majC = 0;
        //majCRatio=-999;
        //rms=-1;
        //simpleMajC=-999;
        smoothMajC = -999;
        //sqJitter=-1;
        //totWidth=-999;

        xTrkVtx = -999;
        yTrkVtx = -999;
        zTrkVtx = -999;
        uTrkVtx = -999;
        vTrkVtx = -999;
        tTrkVtx = -999;
        planeTrkVtx = -1;
        planeTrkBeg = -1;
        planeTrkBegu = -1;
        planeTrkBegv = -1;
        inverseBetaTrk = -999;
        t0Trk = -999;
        chi2TimeTrk = -999;
        ndigitTimeTrk = 0;
        forwardRMSTrk = -999;
        forwardNDOFTrk = -999;

        stripTrkBeg = -1;
        stripTrkBegu = -1;
        stripTrkBegv = -1;
        stripTrkEnd = -1;
        stripTrkEndu = -1;
        stripTrkEndv = -1;
        stripTrkBegIsu = false;
        stripTrkEndIsu = false;
        regionTrkVtx = -1;
        edgeRegionTrkVtx = -1;
        edgeRegionTrkEnd = -1;
        phiTrkVtx = -999;
        phiTrkEnd = -999;

        parallelStripTrkVtx = -999;
        parallelStripTrkVtxNoShift = -999;
        parallelStripTrkEnd = -999;
        stripTrkBegPerpFlag = -10;
        stripTrkEndPerpFlag = -10;
        stripHoveNumTrkVtx = -999;
        stripHoveNumTrkVtxNoShift = -999;
        stripHoveNumTrkEnd = -999;


        xTrkEnd = -999;
        yTrkEnd = -999;
        zTrkEnd = -999;
        uTrkEnd = -999;
        vTrkEnd = -999;
        planeTrkEnd = -1;
        planeTrkEndu = -1;
        planeTrkEndv = -1;

        drTrkFidall = -999;
        dzTrkFidall = -999;
        drTrkFidvtx = -999;
        dzTrkFidvtx = -999;
        drTrkFidend = -999;
        dzTrkFidend = -999;
        traceTrkFidall = -999;
        traceTrkFidvtx = -999;
        traceTrkFidend = -999;

        cosPrTrkVtx = -999;

        //shw variables
        shwExists = false;
        ndigitShw = -1;
        nstripShw = -1;
        nplaneShw = -1;
        shwEnCor = -1;
        shwEnNoCor = -1;
        shwEnMip = -1;
        shwEnLinCCNoCor = -1;
        shwEnLinCCCor = -1;
        shwEnWtCCNoCor = -1;
        shwEnWtCCCor = -1;
        shwEnLinNCNoCor = -1;
        shwEnLinNCCor = -1;
        shwEnWtNCNoCor = -1;
        shwEnWtNCCor = -1;
        shwResolution = -1;
        shwEnkNN = -1;
        shwEnkNNNoCor = -1;
        shwEnReskNN = -1;

        planeShwBeg = -1;
        planeShwEnd = -1;
        planeShwMax = -1;
        xShwVtx = -999;
        yShwVtx = -999;
        zShwVtx = -999;

        // MINOS+ shower variables
        // junting@physics.utexas.edu
        shwEnkNNPlus = -1;
        shwEnkNNPlusNoCor = -1;
        shwRmsT = -1;
        shwPctWidth = -1;
        shwPctLength = -1;
        shwSumHitLength = -1;
        shwSumHitWidth = -1;
        avgShwHitMip = -1;
        maxShwPlaneMip = -1;
        maxShwHitMip = -1;
        shwRmsTAllShws = -1;
        shwPctWidthAllShws = -1;
        shwPctLengthAllShws = -1;
        avgShwHitMipAllShws = -1;
        shwSumHitLengthAllShws = -1;
        shwSumHitWidthAllShws = -1;
        maxShwPlaneMipAllShws = -1;
        maxShwHitMipAllShws = -1;
        nplaneAllShws = -1;
        planeAllShwMax = -1;

        // Initialize tracks and showers
/*
        for (int trkIdx = 1; trkIdx <= 3; ++trkIdx) {
            get_trkExists(this, trkIdx) = false;
            get_trkIndex(this, trkIdx) = -1;
            get_ndigitTrk(this, trkIdx) = -1;
            get_nstripTrk(this, trkIdx) = -1;
            get_trkEnCorRange(this, trkIdx) = -1;
            get_trkEnCorCurv(this, trkIdx) = -1;
            get_trkShwEnNear(this, trkIdx) = -1;
            get_trkShwEnNearDW(this, trkIdx) = -1;
            get_trkMomentumRange(this, trkIdx) = -1;
            get_containedTrk(this, trkIdx) = 0;
            get_trkfitpass(this, trkIdx) = -1;
            get_trkvtxdcosz(this, trkIdx) = -999;
            get_trkvtxdcosy(this, trkIdx) = -999;
            get_trknplane(this, trkIdx) = -999;
            get_charge(this, trkIdx) = 0;
            get_qp(this, trkIdx) = -999;
            get_qp_rangebiased(this, trkIdx) = -999;
            get_sigqp(this, trkIdx) = -1;
            get_qp_sigqp(this, trkIdx) = -999;
            get_chi2(this, trkIdx) = -1;
            get_ndof(this, trkIdx) = 0;
            get_qpFraction(this, trkIdx) = -1;
            get_trkVtxUVDiffPl(this, trkIdx) = -999;
            get_trkLength(this, trkIdx) = -1;
            get_planeTrkNu(this, trkIdx) = -1;
            get_planeTrkNv(this, trkIdx) = -1;
            get_ntrklike(this, trkIdx) = -1;
            get_trkphsigcor(this, trkIdx) = -1;
            get_trkphsigmap(this, trkIdx) = -1;
            get_trkIdMC(this, trkIdx) = 0;
            get_trkds(this, trkIdx) = -1;

            get_trkfitpassSA(this, trkIdx) = -1;
            get_trkvtxdcoszSA(this, trkIdx) = -999;
            get_chargeSA(this, trkIdx) = 0;
            get_qpSA(this, trkIdx) = -999;
            get_sigqpSA(this, trkIdx) = -1;
            get_chi2SA(this, trkIdx) = -1;
            get_ndofSA(this, trkIdx) = 0;
            get_probSA(this, trkIdx) = -1;
            get_xTrkVtxSA(this, trkIdx) = -999;
            get_yTrkVtxSA(this, trkIdx) = -999;
            get_zTrkVtxSA(this, trkIdx) = -999;
            get_uTrkVtxSA(this, trkIdx) = -999;
            get_vTrkVtxSA(this, trkIdx) = -999;

            get_jitter(this, trkIdx) = -1;
            get_jPID(this, trkIdx) = -999;
            get_majC(this, trkIdx) = 0;
            //    get_majCRatio(this, trkIdx) = -999;
            //    get_rms(this, trkIdx) = -1;
            //    get_simpleMajC(this, trkIdx) = -999;
            get_smoothMajC(this, trkIdx) = -999;
            //    get_sqJitter(this, trkIdx) = -1;
            //    get_totWidth(this, trkIdx) = -999;

            get_roID(this, trkIdx) = -999;
            get_knn01TrkActivePlanes(this, trkIdx) = -999;
            get_knn10TrkMeanPh(this, trkIdx) = -999;
            get_knn20LowHighPh(this, trkIdx) = -999;
            get_knn40TrkPhFrac(this, trkIdx) = -999;
            get_roIDNuMuBar(this, trkIdx) = -999;
            get_relativeAngle(this, trkIdx) = -999;

            // jasmine's test variables
            get_jmID(this, trkIdx) = -1;
            get_jmTrackPlane(this, trkIdx) = -999;
            get_jmMeanPh(this, trkIdx) = -999;
            get_jmEndPh(this, trkIdx) = -999;
            get_jmScatteringU(this, trkIdx) = -999;
            get_jmScatteringV(this, trkIdx) = -999;
            get_jmScatteringUV(this, trkIdx) = -999;

            get_xTrkVtx(this, trkIdx) = -999;
            get_yTrkVtx(this, trkIdx) = -999;
            get_zTrkVtx(this, trkIdx) = -999;
            get_uTrkVtx(this, trkIdx) = -999;
            get_vTrkVtx(this, trkIdx) = -999;
            get_tTrkVtx(this, trkIdx) = -999;
            get_planeTrkVtx(this, trkIdx) = -1;
            get_planeTrkBeg(this, trkIdx) = -1;
            get_planeTrkBegu(this, trkIdx) = -1;
            get_planeTrkBegv(this, trkIdx) = -1;
            get_inverseBetaTrk(this, trkIdx) = -999;
            get_t0Trk(this, trkIdx) = -999;
            get_chi2TimeTrk(this, trkIdx) = -999;
            get_ndigitTimeTrk(this, trkIdx) = 0;
            get_forwardRMSTrk(this, trkIdx) = -999;
            get_forwardNDOFTrk(this, trkIdx) = -999;
            get_stripTrkBeg(this, trkIdx) = -1;
            get_stripTrkBegu(this, trkIdx) = -1;
            get_stripTrkBegv(this, trkIdx) = -1;
            get_stripTrkEnd(this, trkIdx) = -1;
            get_stripTrkEndu(this, trkIdx) = -1;
            get_stripTrkEndv(this, trkIdx) = -1;
            get_stripTrkBegIsu(this, trkIdx) = false;
            get_stripTrkEndIsu(this, trkIdx) = false;
            get_regionTrkVtx(this, trkIdx) = -1;
            get_edgeRegionTrkVtx(this, trkIdx) = -1;
            get_edgeRegionTrkEnd(this, trkIdx) = -1;
            get_phiTrkVtx(this, trkIdx) = -999;
            get_phiTrkEnd(this, trkIdx) = -999;

            get_parallelStripTrkVtx(this, trkIdx) = -999;
            get_parallelStripTrkEnd(this, trkIdx) = -999;
            get_stripTrkBegPerpFlag(this, trkIdx) = -10;
            get_stripTrkEndPerpFlag(this, trkIdx) = -10;
            get_stripHoveNumTrkVtx(this, trkIdx) = -999;
            get_stripHoveNumTrkEnd(this, trkIdx) = -999;

            get_xTrkEnd(this, trkIdx) = -999;
            get_yTrkEnd(this, trkIdx) = -999;
            get_zTrkEnd(this, trkIdx) = -999;
            get_uTrkEnd(this, trkIdx) = -999;
            get_vTrkEnd(this, trkIdx) = -999;
            get_planeTrkEnd(this, trkIdx) = -1;
            get_planeTrkEndu(this, trkIdx) = -1;
            get_planeTrkEndv(this, trkIdx) = -1;

            get_drTrkFidall(this, trkIdx) = -999;
            get_dzTrkFidall(this, trkIdx) = -999;
            get_drTrkFidvtx(this, trkIdx) = -999;
            get_dzTrkFidvtx(this, trkIdx) = -999;
            get_drTrkFidend(this, trkIdx) = -999;
            get_dzTrkFidend(this, trkIdx) = -999;
            get_traceTrkFidall(this, trkIdx) = -999;
            get_traceTrkFidvtx(this, trkIdx) = -999;
            get_traceTrkFidend(this, trkIdx) = -999;

            get_cosPrTrkVtx(this, trkIdx) = -999;
        } // end for trkIdx

        for (int shwIdx = 1; shwIdx <= 5; ++shwIdx) {
            get_shwExists(this, shwIdx) = false;
            get_ndigitShw(this, shwIdx) = -1;
            get_nstripShw(this, shwIdx) = -1;
            get_nplaneShw(this, shwIdx) = -1;
            get_shwEnCor(this, shwIdx) = -1;
            get_shwEnNoCor(this, shwIdx) = -1;
            get_shwEnLinCCNoCor(this, shwIdx) = -1;
            get_shwEnLinCCCor(this, shwIdx) = -1;
            get_shwEnWtCCNoCor(this, shwIdx) = -1;
            get_shwEnWtCCCor(this, shwIdx) = -1;
            get_shwEnLinNCNoCor(this, shwIdx) = -1;
            get_shwEnLinNCCor(this, shwIdx) = -1;
            get_shwEnWtNCNoCor(this, shwIdx) = -1;
            get_shwEnWtNCCor(this, shwIdx) = -1;
            get_shwEnMip(this, shwIdx) = -1;
            get_planeShwBeg(this, shwIdx) = -1;
            get_planeShwEnd(this, shwIdx) = -1;
            get_planeShwMax(this, shwIdx) = -1;
            get_xShwVtx(this, shwIdx) = -999;
            get_yShwVtx(this, shwIdx) = -999;
            get_zShwVtx(this, shwIdx) = -999;
        } // end for shwIdx
*/

        ////////////////////////
        //other info calculated
        rEvtVtx = -999;
        rEvtEnd = -999;
        distToEdgeEvtVtx = -999;
        evtVtxUVDiffPl = -999;

        rTrkVtx = -999;
        rTrkEnd = -999;
        sigqp_qp = -999;
        chi2PerNdof = 999;
        prob = -1;

        containmentFlag = -1;
        containmentFlagCC0093Std = -1;
        containmentFlagCC0250Std = -1;
        containmentFlagPitt = -1;
        usedRange = 0;
        usedCurv = 0;


        /////////
        //weights
        /////////
        rw = 1;//default of no reweight
        fluxErr = -1;
        rwActual = 1;//default of no reweight
        generatorWeight = 1;//no weight
        detectorWeight = 1;//no weight
        coilCorrWeight = 1;//no weight
        rockxsecCorrWeight = 1;//no weight
        anaWeightCC2010 = 1; // Weight currently unknown

        trkEnWeight = 1;//no weight
        shwEnWeight = 1;//no weight
        beamWeight = 1;//no weight
        detectorWeightNMB = 1;//no weight
        detectorWeightNM = 1;//no weight

        //energies with and without weights
        energyRw = -1;
        energyNoRw = -1;
        trkEnRw = -1;
        trkEnNoRw = -1;
        shwEnRw = -1;
        shwEnNoRw = -1;

        //pids
        dpID = -999;
        abID = -999;
        roID = -999;
        knn01TrkActivePlanes = -999;
        knn10TrkMeanPh = -999;
        knn20LowHighPh = -999;
        knn40TrkPhFrac = -999;
        roIDNuMuBar = -999;
        relativeAngle = -999;
        poID = -999;
        poIDKin = -999;
        roIDPlus = -999;  // roID for MINOS+
        // junting@physics.utexas.edu
        // jasmine ID
        jmID = -1;
        jmTrackPlane = -999;
        jmMeanPh = -999;
        jmEndPh = -999;
        jmScatteringU = -999;
        jmScatteringV = -999;
        jmScatteringUV = -999;
        jmEventknnID = -999;
        jmEventknn208 = -999;
        jmEventknn207 = -999;
        jmEventknn206 = -999;
        jmEventknn205 = -999;
        jmEventknn204 = -999;


        ////////////////////////////
        // NC variables
        ////////////////////////////
        // event
        closeTimeDeltaZ = -1;
        edgeActivityStrips = -1;
        edgeActivityPH = -1;
        oppEdgeStrips = -1;
        oppEdgePH = -1;
        vtxMetersToCoilEvt = -1;
        vtxMetersToCloseEdgeEvt = -1;
        minTimeSeparation = -1;
        slicePHFraction = -1;
        maxConsecutivePlanes = -1;

        // shower
        transverseRMSU = -1;
        transverseRMSV = -1;
        // track
        dtdz = -1;
        endMetersToCloseEdge = -1;
        vtxMetersToCloseEdgeTrk = -1;
        vtxMetersToCoilTrk = -1;
        traceEndZ = -1;

        //beam variables
        pot = -1;
        potDB = -1;
        potSinceLastEvt = 0;    //set to zero in case it's used to count
        potSinceLastEvtGood = 0;//rather than being set to a value
        potSinceLastEvtBad = 0;
        potSinceLastEvtDB = 0;
        potSinceLastEvtGoodDB = 0;
        potSinceLastEvtBadDB = 0;
        nBatch = 0;
        potListSinceLastEvtGood.clear();

        runPeriod = -1;
        hornIsReverse = false;
        beamTypeDB = 0;
        beamType = 0;
        intensity = -1;
        hornCur = -999999;//might want 0 horn current
        goodBeam = false;
        goodBeamSntp = false;

        //////////////////
        //truth variables
        //////////////////

        //////////////////////////////////////////////////////////////////
        //EVERYTIME A TRUTH VARIABLE IS ADDED TO THIS CLASS IT MUST
        //ALSO BE ADDED TO NuMCEvent
        ///////////////////////////////////////////////////////////////////

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
        mcTrk = -1;
        mcShw = -1;
        mcEvt = -1;

        mcTrk1 = -1;
        mcTrk2 = -1;
        mcTrk3 = -1;

        mcShw1 = -1;
        mcShw2 = -1;
        mcShw3 = -1;
        mcShw4 = -1;
        mcShw5 = -1;

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

        InukeNwts = -999;
        InukePiCExchgP = -999;   //0
        InukePiCExchgN = -999;   //1
        InukePiEScatP = -999;  //2
        InukePiEScatN = -999;  //3
        InukePiInEScatP = -999;  //4
        InukePiInEScatN = -999;  //5
        InukePiAbsorbP = -999;   //6
        InukePiAbsorbN = -999;   //7
        InukePi2PiP = -999;      //8
        InukePi2PiN = -999;      //9
        InukeNknockP = -999;     //10
        InukeNknockN = -999;     //11
        InukeNNPiP = -999;       //12
        InukeNNPiN = -999;       //13
        InukeFormTP = -999;      //14
        InukeFormTN = -999;      //15
        InukePiXsecP = -999;     //16
        InukePiXsecN = -999;     //17
        InukeNXsecP = -999;      //18
        InukeNXsecN = -999;      //19
        InukeNucrad = -999;
        InukeWrad = -999;

        //////////////////////////////////////////////////////////////////
        //EVERYTIME A TRUTH VARIABLE IS ADDED TO THIS CLASS IT MUST
        //ALSO BE ADDED TO NuMCEvent
        ///////////////////////////////////////////////////////////////////

        ////////////////////////////
        //program control variables
        ////////////////////////////
        anaVersion = 0;
        releaseType = -1;//the value of Conventions/ReleaseType::kUnknown
        recoVersion = -1;//the value of Conventions/ReleaseType::kUnknown
        mcVersion = -1;//the value of Conventions/ReleaseType::kUnknown
        reweightVersion = 0;
        //beamType=0;

        useGeneratorReweight = 1;//default is to do the generatorReweighting
        sGeneratorConfigName = "Unknown";
        generatorConfigNo = -1;

        //set all these to true by default, so they are used and cut on
        useDBForDataQuality = true;
        useDBForSpillTiming = true;
        useDBForBeamInfo = true;

        cutOnDataQuality = true;
        cutOnSpillTiming = true;
        cutOnBeamInfo = true;

        applyEnergyShifts = false;
        applyBeamWeight = true;
        apply1SigmaWeight = false;
        applyDetectorWeight = false;
        applyGeneratorWeight = false;

        calcMajCurv = true;
        calcRoID = true;
        calcJmID = true;
    }

//......................................................................
} // end namespace
