############################################################################
# This code is an adaptation of SeitzModel.c in python. It needs the REFPROP
# wrapper on github (ctREFPROP) and the REFPROP library (.dll or .so)
# as well as the fluids files
#
# By: Mathieu Laurin
#
# v1.0 Initial code 24/09/19 ML
# v1.1 P and T can now be scalar or arrays 27/09/19 ML
# v1.2 Mixtures can now be defined 12/05/20 AER
#
# INPUTS:
#               PIn   -   Pressure (psia) of the superheated liquid.  Have to be a
#                          list or a scalar. If P is also a list of many values, they
#                          must have the same number of elements.
#               TIn   -   Temperature (C) of the superheated liquid. Hate to be a
#                          list or a scalar. If T is also a list of many values, they
#                          must have the same number of elements.
#           Fluid   -   Name of the fluid or array of fluid names.  In the
#                        case of an array, input a string with elements
#                        delimited by vertical bars '|'. Any names used
#                        should match the stem of the .FLD file (case
#                        insensitive) located in the FLUIDS folder.
#                        Default is 'C3F8'.
#
# OUTPUTS:  The output contains the fields described below.
#           Except for the 'units' field, every field has the same
#           dimensions as P, or if P is scalar and T is not, the same
#           dimensions as T.
#
#   Thermodynamic Fields
#
#               P   -   (psia) Pressure of the superheated liquid
#               T   -   (C) Temperature of the superheated liquid
#               Q   -   (keV) Minimum heat input to create critical bubble - allowing
#                       discontinuity in chemical potential at bubble surface.
#            Q_ER   -   (keV) Minimum heat input to create critical bubble assuming mass
#                       transport across the bubble surface.
#              Rc   -   (nm) Radius of critical bubble
#            Pvap   -   (psia) Vapor pressure of saturated fluid at
#                         temperature T
#            Pbub   -   (psia) Pressure inside hot-spike bubble (vapor pressure of
#                         superheated fluid at T and P, and liquid mixture composition)
#          DeltaH   -   (J/kg) Change in specific enthalpy between superheated liquid and
#                         bubble vapor states assuming the same composition in each
#          DeltaS   -   (J/kg-K) Change in specific entropy between
#                         superheated liquid and bubble vapor states
#                         (should equal DeltaH/T )
#           Gibbs   -   (J/kg) Specific Gibbs free energy of superheated
#                         liquid (see NIST documentation for reference
#                         point) -- should equal Gibbs free energy of
#                         bubble vapor
#           Rho_l   -   (g/cc) Density of superheated liquid
#           Rho_b   -   (g/cc) Density of bubble vapor
#           Sigma   -   (N/m) Surface tension at temperature T
#        dSigmadT   -   (N/m-K) Derivative of surface tension versus
#                         temperature at temperature T
#           ERrej   -   (keV^-1) Baxter-model electron recoil rejection factor (placeholder)
#           Pb_ER   -   (psia) Pressure inside adiabatically bubble (vapor pressure
#                        of superheated fluid at T and P, and equilibrium mixture
#                        composition) (placeholder)
#           Rc_ER   -   (nm) Radius of critical ER bubble (= Rc for pure fluids) (placeholder)
#        Rho_b_ER   -   (g/cc) Density of ER bubble (= Rho_b for pure fluids) (placeholder)
#         Xbub_ER   -   (mol/mol) Critical bubble molar gas composition ratio.
#            Xbub   -   (mol/mol) Critical bubble molar gas composition ratio.  Should be same
#                       as the liquid composition ratio.
#       DeltaH_ER   -   (J/kg) Change in specific enthalpy between superheated liquid and bubble
#                       vapor states at the equilibrium composition of the critical bubble.
#            Eion   -   (keV) Energy for bubble formation by cavitation.
#            ER_x   -   (GeV cm^2/g) Baxter parameter for ER rejection.
#
#   Consisitency Checks
#
#           P_err   -   (psia) Preadback - P, where Preadback is the
#                         pressure corresponding to T and Rho_l
#       Rho_b_err   -   (g/cc) Size of the last step in the iterative
#                         process to find Rho_b
#    Rho_b_ER_err   -   (g/cc) Size of the last step in the iterative
#                         process to find Rho_b_ER
#           G_err   -   (J/kg) Gibbs free energy of the superheated liquid
#                         minus that of the bubble vapor.
#       Sigma_err   -   (N/m) Difference between surface tension given by
#                         the SURTEN and SURFT subroutines (should be 0)
#    dSigmadT_err   -   (N/m) Difference between surface tension derivative
#                         given by the SURTEN and SURFT subroutines (should
#                         be 0)
#          err_ID   -   Error code given by REFPROP code.  Integer part is
#                         the error code, fractional part identifies which
#                         REFPROP call produced the error
# interp_iterations -   Number of iterations used to determind Rho_b
#                         (matching the bubble and liquid Gibbs energies)
#           units   -   Struct containing all of the above fields, each
#                         containing a string giving the units of that
#                         field.
############################################################################

#System import
import os
import math
from collections import namedtuple
import numpy as np

# REFPROP wrapper
import ctREFPROP.ctREFPROP as ct

class SeitzModeloutput_tuple(namedtuple('SeitzModeloutput', ["P", "T", "Q", "Q_ER", "Rc", "Pvap", "Pbub", "DeltaH", "DeltaS", "Gibbs", "Rho_l", "Rho_b", "Sigma", "dSigmadT", "P_err", "Rho_b_err", "G_err", "Sigma_err", "dSigmadT_err", "errID", "interp_iterations", "Lt", "DeltaE_spike", "DeltaW_spike", "Tspike", "Pspike", "Q0_h", "Q1_h", "Q0_s", "Q1_s", "Q0_ds", "Q1_ds", "Q0_w", "Q1_w", "Q0", "Q1", "W0_ext", "W1_ext", "Units", "ERrej", "Pb_ER", "Rc_ER", "Rho_b_ER", "Xbub_ER", "Xbub", "DeltaH_ER", "Rho_b_ER_err", "Q0_h_ER", "Q0_ER", "Q1_ER", "Eion", "ER_x"])):
    def __repr__(self):
        string = ""
        for prop in self._fields:
            string += "\t{0:>20s} : {1:<20s}\n".format(prop, str(getattr(self, prop)))
        return string

def SeitzModel(PIn, TIn, Fluid = "R218", Mix_ratio = [1]):

    # CONSTANTS! must not be overwritten!
    C_maxtrials = 1000
    absolute_zero_C = -273.15
    kPa_per_psi = 6.89475729317
    psi_per_kPa = 1 / kPa_per_psi
    keV_per_J = 6.24150636309e15
    Na = 6.0221415E23

    # Units of all the variables
    Fields = ["P", "T", "Q", "Q_ER", "Rc", "Pvap", "Pbub", "DeltaH", "DeltaS", "Gibbs", "Rho_l", "Rho_b", "Sigma", "dSigmadT", "P_err", "Rho_b_err", "G_err", "Sigma_err", "dSigmadT_err", "errID", "interp_iterations", "Lt", "DeltaE_spike", "DeltaW_spike", "Tspike", "Pspike", "Q0_h", "Q1_h", "Q0_s", "Q1_s", "Q0_ds", "Q1_ds", "Q0_w", "Q1_w", "Q0", "Q1", "W0_ext", "W1_ext", "ERrej", "Pb_ER", "Rc_ER", "Rho_b_ER", "Xbub_ER", "Xbub", "DeltaH_ER", "Rho_b_ER_err", "Q0_h_ER", "Q0_ER", "Q1_ER", "Eion", "ER_x"]
    Units = dict(zip(Fields, ("psia", "C", "keV", "keV", "nm", "psia", "psia", "J/kg", "J/kg-K", "J/kg", "g/cc", "g/cc", "N/m", "N/m-K", "psia", "g/cc", "J/kg", "N/m", "N/m-K", "ierr.SeitzStep", "count", "nm", "J/kg", "J/kg", "C", "psia", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "1/keV", "psia", "nm", "g/cc", "mol/mol", "mol/mol", "J/kg", "g/cc", "keV", "keV", "keV", "keV", "MeV cm^2/g")))

    # Open the library in the same folder as this module
    RPPath = os.path.dirname(os.path.abspath(__file__))
    RP = ct.REFPROPFunctionLibrary(RPPath)
    RP.SETPATHdll(RPPath)

    # Check if the inputs are the right type
    if not isinstance(PIn, (list, int, float, np.ndarray)) or not isinstance(TIn, (list, int, float, np.ndarray)) or not isinstance(Fluid, str):
        print("Incorrect type passed as argument.")
        return

    # Check if PIn contains something
    if isinstance(PIn, (list, np.ndarray)):
        if len(PIn) == 0:
            print("Empty pressure list")
            return

    # Check if TIn contains something
    if isinstance(TIn, (list, np.ndarray)):
        if len(TIn) == 0:
            print("Empty temperature list")
            return

    # Check if PIn is the same size as TIn if both are multivalued arrays
    if isinstance(PIn, (list, np.ndarray)) and isinstance(TIn, (list, np.ndarray)):
        if  len(PIn) > 1 and len(TIn) > 1 and len(PIn) != len(TIn):
            print("Pressure and temperature lists have different lengths")
            return

    # Make sure Fluid is uppercase and patch for C3F8 which doesn't exist anymore
    Fluid = Fluid.upper()
    Fluid = Fluid.replace('C3F8','R218')

    # Check if each element of Fluid exists
    FluidPath = os.path.join(RPPath, "FLUIDS")
    for subFluid in Fluid.split('|'):
        if not os.path.isfile(os.path.join(FluidPath, subFluid + ".FLD")):
            print("Fluid %s does not exist" % subFluid)
            return

    # Open Fluid
    ierr, herr = RP.SETUPdll(len(Fluid.split('|')), Fluid, os.path.join(FluidPath, "HMX.BNC"), "DEF")
    if ierr != 0:
    # Print and continue if 'mixture parameters not found' error occurs.
        if ierr == 4294967179: # warning code -117 in uint
            print(herr)
        else:
            return ierr, herr

    # Get Fluid properties
    #wmm = RP.WMOLdll(Mix_ratio)
    #one_over_wmm = 1 / wmm

    # Transform scalars into arrays and find n_PT
    if isinstance(PIn, (int, float)) and isinstance(TIn, (int, float)):
        n_PT = 1
        T = [TIn]
        P = [PIn]
    elif isinstance(PIn, (int, float)):
        n_PT = len(TIn)
        P = [PIn] * n_PT
        T = TIn
    elif isinstance(TIn, (int, float)):
        n_PT = len(PIn)
        P = PIn
        T = [TIn] * n_PT
    else:
        if len(PIn) > len(TIn):
            n_PT = len(PIn)
            P = PIn
            T = [TIn[0]] * n_PT
        elif len(TIn) > len(PIn):
            n_PT = len(TIn)
            P = [PIn[0]] * n_PT
            T = TIn
        else:
            n_PT = len(PIn)
            P = PIn
            T = TIn

    #Allocate array space
    Q = [None] * n_PT
    Q_ER = [None] * n_PT
    Rc = [None] * n_PT
    Pvap = [None] * n_PT
    Pbub = [None] * n_PT
    DeltaH = [None] * n_PT
    DeltaS = [None] * n_PT
    Gibbs = [None] * n_PT
    Rho_l = [None] * n_PT
    Rho_b = [None] * n_PT
    Sigma = [None] * n_PT
    dSigmadT = [None] * n_PT
    P_err = [None] * n_PT
    Rho_b_err = [None] * n_PT
    G_err = [None] * n_PT
    Sigma_err = [None] * n_PT
    dSigmadT_err = [None] * n_PT
    errID = [None] * n_PT
    interp_iterations = [None] * n_PT
    Lt = [None] * n_PT
    DeltaE_spike = [None] * n_PT
    DeltaW_spike = [None] * n_PT
    Tspike = [None] * n_PT
    Pspike = [None] * n_PT
    Q0_h = [None] * n_PT
    Q1_h = [None] * n_PT
    Q0_s = [None] * n_PT
    Q1_s = [None] * n_PT
    Q0_ds = [None] * n_PT
    Q1_ds = [None] * n_PT
    Q0_w = [None] * n_PT
    Q1_w = [None] * n_PT
    Q0 = [None] * n_PT
    Q1 = [None] * n_PT
    W0_ext = [None] * n_PT
    W1_ext = [None] * n_PT
    ERrej = [None] * n_PT
    Pb_ER = [None] * n_PT
    Rc_ER = [None] * n_PT
    Rho_b_ER = [None] * n_PT
    Xbub_ER = [None] * n_PT
    Xbub = [None] * n_PT
    DeltaH_ER = [None] * n_PT
    Rho_b_ER_err = [None] * n_PT
    Q0_h_ER = [None] * n_PT
    Q0_ER = [None] * n_PT
    Q1_ER = [None] * n_PT
    ER_x  = [None] * n_PT
    Eion = [None] * n_PT

    # Loop over PT values
    for i_PT in range(n_PT):
        t = T[i_PT] - absolute_zero_C

        # Find saturation values at this temperature
        pv, rhol, rhov, xliq, xvap, ierr, herr = RP.SATTdll(t, Mix_ratio, 1)
        #pv, rholv, rhovv, xliqv, xvapv, ierr, herr = RP.SATTdll(t, Mix_ratio, 2)
        if ierr != 0:
            errID[i_PT] = ierr + 0.1;
            continue

        Pvap[i_PT] = pv * psi_per_kPa

        # Find liquid density
        pl = P[i_PT] * kPa_per_psi

        rhol, ierr, herr = RP.TPRHOdll(t, pl, Mix_ratio, -1, 1, rhol)
        if ierr != 0:
            errID[i_PT] = ierr + 0.2
            continue

        Rho_l[i_PT] = rhol * RP.WMOLdll(Mix_ratio) * 0.001

        # Calculate intermolecular spacing then tolman length
        deltar = pow(1.e24 / (rhol * Na), 1/3) 
        Lt[i_PT] = deltar * pow(2., 1/6) / 3.

        def THERM_superheat(tT, rholT, xlT, rhovT, xvT, xfix = False):
            # Find the bubble density, composition, and other thermodynamic parameters
            # tT - temperature, rholT - liquid density, xlT - liquid composition
            # rhovT - guessed bubble density, xvT - guessed bubble composition,
            # xfix - Boolean to fix the bubble composition (vary density only)
            # Output -
            # bubble density, bubble composition, last step size, number of iterations

            xvT = np.array(xvT)
            xlT = np.array(xlT)
            uT = RP.CHEMPOTdll(tT,rholT,xlT)[0][:len(xvT)]
            utarget = uT
            gtarget = np.sum(uT * xlT)

            # Find gas pressure that matches this gibbs potential
            r1 = rhovT * xvT
            u1 = RP.CHEMPOTdll(tT, rhovT, xvT)[0][:len(xvT)]
            g1 = np.sum(u1 * xvT)

            rho2, ierr, herr = RP.TPRHOdll(tT, pl, xvT, 2, 1, rhovT)
            r2 = rho2 * xvT
            if ierr != 0:
                errID[i_PT] = ierr + 0.3

            u2 = RP.CHEMPOTdll(tT, rho2, xvT)[0][:len(xvT)]
            g2 = np.sum(u2 * xvT)

            if gtarget > (2*g1-g2) or gtarget < (2*g2-g1): 
                # AER 20200525 - Modified error criterion for mixtures.  For pure fluids, g2 < gtarget < g1.
                errID[i_PT] = 1.4

            i_trial = 0
            rshift = r1
            rguess = r1
            while (i_trial < C_maxtrials):
                if xfix:
                    rguess = r1 + ((gtarget - g1) * (r2 - r1) / (g2 - g1))
                else:
                    rguess = r1 + np.nan_to_num(np.subtract(utarget, u1) * np.subtract(r2, r1) / np.subtract(u2, u1))
                unew = RP.CHEMPOTdll(tT, np.sum(rguess), rguess/np.sum(rguess))[0][:len(xvT)]
                gnew = np.sum(unew * rguess) / np.sum(rguess)
                if utarget > unew:
                    u2 = unew
                    rshift = rguess - r2
                    r2 = rguess
                    g2 = np.sum(u2 * rguess) / np.sum(rguess)
                else:
                    u1 = unew
                    rshift = r1 - rguess
                    r1 = rguess
                    g1 = np.sum(u1 * rguess) / np.sum(rguess)
                i_trial += 1
                
                if xfix:
                    break_condition = abs((gtarget - gnew) / gtarget) < 1.0e-8
                else:
                    break_condition = np.amax(np.absolute(np.subtract(utarget, unew)/utarget)) < 1.0E-5
                if (np.amax(np.nan_to_num(rshift/rguess)) < 1.0E-8) and break_condition: break 
    
            if i_trial == C_maxtrials:
                #print(np.amax(np.absolute(np.subtract(utarget, unew)/utarget)), utarget, unew, xfix)
                errID[i_PT] = 2.4

            if xfix:
                rT = r1 + ((gtarget - g1) * (r2 - r1) / (g2 - g1))
            else:
                rT = r1 + np.nan_to_num(np.subtract(utarget, u1) * np.subtract(r2, r1) / np.subtract(u2, u1))
            return np.sum(rT), rT/np.sum(rT), np.sum(rshift), i_trial

        rhob, Xbub[i_PT], rhoerr, niter = THERM_superheat(t, rhol, Mix_ratio, rhov, Mix_ratio, True)
        #if errID[i_PT]: continue

        interp_iterations[i_PT] = [niter]
        Rho_b[i_PT] = rhob * RP.WMOLdll(Xbub[i_PT]) * 0.001
        Rho_b_err[i_PT] = rhoerr * RP.WMOLdll(Xbub[i_PT]) * 0.001

        # Find Thermodynamic properties of the superheated liquid
        preadback, el, hl, sl, cvl, cpl, wl, hjt = RP.THERMdll(t, rhol, Mix_ratio)
        if errID[i_PT]: continue
        P_err[i_PT] = [(preadback - pl) * psi_per_kPa,]
        A, Gl = RP.AGdll(t, rhol, Mix_ratio)

        # Find Thermodynamic properties of the bubble vapor
        preadback, ev, hv, sv, cvv, cpv, wv, hjt = RP.THERMdll(t, rhob, Xbub[i_PT])
        A, Gv = RP.AGdll(t, rhob, Xbub[i_PT])
        one_over_wmm = 1 / RP.WMOLdll(Xbub[i_PT])
        DeltaH[i_PT] = (hv - hl) * 1.0E3 * one_over_wmm
        DeltaS[i_PT] = (sv - sl) * 1.0E3 * one_over_wmm
        G_err[i_PT] = (Gv - Gl) * 1.0E3 * one_over_wmm
        Gibbs[i_PT] = Gv * 1.0E3 * one_over_wmm
        Pbub[i_PT] = preadback * psi_per_kPa

        # Recalculate critical bubble properties for adiabatic growth (may be different from the sudden growth model for mixtures)

        rhob, Xbub_ER[i_PT], rhoerr, niter = THERM_superheat(t, rhol, Mix_ratio, rhob, Mix_ratio)
        if errID[i_PT]: continue

        preadback, eb, hb, sb, cvb, cpb, wb, hjt = RP.THERMdll(t, rhob, Xbub[i_PT])
        Pb_ER[i_PT] = preadback * psi_per_kPa
        
        interp_iterations[i_PT].append(niter)
        
        ER_dmix = (Xbub[i_PT] - 0.001 * Xbub_ER[i_PT]) / 0.999
        rhodl, ierr, herr = RP.TPRHOdll(t, pl, ER_dmix, -1, 1, rhol)
        preadback, edl, hdl, sdl, cvdl, cpdl, wdl, hjt = RP.THERMdll(t, rhodl, ER_dmix)
        P_err[i_PT].append((preadback - pl) * psi_per_kPa)

        DeltaH_ER[i_PT] = (hb + 999. * hdl - 1000. * hl) * 1.0E3 / RP.WMOLdll(Xbub_ER[i_PT])
        
        Rho_b_ER[i_PT] = rhob * RP.WMOLdll(Xbub_ER[i_PT]) * 0.001
        Rho_b_ER_err[i_PT] = rhoerr * RP.WMOLdll(Xbub_ER[i_PT]) * 0.001

        #Draw adiabat connecting bubble vapor to liquid density
        preadback, e, h, s, cv, cp, w, hjt = RP.THERMdll(t, rhov, xvap)
        DeltaW_spike[i_PT] = e # temp value
        tspike, ierr, herr = RP.DSFL1dll(rhol, s, xvap)
        pspike, e, h, s, cv, cp, w, hjt = RP.THERMdll(tspike, rhol, xvap)
        wmmvap = RP.WMOLdll(xvap)
        DeltaE_spike[i_PT] = (e - el) * 1.0e3 / wmmvap
        Tspike[i_PT] = tspike + absolute_zero_C
        Pspike[i_PT] = pspike * psi_per_kPa
        DeltaW_spike[i_PT] = (e - DeltaW_spike[i_PT]) * 1.0e3 / wmmvap

        # Find Surface Tension
        sigma, ierr, herr = RP.SURFTdll(t, -1, Mix_ratio)

        Sigma[i_PT] = sigma;

        sigma, ierr, herr = RP.STNdll(t, rhol, rhob, Mix_ratio, Xbub_ER[i_PT])
        if ierr != 0:
            errID[i_PT] = ierr + 0.5
            continue

        if ierr != 0:
            errID[i_PT] = ierr + 0.55
            continue

        Sigma_err[i_PT] = sigma - Sigma[i_PT]

        # Find Critical Radius
        Rc[i_PT] = (2.0 * Sigma[i_PT] / ((Pbub[i_PT] - P[i_PT]) * kPa_per_psi)) * 1.0E6
        Rc_ER[i_PT] = (2.0 * Sigma[i_PT] / ((Pb_ER[i_PT] - P[i_PT]) * kPa_per_psi)) * 1.0E6
        # Find dSigmadT
        tlow = t - 0.0005
        thigh = t + 0.0005

        siglow, ierr, herr = RP.SURFTdll(tlow, -1, Mix_ratio)
        if ierr != 0:
            errID[i_PT] = ierr + 0.7
            continue

        sighigh, ierr, herr = RP.SURFTdll(thigh, -1, Mix_ratio)
        if ierr != 0:
            errID[i_PT] = ierr + 0.75
            continue

        dSigmadT[i_PT] = (sighigh - siglow) / (thigh - tlow)

        siglow, ierr, herr = RP.STNdll(tlow, rhol, rhob, Mix_ratio, Xbub_ER[i_PT])
        if ierr != 0:
            errID[i_PT] = ierr + 0.6
            continue

        sighigh, ierr, herr = RP.STNdll(thigh, rhol, rhob, Mix_ratio, Xbub_ER[i_PT])
        if ierr != 0:
            errID[i_PT] = ierr + 0.65
            continue

        dSigmadT_err[i_PT] = ((sighigh - siglow) / (thigh - tlow)) - dSigmadT[i_PT]
        
        # Calculate Q
        V = (4/3) * math.pi * Rc[i_PT] * Rc[i_PT] * Rc[i_PT] * 1.0E-24; # L
        S = 4 * math.pi * Rc[i_PT] * Rc[i_PT] * 1.0E-18 # m^2

        Q0_h[i_PT] = Rho_b[i_PT] * V * DeltaH[i_PT] * keV_per_J
        Q1_h[i_PT] = 6 * Q0_h[i_PT]
        Q0_s[i_PT] = S * Sigma[i_PT] * keV_per_J
        Q1_s[i_PT] = 6 * Q0_s[i_PT]
        Q0_ds[i_PT] = -S * t * dSigmadT[i_PT] * keV_per_J
        Q1_ds[i_PT] = 4 * Q0_ds[i_PT] - 0.5 * Q1_h[i_PT]
        Q0_w[i_PT] = -2 * Q0_s[i_PT] / 3
        Q1_w[i_PT] = 6 * Q0_w[i_PT]
        Q0[i_PT] = Q0_h[i_PT] + Q0_s[i_PT] + Q0_ds[i_PT] + Q0_w[i_PT]
        Q1[i_PT] = Q1_h[i_PT] + Q1_s[i_PT] + Q1_ds[i_PT] + Q1_w[i_PT]

        Q[i_PT] = Q0[i_PT] + (Q1[i_PT] * Lt[i_PT] / Rc[i_PT])

        W0_ext[i_PT] = V * P[i_PT] * ((Rho_l[i_PT] - Rho_b[i_PT]) / Rho_l[i_PT]) * kPa_per_psi * keV_per_J
        W1_ext[i_PT] = 3.0 * W0_ext[i_PT]

        # Calculate Q_ER
        V = (4/3) * math.pi * Rc_ER[i_PT] * Rc_ER[i_PT] * Rc_ER[i_PT] * 1.0E-24; # L
        S = 4 * math.pi * Rc_ER[i_PT] * Rc_ER[i_PT] * 1.0E-18 # m^2

        Q0_h_ER[i_PT] = Rho_b_ER[i_PT] * V * DeltaH_ER[i_PT] * keV_per_J
        Q1_h_ER = 6 * Q0_h_ER[i_PT]
        Q0_s_ER = S * Sigma[i_PT] * keV_per_J
        Q1_s_ER = 6 * Q0_s_ER
        Q0_ds_ER = -S * t * dSigmadT[i_PT] * keV_per_J
        Q1_ds_ER = 4 * Q0_ds_ER - 0.5 * Q1_h_ER
        Q0_w_ER = -2 * Q0_s_ER / 3
        Q1_w_ER = 6 * Q0_w_ER
        Q0_ER[i_PT] = Q0_h_ER[i_PT] + Q0_s_ER + Q0_ds_ER + Q0_w_ER
        Q1_ER[i_PT] = Q1_h_ER + Q1_s_ER + Q1_ds_ER + Q1_w_ER

        Q_ER[i_PT] = Q0_ER[i_PT] + (Q1_ER[i_PT] * Lt[i_PT] / Rc_ER[i_PT])

        Eion[i_PT] = Q0_s_ER + Q0_ds_ER + V * pl * keV_per_J
        ER_x[i_PT] = Eion[i_PT] / Rho_l[i_PT] / Rc_ER[i_PT] * (rhol/rhob)**(1./3.) * 1.e4
        ERrej[i_PT] = 17 * math.exp(-ER_x[i_PT]/37.)

#    SeitzModeloutput_tuple = namedtuple('SeitzModeloutput', ["P", "T", "Q", "Rc", "Pvap", "Pbub", "DeltaH", "DeltaS", "Gibbs", "Rho_l", "Rho_b", "Sigma", "dSigmadT", "P_err", "Rho_b_err", "G_err", "Sigma_err", "dSigmadT_err", "errID", "interp_iterations", "Lt", "DeltaE_spike", "DeltaW_spike", "Tspike", "Pspike", "Q0_h", "Q1_h", "Q0_s", "Q1_s", "Q0_ds", "Q1_ds", "Q0_w", "Q1_w", "Q0", "Q1", "W0_ext", "W1_ext", "Units"])

    ## kludge to return single values instead of length-1 arrays
    if len(P)>1: return SeitzModeloutput_tuple(P, T, Q, Q_ER, Rc, Pvap, Pbub, DeltaH, DeltaS, Gibbs, Rho_l, Rho_b, Sigma, dSigmadT, P_err, Rho_b_err, G_err, Sigma_err, dSigmadT_err, errID, interp_iterations, Lt, DeltaE_spike, DeltaW_spike, Tspike, Pspike, Q0_h, Q1_h, Q0_s, Q1_s, Q0_ds, Q1_ds, Q0_w, Q1_w, Q0, Q1, W0_ext, W1_ext, Units, ERrej, Pb_ER, Rc_ER, Rho_b_ER, Xbub_ER, Xbub, DeltaH_ER, Rho_b_ER_err, Q0_h_ER, Q0_ER, Q1_ER, Eion, ER_x )
    else: return SeitzModeloutput_tuple(P[0], T[0], Q[0], Q_ER[0], Rc[0], Pvap[0], Pbub[0], DeltaH[0], DeltaS[0], Gibbs[0], Rho_l[0], Rho_b[0], Sigma[0], dSigmadT[0], P_err[0], Rho_b_err[0], G_err[0], Sigma_err[0], dSigmadT_err[0], errID[0], interp_iterations[0], Lt[0], DeltaE_spike[0], DeltaW_spike[0], Tspike[0], Pspike[0], Q0_h[0], Q1_h[0], Q0_s[0], Q1_s[0], Q0_ds[0], Q1_ds[0], Q0_w[0], Q1_w[0], Q0[0], Q1[0], W0_ext[0], W1_ext[0], Units, ERrej[0], Pb_ER[0], Rc_ER[0], Rho_b_ER[0], Xbub_ER[0], Xbub[0], DeltaH_ER[0], Rho_b_ER_err[0], Q0_h_ER[0], Q0_ER[0], Q1_ER[0], Eion[0], ER_x[0])
