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
# v1.3 Disabled a print warning and put it into the err_ID
#
# INPUTS:
#               PIn   -   Pressure (psia) of the superheated liquid.  Have to be a
#                          list or a scalar. If P is also a list of many values, they
#                          must have the same number of elements.
#               TIn   -   Temperature (C) of the superheated liquid. Hate to be a
#                          list or a scalar. If T is also a list of many values, they
#                          must have the same number of elements.
#               Fluid -   Name of the fluid or array of fluid names.  In the
#                         case of an array, input a string with elements
#                         delimited by vertical bars '|' (C3F8|C4F10). Any names used
#                         should match the stem of the .FLD file (case
#                         insensitive) located in the FLUIDS folder.
#                         Default is 'C3F8'.
#               Mixture - List containing the mixture ratio ([0.4, 0.6]). The number of
#                         The number of elements must be equal to the number of fluids
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
#               Q   -   (keV) Minimum heat input to create critical bubble
#              Rc   -   (nm) Radius of critical bubble
#            Pvap   -   (psia) Vapor pressure of saturated fluid at
#                         temperature T
#            Pbub   -   (psia) Pressure inside hot-spike bubble (vapor pressure of
#                         superheated fluid at T and P, and liquid mixture composition)
#          DeltaH   -   (J/kg) Change in specific enthalpy between
#                         superheated liquid and bubble vapor states
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
#            PbER   -   (psia) Pressure inside adiabatically bubble (vapor pressure
#                        of superheated fluid at T and P, and equilibrium mixture
#                        composition) (placeholder)
#            RcER   -   (nm) Radius of critical ER bubble (= Rc for pure fluids) (placeholder)
#         Rho_bER   -   (g/cc) Density of ER bubble (= Rho_bER for pure fluids) (placeholder)
#          mix_ER   -   Critical bubble molar gas composition ratio (placeholder)
#
#   Consisitency Checks
#
#           P_err   -   (psia) Preadback - P, where Preadback is the
#                         pressure corresponding to T and Rho_l
#       Rho_b_err   -   (g/cc) Size of the last step in the iterative
#                         process to find Rho_b
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

class SeitzModeloutput_tuple(namedtuple('SeitzModeloutput', ["P", "T", "Q", "Rc", "Pvap", "Pbub", "DeltaH", "DeltaS", "Gibbs", "Rho_l", "Rho_b", "Sigma", "dSigmadT", "P_err", "Rho_b_err", "G_err", "Sigma_err", "dSigmadT_err", "errID", "interp_iterations", "Lt", "DeltaE_spike", "DeltaW_spike", "Tspike", "Pspike", "Q0_h", "Q1_h", "Q0_s", "Q1_s", "Q0_ds", "Q1_ds", "Q0_w", "Q1_w", "Q0", "Q1", "W0_ext", "W1_ext", "Units", "ERrej", "PbER", "RcER", "Rho_bER", "mix_ER"])):
    def __repr__(self):
        string = ""
        for prop in self._fields:
            string += "\t{0:>20s} : {1:<20s}\n".format(prop, str(getattr(self, prop)))
        return string

class SeitzModelOutput(object):
    """
    The class which contains the data from :func:`~SeitzModel`.

    Each property calculated by `SeitzModel` is accessed as attributes of this
    class, e.g. SeitzModelOutput.Q.

    Parameters
    ----------
    fields : list of str
        The names of the properties.
    data : list of float or list of array of float
        The values of the properties calculated by `SeitzModel`.  Must be same
        length as *fields*.
    units : list of str
        The units of each property.   Must be same length as *fields*.
    """

    def __init__(self, fields, data, units):
        assert( len(fields) == len(data) == len(units) )
        for field, datum in zip( fields, data ):
            if field=="Units": continue
            if type(datum)==list: setattr(self, field, np.array(datum))
            else: setattr(self, field, datum)
        setattr(self, "Units", units)

    def __repr__(self):
        ## A super kludged way of imitating the old MATLAB-style printing
        string = ""
        for prop in vars(self):
            if getattr(self, prop) is not None and prop!="Units":
                if type(getattr(self, prop)) is not np.ndarray:
                    string += "\t{0:>20s} : {1:<15.5f} {2:s}\n".format(prop, getattr(self, prop), self.Units[prop])
                else:
                    string += "\t{0:>20s} : [".format(prop)
                    for i in range(len(getattr(self, prop))):
                        if getattr(self, prop)[i] is not None:
                            string += "{0:.3f}, ".format(getattr(self, prop)[i])
                        else:
                            string += "None, "
                    string = string[:-2] + "] {0:s}\n".format(self.Units[prop])
            else:
                string += "\t{0:>20s} : {1:<15s}\n".format(prop, str(getattr(self, prop)))
        return string

def SeitzModel(PIn, TIn, Fluid = "R218", Mix_ratio = [1]):
    """
    Compute fluid properties, including Seitz threshold, using REFPROP.

    Parameters
    ----------
    PIn : float or list of float
        The pressure of the fluid, in psia, at which the fluid properties are
        computed.
    TIn : float or list of float
        The temperature of the fluid, in ℃, at which the fluid properties are
        computed.
    Fluid : str, default="R218"
        The fluid for which the Seitz threshold will be calculated. Default is
        R218, a.k.a. C_3F_8.  Mixtures may be used by including the name of
        each substance separated by vertical bars (:code:`|`), and using the
        *Mix_ratio* parameter.
    Mix_ratio : list of float, default=[1]
        The mixture ratio for each substance included in the *Fluid* parameter.
        Must sum to 1.

    Returns
    -------
    :class:`~SeitzModelOutput`
        The object which contains the computed data.  Data are accessible as
        attributes of the object.

    Notes
    -----
    *PIn* and *TIn* may each be passed as floats or lists.

    If both are lists and the same length, the properties are calculated for
    each P/T pair, i.e. for (PIn[0], TIn[0]), (PIn[1], TIn[1]), etc. and the
    attributes of :class:`~SeitsModelOutput` are numpy arrays.

    If one is a float while the other is a list, the properties are calculated
    as though the float were a list of the same length of the one passed as a
    list, with all values the same as the float.

    If they are boths lists but the lengths differ, a ValueError is raised.

    Raises
    ------
    ValueError
        If *PIn* and *TIn* are lists or arrays of different lengths.

    Examples
    --------
    >>> # Print the Seitz energy threshold
    >>> Seitz = SeitzModel(30, 15)
    >>> Seitz.Q 
    2.726618642546287

    >>> # Two temperatures, one pressure
    >>> Seitz = SeitzModel(30, [14, 15])
    >>> Seitz.Q
    array([3.15461861, 2.72661864])

    >>> # Two different thresholds
    >>> Seitz = SeitzModel([25, 30], [14, 15])
    >>> Seitz.Q
    array([3.02997352, 2.72661864])

    >>> # Mixture of substances: 95% C3F8, 5% C4F10
    >>> Seitz = SeitzModel(30, 16, "c3f8|c4f10", [0.95, .05])
    >>> Seitz.Q
    2.9517274123683976
    """

    # CONSTANTS! must not be overwritten!
    C_maxtrials = 1000
    absolute_zero_C = -273.15
    kPa_per_psi = 6.89475729317
    psi_per_kPa = 1 / kPa_per_psi
    keV_per_J = 6.24150636309e15
    Na = 6.0221415E23

    # Units of all the variables
    Fields = ["P", "T", "Q", "Rc", "Pvap", "Pbub", "DeltaH", "DeltaS", "Gibbs", "Rho_l", "Rho_b", "Sigma", "dSigmadT", "P_err", "Rho_b_err", "G_err", "Sigma_err", "dSigmadT_err", "errID", "interp_iterations", "Lt", "DeltaE_spike", "DeltaW_spike", "Tspike", "Pspike", "Q0_h", "Q1_h", "Q0_s", "Q1_s", "Q0_ds", "Q1_ds", "Q0_w", "Q1_w", "Q0", "Q1", "W0_ext", "W1_ext", "Units", "ERrej", "PbER", "RcER", "Rho_bER", "mix_ER"]
    Units = dict(zip(Fields, ("psia", "C", "keV", "nm", "psia", "psia", "J/kg", "J/kg-K", "J/kg", "g/cc", "g/cc", "N/m", "N/m-K", "psia", "g/cc", "J/kg", "N/m", "N/m-K", "ierr.SeitzStep", "count", "nm", "J/kg", "J/kg", "C", "psia", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "keV", "", "?", "?", "?", "?", "?", "?")))

    # Open the library in the same folder as this module
    RPPath = os.path.dirname( os.path.join( os.path.abspath(__file__)) )
    RP = ct.REFPROPFunctionLibrary(RPPath)
    RP.SETPATHdll( RPPath )

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
            raise ValueError("Pressure and temperature lists have different lengths")
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
    # Fixed the warning number and disabled the print since it can be overwelming in loops
        if ierr == -117:
            #print(herr)
            pass
        else:
            return ierr, herr

    # Get Fluid properties
    wmm = RP.WMOLdll(Mix_ratio)
    #wmm, Ttrp, Tnbpt, Tc, Pc, Dc, Zc, acf, dip, Rgas = RP.INFOdll(Mix_ratio)
    one_over_wmm = 1 / wmm

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
    PbER = [None] * n_PT
    RcER = [None] * n_PT
    Rho_bER = [None] * n_PT
    mix_ER = [None] * n_PT

    # Loop over PT values
    for i_PT in range(n_PT):
        # If the mixture warning shows up, it will be put in errID
        if ierr == -117:
            errID[i_PT] = ierr + 0.05;

        t = T[i_PT] - absolute_zero_C

        # Find saturation values at this temperature
        p, rhol, rhov, xliq, xvap, ierr, herr = RP.SATTdll(t, Mix_ratio, 1)
        pv, rholv, rhovv, xliqv, xvapv, ierr, herr = RP.SATTdll(t, Mix_ratio, 2)
        if ierr != 0:
            errID[i_PT] = ierr + 0.1;
            continue

        Pvap[i_PT] = p * psi_per_kPa

        # Find liquid density
        p = P[i_PT] * kPa_per_psi

        rho, ierr, herr = RP.TPRHOdll(t, p, Mix_ratio, -1, 1, rhol)
        if ierr != 0:
            errID[i_PT] = ierr + 0.2
            continue

        Rho_l[i_PT] = rho * wmm * 0.001
        rhol = rho

        Lt[i_PT] = (2/3) * pow(7.5e23 / (rho * Na * math.sqrt(18)), 1/3)

        # Find Thermodynamic properties of the superheated liquid
        preadback, e, h, s, cv, cp, w, hjt = RP.THERMdll(t, rhol, Mix_ratio)
        A, G = RP.AGdll(t, rhol, Mix_ratio)
        DeltaH[i_PT] = h # temp value
        DeltaS[i_PT] = s # temp value
        DeltaE_spike[i_PT] = e # temp value
        gtarget = G
        P_err[i_PT] = (preadback - p) * psi_per_kPa

        # Find gas pressure that matches this gibbs potential
        r1 = rhovv # Use gas density using same composition as liquid
        r2 = rhovv
        A, g1 = RP.AGdll(t, r1, Mix_ratio)

        r2, ierr, herr = RP.TPRHOdll(t, p, Mix_ratio, 2, 1, r2)
        if ierr != 0:
            errID[i_PT] = ierr + 0.3
            continue

        A, g2 = RP.AGdll(t, r2, Mix_ratio)

        if gtarget > (2*g1-g2) or gtarget < (2*g2-g1): 
            # AER 20200525 - Modified error criterion for mixtures.  For pure fluids, g2 < gtarget < g1.
            errID[i_PT] = 1.4
            continue

        i_trial = 0
        rshift = r1
        rguess = r1
        while (i_trial < C_maxtrials) and ((((rshift/rguess) > 1.0E-8)) or (((gtarget - gnew)/gtarget) < -1.0E-8) or (((gtarget - gnew)/gtarget) > 1.0E-8)):
            rguess = r1 + ((gtarget - g1) * (r2 - r1) / (g2 - g1))
            A, gnew = RP.AGdll(t, rguess, Mix_ratio)
            if gtarget > gnew:
                g2 = gnew
                rshift = rguess - r2
                r2 = rguess
            else:
                g1 = gnew
                rshift = r1 - rguess
                r1 = rguess
            i_trial += 1

        interp_iterations[i_PT] = i_trial
        if i_trial == C_maxtrials:
            errID[i_PT] = 2.4
            continue

        rhov = r1 + ((gtarget - g1) * (r2 - r1) / (g2 - g1))
        Rho_b[i_PT] = rhov * wmm * 0.001
        Rho_b_err[i_PT] = rshift * wmm * 0.001

        # Find Thermodynamic properties of the bubble vapor
        preadback, e, h, s, cv, cp, w, hjt = RP.THERMdll(t, rhov, Mix_ratio)
        A, G = RP.AGdll(t, rhov, Mix_ratio)
        DeltaW_spike[i_PT] = e # temp value
        DeltaH[i_PT] = (h - DeltaH[i_PT]) * 1.0E3 * one_over_wmm
        DeltaS[i_PT] = (s - DeltaS[i_PT]) * 1.0E3 * one_over_wmm
        G_err[i_PT] = (G - gtarget) * 1.0E3 * one_over_wmm
        Gibbs[i_PT] = G * 1.0E3 * one_over_wmm
        Pbub[i_PT] = preadback * psi_per_kPa

        #Draw adiabat connecting bubble vapor to liquid density
        tspike, ierr, herr = RP.DSFL1dll(rhol, s, Mix_ratio)
        pspike, e, h, s, cv, cp, w, hjt = RP.THERMdll(tspike, rhol, Mix_ratio)
        DeltaE_spike[i_PT] = (e - DeltaE_spike[i_PT]) * 1.0e3 * one_over_wmm
        Tspike[i_PT] = tspike + absolute_zero_C
        Pspike[i_PT] = pspike * psi_per_kPa
        DeltaW_spike[i_PT] = (e - DeltaW_spike[i_PT]) * 1.0e3 * one_over_wmm

        # Find Surface Tension
        if len(Mix_ratio) == 1: # STNdll doesn't work with mixtures.
            sigma, ierr, herr = RP.STNdll(t, rhol, rhovv, Mix_ratio, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.5
                continue
        else: 
            sigma, ierr, herr = RP.SURFTdll(t, -1, Mix_ratio)

        Sigma[i_PT] = sigma;
        #if len(Mix_ratio) == 1:
        sigma, ierr, herr = RP.SURFTdll(t, -1, Mix_ratio)
        #else:
            #rhols, sigma, ierr, herr = RP.SURFTdll(t, -1, Mix_ratio)
        if ierr != 0:
            errID[i_PT] = ierr + 0.55
            continue

        Sigma_err[i_PT] = sigma - Sigma[i_PT]

        # Find Critical Radius
        Rc[i_PT] = (2.0 * Sigma[i_PT] / ((Pbub[i_PT] - P[i_PT]) * kPa_per_psi)) * 1.0E6
        RcER[i_PT] = (2.0 * Sigma[i_PT] / ((Pbub[i_PT] - P[i_PT]) * kPa_per_psi)) * 1.0E6
        # Find dSigmadT
        tlow = t - 0.0005
        thigh = t + 0.0005

        if len(Mix_ratio) == 1:
            siglow, ierr, herr = RP.STNdll(tlow, rhol, rhovv, Mix_ratio, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.6
                continue

            sighigh, ierr, herr = RP.STNdll(thigh, rhol, rhovv, Mix_ratio, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.65
                continue

            dSigmadT[i_PT] = (sighigh - siglow) / (thigh - tlow)

            siglow, ierr, herr = RP.SURFTdll(tlow, -1, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.7
                continue

            sighigh, ierr, herr = RP.SURFTdll(thigh, -1, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.75
                continue

            dSigmadT_err[i_PT] = ((sighigh - siglow) / (thigh - tlow)) - dSigmadT[i_PT]
        else:
            siglow, ierr, herr = RP.SURFTdll(tlow, -1, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.7
                continue

            sighigh, ierr, herr = RP.SURFTdll(thigh, -1, Mix_ratio)
            if ierr != 0:
                errID[i_PT] = ierr + 0.75
                continue

            dSigmadT[i_PT] = ((sighigh - siglow) / (thigh - tlow))

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

#    SeitzModeloutput_tuple = namedtuple('SeitzModeloutput', ["P", "T", "Q", "Rc", "Pvap", "Pbub", "DeltaH", "DeltaS", "Gibbs", "Rho_l", "Rho_b", "Sigma", "dSigmadT", "P_err", "Rho_b_err", "G_err", "Sigma_err", "dSigmadT_err", "errID", "interp_iterations", "Lt", "DeltaE_spike", "DeltaW_spike", "Tspike", "Pspike", "Q0_h", "Q1_h", "Q0_s", "Q1_s", "Q0_ds", "Q1_ds", "Q0_w", "Q1_w", "Q0", "Q1", "W0_ext", "W1_ext", "Units"])

    ## kludge to return single values instead of length-1 arrays
#    if len(P)>1: return SeitzModeloutput_tuple(P, T, Q, Rc, Pvap, Pbub, DeltaH, DeltaS, Gibbs, Rho_l, Rho_b, Sigma, dSigmadT, P_err, Rho_b_err, G_err, Sigma_err, dSigmadT_err, errID, interp_iterations, Lt, DeltaE_spike, DeltaW_spike, Tspike, Pspike, Q0_h, Q1_h, Q0_s, Q1_s, Q0_ds, Q1_ds, Q0_w, Q1_w, Q0, Q1, W0_ext, W1_ext, Units, ERrej, PbER, RcER, Rho_bER, mix_ER)
#    else: return SeitzModeloutput_tuple(P[0], T[0], Q[0], Rc[0], Pvap[0], Pbub[0], DeltaH[0], DeltaS[0], Gibbs[0], Rho_l[0], Rho_b[0], Sigma[0], dSigmadT[0], P_err[0], Rho_b_err[0], G_err[0], Sigma_err[0], dSigmadT_err[0], errID[0], interp_iterations[0], Lt[0], DeltaE_spike[0], DeltaW_spike[0], Tspike[0], Pspike[0], Q0_h[0], Q1_h[0], Q0_s[0], Q1_s[0], Q0_ds[0], Q1_ds[0], Q0_w[0], Q1_w[0], Q0[0], Q1[0], W0_ext[0], W1_ext[0], Units, ERrej[0], PbER[0], RcER[0], Rho_bER[0], mix_ER[0])
    if len(P)>1: return SeitzModelOutput(Fields, [P, T, Q, Rc, Pvap, Pbub, DeltaH, DeltaS, Gibbs, Rho_l, Rho_b, Sigma, dSigmadT, P_err, Rho_b_err, G_err, Sigma_err, dSigmadT_err, errID, interp_iterations, Lt, DeltaE_spike, DeltaW_spike, Tspike, Pspike, Q0_h, Q1_h, Q0_s, Q1_s, Q0_ds, Q1_ds, Q0_w, Q1_w, Q0, Q1, W0_ext, W1_ext, Units, ERrej, PbER, RcER, Rho_bER, mix_ER], Units)
    else: return SeitzModelOutput(Fields, [P[0], T[0], Q[0], Rc[0], Pvap[0], Pbub[0], DeltaH[0], DeltaS[0], Gibbs[0], Rho_l[0], Rho_b[0], Sigma[0], dSigmadT[0], P_err[0], Rho_b_err[0], G_err[0], Sigma_err[0], dSigmadT_err[0], errID[0], interp_iterations[0], Lt[0], DeltaE_spike[0], DeltaW_spike[0], Tspike[0], Pspike[0], Q0_h[0], Q1_h[0], Q0_s[0], Q1_s[0], Q0_ds[0], Q1_ds[0], Q0_w[0], Q1_w[0], Q0[0], Q1[0], W0_ext[0], W1_ext[0], Units, ERrej[0], PbER[0], RcER[0], Rho_bER[0], mix_ER[0]], Units)
