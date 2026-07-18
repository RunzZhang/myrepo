from scipy.interpolate import interp1d
import numpy as np
"""G4
-> Found Active Target Material: GXe
   Energy [MeV] Photoelectric [cm2/g]      Compton [cm2/g]
---------------------------------------------------------
         0.0010            8.9160e+03            5.1529e-03
         0.0100            1.6616e+02            5.5061e-02
         0.1000            1.7872e+00            1.0773e-01
         1.0000            4.3402e-03            5.2718e-02
---------------------------------------------------------
 -> Found Active Target Material: LAr
   Energy [MeV] Photoelectric [cm2/g]      Compton [cm2/g]
---------------------------------------------------------
         0.0010            3.1799e+03            1.1896e-02
         0.0100            6.1636e+01            1.0162e-01
         0.1000            5.5435e-02            1.2794e-01
         1.0000            8.6108e-05            5.7715e-02

NIST: Argon
Photon    Incoher. Photoel.
Energy    Scatter. Absorb.

1.000E-03 7.079E-03 3.181E+03
1.500E-03 1.420E-02 1.102E+03
2.000E-03 2.204E-02 5.094E+02
3.000E-03 3.717E-02 1.682E+02
3.203E-03 4.001E-02 1.404E+02
3.203E-03 4.001E-02 1.273E+03
4.000E-03 5.023E-02 7.554E+02
5.000E-03 6.099E-02 4.210E+02
6.000E-03 6.971E-02 2.581E+02
8.000E-03 8.293E-02 1.170E+02
1.000E-02 9.288E-02 6.232E+01
1.500E-02 1.105E-01 1.927E+01
2.000E-02 1.213E-01 8.207E+00
3.000E-02 1.316E-01 2.403E+00
4.000E-02 1.351E-01 9.909E-01
5.000E-02 1.359E-01 4.946E-01
6.000E-02 1.353E-01 2.793E-01
8.000E-02 1.322E-01 1.127E-01
1.000E-01 1.280E-01 5.564E-02
1.500E-01 1.176E-01 1.544E-02
2.000E-01 1.087E-01 6.280E-03
3.000E-01 9.518E-02 1.832E-03
4.000E-01 8.554E-02 7.979E-04
5.000E-01 7.822E-02 4.345E-04
6.000E-01 7.244E-02 2.721E-04
8.000E-01 6.369E-02 1.381E-04
1.000E+00 5.730E-02 8.587E-05

Xenon

Photon    Incoher. Photoel.
Energy    Scatter. Absorb.

1.000E-03 4.417E-03 9.403E+03
1.072E-03 4.969E-03 8.133E+03
1.149E-03 5.582E-03 7.032E+03
1.149E-03 5.582E-03 7.334E+03
1.500E-03 8.522E-03 4.077E+03
2.000E-03 1.285E-02 2.080E+03
3.000E-03 2.123E-02 7.715E+02
4.000E-03 2.874E-02 3.730E+02
4.782E-03 3.383E-02 2.357E+02
4.782E-03 3.383E-02 6.890E+02
5.000E-03 3.515E-02 6.344E+02
5.104E-03 3.577E-02 5.995E+02
5.104E-03 3.577E-02 8.133E+02
5.275E-03 3.677E-02 7.515E+02
5.453E-03 3.777E-02 6.945E+02
5.453E-03 3.777E-02 8.018E+02
6.000E-03 4.071E-02 6.330E+02
8.000E-03 5.027E-02 2.998E+02
1.000E-02 5.848E-02 1.662E+02
1.500E-02 7.440E-02 5.550E+01
2.000E-02 8.486E-02 2.510E+01
3.000E-02 9.729E-02 8.078E+00
3.456E-02 1.008E-01 5.417E+00
3.456E-02 1.008E-01 3.244E+01
4.000E-02 1.038E-01 2.211E+01
5.000E-02 1.073E-01 1.227E+01
6.000E-02 1.091E-01 7.454E+00
8.000E-02 1.096E-01 3.363E+00
1.000E-01 1.081E-01 1.793E+00
1.500E-01 1.019E-01 5.651E-01
2.000E-01 9.555E-02 2.490E-01
3.000E-01 8.495E-02 8.009E-02
4.000E-01 7.688E-02 3.699E-02
5.000E-01 7.064E-02 2.088E-02
6.000E-01 6.559E-02 1.338E-02
8.000E-01 5.784E-02 6.940E-03
1.000E+00 5.211E-02 4.335E-03
"""

import io
import matplotlib.pyplot as plt
import pandas as pd

# ------------------------------------------------------------------
# 1. RAW DATA INPUT STRING
# ------------------------------------------------------------------
# Geant4 discrete extracted points
g4_xe_data = """
Energy_MeV Photoelectric_cm2_g Compton_cm2_g
0.0010            8.9160e+03            5.1529e-03
0.0100            1.6616e+02            5.5061e-02
0.1000            1.7872e+00            1.0773e-01
1.0000            4.3402e-03            5.2718e-02
"""

g4_ar_data = """
Energy_MeV Photoelectric_cm2_g Compton_cm2_g
0.0010            3.1799e+03            1.1896e-02
0.0100            6.1636e+01            1.0162e-01
0.1000            5.5435e-02            1.2794e-01
1.0000            8.6108e-05            5.7715e-02
"""

# NIST continuous database lines
nist_ar_data = """
Energy_MeV Compton_cm2_g Photoelectric_cm2_g
1.000E-03 7.079E-03 3.181E+03
1.500E-03 1.420E-02 1.102E+03
2.000E-03 2.204E-02 5.094E+02
3.000E-03 3.717E-02 1.682E+02
3.203E-03 4.001E-02 1.404E+02
3.203E-03 4.001E-02 1.273E+03
4.000E-03 5.023E-02 7.554E+02
5.000E-03 6.099E-02 4.210E+02
6.000E-03 6.971E-02 2.581E+02
8.000E-03 8.293E-02 1.170E+02
1.000E-02 9.288E-02 6.232E+01
1.500E-02 1.105E-01 1.927E+01
2.000E-02 1.213E-01 8.207E+00
3.000E-02 1.316E-01 2.403E+00
4.000E-02 1.351E-01 9.909E-01
5.000E-02 1.359E-01 4.946E-01
6.000E-02 1.353E-01 2.793E-01
8.000E-02 1.322E-01 1.127E-01
1.000E-01 1.280E-01 5.564E-02
1.500E-01 1.176E-01 1.544E-02
2.000E-01 1.087E-01 6.280E-03
3.000E-01 9.518E-02 1.832E-03
4.000E-01 8.554E-02 7.979E-04
5.000E-01 7.822E-02 4.345E-04
6.000E-01 7.244E-02 2.721E-04
8.000E-01 6.369E-02 1.381E-04
1.000E+00 5.730E-02 8.587E-05
"""

nist_xe_data = """
Energy_MeV Compton_cm2_g Photoelectric_cm2_g
1.000E-03 4.417E-03 9.403E+03 
1.072E-03 4.969E-03 8.133E+03 
1.149E-03 5.582E-03 7.032E+03 
1.149E-03 5.582E-03 7.334E+03 
1.500E-03 8.522E-03 4.077E+03 
2.000E-03 1.285E-02 2.080E+03 
3.000E-03 2.123E-02 7.715E+02 
4.000E-03 2.874E-02 3.730E+02 
4.782E-03 3.383E-02 2.357E+02 
4.782E-03 3.383E-02 6.890E+02 
5.000E-03 3.515E-02 6.344E+02 
5.104E-03 3.577E-02 5.995E+02 
5.104E-03 3.577E-02 8.133E+02 
5.275E-03 3.677E-02 7.515E+02 
5.453E-03 3.777E-02 6.945E+02 
5.453E-03 3.777E-02 8.018E+02 
6.000E-03 4.071E-02 6.330E+02 
8.000E-03 5.027E-02 2.998E+02 
1.000E-02 5.848E-02 1.662E+02 
1.500E-02 7.440E-02 5.550E+01 
2.000E-02 8.486E-02 2.510E+01 
3.000E-02 9.729E-02 8.078E+00 
3.456E-02 1.008E-01 5.417E+00 
3.456E-02 1.008E-01 3.244E+01 
4.000E-02 1.038E-01 2.211E+01 
5.000E-02 1.073E-01 1.227E+01 
6.000E-02 1.091E-01 7.454E+00 
8.000E-02 1.096E-01 3.363E+00 
1.000E-01 1.081E-01 1.793E+00 
1.500E-01 1.019E-01 5.651E-01 
2.000E-01 9.555E-02 2.490E-01 
3.000E-01 8.495E-02 8.009E-02 
4.000E-01 7.688E-02 3.699E-02 
5.000E-01 7.064E-02 2.088E-02 
6.000E-01 6.559E-02 1.338E-02 
8.000E-01 5.784E-02 6.940E-03 
1.000E+00 5.211E-02 4.335E-03 
"""

# ------------------------------------------------------------------
# 2. PARSE STRINGS INTO PANDAS DATAFRAMES
# ------------------------------------------------------------------
df_g4_xe = pd.read_csv(io.StringIO(g4_xe_data.strip()), sep=r"\s+")
df_g4_ar = pd.read_csv(io.StringIO(g4_ar_data.strip()), sep=r"\s+")
df_nist_ar = pd.read_csv(io.StringIO(nist_ar_data.strip()), sep=r"\s+")
df_nist_xe = pd.read_csv(io.StringIO(nist_xe_data.strip()), sep=r"\s+")

# ------------------------------------------------------------------
# 3. PLOTTING INITIALIZATION
# ------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# Styling parameters
linewidth = 1.8
markersize = 7

# --- GRAPH 1: ARGON ---
# NIST reference continuous lines
ax1.plot(
    df_nist_ar["Energy_MeV"],
    df_nist_ar["Photoelectric_cm2_g"],
    label="NIST Photoelectric",
    color="darkred",
    linestyle="-",
    linewidth=linewidth,
)
ax1.plot(
    df_nist_ar["Energy_MeV"],
    df_nist_ar["Compton_cm2_g"],
    label="NIST no coherent",
    color="navy",
    linestyle="-",
    linewidth=linewidth,
)

# Geant4 extracted discrete markers
ax1.plot(
    df_g4_ar["Energy_MeV"],
    df_g4_ar["Photoelectric_cm2_g"],
    label="Geant4 Photoelectric (LAr)",
    color="orange",
    marker="o",
    linestyle="None",
    markersize=markersize,
    markeredgecolor="black",
)
ax1.plot(
    df_g4_ar["Energy_MeV"],
    df_g4_ar["Compton_cm2_g"],
    label="Geant4 Compton (LAr)",
    color="cyan",
    marker="s",
    linestyle="None",
    markersize=markersize,
    markeredgecolor="black",
)

ax1.set_title("Photon Cross Sections in Argon ($Z=18$)", fontsize=14, weight="bold")
ax1.set_xlabel("Photon Energy [MeV]", fontsize=12)
ax1.set_ylabel("Mass Cross Section [$\mathrm{cm}^2/\mathrm{g}$]", fontsize=12)
ax1.set_xscale("log")
ax1.set_yscale("log")
# ax1.grid(True, which="both", linestyle="--", alpha=0.5)
ax1.legend(fontsize=10, loc="lower left")

# --- GRAPH 2: XENON ---
# NIST reference continuous lines
ax2.plot(
    df_nist_xe["Energy_MeV"],
    df_nist_xe["Photoelectric_cm2_g"],
    label="NIST Photoelectric",
    color="darkred",
    linestyle="-",
    linewidth=linewidth,
)
ax2.plot(
    df_nist_xe["Energy_MeV"],
    df_nist_xe["Compton_cm2_g"],
    label="NIST no coherent",
    color="navy",
    linestyle="-",
    linewidth=linewidth,
)

# Geant4 extracted discrete markers
ax2.plot(
    df_g4_xe["Energy_MeV"],
    df_g4_xe["Photoelectric_cm2_g"],
    label="Geant4 Photoelectric (GXe)",
    color="orange",
    marker="o",
    linestyle="None",
    markersize=markersize,
    markeredgecolor="black",
)
ax2.plot(
    df_g4_xe["Energy_MeV"],
    df_g4_xe["Compton_cm2_g"],
    label="Geant4 Compton (GXe)",
    color="cyan",
    marker="s",
    linestyle="None",
    markersize=markersize,
    markeredgecolor="black",
)

ax2.set_title("Photon Cross Sections in Xenon ($Z=54$)", fontsize=14, weight="bold")
ax2.set_xlabel("Photon Energy [MeV]", fontsize=12)
ax2.set_xscale("log")
# ax2.grid(True, which="both", linestyle="--", alpha=0.5)
ax2.legend(fontsize=10, loc="lower left")

# Final polishing layout adjustments
plt.tight_layout()
# plt.show()


def _get_log_interpolator(df, energy_col, cross_section_col):
    """Generates a log-log interpolator for continuous cross-section parsing."""
    return interp1d(
        np.log10(df[energy_col]),
        np.log10(df[cross_section_col]),
        kind="linear", fill_value="extrapolate"
    )


# ------------------------------------------------------------------
# 3. EXPORTABLE PROBABILITY CALCULATORS
# ------------------------------------------------------------------
def calculate_doped_compton_probabilities(energy_mev, mass_fraction_xe, mass_fraction_ar):
    """
    Calculates the relative interaction probabilities for Compton scattering
    on Argon vs Xenon within a doped mixture given target mass fractions.
    """
    total_mass = mass_fraction_xe + mass_fraction_ar
    w_xe = mass_fraction_xe / total_mass
    w_ar = mass_fraction_ar / total_mass

    interp_comp_ar = _get_log_interpolator(df_nist_ar, "Energy_MeV", "Compton_cm2_g")
    interp_comp_xe = _get_log_interpolator(df_nist_xe, "Energy_MeV", "Compton_cm2_g")

    log_E = np.log10(energy_mev)
    sigma_mass_ar = 10 ** interp_comp_ar(log_E)
    sigma_mass_xe = 10 ** interp_comp_xe(log_E)

    total_macroscopic_compton = (w_ar * sigma_mass_ar) + (w_xe * sigma_mass_xe)

    return {
        "Energy_MeV": energy_mev,
        "Ar_Mass_CrossSection_cm2_g": sigma_mass_ar,
        "Xe_Mass_CrossSection_cm2_g": sigma_mass_xe,
        "Ar_Interaction_Probability": (w_ar * sigma_mass_ar) / total_macroscopic_compton,
        "Xe_Interaction_Probability": (w_xe * sigma_mass_xe) / total_macroscopic_compton
    }


def calculate_doped_photoelectric_probabilities(energy_mev, mass_fraction_xe, mass_fraction_ar):
    """
    Calculates the relative interaction probabilities for the Photoelectric Effect
    on Argon vs Xenon within a doped mixture given target mass fractions.
    """
    total_mass = mass_fraction_xe + mass_fraction_ar
    w_xe = mass_fraction_xe / total_mass
    w_ar = mass_fraction_ar / total_mass

    interp_phot_ar = _get_log_interpolator(df_nist_ar, "Energy_MeV", "Photoelectric_cm2_g")
    interp_phot_xe = _get_log_interpolator(df_nist_xe, "Energy_MeV", "Photoelectric_cm2_g")

    log_E = np.log10(energy_mev)
    sigma_mass_ar = 10 ** interp_phot_ar(log_E)
    sigma_mass_xe = 10 ** interp_phot_xe(log_E)

    total_macroscopic_photoelectric = (w_ar * sigma_mass_ar) + (w_xe * sigma_mass_xe)

    return {
        "Energy_MeV": energy_mev,
        "Ar_Mass_CrossSection_cm2_g": sigma_mass_ar,
        "Xe_Mass_CrossSection_cm2_g": sigma_mass_xe,
        "Ar_Interaction_Probability": (w_ar * sigma_mass_ar) / total_macroscopic_photoelectric,
        "Xe_Interaction_Probability": (w_xe * sigma_mass_xe) / total_macroscopic_photoelectric
    }


# ------------------------------------------------------------------
# 4. STANDALONE TESTING DIAGNOSTICS
# ------------------------------------------------------------------
if __name__ == "__main__":
    gamma_energy = 0.040  # 40 keV (near the Xenon K-edge region)
    mass_xe = 0.01  # 1% Xenon dopant
    mass_ar = 0.99  # 99% Argon base

    print(f"--- Standalone Verification Analysis at {gamma_energy * 1000:.1f} keV ---")
    print(f"Mixture: {mass_ar * 100:.1f}% Argon, {mass_xe * 100:.1f}% Xenon by mass\n")

    # Run Compton Check
    comp_res = calculate_doped_compton_probabilities(gamma_energy, mass_xe, mass_ar)
    print(f"[COMPTON SCATTERING]")
    print(f"  Ar Cross Section: {comp_res['Ar_Mass_CrossSection_cm2_g']:.4f} cm2/g")
    print(f"  Xe Cross Section: {comp_res['Xe_Mass_CrossSection_cm2_g']:.4f} cm2/g")
    print(f"  ==> Scatter off Ar: {comp_res['Ar_Interaction_Probability'] * 100:.2f}%")
    print(f"  ==> Scatter off Xe: {comp_res['Xe_Interaction_Probability'] * 100:.2f}%\n")

    # Run Photoelectric Check
    phot_res = calculate_doped_photoelectric_probabilities(gamma_energy, mass_xe, mass_ar)
    print(f"[PHOTOELECTRIC ABSORPTION]")
    print(f"  Ar Cross Section: {phot_res['Ar_Mass_CrossSection_cm2_g']:.4f} cm2/g")
    print(f"  Xe Cross Section: {phot_res['Xe_Mass_CrossSection_cm2_g']:.4f} cm2/g")
    print(f"  ==> Absorb by Ar: {phot_res['Ar_Interaction_Probability'] * 100:.2f}%")
    print(f"  ==> Absorb by Xe: {phot_res['Xe_Interaction_Probability'] * 100:.2f}%")