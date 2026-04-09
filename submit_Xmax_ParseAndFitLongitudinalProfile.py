#!/usr/bin/env python3
#
# File that parses the longitudinal profiles from the CORSIKA .long files
# and fits the profiles to extract Xmax, R, and L parameters.
# Also extracts the muon and EM numbers at ground level for crosschecks with
# particle data block file (read from corsikaReader.cpp).
#
# Usage:
# python3 ParseAndFitLongitudinalProfile_Auger.py <InputLongFile> --zen <ZenithAngleInRadians> [--removeFinal20gcm2]
#

import numpy as np
from scipy.optimize import curve_fit
import argparse
import os 
import glob

parser = argparse.ArgumentParser()
# parser.add_argument("input", type=str, help="Input data files.")
# parser.add_argument("--zen", required=True, type=float, help="Zenith angle of the shower (rad)")

parser.add_argument("--longbase", required=True, type=str,
                    help="Base directory of .long files")
parser.add_argument("--zenfile", required=True, type=str,
                    help=".npz file containing zenith info")

parser.add_argument("--outfile", required=True, type=str,
                    help="Output .dat file")

parser.add_argument("--removeFinal20gcm2", action="store_true", help="Remove final 20 g/cm2 from longitudinal profile before fitting.")

args = parser.parse_args()

# def load_zeniths(zenfile):
#     zenith_dict = {}

#     with open(zenfile) as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue

#             cols = [c.strip() for c in line.split()]

#             run = int(cols[0])                 # RunNumber
#             zenith = float(cols[3])       # column 3 = Zenith

#             zenith_dict[run] = zenith

#     return zenith_dict
def load_zeniths(zenfile):

    data = np.load(zenfile)
    
    # Convert arrays to native Python lists first
    runs = data['run'].tolist()
    zeniths = data['zenith'].tolist()

    # Zip them into a dictionary
    zenith_dict = dict(zip(runs, zeniths))

    return zenith_dict

long_files = glob.glob(
    os.path.join(args.longbase, "**", "DAT*.long"),
    recursive=True
)

def get_run_from_longfile(longfile):
    # .../000123/DAT000123.long → 000123
    return int(os.path.basename(longfile)[3:9])



def LongFileParser(filename):
    depths = []
    positrons = []
    electrons = []
    muPlus = []
    muMinus = []
    chargedParticles = []

    file = open(filename, "r")
    nLines = 0
    xmax = -1
    Rcorsika = -1
    Lcorsika = -1
    energyDeposit = 0

    # x0 = np.nan
    # lambdaApprox = np.nan
    # Rcorsika = np.nan
    # Lcorsika = np.nan

    for iline, line in enumerate(file):
        if iline == 0:
            nLines = line.split()[3]
            continue

        if iline <= 2:
            continue  # Skip header

        cols = line.split()

        if "ENERGY" in cols and "DEPOSIT" in cols:
            energyDeposit += 1

        if "PARAMETERS" in cols:
            xmax = float(cols[4])
            x0 = float(cols[3])
            lambdaApprox = float(cols[5]) # Drop 1st + 2nd order corrections

            x0Prime = x0 - xmax
            Rcorsika = float(np.sqrt(lambdaApprox / abs(x0Prime))) # Shape parameter
            Lcorsika = float(np.sqrt(abs(x0Prime * lambdaApprox))) # Characteristic width

        if len(cols) != 10 or "FIT" in cols:
            continue

        if energyDeposit > 0:
            continue  # Skip all energy deposit lines in .long file

        depths.append(float(cols[0]))
        positrons.append(float(cols[2]))
        electrons.append(float(cols[3]))
        muPlus.append(float(cols[4]))
        muMinus.append(float(cols[5]))
        chargedParticles.append(float(cols[7]))

    file.close()

    return depths, positrons, electrons, muPlus, muMinus, chargedParticles, xmax, Rcorsika, Lcorsika, x0, lambdaApprox


def FindGroundIndex(ground, depths):
    absDepthToGround = []

    for i in range(len(depths)):
        absDepthToGround.append(abs(depths[i] - ground))

    minDepthToGround = min(absDepthToGround)
    indGround = absDepthToGround.index(minDepthToGround)

    return indGround



def GHFunction(X, Nmax, Xmax, X0, lamb):
    return Nmax * ((X - X0) / (Xmax - X0)) ** ((Xmax - X0) / lamb) * np.exp((Xmax - X) / lamb)

def GHFunctionWithABS(X, Nmax, Xmax, X0, lamb):
    return Nmax * (abs(X - X0) / (Xmax - X0)) ** ((Xmax - X0) / lamb) * np.exp((Xmax - X) / lamb)


def AndringaFunction(X, Xmax, R, L):
    return (1 + (R * (X - Xmax) / L)) ** (1 / (R**2)) * np.exp(-1. * (X - Xmax) / (L * R))

def AndringaFunctionWithABS(X, Xmax, R, L):
    return abs(1 + (R * (X - Xmax) / L)) ** (1 / (R**2)) * np.exp(-1. * (X - Xmax) / (L * R))


def FitLongitudinalProfile(depths, chargedParticles, NmaxGuess, XmaxGuess, X0Guess, lambGuess, shift=False, absoluteValue=False):

    if (shift == True) and (absoluteValue == False):
        # Shift all depths by 100 g/cm^2, makes fits more robust against large X0 values
        depthArray = np.asarray(depths) + 100.0
    else:
        depthArray = np.asarray(depths)

    particleArray = np.asarray(chargedParticles)

    # Assume Poissonian-like relative uncertainty for N, i.e. sqrt(N) / N = 1 / sqrt(N)
    # Only used as an uncertainty estimate in weighting the fit, applying emphasis near Xmax
    # Not to be used as a true, physical uncertainty, on the fluctuations in N
    uncerts = 1.0 / np.sqrt(particleArray)

    # Insert a try-except statement for poor fits...
    try:
        if (shift == True) and (absoluteValue == False):
            popt, pcov = curve_fit(GHFunction, depthArray, particleArray, p0=[NmaxGuess, XmaxGuess+100.0, X0Guess, lambGuess], sigma=uncerts)
        elif (absoluteValue == True) and (shift == False):
            popt, pcov = curve_fit(GHFunctionWithABS, depthArray, particleArray, p0=[NmaxGuess, XmaxGuess, X0Guess, lambGuess], sigma=uncerts)
        else:
            popt, pcov = curve_fit(GHFunction, depthArray, particleArray, p0=[NmaxGuess, XmaxGuess, X0Guess, lambGuess], sigma=uncerts)
    except RuntimeError:
        print('RuntimeError')
        return np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf

    perr = np.sqrt(np.diag(pcov))

    NmaxFit = popt[0]
    XmaxFit = popt[1]
    X0Fit = popt[2]
    lambFit = popt[3]

    NmaxSigma = perr[0]
    XmaxSigma = perr[1]
    X0Sigma = perr[2]
    lambSigma = perr[3]

    X0Prime = abs(X0Fit - XmaxFit)

    RFit = np.sqrt(lambFit / X0Prime)
    LFit = np.sqrt(X0Prime * lambFit)

    # Squared terms found from error prop.
    sigmaRterm1 = (lambSigma ** 2) / (4 * lambFit * X0Prime)
    sigmaRterm2 = (lambFit * (X0Sigma ** 2 + XmaxSigma ** 2)) / (4 * X0Prime ** 3)

    RFitSigma = np.sqrt(sigmaRterm1 + sigmaRterm2)

    # Squared terms found from error prop.
    sigmaLterm1 = (X0Prime * lambSigma ** 2) / (4 * lambFit)
    sigmaLterm2 = (lambFit * (X0Sigma ** 2 + XmaxSigma ** 2)) / (4 * X0Prime)

    LFitSigma = np.sqrt(sigmaLterm1 + sigmaLterm2)

    if (shift == True) and (absoluteValue == False):
        XmaxFit = popt[1] - 100.0 # Shift Xmax from fit back to real xmax
        X0Fit = popt[2] - 100.0 # Shift X0 value from fit back to true X0 for observed profile    

    return RFit, RFitSigma, LFit, LFitSigma, XmaxFit, XmaxSigma, *popt


def FitLongitudinalProfileAndringa(depths, NprimeArray, XmaxGuess, RGuess, LGuess, shift=False, absoluteValue=False):

    if (shift == True) and (absoluteValue == False):
        # Shift all depths by 100 g/cm^2, makes fits more robust against large X0 values
        depthArray = np.asarray(depths) + 100.0
    else:
        depthArray = np.asarray(depths)

    # Assume Poissonian-like relative uncertainty for N, i.e. sqrt(N) / N = 1 / sqrt(N)
    # Only used as an uncertainty estimate in weighting the fit, applying emphasis near Xmax
    # Not to be used as a true, physical uncertainty, on the fluctuations in N
    uncerts = 1.0 / np.sqrt(NprimeArray)

    # Insert a try-except statement for poor fits...
    try:
        if (shift == True) and (absoluteValue == False):
            popt, pcov = curve_fit(AndringaFunction, depthArray, NprimeArray, p0=[XmaxGuess+100.0, RGuess, LGuess], sigma=uncerts)
        elif (absoluteValue == True) and (shift == False):
            popt, pcov = curve_fit(AndringaFunctionWithABS, depthArray, NprimeArray, p0=[XmaxGuess, RGuess, LGuess], sigma=uncerts)
        else:
            popt, pcov = curve_fit(AndringaFunction, depthArray, NprimeArray, p0=[XmaxGuess, RGuess, LGuess], sigma=uncerts)
    except RuntimeError:
        return np.inf, np.inf, np.inf, np.inf, np.inf, np.inf

    perr = np.sqrt(np.diag(pcov))

    XmaxAndringaFit = popt[0]
    RAndringaFit = popt[1]
    LAndringaFit = popt[2]

    XmaxAndringaSigma = perr[0]
    RAndringaSigma = perr[1]
    LAndringaSigma = perr[2]

    if (shift == True) and (absoluteValue == False):
        XmaxAndringaFit = popt[0] - 100.0 # Shift Xmax from fit back to real xmax

    return XmaxAndringaFit, XmaxAndringaSigma, RAndringaFit, RAndringaSigma, LAndringaFit, LAndringaSigma, *popt

def remove_zeros(listToUpdate, pairedList):
    for i in reversed(range(len(listToUpdate))):
        if listToUpdate[i] == 0:
            del listToUpdate[i]
            del pairedList[i]
    return listToUpdate, pairedList



# depths, positrons, electrons, muPlus, muMinus, chargedParticles, xmax, Rcorsika, Lcorsika, x0, lambdaApprox = LongFileParser(args.input)

# totalEM = np.array(positrons) + np.array(electrons)
# ixmax = np.argmax(totalEM)

# if xmax > 1700:  # Sometimes corsika fits of the profile fail
#     xmax = depths[ixmax]
#     Rcorsika = -1
#     Lcorsika = -1

# # TO-DO: Insert a keyword that defines the observation level from either a settings file or command line argument
# ground = 870 / np.cos(args.zen)


zenith_dict = load_zeniths(args.zenfile)
results = {}


###########################  MAIN LOOP ###################################

for longfile in long_files:

     # skip empty files
    if os.path.getsize(longfile) == 0:
        print(f"WARNING: skipping empty file {longfile}")
        continue

    run = get_run_from_longfile(longfile)
    if run not in zenith_dict:
        print(f"WARNING: no zenith for run {run}, skipping")
        continue

    try:
        depths, positrons, electrons, muPlus, muMinus, chargedParticles, \
            xmax, Rcorsika, Lcorsika, x0, lambdaApprox = LongFileParser(longfile)
    except Exception as e:
        print(f"WARNING: skipping corrupted file {longfile} ({e})")
        continue

    zen = zenith_dict[run]
    
    totalEM = np.array(positrons) + np.array(electrons)
    ixmax = np.argmax(totalEM)


    ground = 870 / np.cos(zen)

    indGround = FindGroundIndex(ground, depths)
    
    muPlusGround = muPlus[indGround]
    muMinusGround = muMinus[indGround]
    
    positronsGround = positrons[indGround]
    electronsGround = electrons[indGround]
    
    # Remove final 20g/cm2 from fit because they are not physical
    # My understanding is some part of the shower front reaches ground which causes dip in particle numbers...
    if args.removeFinal20gcm2 == True:
        depthSpacing = depths[1] - depths[0]
        numPointsToRemove = int(20.0 / depthSpacing)
        del depths[-numPointsToRemove:]
        del positrons[-numPointsToRemove:]
        del electrons[-numPointsToRemove:]
        del muPlus[-numPointsToRemove:]
        del muMinus[-numPointsToRemove:]
        del chargedParticles[-numPointsToRemove:]
    
    Nmax = np.max(totalEM)
    NmaxGuess = Nmax
    # Prevent errors if xmax is near ground level (i.e. ixmax is last index)
    if ixmax >= len(depths):
        XmaxGuess = xmax
    else:
        XmaxGuess = depths[ixmax]
    X0Guess = 0
    lambGuess = 80.0
    
    totalEMUpdate = np.array(positrons) + np.array(electrons)
    totalEMList = totalEMUpdate.tolist()
    
    # Put removal of zeros after all index calls are done, to ensure all indices are the same
    totalEMList, depths = remove_zeros(totalEMList, depths)
    
    totalEMAfterCuts = np.array(totalEMList)
    
    NPrimeValsCORSIKA = totalEMAfterCuts / Nmax
    RGuess = np.sqrt(lambGuess / abs(X0Guess - XmaxGuess))
    LGuess = np.sqrt(abs(X0Guess - XmaxGuess) * lambGuess)
    
    RFitShift, RFitSigmaShift, LFitShift, LFitSigmaShift, XmaxFitShift, XmaxSigmaShift, NmaxFit, XmaxFit, X0Fit, lambFit  = FitLongitudinalProfile(depths, totalEMList, NmaxGuess, XmaxGuess, X0Guess, lambGuess, shift=False, absoluteValue=False)
    # XmaxAndringaFit, XmaxAndringaSigma, RAndringaFit, RAndringaSigma, LAndringaFit, LAndringaSigma,  Xmax, R, L = FitLongitudinalProfileAndringa(depths, NPrimeValsCORSIKA, XmaxGuess, RGuess, LGuess, shift=False, absoluteValue=False) #parameterized


    # results[run] = (
    #     run, 
    #     xmax,
    #     round(muPlusGround + muMinusGround),
    #     totalEM[ixmax],
    #     round(positronsGround + electronsGround),
    #     Rcorsika,
    #     Lcorsika,
    #     RFitShift,
    #     RFitSigmaShift,
    #     LFitShift,
    #     LFitSigmaShift,
    #     XmaxFitShift,
    #     XmaxSigmaShift,
    #     RAndringaFit,
    #     RAndringaSigma,
    #     LAndringaFit,
    #     LAndringaSigma,
    #     XmaxAndringaFit,
    #     XmaxAndringaSigma
    # )

    results[run] = (
        run, 
        xmax,
        round(muPlusGround + muMinusGround),
        totalEM[ixmax],
        round(positronsGround + electronsGround),
        NmaxFit, XmaxFit, X0Fit, lambFit
    )

    # results[run] = (
    #     run, 
    #     xmax,
    #     round(muPlusGround + muMinusGround),
    #     totalEM[ixmax],
    #     round(positronsGround + electronsGround),
    #     Xmax, R, L
    # )

output_dir = os.path.dirname(args.outfile)
if output_dir and not os.path.exists(output_dir):
    print(f"Creating output directory: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

with open(args.outfile, "w") as f:
    # Header
    # f.write(
    #     "# run xmax mupm_ground totalEM_Xmax epm_ground "
    #     "Rcorsika Lcorsika "
    #     "RFitShift RFitSigmaShift "
    #     "LFitShift LFitSigmaShift "
    #     "XmaxFitShift XmaxSigmaShift "
    #     "RAndringaFit RAndringaSigma "
    #     "LAndringaFit LAndringaSigma "
    #     "XmaxAndringaFit XmaxAndringaSigma\n"
    # )

    f.write(
        "# run xmax mupm_ground totalEM_Xmax epm_ground "
        "NmaxFit, XmaxFit, X0Fit, lambFits \n"
    )

    # Write runs in order 0 → 199
    for run in range(200):
        if run not in results:
            continue   # skip missing runs safely

        line = " ".join(f"{val:.6g}" for val in results[run])
        f.write(line + "\n")
print(f"File write successful: {os.path.exists(args.outfile)}")

