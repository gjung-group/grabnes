"""Cases of model_survey.py: BASES (input texts), CASES (name -> (base, changes[, structure])), STRUCTURES.

A case named like its base is the base itself; it must come before the cases derived from it.
Every switch of the solver appears in at least one case. The companion values are small,
arbitrary numbers chosen so that the switch has a visible effect; they are not physical
recommendations.
"""

import math
import os

T, F = ".true.", ".false."

COMMON = """
SuperCell 1
Neigh.fastNNnotsquare .true.
Neigh.LayerNeighbors 100
Neigh.LayerDistFactor 6.2
TB.NeighLevels 5
"""

BILAYER = """
TypeOfBL Koshino
MoireCellParameters 3 2 2 3
twoLayers .true.
twoLayersZ1 18.33
twoLayersZ2 21.67
InterlayerDistance 3.34
CellHeight 40.0
""" + COMMON

# effective graphene/hBN moire model on an 11 x 11 graphene cell (one moire period)
EFFECTIVE = """
TypeOfSystem Graphene
CellSize 11
MoirePotential .true.
MoireJeil .true.
MoirePotC0   0-0.01013
MoirePotCz   0-0.00901
MoirePotCab   0.01134
MoirePotPhi0  1.510233401750693
MoirePotPhiz  0.147131255943122
MoirePotPhiab  0.342084533390889
Latticepercent 0-0.0909090909090909
""" + COMMON

BASES = {
    "graphene": "TypeOfSystem Graphene\nCellSize 6\n" + COMMON,
    "hbn": "TypeOfSystem BoronNitride\nCellSize 6\n" + COMMON,
    "tbg": "TypeOfSystem TwistedBilayerBasedOnMoireCell\n" + BILAYER,
    "eff": EFFECTIVE,
}

CASES = {}
STRUCTURES = {}


def case(name, base, structure=None, **changes):
    ch = {k.replace("__", "."): v for k, v in changes.items()}
    CASES[name] = (base, ch, structure) if structure else (base, ch)


def _with(template, changes):
    return template + "".join(f"{k.replace('__', '.')} {v}\n" for k, v in changes.items())


def sys_case(name, typeofsystem, template=BILAYER, structure=None, **changes):
    """A geometry of its own: registered as a base (its changes included) and as its own case."""
    BASES[name] = _with(f"TypeOfSystem {typeofsystem}\n" + template, changes)
    case(name, name, structure)


def toggles(base, flags, structure=None, **extra):
    """One case per switch: the switch set against its default, plus the companion values."""
    for flag in flags:
        value = T
        if flag.endswith("=F"):
            flag, value = flag[:-2], F
        ch = dict(extra)
        ch[flag.replace(".", "__")] = value
        case(f"{base}+{flag}", base, structure, **ch)


# ------------------------------------------------------------------ structures for ReadXYZ
A_G = 2.46
Z0, DZ = 10.0, 3.35


def stack(layers, mn=(3, 2), cell_height=None, relax=False):
    """Write generate.xyz, layerIndex.dat and sublatticesSorted.dat for a stack of honeycomb layers.

    layers: list of (kind, twisted, shift) with kind 'C' or 'BN'; a twisted layer is rotated to the
    other commensurate orientation of the (m, n) cell; shift 1 moves the layer by one bond (Bernal).
    With relax the atoms of generate.xyz are displaced by a smooth, cell-periodic field (a stand-in for a
    relaxed structure) while generateInit.xyz keeps the rigid positions. The auxiliary tables that some
    switches read (interlayerDistances.dat, displacements.txt) are always written."""
    m, n = mn
    h = cell_height or (2 * Z0 + DZ * (len(layers) - 1))

    def lattice(mm, nn):
        a1 = (A_G, 0.0)
        a2 = (A_G / 2, A_G * math.sqrt(3) / 2)
        L = (mm * a1[0] + nn * a2[0], mm * a1[1] + nn * a2[1])
        return a1, a2, math.atan2(L[1], L[0])

    a1, a2, ang0 = lattice(m, n)
    L1 = (m * a1[0] + n * a2[0], m * a1[1] + n * a2[1])
    L2 = (L1[0] * 0.5 - L1[1] * math.sqrt(3) / 2, L1[0] * math.sqrt(3) / 2 + L1[1] * 0.5)
    det = L1[0] * L2[1] - L1[1] * L2[0]
    atoms = []
    for il, (kind, twisted, shift) in enumerate(layers):
        _, _, ang = lattice(n, m) if twisted else (None, None, ang0)
        rot = ang0 - ang
        c, s = math.cos(rot), math.sin(rot)
        b1 = (c * a1[0] - s * a1[1], s * a1[0] + c * a1[1])
        b2 = (c * a2[0] - s * a2[1], s * a2[0] + c * a2[1])
        basis = [(0.0, 0.0), ((b1[0] + b2[0]) / 3, (b1[1] + b2[1]) / 3)]
        sh = ((b1[0] + b2[0]) / 3 * shift, (b1[1] + b2[1]) / 3 * shift)
        seen = set()
        rng = 3 * (m + n)
        for i in range(-rng, rng + 1):
            for j in range(-rng, rng + 1):
                for ib, b in enumerate(basis):
                    x = i * b1[0] + j * b2[0] + b[0] + sh[0]
                    y = i * b1[1] + j * b2[1] + b[1] + sh[1]
                    f1 = (x * L2[1] - y * L2[0]) / det
                    f2 = (-x * L1[1] + y * L1[0]) / det
                    f1r, f2r = round(f1 % 1.0, 6) % 1.0, round(f2 % 1.0, 6) % 1.0
                    if (f1r, f2r) in seen or not (-1e-7 <= f1 < 1 - 1e-7 and -1e-7 <= f2 < 1 - 1e-7):
                        continue
                    seen.add((f1r, f2r))
                    species = ib + 1 if kind == "C" else ib + 3
                    atoms.append((kind if kind == "C" else ("B" if ib == 0 else "N"), x, y, Z0 + DZ * il, il + 1, species,
                                  f1 % 1.0, f2 % 1.0))
    expect = 2 * (m * m + m * n + n * n) * len(layers)
    assert len(atoms) == expect, (len(atoms), expect)

    def moved(a):
        """Smooth periodic displacement: 0.03 A in plane, 0.15 A out of plane."""
        el, x, y, z, layer, species, f1, f2 = a
        t1, t2 = 2 * math.pi * f1, 2 * math.pi * f2
        sign = 1 if layer % 2 else -1
        return (x + 0.03 * math.sin(t1) * sign, y + 0.03 * math.cos(t2) * sign, z + 0.15 * math.cos(t1 + t2) * sign)

    def write(d):
        header = f"{L1[0]:.10f} {L1[1]:.10f} 0.0\n{L2[0]:.10f} {L2[1]:.10f} 0.0\n0.0 0.0 {h:.4f}\n{len(atoms)}\n"
        with open(os.path.join(d, "generateInit.xyz"), "w") as f:
            f.write(header + "".join(f"{a[0]} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in atoms))
        with open(os.path.join(d, "generate.xyz"), "w") as f:
            f.write(header)
            for a in atoms:
                x, y, z = moved(a) if relax else a[1:4]
                f.write(f"{a[0]} {x:.10f} {y:.10f} {z:.10f}\n")
        with open(os.path.join(d, "layerIndex.dat"), "w") as f:
            f.write("".join(f"{a[4]}\n" for a in atoms))
        with open(os.path.join(d, "sublatticesSorted.dat"), "w") as f:
            f.write("".join(f"{i + 1} {a[5]}\n" for i, a in enumerate(atoms)))
        with open(os.path.join(d, "interlayerDistances.dat"), "w") as f:
            f.write("".join(f"{DZ + 0.05 * math.sin(2 * math.pi * a[6]):.8f}\n" for a in atoms))
        # stacking displacement of every atom with respect to the neighbouring layer(s): x, y, distance, twice
        # (the encapsulated-trilayer reader takes six columns, the others the first three)
        with open(os.path.join(d, "displacements.txt"), "w") as f:
            for a in atoms:
                dx, dy = 0.8 * math.cos(2 * math.pi * a[6]), 0.8 * math.sin(2 * math.pi * a[7])
                f.write(f"{dx:.8f} {dy:.8f} {DZ:.8f} {-dy:.8f} {dx:.8f} {DZ:.8f}\n")
    return write


C, CT, CB, BN, BNT = ("C", False, 0), ("C", True, 0), ("C", False, 1), ("BN", False, 0), ("BN", True, 0)
STRUCTURES.update({
    "x1": stack([C]),
    "x2": stack([C, CT]),
    "x2gbn": stack([C, BNT]),
    "x2bnbn": stack([BN, BNT]),
    "x3": stack([C, CT, C]),
    "x3enc": stack([BNT, C, BNT]),
    "x4": stack([C, CT, CT, C]),
    "x4b": stack([C, CB, CT, CT]),
    "x4enc": stack([BNT, C, CB, BNT]),
    "x5": stack([C, CT, C, CT, C]),
    "x6": stack([C, CT] * 3),
    "x7": stack([C, CT, C, CT, C, CT, C]),
    "x8": stack([C, CT] * 4),
    "x10": stack([C, CT] * 5),
    "x20": stack([C, CT] * 10),
    "x2r": stack([C, CT], relax=True),
    "x2gbn_r": stack([C, BNT], relax=True),
    "x3enc_r": stack([BNT, C, BNT], relax=True),
})

XYZ = """
TypeOfSystem ReadXYZ
TypeOfBL Koshino
readLayerIndex .true.
""" + COMMON


def xyz_case(name, structure, **changes):
    BASES[name] = _with(XYZ, changes)
    case(name, name, structure)


# ================================================================== A. geometries (TypeOfSystem)
case("graphene", "graphene")
case("hbn", "hbn")
case("tbg", "tbg")
case("eff", "eff")
sys_case("sys_TwistedBilayer", "TwistedBilayer", CellSize="6")
sys_case("sys_TwistedBilayerRectangular", "TwistedBilayerBasedOnMoireCellRectangular")
sys_case("sys_Trilayer", "TrilayerBasedOnMoireCell")
sys_case("sys_MoireEncapsulatedBilayer", "MoireEncapsulatedBilayer", CellSize="6")
sys_case("sys_MoireEncapsulatedBilayerMC", "MoireEncapsulatedBilayerBasedOnMoireCell")
sys_case("sys_Graphene_Over_BN", "Graphene_Over_BN")
sys_case("sys_BilayerGraphene_Over_BN", "BilayerGraphene_Over_BN")
sys_case("sys_BLtoSLYoungju", "BLtoSLYoungju", template="CellSize 6\n" + COMMON)
sys_case("sys_Ribbons_zigzag", "Ribbons", template="RibbonType Zigzag\nCellSize 6\n" + COMMON)
sys_case("sys_Ribbons_armchair", "Ribbons", template="RibbonType Armchair\nCellSize 6\n" + COMMON)
sys_case("sys_Hybrid", "Hybrid", template="CellSize 6\n" + COMMON)
case("graphene+SuperCell2", "graphene", SuperCell="2")
case("graphene+SuperCellAsymmetric", "graphene", SuperCellAsymmetric=T, SuperCellX="2", SuperCellY="1")
case("tbg+basedOnMoireCellParameters", "tbg", basedOnMoireCellParameters=T)

# ================================================================== B. layer stacks read from a file
xyz_case("xyz1", "x1", oneLayer=T)
xyz_case("xyz2", "x2", twoLayers=T)
xyz_case("xyz2_gbn", "x2gbn", GBNtwoLayers=T, **{"TB.Hopping": "3.5"})
toggles("xyz2_gbn", ["GBNtwoLayersF2G2s", "GBNOffDiag", "GBNuseHarmonicApprox", "dontUseInplaneMoire", "MoirePotential",
                     "realStrain", "shellsFromRigidPositions", "renormalizeCoupling", "useBNGKaxiras", "useBNGSrivani"],
        "x2gbn", couplingFactor="0.5")
xyz_case("xyz2_bnbn", "x2bnbn", BNBNtwoLayers=T)
toggles("xyz2_bnbn", ["BNBNDiag"], "x2bnbn")
xyz_case("xyz3", "x3", threeLayers=T, middleTwist=T)
toggles("xyz3", ["middleTwist=F", "threeLayerShort", "MoireTrilayer", "TrilayerFanZhang"], "x3")
case("xyz3+bilayerF2G2", "xyz3", "x3", middleTwist=F, forceBilayerF2G2Intralayer=T)
xyz_case("xyz3_enc", "x3enc", encapsulatedThreeLayers=T)
toggles("xyz3_enc", ["GBNOffDiag", "removeTopMoireInL2"], "x3enc")
xyz_case("xyz4", "x4", fourLayers=T, forceBilayerF2G2Intralayer=T)
case("xyz4+forceBilayerF2G2Intralayer=F", "xyz4", "x4", forceBilayerF2G2Intralayer=F)
xyz_case("xyz4_sandwiched", "x4", fourLayersSandwiched=T, middleTwist=T)
case("xyz4_sandwiched+middleTwist=F", "xyz4_sandwiched", "x4", middleTwist=F)
case("xyz4_sandwiched+bilayerF2G2", "xyz4_sandwiched", "x4", middleTwist=F, forceBilayerF2G2Intralayer=T)
toggles("xyz4_sandwiched", ["fourLayerOnsiteShifts", "differentCouplings", "MoireBilayerElectricField", "helicalTwistedMBM",
                            "deactivateInterlayer12", "deactivateInterlayer23", "deactivateInterlayer34"], "x4",
        fourLayerShift1="0.01", fourLayerShift2="0.02", fourLayerShift3="0-0.02", fourLayerShift4="0-0.01",
        renormalizeCoupling=T, couplingFactor="0.9", couplingFactor2="0.5", MoireBilayerElectricShift="0.05")
xyz_case("xyz4_helical", "x4b", helicalTwistedMBM=T, fourLayersSandwiched=T, middleTwist=T)
toggles("xyz4_helical", ["helicalTwistedMBM_CDW"], "x4b", CDWAmplitude="0.02")
case("xyz4_helical+CDWUseMassTerm", "xyz4_helical", "x4b", helicalTwistedMBM_CDW=T, CDWUseMassTerm=T, CDWAmplitude="0.02")
xyz_case("xyz4_enc", "x4enc", encapsulatedFourLayers=T)
toggles("xyz4_enc", ["encapsulatedFourLayersF2G2", "GBNOffDiag"], "x4enc")
xyz_case("xyz4_t2GBN", "x4enc", t2GBN=T)
xyz_case("xyz4_BNt2GBN", "x4enc", BNt2GBN=T)
xyz_case("xyz4_t3GwithBN", "x4enc", t3GwithBN=T)
xyz_case("xyz3_t2BG", "x3", t2BG=T)
xyz_case("xyz4_t3BG", "x4", t3BG=T)
xyz_case("xyz5", "x5", fiveLayersSandwiched=T, middleTwist=T)
xyz_case("xyz5_enc", "x5", encapsulatedFiveLayers=T)
xyz_case("xyz6", "x6", sixLayersSandwiched=T, middleTwist=T)
xyz_case("xyz6_enc", "x6", encapsulatedSixLayers=T)
xyz_case("xyz7", "x7", sevenLayersSandwiched=T, middleTwist=T)
xyz_case("xyz7_enc", "x7", encapsulatedSevenLayers=T)
xyz_case("xyz8", "x8", eightLayersSandwiched=T, middleTwist=T)
xyz_case("xyz10", "x10", tenLayersSandwiched=T, middleTwist=T)
xyz_case("xyz20", "x20", twentyLayersSandwiched=T, middleTwist=T)
toggles("xyz2", ["useSublatticeFile=F", "readLayerIndex=F", "readRigidXYZ", "readInterlayerDistances", "tBGuseDisplacementFile",
                 "GBNuseDisplacementFile", "invertDisplacements", "singleLayerXYZ", "BernalReadXYZ", "AtomsOrderDeactivated=F",
                 "useLayerSpecificOnsiteEnergyTerms", "shellsFromRigidPositions"], "x2",
        twoLayersZ1="10.0", twoLayersZ2="13.35")
toggles("tbg", ["Bulk", "BulkSmall", "nonBulkSmall", "twistedBLAddShift", "BernalShift", "bridgeShift",
                "createBLDomainBoundary"])
case("tbg+BLDomainBoundaryTypeArmAA", "tbg", createBLDomainBoundary=T, BLDomainBoundaryTypeArmAA=T)
case("tbg+BLDomainBoundaryTypeArmSP", "tbg", createBLDomainBoundary=T, BLDomainBoundaryTypeArmSP=T)
case("sys_Trilayer+TrilayerAddShift", "sys_Trilayer", TrilayerAddShift=T, TrilayerShiftFactor="1")

# ================================================================== C. interlayer models (TypeOfBL)
for bl in ["Jeil", "BLKaxiras", "BLSrivani", "Srivani", "Mayou", "HTC", "None"]:
    case(f"tbg+TypeOfBL_{bl}", "tbg", TypeOfBL=bl)
case("graphene+TypeOfSL_HTC", "graphene", TypeOfSL="HTC")
toggles("tbg", ["KoshinoSR", "corrugatedInterlayerTwoCenter", "onlyvppsigma", "renormalizeCoupling", "renormalizeHoppings",
                "deactivateInterlayer", "deactivateInterlayerBG", "deactivateInterlayerTwisted", "setnLayersToZero",
                "addSecondLayerInteractions"], couplingFactor="0.5")
for bl in ["Jeil", "BLKaxiras", "BLSrivani"]:
    for flag in ["BilayerOneParameter", "BilayerThreeParameters", "useOnlyVAB", "onlyV0", "switchV3Sign", "deactivateV3",
                 "deactivateV6", "oldParameterSet", "newFittingFunctions", "sublatticeDependent", "sublatticeIndependent",
                 "useTheta", "useThetaIJ", "findThetasGeometrically", "addExponentialDecayForDihedral", "oppositedxdy",
                 "addPressureDependence=F", "changeLatticeParameterForSrivaniModel", "useBNGKaxiras", "useBNGSrivani"]:
        value = F if flag.endswith("=F") else T
        case(f"tbg+{bl}+{flag}", "tbg", **{"TypeOfBL": bl, flag.replace("=F", ""): value})
xyz_case("xyz6_t3BG", "x6", t3BG=T)
toggles("xyz6_t3BG", ["deactivateInterlayert3BG1to2", "deactivateInterlayert3BG2to3", "deactivateInterlayert3BG3to4",
                      "deactivateInterlayert3BG4to5", "deactivateInterlayert3BG5to6"], "x6")
toggles("xyz4_t2GBN", ["deactivateInterlayert2GBN1to2", "deactivateInterlayert2GBN2to3"], "x4enc")
toggles("xyz3_t2BG", ["deactivateInterlayert2BG2to3"], "x3")

# ================================================================== D. intralayer models
for n in range(1, 9):
    case(f"graphene+NeighLevels{n}", "graphene", **{"TB.NeighLevels": str(n)})
toggles("graphene", ["F2G2Model=F", "KoshinoIntralayer", "MayouIntralayer", "useOldGrapheneF2G2", "removeF2G2Flag",
                     "forceBilayerF2G2Intralayer", "realStrain", "strainedMoire", "MoireStrain", "RandomStrain",
                     "realisticBubbles", "BigBubble", "manyBubbles", "bubbleInPlaneStrain", "bubbleGaussian",
                     "deactivateIntrasublattice", "deactivateIntersublattice", "deactivateIntraSublatticeForC",
                     "deactivateASubLattice", "deactivateBSubLattice", "Neigh.CutAtNN3"])
case("graphene+F2G2off+NeighLevels8", "graphene", F2G2Model=F, **{"TB.NeighLevels": "8"})
case("graphene+realStrain+onlyFirstNeighbor", "graphene", realStrain=T, onlyFirstNeighborRealStrain=T,
     realStrainReferenceLatticeConstant="2.44")
case("graphene+realStrain+reference", "graphene", realStrain=T, realStrainReferenceLatticeConstant="2.44")
case("graphene+realStrain+beta", "graphene", realStrain=T, realStrainReferenceLatticeConstant="2.44", realStrainBeta="2.0")
case("graphene+realStrain+periodicStrain", "graphene", realStrain=T, periodicStrain=T)
toggles("tbg", ["KoshinoIntralayer", "MayouIntralayer", "F2G2Model=F", "realStrain", "forceBilayerF2G2Intralayer"])
case("graphene+IntralayerRadius", "graphene", **{"Neigh.IntralayerRadius": "5.1134"})

# ================================================================== E. effective moire models
toggles("eff", ["MoireOffDiag", "MoireOnlyH0", "MoireOnlyHZ", "MoireNoH0AndHZ", "MoireH0AndHZ", "MoireSymmetricPot",
                "switchHzjj", "MoireKekule", "MoireTrilayer", "TrilayerFanZhang", "distanceDependentEffectiveModel",
                "addDisplacements", "sublatticeBasis", "MoireBilayerElectricField", "MoireBLDeactivateUpperLayer",
                "MoiretDBLDeactivateUpperLayers", "MoireTwisted", "useLayerSpecificOnsiteEnergyTerms",
                "addSublatticeMassterm", "onlyBottomLayerMassTerm", "addOnsiteEnergyShift"], MoireTwistAngle="0.5")
case("eff+MoireOffDiagMidpoint", "eff", MoireOffDiag=T, MoireOffDiagMidpoint=T)
case("eff+MoireSachs", "eff", MoireJeil=F, MoireSachs=T)
case("eff+MoireAddSecondMoire", "eff", MoireAddSecondMoire=T, MoireOffDiag=T)
case("eff+MoireAddSecondMoire+Midpoint", "eff", MoireAddSecondMoire=T, MoireOffDiag=T, MoireOffDiagMidpoint=T)
case("eff+MoireAddSecondMoire+Twisted2", "eff", MoireAddSecondMoire=T, MoireTwisted2=T, MoireTwistAngle2="0.5")
case("eff+MoireSecondMoireRotateFirst_off", "eff", MoireAddSecondMoire=T, MoireSecondMoireRotateFirst=F)
case("eff+MoireBilayerElectricFieldInvert", "eff", MoireBilayerElectricField=T, MoireBilayerElectricFieldInvert=T)
toggles("tbg", ["tBGDiag", "tBGDiagPRB", "tBGOffDiag", "tBGSwitchDxDy", "MoireBilayerElectricField", "addSublatticeMassterm",
                "onlyBottomLayerMassTerm", "MoirePotential"])
case("tbg+tBGOffDiagPRB", "tbg", tBGOffDiag=T, tBGOffDiagPRB=T)
toggles("sys_MoireEncapsulatedBilayerMC", ["MoirePotential", "MoireEncapsulatedBilayer", "MoireBLDeactivateUpperLayer",
                                           "MoiretDBLDeactivateUpperLayers", "MoireBilayerElectricField"], MoireJeil=T)
case("sys_MoireEncapsulatedBilayerMC+MoireOffDiag", "sys_MoireEncapsulatedBilayerMC", MoirePotential=T, MoireJeil=T,
     MoireOffDiag=T)

# ================================================================== F. artificial potentials and disorder
case("graphene+moireCDW", "graphene", moireCDW=T, **{"moireCDW.Amplitude": "0.02", "moireCDW.Denominator": "1",
     "&q": "&begin moireCDW.Qvectors 1\n1 0\n&end moireCDW.Qvectors"})
case("graphene+moireCDW_mass", "graphene", moireCDW=T, **{"moireCDW.MassAmplitude": "0.02", "moireCDW.Denominator": "1",
     "&q": "&begin moireCDW.Qvectors 2\n1 0\n0 1\n&end moireCDW.Qvectors"})
case("graphene+moireCDW_incommensurate", "graphene", moireCDW=T, **{"moireCDW.Amplitude": "0.02",
     "moireCDW.Denominator": "5", "&q": "&begin moireCDW.Qvectors 1\n1 0\n&end moireCDW.Qvectors"})
case("graphene+helicalTwistedMBM_CDW", "graphene", helicalTwistedMBM_CDW=T, CDWAmplitude="0.02")
toggles("graphene", ["sinusModulation", "sinusModulationUsingPeriod", "sinusModulationUsingPeriodYDirection",
                     "cosinusModulationUsingPeriod", "cosinusModulationUsingPeriodYDirection", "SquareFunction",
                     "SquareFunction2", "SquareChecker2219", "TwoDimensional", "TwoDimension", "ArmChairShape", "AddZTerm",
                     "Zterm1D", "Zterm1DKink", "PNP", "PNPKink", "Bubbles", "printBubble", "Anderson", "GaussDisorder",
                     "deltaDisorder", "SublatticeDisorder", "addSublatticeMassterm", "addOnsiteEnergyShift"],
        sinusModulationPeriod="7.38", cosinusModulationPeriod="7.38", sinusFactor="0.01", cosinusFactor="0.01",
        AndersonAmp="0.1", onsiteShift="0.05", checkerDivider="2")
case("graphene+sinusModulationAddMassTerm", "graphene", sinusModulationUsingPeriod=T, sinusModulationAddMassTerm=T,
     sinusModulationPeriod="7.38")
case("graphene+cosinusModulationAddMassTerm", "graphene", cosinusModulationUsingPeriod=T, cosinusModulationAddMassTerm=T,
     cosinusModulationPeriod="7.38")

# ================================================================== G. spin, topology, field
toggles("graphene", ["ZeemanTerm", "PseudoZeemanTerm", "SpinPolarized", "IntrinsicSOCterm", "IsingSOCterm", "RashbaSOCterm",
                     "PIASOCterm", "FrankMagneticField", "EnableSCF=F"],
        LambdaI="0.01", LambdaIsing="0.01", LambdaR="0.01", LambdaPIA="0.01")
case("graphene+SOCLayerControl", "graphene", IntrinsicSOCterm=T, LambdaI="0.01", SOCLayerControl=T)
case("graphene+HaldaneNNN", "graphene", HaldaneNNN=T, HaldaneT2="0.05")
for flag in ["HaldaneSpecifyFlux", "HaldaneSpecifyPhase", "HaldaneSpecifyRange", "paperOrientation", "HaldaneOppositePhase",
             "HaldaneBothLayers", "HaldaneSetFluxQ", "HaldaneLayerControl"]:
    case(f"graphene+HaldaneNNN+{flag}", "graphene", HaldaneNNN=T, HaldaneT2="0.05", HaldanePhase="0.3", **{flag: T})
case("graphene+MagField1", "graphene", **{"MagField.Integer": "1"})
case("tbg+MagField1", "tbg", **{"MagField.Integer": "1"})
case("tbg+HaldaneBothLayers", "tbg", HaldaneNNN=T, HaldaneT2="0.05", HaldaneBothLayers=T)
toggles("tbg", ["IsingSOCterm", "RashbaSOCterm"], LambdaIsing="0.01", LambdaR="0.01")

# ================================================================== I. neighbour searches and table files
# NeighList does not write the translation table, so it is compared through the eigenvalues: the same
# bands as the default search are required (identical fingerprints, hence INERT)
CASES["graphene@bands"] = ("graphene", {"@mode": "bands", "TB.NeighLevels": "3"})
CASES["graphene@bands+NeighList"] = ("graphene", {"@mode": "bands", "TB.NeighLevels": "3", "Neigh.fastNNnotsquare": F,
                                                    "Neigh.LayerNeighbors": None, "Neigh.LayerDistFactor": None})
# with an interlayer search requested on a monolayer, NeighList builds a different Hamiltonian
CASES["graphene@bands+NeighList+LayerNeighbors"] = ("graphene", {"@mode": "bands", "TB.NeighLevels": "3",
                                                                   "Neigh.fastNNnotsquare": F})
case("graphene+fastNN", "graphene", **{"Neigh.fastNN": T})
case("graphene+fastNNnotsquareNotRectangle", "graphene", **{"Neigh.fastNNnotsquare": F, "Neigh.fastNNnotsquareNotRectangle": T})

# ------------------------------------------------------------------ second-round variants
case("tbg+TypeOfBL_Jeil+NeighLevels1", "tbg", TypeOfBL="Jeil", **{"TB.NeighLevels": "1"})
case("tbg+TypeOfBL_BLKaxiras+NeighLevels1", "tbg", TypeOfBL="BLKaxiras", **{"TB.NeighLevels": "1"})
case("sys_Trilayer+TypeOfBL_Jeil", "sys_Trilayer", TypeOfBL="Jeil")
case("sys_MoireEncapsulatedBilayerMC+TypeOfBL_Jeil", "sys_MoireEncapsulatedBilayerMC", TypeOfBL="Jeil")
case("sys_MoireEncapsulatedBilayerMC+TypeOfBL_BLKaxiras", "sys_MoireEncapsulatedBilayerMC", TypeOfBL="BLKaxiras")

# ------------------------------------------------------------------ displaced ("relaxed") structures and tables
xyz_case("xyz2_relaxed", "x2r", twoLayers=T)
toggles("xyz2_relaxed", ["realStrain", "readRigidXYZ", "readInterlayerDistances", "tBGuseDisplacementFile",
                         "corrugatedInterlayerTwoCenter", "KoshinoSR"], "x2r", realStrainReferenceLatticeConstant="2.44")
case("xyz2_relaxed+tBGOffDiag+file", "xyz2_relaxed", "x2r", tBGOffDiag=T, tBGuseDisplacementFile=T)
case("xyz2_relaxed+tBGDiag+file", "xyz2_relaxed", "x2r", tBGDiag=T, tBGuseDisplacementFile=T)
case("xyz2_relaxed+tBGOffDiag+invert", "xyz2_relaxed", "x2r", tBGOffDiag=T, tBGuseDisplacementFile=T, invertDisplacements=T)
xyz_case("xyz2_gbn_relaxed", "x2gbn_r", GBNtwoLayers=T, **{"TB.Hopping": "3.5"})
toggles("xyz2_gbn_relaxed", ["realStrain", "shellsFromRigidPositions", "GBNOffDiag", "GBNuseDisplacementFile",
                             "readInterlayerDistances", "renormalizeCoupling"], "x2gbn_r",
        realStrainReferenceLatticeConstant="2.44", realStrainReferenceLatticeConstantBN="2.44", couplingFactor="0.5")
case("xyz2_gbn_relaxed+realStrain+shellsFromRigid", "xyz2_gbn_relaxed", "x2gbn_r", realStrain=T, shellsFromRigidPositions=T,
     realStrainReferenceLatticeConstant="2.44", realStrainReferenceLatticeConstantBN="2.44")
case("xyz2_gbn_relaxed+GBNOffDiag+file", "xyz2_gbn_relaxed", "x2gbn_r", GBNOffDiag=T, GBNuseDisplacementFile=T)
case("xyz2_gbn_relaxed+GBNOffDiag+file+harmonic", "xyz2_gbn_relaxed", "x2gbn_r", GBNOffDiag=T, GBNuseDisplacementFile=T,
     GBNuseHarmonicApprox=T)
case("xyz2_gbn_relaxed+IntralayerRadius", "xyz2_gbn_relaxed", "x2gbn_r", **{"Neigh.IntralayerRadius": "5.1134"})
xyz_case("xyz3_enc_relaxed", "x3enc_r", encapsulatedThreeLayers=T)
case("xyz3_enc_relaxed+GBNOffDiag+file", "xyz3_enc_relaxed", "x3enc_r", GBNOffDiag=T, GBNuseDisplacementFile=T)

# ------------------------------------------------------------------ terms seen only in the eigenvalues
# (spin and spin-orbit terms are added when the Hamiltonian matrix is built, not in the hopping table)
BASES["graphene_small"] = "TypeOfSystem Graphene\nCellSize 2\n" + COMMON
BASES["tbg_small"] = BASES["tbg"]
for b in ("graphene_small", "tbg_small"):
    CASES[f"{b}@bands"] = (b, {"@mode": "bands"})
    for flag, extra in [("ZeemanTerm", {}), ("PseudoZeemanTerm", {}), ("SpinPolarized", {}),
                        ("IntrinsicSOCterm", {"LambdaI": "0.01"}), ("IsingSOCterm", {"LambdaIsing": "0.01"}),
                        ("RashbaSOCterm", {"LambdaR": "0.01"}), ("PIASOCterm", {"LambdaPIA": "0.01"}),
                        ("EnableSCF", None), ("HaldaneNNN", {"HaldaneT2": "0.05"}), ("MagField.Integer", None)]:
        ch = {"@mode": "bands"}
        if flag == "EnableSCF":
            ch.update(SpinPolarized=T, EnableSCF=F)
        elif flag == "MagField.Integer":
            ch["MagField.Integer"] = "1"
        else:
            ch[flag] = T
            ch.update(extra)
        CASES[f"{b}@bands+{flag}"] = (b, ch)
    CASES[f"{b}@bands+SOCLayerControl"] = (b, {"@mode": "bands", "IsingSOCterm": T, "LambdaIsing": "0.01",
                                              "SOCLayerControl": T, "SOCLayers": "1"})
    CASES[f"{b}@bands+Zeeman+Spin-1"] = (b, {"@mode": "bands", "ZeemanTerm": T, "Spin": "0-1"})

# ------------------------------------------------------------------ switches placed where they act
# (found by trying every switch that was inert on its first base on every other base)
case('xyz1+SuperCellAsymmetric', 'xyz1', 'x1', **{'SuperCellAsymmetric': '.true.', 'SuperCellX': '2', 'SuperCellY': '1'})
case('xyz3+basedOnMoireCellParameters', 'xyz3', 'x3', **{'basedOnMoireCellParameters': '.true.'})
case('tbg+useBNGKaxiras', 'tbg', None, **{'couplingFactor': '0.5', 'useBNGKaxiras': '.true.'})
case('xyz2_gbn+threeLayerShort', 'xyz2_gbn', 'x2gbn', **{'threeLayerShort': '.true.'})
case('tbg+useSublatticeFile', 'tbg', None, **{'twoLayersZ1': '10.0', 'twoLayersZ2': '13.35', 'useSublatticeFile': '.false.'})
case('tbg+readLayerIndex', 'tbg', None, **{'twoLayersZ1': '10.0', 'twoLayersZ2': '13.35', 'readLayerIndex': '.false.'})
case('tbg+readRigidXYZ', 'tbg', None, **{'twoLayersZ1': '10.0', 'twoLayersZ2': '13.35', 'readRigidXYZ': '.true.'})
case('tbg+invertDisplacements', 'tbg', None, **{'twoLayersZ1': '10.0', 'twoLayersZ2': '13.35', 'invertDisplacements': '.true.'})
case('tbg+useLayerSpecificOnsiteEnergyTerms', 'tbg', None, **{'twoLayersZ1': '10.0', 'twoLayersZ2': '13.35', 'useLayerSpecificOnsiteEnergyTerms': '.true.'})
case('eff+TypeOfBL_HTC', 'eff', None, **{'TypeOfBL': 'HTC'})
case('eff+TypeOfSL_HTC', 'eff', None, **{'TypeOfSL': 'HTC'})
case('xyz3+corrugatedInterlayerTwoCenter', 'xyz3', 'x3', **{'couplingFactor': '0.5', 'corrugatedInterlayerTwoCenter': '.true.'})
case('xyz3_enc+deactivateInterlayerBG', 'xyz3_enc', 'x3enc', **{'couplingFactor': '0.5', 'deactivateInterlayerBG': '.true.'})
case('xyz3+deactivateInterlayerTwisted', 'xyz3', 'x3', **{'couplingFactor': '0.5', 'deactivateInterlayerTwisted': '.true.'})
case('xyz3+addSecondLayerInteractions', 'xyz3', 'x3', **{'couplingFactor': '0.5', 'addSecondLayerInteractions': '.true.'})
case('tbg+useOldGrapheneF2G2', 'tbg', None, **{'useOldGrapheneF2G2': '.true.'})
case('sys_MoireEncapsulatedBilayer+removeF2G2Flag', 'sys_MoireEncapsulatedBilayer', None, **{'removeF2G2Flag': '.true.'})
case('sys_Graphene_Over_BN+realStrain', 'sys_Graphene_Over_BN', None, **{'realStrain': '.true.'})
case('tbg+strainedMoire', 'tbg', None, **{'strainedMoire': '.true.'})
case('tbg+deactivateIntrasublattice', 'tbg', None, **{'deactivateIntrasublattice': '.true.'})
case('tbg+deactivateIntersublattice', 'tbg', None, **{'deactivateIntersublattice': '.true.'})
case('xyz3+deactivateIntraSublatticeForC', 'xyz3', 'x3', **{'deactivateIntraSublatticeForC': '.true.'})
case('xyz3+distanceDependentEffectiveModel', 'xyz3', 'x3', **{'MoireTwistAngle': '0.5', 'distanceDependentEffectiveModel': '.true.'})
case('xyz3+MoireBLDeactivateUpperLayer', 'xyz3', 'x3', **{'MoireTwistAngle': '0.5', 'MoireBLDeactivateUpperLayer': '.true.'})
case('xyz3+MoiretDBLDeactivateUpperLayers', 'xyz3', 'x3', **{'MoireTwistAngle': '0.5', 'MoiretDBLDeactivateUpperLayers': '.true.'})
case('xyz3+useLayerSpecificOnsiteEnergyTerms', 'xyz3', 'x3', **{'MoireTwistAngle': '0.5', 'useLayerSpecificOnsiteEnergyTerms': '.true.'})
case('eff+tBGDiag', 'eff', None, **{'tBGDiag': '.true.'})
case('eff+TypeOfBL_BLKaxiras', 'eff', None, **{'TypeOfBL': 'BLKaxiras'})

# ------------------------------------------------------------------ defaults
# middleTwist is the default for the multilayer stacks (value None removes the key from the input)
case("xyz3+middleTwist_default", "xyz3", "x3", middleTwist=None)
case("xyz4_sandwiched+middleTwist_default", "xyz4_sandwiched", "x4", middleTwist=None)
case("xyz4_sandwiched+bilayerF2G2_default", "xyz4_sandwiched", "x4", middleTwist=None, forceBilayerF2G2Intralayer=T)

# ------------------------------------------------------------------ results that depend on the compiler
# With identical sources and inputs these cases give a different Hamiltonian (or end differently) with the
# checked GNU build, the optimised GNU build and the Intel build: they use variables that are never set
# in this configuration, or random numbers. Their fingerprint is not compared; they are listed so that
# the defect stays visible until each one is repaired or refused. Found 2026-10-09; 51 before the
# "not set" markers of ham.F90, 33 with them. With the Intel compiler a marked value does not always
# propagate (its default floating-point model may drop a NaN), which is why some cases are refused by the
# GNU builds and run with Intel.
COMPILER_DEPENDENT = {
    'eff+MoireAddSecondMoire',
    'eff+MoireAddSecondMoire+Midpoint',
    'eff+MoireAddSecondMoire+Twisted2',
    'eff+MoireSecondMoireRotateFirst_off',
    'eff+MoireTrilayer',
    'graphene+Anderson',
    'graphene+Bubbles',
    'graphene+GaussDisorder',
    'graphene+MoireStrain',
    'graphene+PNP',
    'graphene+PNPKink',
    'graphene+SquareChecker2219',
    'graphene+SquareFunction',
    'graphene+SquareFunction2',
    'graphene+SublatticeDisorder',
    'graphene+Zterm1D',
    'graphene+Zterm1DKink',
    'graphene+deltaDisorder',
    'graphene+realisticBubbles',
    'graphene+sinusModulation',
    'sys_MoireEncapsulatedBilayer',
    'sys_MoireEncapsulatedBilayer+removeF2G2Flag',
    'sys_Trilayer+TrilayerAddShift',
    'tbg+BLKaxiras+deactivateV6',
    'tbg+BLKaxiras+findThetasGeometrically',
    'tbg+BLKaxiras+newFittingFunctions',
    'tbg+BLKaxiras+oldParameterSet',
    'tbg+BLKaxiras+oppositedxdy',
    'tbg+BLKaxiras+sublatticeDependent',
    'tbg+BLKaxiras+sublatticeIndependent',
    'tbg+BLKaxiras+switchV3Sign',
    'tbg+TypeOfBL_BLKaxiras+NeighLevels1',
    'tbg+twistedBLAddShift',
}
