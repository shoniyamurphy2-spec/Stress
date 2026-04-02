"""
Stress Concentration & Safety Factor Estimator
for Load-Bearing Mechanical Components
"""

import math
from dataclasses import dataclass, field
from typing import Optional
from enum import Enum


# ─────────────────────────────────────────────
# Enumerations
# ─────────────────────────────────────────────

class LoadType(Enum):
    AXIAL       = "Axial (Tension/Compression)"
    BENDING     = "Bending"
    TORSION     = "Torsion"
    COMBINED    = "Combined (Axial + Bending)"


class FeatureType(Enum):
    CIRCULAR_HOLE      = "Circular Hole in Plate"
    SHOULDER_FILLET    = "Shoulder Fillet"
    GROOVE             = "U-Groove / Circumferential Groove"
    NOTCH              = "V-Notch"
    KEYWAY             = "Keyway"
    PRESS_FIT          = "Press Fit"


class FailureCriterion(Enum):
    VON_MISES    = "Von Mises (Distortion Energy)"
    TRESCA       = "Tresca (Max Shear Stress)"
    MAX_NORMAL   = "Maximum Normal Stress"


# ─────────────────────────────────────────────
# Data Classes
# ─────────────────────────────────────────────

@dataclass
class Material:
    name: str
    yield_strength: float        # MPa
    ultimate_strength: float     # MPa
    endurance_limit: float       # MPa  (approx 0.5 * Su for steels)
    elastic_modulus: float       # GPa
    poisson_ratio: float = 0.3

    @classmethod
    def from_library(cls, key: str) -> "Material":
        return MATERIAL_LIBRARY[key]


@dataclass
class GeometryFeature:
    feature_type: FeatureType
    # Shared dimensions (mm)
    width:  Optional[float] = None   # plate width or shaft diameter D
    height: Optional[float] = None   # plate height or reduced diameter d
    hole_diameter: Optional[float] = None
    fillet_radius: Optional[float] = None
    groove_depth:  Optional[float] = None
    groove_radius: Optional[float] = None


@dataclass
class LoadCase:
    load_type: LoadType
    axial_force:    float = 0.0   # N
    bending_moment: float = 0.0   # N·mm
    torque:         float = 0.0   # N·mm
    cross_section_area: float = 0.0   # mm²
    section_modulus:    float = 0.0   # mm³
    polar_section_mod:  float = 0.0   # mm³


@dataclass
class AnalysisResult:
    nominal_stress:   float   # MPa
    Kt:               float   # theoretical stress concentration factor
    peak_stress:      float   # MPa  (Kt × σ_nom)
    safety_factor_yield:   float
    safety_factor_fatigue: float
    safety_factor_ultimate: float
    governing_criterion: str
    warnings: list[str] = field(default_factory=list)
    details: dict = field(default_factory=dict)


# ─────────────────────────────────────────────
# Material Library
# ─────────────────────────────────────────────

MATERIAL_LIBRARY: dict[str, Material] = {
    "1020_steel": Material(
        name="AISI 1020 Steel (HR)",
        yield_strength=210, ultimate_strength=380,
        endurance_limit=190, elastic_modulus=200),
    "4140_steel": Material(
        name="AISI 4140 Steel (Q&T 600°F)",
        yield_strength=1570, ultimate_strength=1720,
        endurance_limit=860, elastic_modulus=207),
    "304_ss": Material(
        name="304 Stainless Steel (Annealed)",
        yield_strength=215, ultimate_strength=505,
        endurance_limit=240, elastic_modulus=193),
    "al6061_t6": Material(
        name="Aluminum 6061-T6",
        yield_strength=276, ultimate_strength=310,
        endurance_limit=97, elastic_modulus=69,
        poisson_ratio=0.33),
    "al7075_t6": Material(
        name="Aluminum 7075-T6",
        yield_strength=503, ultimate_strength=572,
        endurance_limit=160, elastic_modulus=72,
        poisson_ratio=0.33),
    "ti6al4v": Material(
        name="Titanium Ti-6Al-4V",
        yield_strength=880, ultimate_strength=950,
        endurance_limit=515, elastic_modulus=114,
        poisson_ratio=0.34),
    "cast_iron_g3000": Material(
        name="Gray Cast Iron G3000",
        yield_strength=207, ultimate_strength=214,
        endurance_limit=100, elastic_modulus=100),
}


# ─────────────────────────────────────────────
# Stress Concentration Factor (Kt) Calculators
# Using Pilkey & Pilkey curve-fit formulas
# ─────────────────────────────────────────────

def kt_circular_hole_axial(w: float, d: float) -> float:
    """
    Infinite plate with circular hole under axial tension.
    w = plate width, d = hole diameter.
    Uses Pilkey (2008) finite-width correction.
    """
    ratio = d / w
    if ratio >= 1.0:
        raise ValueError("Hole diameter must be less than plate width.")
    # Heywood's formula (finite width)
    Kt = (3 - 3.13 * ratio + 3.66 * ratio**2 - 1.53 * ratio**3)
    return max(Kt, 1.0)


def kt_shoulder_fillet_axial(D: float, d: float, r: float) -> float:
    """
    Stepped shaft / plate under axial load.
    D = large diameter, d = small diameter, r = fillet radius.
    Neuber / Pilkey curve fit.
    """
    if r <= 0:
        raise ValueError("Fillet radius must be > 0.")
    t = (D - d) / 2          # step height
    ratio_r_d = r / d
    ratio_t_r = t / r
    # Simplified Pilkey axial shoulder fillet
    C1 =  0.926 + 1.157 * math.sqrt(t / d) - 0.099 * (t / d)
    C2 = -0.012 - 1.139 * math.sqrt(t / d) + 0.217 * (t / d)
    C3 = -0.302 + 3.061 * math.sqrt(t / d) - 0.739 * (t / d)
    C4 =  0.387 - 3.070 * math.sqrt(t / d) + 0.610 * (t / d)
    Kt = C1 + C2 * (2*t/d) + C3 * (2*t/d)**2 + C4 * (2*t/d)**3
    # Adjust with r/d
    correction = 1.0 / (ratio_r_d ** 0.3)
    return max(Kt * correction / 3.5, 1.0)   # empirical normalisation


def kt_shoulder_fillet_bending(D: float, d: float, r: float) -> float:
    """Shoulder fillet under bending (Pilkey 2008 Table 2-8 curve fit)."""
    if r <= 0:
        raise ValueError("Fillet radius must be > 0.")
    t = (D - d) / 2
    ratio_r_d = r / d
    A = -0.291 + 1.025 * (t / r)**0.5 - 0.025 * (t / r)
    Kt = 1 + A * (2 * t / r) ** 0.5
    return max(Kt, 1.0)


def kt_shoulder_fillet_torsion(D: float, d: float, r: float) -> float:
    """Shoulder fillet under torsion."""
    if r <= 0:
        raise ValueError("Fillet radius must be > 0.")
    t = (D - d) / 2
    ratio_r_d = r / d
    Kt = 1 + 0.6 * (t / r) ** 0.5
    return max(Kt, 1.0)


def kt_u_groove(D: float, d: float, r: float) -> float:
    """Circumferential U-groove in a round bar (axial/bending)."""
    if r <= 0:
        raise ValueError("Groove root radius must be > 0.")
    t = (D - d) / 2
    a = t / r
    Kt = 1 + 0.78 * a ** 0.5 + 0.20 * a
    return max(Kt, 1.0)


def kt_keyway(d: float) -> tuple[float, float]:
    """
    Keyway in a shaft.
    Returns (Kt_bending, Kt_torsion) using Pilkey / Peterson values.
    """
    # Standard ANSI keyway proportions: width ≈ d/4, depth ≈ d/8
    Kt_bending = 2.14
    Kt_torsion = 3.0
    return Kt_bending, Kt_torsion


def kt_press_fit(d: float) -> float:
    """Press fit hub on shaft – stress concentration at fit edge."""
    # Peterson recommends Kt ~ 1.9–3.0; use 2.0 as conservative baseline
    return 2.0


# ─────────────────────────────────────────────
# Nominal Stress Calculators
# ─────────────────────────────────────────────

def nominal_stress_axial(force: float, area: float) -> float:
    """σ_nom = F / A  (MPa)"""
    if area <= 0:
        raise ValueError("Cross-sectional area must be > 0.")
    return force / area


def nominal_stress_bending(moment: float, section_mod: float) -> float:
    """σ_nom = M / Z  (MPa)"""
    if section_mod <= 0:
        raise ValueError("Section modulus must be > 0.")
    return moment / section_mod


def nominal_shear_torsion(torque: float, polar_mod: float) -> float:
    """τ_nom = T / Zp  (MPa)"""
    if polar_mod <= 0:
        raise ValueError("Polar section modulus must be > 0.")
    return torque / polar_mod


# ─────────────────────────────────────────────
# Section Property Helpers
# ─────────────────────────────────────────────

def circular_section(d: float) -> tuple[float, float, float]:
    """Returns (area mm², section_modulus mm³, polar_section_mod mm³)."""
    A  = math.pi * d**2 / 4
    Z  = math.pi * d**3 / 32
    Zp = math.pi * d**3 / 16
    return A, Z, Zp


def rectangular_section(b: float, h: float) -> tuple[float, float]:
    """Returns (area mm², section_modulus mm³) for rectangular plate."""
    A = b * h
    Z = b * h**2 / 6
    return A, Z


# ─────────────────────────────────────────────
# Main Analysis Engine
# ─────────────────────────────────────────────

def analyse(
    feature: GeometryFeature,
    load: LoadCase,
    material: Material,
    criterion: FailureCriterion = FailureCriterion.VON_MISES,
) -> AnalysisResult:
    """
    Full stress concentration + safety factor analysis.
    """
    warnings: list[str] = []
    details: dict = {}

    # ── 1. Compute Kt ──────────────────────────────
    ft = feature.feature_type

    if ft == FeatureType.CIRCULAR_HOLE:
        if not (feature.width and feature.hole_diameter):
            raise ValueError("Circular hole requires width and hole_diameter.")
        Kt = kt_circular_hole_axial(feature.width, feature.hole_diameter)
        net_area = (feature.width - feature.hole_diameter) * (feature.height or feature.width)
        load.cross_section_area = net_area

    elif ft == FeatureType.SHOULDER_FILLET:
        D = feature.width
        d = feature.height
        r = feature.fillet_radius
        if not (D and d and r):
            raise ValueError("Shoulder fillet requires width(D), height(d), fillet_radius.")
        if load.load_type == LoadType.TORSION:
            Kt = kt_shoulder_fillet_torsion(D, d, r)
        elif load.load_type == LoadType.BENDING:
            Kt = kt_shoulder_fillet_bending(D, d, r)
        else:
            Kt = kt_shoulder_fillet_axial(D, d, r)
        # Update section props to reduced section
        A, Z, Zp = circular_section(d)
        if load.cross_section_area == 0:  load.cross_section_area = A
        if load.section_modulus    == 0:  load.section_modulus    = Z
        if load.polar_section_mod  == 0:  load.polar_section_mod  = Zp

    elif ft == FeatureType.GROOVE:
        D = feature.width
        d = feature.height
        r = feature.groove_radius
        if not (D and d and r):
            raise ValueError("Groove requires width(D), height(d), groove_radius.")
        Kt = kt_u_groove(D, d, r)
        A, Z, Zp = circular_section(d)
        if load.cross_section_area == 0:  load.cross_section_area = A
        if load.section_modulus    == 0:  load.section_modulus    = Z
        if load.polar_section_mod  == 0:  load.polar_section_mod  = Zp

    elif ft == FeatureType.KEYWAY:
        d = feature.width
        if not d:
            raise ValueError("Keyway requires shaft diameter (width).")
        Kt_b, Kt_t = kt_keyway(d)
        Kt = Kt_t if load.load_type == LoadType.TORSION else Kt_b
        A, Z, Zp = circular_section(d)
        if load.cross_section_area == 0:  load.cross_section_area = A
        if load.section_modulus    == 0:  load.section_modulus    = Z
        if load.polar_section_mod  == 0:  load.polar_section_mod  = Zp

    elif ft == FeatureType.PRESS_FIT:
        d = feature.width
        if not d:
            raise ValueError("Press fit requires shaft diameter (width).")
        Kt = kt_press_fit(d)
        A, Z, Zp = circular_section(d)
        if load.cross_section_area == 0:  load.cross_section_area = A
        if load.section_modulus    == 0:  load.section_modulus    = Z
        if load.polar_section_mod  == 0:  load.polar_section_mod  = Zp

    else:
        raise ValueError(f"Unsupported feature type: {ft}")

    # ── 2. Nominal stress ──────────────────────────
    sigma_nom = 0.0
    tau_nom   = 0.0

    lt = load.load_type
    if lt in (LoadType.AXIAL, LoadType.COMBINED):
        sigma_nom += nominal_stress_axial(load.axial_force, load.cross_section_area)
    if lt in (LoadType.BENDING, LoadType.COMBINED):
        sigma_nom += nominal_stress_bending(load.bending_moment, load.section_modulus)
    if lt == LoadType.TORSION:
        tau_nom = nominal_shear_torsion(load.torque, load.polar_section_mod)

    # ── 3. Peak (local) stress ─────────────────────
    sigma_peak = Kt * sigma_nom
    tau_peak   = Kt * tau_nom

    # ── 4. Equivalent stress (failure criterion) ───
    if criterion == FailureCriterion.VON_MISES:
        sigma_eq = math.sqrt(sigma_peak**2 + 3 * tau_peak**2)
        crit_name = "Von Mises"
    elif criterion == FailureCriterion.TRESCA:
        sigma_eq = math.sqrt(sigma_peak**2 + 4 * tau_peak**2)
        crit_name = "Tresca"
    else:
        sigma_eq = sigma_peak + tau_peak   # conservative bound
        crit_name = "Max Normal"

    if sigma_eq == 0:
        raise ValueError("Equivalent stress is zero – check inputs.")

    # ── 5. Safety factors ─────────────────────────
    Sy  = material.yield_strength
    Su  = material.ultimate_strength
    Se  = material.endurance_limit

    SF_yield   = Sy  / sigma_eq
    SF_fatigue = Se  / sigma_eq
    SF_ult     = Su  / sigma_eq

    # ── 6. Warnings ───────────────────────────────
    if Kt > 4.0:
        warnings.append(f"Very high Kt = {Kt:.2f} – consider redesign (larger fillet, groove radius, etc.)")
    if SF_yield < 1.0:
        warnings.append("⚠  YIELDING EXPECTED under current loading!")
    elif SF_yield < 1.5:
        warnings.append("Low static safety factor – consider increasing to ≥ 1.5 for ductile materials.")
    if SF_fatigue < 1.0:
        warnings.append("⚠  FATIGUE FAILURE LIKELY under cyclic loading!")
    elif SF_fatigue < 2.0:
        warnings.append("Low fatigue safety factor – standard design targets ≥ 2.0.")

    details = {
        "nominal_normal_stress_MPa": round(sigma_nom, 3),
        "nominal_shear_stress_MPa":  round(tau_nom, 3),
        "peak_normal_stress_MPa":    round(sigma_peak, 3),
        "peak_shear_stress_MPa":     round(tau_peak, 3),
        "equivalent_stress_MPa":     round(sigma_eq, 3),
        "failure_criterion":         crit_name,
        "Kt":                        round(Kt, 4),
        "material":                  material.name,
        "Sy_MPa":                    Sy,
        "Su_MPa":                    Su,
        "Se_MPa":                    Se,
    }

    return AnalysisResult(
        nominal_stress=round(sigma_nom if sigma_nom else tau_nom, 3),
        Kt=round(Kt, 4),
        peak_stress=round(sigma_eq, 3),
        safety_factor_yield=round(SF_yield, 3),
        safety_factor_fatigue=round(SF_fatigue, 3),
        safety_factor_ultimate=round(SF_ult, 3),
        governing_criterion=crit_name,
        warnings=warnings,
        details=details,
    )


# ─────────────────────────────────────────────
# Pretty-print report
# ─────────────────────────────────────────────

def print_report(result: AnalysisResult) -> None:
    line = "─" * 58
    print(f"\n{'═'*58}")
    print(f"  STRESS CONCENTRATION & SAFETY FACTOR REPORT")
    print(f"{'═'*58}")
    d = result.details
    print(f"  Material             : {d['material']}")
    print(f"  Failure Criterion    : {d['failure_criterion']}")
    print(line)
    print(f"  Nominal Stress       : {d['nominal_normal_stress_MPa']:>10.3f}  MPa (normal)")
    if d['nominal_shear_stress_MPa']:
        print(f"  Nominal Shear Stress : {d['nominal_shear_stress_MPa']:>10.3f}  MPa")
    print(f"  Stress Conc. Kt      : {result.Kt:>10.4f}")
    print(f"  Peak Normal Stress   : {d['peak_normal_stress_MPa']:>10.3f}  MPa")
    if d['peak_shear_stress_MPa']:
        print(f"  Peak Shear Stress    : {d['peak_shear_stress_MPa']:>10.3f}  MPa")
    print(f"  Equivalent Stress    : {d['equivalent_stress_MPa']:>10.3f}  MPa")
    print(line)
    print(f"  Sy  (Yield Strength) : {d['Sy_MPa']:>10.1f}  MPa")
    print(f"  Su  (UTS)            : {d['Su_MPa']:>10.1f}  MPa")
    print(f"  Se  (Endurance Lim.) : {d['Se_MPa']:>10.1f}  MPa")
    print(line)
    print(f"  Safety Factor (Yield)   : {result.safety_factor_yield:>8.3f}")
    print(f"  Safety Factor (Fatigue) : {result.safety_factor_fatigue:>8.3f}")
    print(f"  Safety Factor (UTS)     : {result.safety_factor_ultimate:>8.3f}")
    print(line)
    if result.warnings:
        print("  WARNINGS:")
        for w in result.warnings:
            print(f"    • {w}")
    else:
        print("  ✓  No critical warnings.")
    print(f"{'═'*58}\n")


# ─────────────────────────────────────────────
# CLI Demo – three worked examples
# ─────────────────────────────────────────────

def demo():
    print("\n" + "★"*58)
    print("  Stress Concentration Estimator – Worked Examples")
    print("★"*58)

    # ── Example 1: Shoulder fillet on 4140 shaft under bending ──
    print("\n[Example 1] AISI 4140 Shaft – Shoulder Fillet – Bending")
    mat1 = Material.from_library("4140_steel")
    geo1 = GeometryFeature(
        feature_type=FeatureType.SHOULDER_FILLET,
        width=50.0,           # D = 50 mm (large diameter)
        height=38.0,          # d = 38 mm (reduced diameter)
        fillet_radius=3.0,    # r = 3 mm
    )
    _, Z1, _ = circular_section(38.0)
    lc1 = LoadCase(
        load_type=LoadType.BENDING,
        bending_moment=1_200_000,   # 1200 N·m = 1,200,000 N·mm
        section_modulus=Z1,
    )
    r1 = analyse(geo1, lc1, mat1, FailureCriterion.VON_MISES)
    print_report(r1)

    # ── Example 2: Al 6061-T6 plate with circular hole under axial ──
    print("[Example 2] Al 6061-T6 – Circular Hole – Axial Tension")
    mat2 = Material.from_library("al6061_t6")
    geo2 = GeometryFeature(
        feature_type=FeatureType.CIRCULAR_HOLE,
        width=80.0,           # plate width 80 mm
        height=10.0,          # plate thickness 10 mm
        hole_diameter=20.0,   # hole ø 20 mm
    )
    net_area = (80 - 20) * 10   # net cross-section
    lc2 = LoadCase(
        load_type=LoadType.AXIAL,
        axial_force=45_000,         # 45 kN
        cross_section_area=net_area,
    )
    r2 = analyse(geo2, lc2, mat2, FailureCriterion.VON_MISES)
    print_report(r2)

    # ── Example 3: Steel shaft with keyway under torsion ──
    print("[Example 3] AISI 1020 Shaft – Keyway – Torsion")
    mat3 = Material.from_library("1020_steel")
    d3 = 40.0
    _, _, Zp3 = circular_section(d3)
    geo3 = GeometryFeature(
        feature_type=FeatureType.KEYWAY,
        width=d3,
    )
    lc3 = LoadCase(
        load_type=LoadType.TORSION,
        torque=850_000,           # 850 N·m
        polar_section_mod=Zp3,
    )
    r3 = analyse(geo3, lc3, mat3, FailureCriterion.TRESCA)
    print_report(r3)


if __name__ == "__main__":
    demo()
