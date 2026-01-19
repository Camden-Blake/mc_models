from dataclasses import dataclass

import mcnpy as mp
from mcnpy.elements import *

UNIT = 'ATOM'
LIBRARY = '01C'

H1_a1   = mp.Nuclide(name = 'H1',   fraction = 1, unit = UNIT, library = LIBRARY)
H1_a2   = mp.Nuclide(name = 'H1',   fraction = 2, unit = UNIT, library = LIBRARY)
H2_a2   = mp.Nuclide(name = 'H2',   fraction = 2, unit = UNIT, library = LIBRARY)
BE9_a1  = mp.Nuclide(name = 'BE9',  fraction = 1, unit = UNIT, library = LIBRARY)
C12_a1  = mp.Nuclide(name = 'C12',  fraction = 1, unit = UNIT, library = LIBRARY)
O16_a1  = mp.Nuclide(name = 'O16',  fraction = 1, unit = UNIT, library = LIBRARY)
SI28_a1 = mp.Nuclide(name = 'SI28', fraction = 1, unit = UNIT, library = LIBRARY)
Y89_a1  = mp.Nuclide(name = 'Y89',  fraction = 1, unit = UNIT, library = LIBRARY)
ZR91_a1 = mp.Nuclide(name = 'ZR91', fraction = 1, unit = UNIT, library = LIBRARY)

@dataclass
class TSLMaterial:
    """
    Holds a list of Isotopes that comprise the material and the identifier to add Sab treatment.
    """
    name: str
    mcnp_placeholder_tsl: str
    tsl_identifier: str
    otf_identifier: str
    nuclides: list[mp.Nuclide]
    density: float = 1.0

LWTR = TSLMaterial(
    name = 'LWTR',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'h-h2o.40t',
    otf_identifier = 'tsl-HinH2O_OTF',
    nuclides = [H1_a2, O16_a1],
    density = 1,
)
    
HWTR = TSLMaterial(
    name = 'HWTR',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'd-d2o.40t',
    otf_identifier = 'tsl-DinD2O_OTF',
    nuclides = [H2_a2, O16_a1],
    density = 1.107,
)
    
GRPH = TSLMaterial(
    name = 'GRPH',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'grph.40t',
    otf_identifier = 'tsl-crystalline-graphite_OTF',
    nuclides = [C12_a1],
    density = 1.7,
)
    
GRPH10P = TSLMaterial(
    name = 'GRPH10P',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'grph10.40t',
    otf_identifier = 'tsl-reactor-graphite-10P_OTF',
    nuclides = [C12_a1],
    density = 1.7,
)
    
GRPH30P = TSLMaterial(
    name = 'GRPH30P',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'grph30.40t',
    otf_identifier = 'tsl-reactor-graphite-30P_OTF',
    nuclides = [C12_a1],
    density = 1.7,
)
    
HinYH = TSLMaterial(
    name = 'HinYH',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'h-yh2.40t',
    otf_identifier = 'tsl-HinYH2_OTF',
    nuclides = [H1_a2, Y89_a1],
    density = 4.2,
)
    
HinZrH = TSLMaterial(
    name = 'HinZrH',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'h-zrh.40t',
    otf_identifier = 'tsl-HinZrH_OTF',
    nuclides = [H1_a2, ZR91_a1],
    density = 5.7,
)
    
CinSiC = TSLMaterial(
    name = 'CinSiC',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'c-sic.40t',
    otf_identifier = 'tsl-CinSiC_OTF',
    nuclides = [C12_a1, SI28_a1],
    density = 3.16,
)
    
SiinSiC = TSLMaterial(
    name = 'SiinSiC',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'si-sic.40t',
    otf_identifier = 'tsl-SiinSiC_OTF',
    nuclides = [C12_a1, SI28_a1],
    density = 3.16,
)
    
BeMetal = TSLMaterial(
    name = 'BeMetal',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'be-met.40t',
    otf_identifier = 'tsl-Be-metal_OTF',
    nuclides = [BE9_a1],
    density = 1.85,
)
    
BeinBeO = TSLMaterial(
    name = 'BeinBeO',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'be-beo.40t',
    otf_identifier = 'tsl-BeinBeO_OTF',
    nuclides = [BE9_a1, O16_a1],
    density = 3,
)
    
OinBeO = TSLMaterial(
    name = 'OinBeO',
    mcnp_placeholder_tsl = 'LWTR',
    tsl_identifier = 'o-beo.40t',
    otf_identifier = 'tsl-OinBeO_OTF',
    nuclides = [BE9_a1, O16_a1],
    density = 3,
)

MATERIALS = [
    LWTR,
    HWTR,
    GRPH,
    GRPH10P,
    GRPH30P,
    HinYH,
    HinZrH,
    SiinSiC,
    CinSiC,
    BeMetal,
    BeinBeO,
    OinBeO,
]

MATERIAL_NAMES = [material.name for material in MATERIALS]