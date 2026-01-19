from dataclasses import dataclass

@dataclass
class Isotope:
    """
    Holds the string of an isotope with its relative atom abundance.
    """
    identifier: str
    abundance: float

@dataclass
class TSLMaterial:
    """
    Holds a list of Isotopes that comprise the material and the identifier to add Sab treatment.
    """
    name: str
    tsl_identifier: str
    otf_identifier: str
    composition: list[Isotope]
    density: float = 1.0

LWTR = TSLMaterial(
    name = 'LWTR',
    tsl_identifier = 'c_H_in_H2O',
    otf_identifier = 'tsl-HinH2O_OTF',
    composition = [
        Isotope(
            identifier = "H1",
            abundance = 2
        ),
        Isotope(
            identifier = "O16",
            abundance = 1
        ),
    ],
    density = 1
)
    
HWTR = TSLMaterial(
    name = 'HWTR',
    tsl_identifier = 'c_D_in_D2O',
    otf_identifier = 'tsl-DinD2O_OTF',
    composition = [
        Isotope(
            identifier = "H2",
            abundance = 2
        ),
        Isotope(
            identifier = "O16",
            abundance = 1
        ),
    ],
    density = 1.107
)
    
GRPH = TSLMaterial(
    name = 'GRPH',
    tsl_identifier = 'c_Graphite',
    otf_identifier = 'tsl-crystalline-graphite_OTF',
    composition = [
        Isotope(
            identifier = "C12",
            abundance = 1
        ),
    ],
    density = 1.7
)
    
GRPH10P = TSLMaterial(
    name = 'GRPH10P',
    tsl_identifier = 'c_Graphite_10p',
    otf_identifier = 'tsl-reactor-graphite-10P_OTF',
    composition = [
        Isotope(
            identifier = "C12",
            abundance = 1
        ),
    ],
    density = 1.7
)
    
GRPH30P = TSLMaterial(
    name = 'GRPH30P',
    tsl_identifier = 'c_Graphite_30p',
    otf_identifier = 'tsl-reactor-graphite-30P_OTF',
    composition = [
        Isotope(
            identifier = "C12",
            abundance = 1
        ),
    ],
    density = 1.7
)
    
HinYH = TSLMaterial(
    name = 'HinYH',
    tsl_identifier = 'c_H_in_YH2',
    otf_identifier = 'tsl-HinYH2_OTF',
    composition = [
        Isotope(
            identifier = "H1",
            abundance = 2
        ),
        Isotope(
            identifier = "Y89",
            abundance = 1
        ),
    ],
    density = 4.2
)
    
HinZrH = TSLMaterial(
    name = 'HinZrH',
    tsl_identifier = 'c_H_in_ZrH',
    otf_identifier = 'tsl-HinZrH_OTF',
    composition = [
        Isotope(
            identifier = "H1",
            abundance = 2
        ),
        Isotope(
            identifier = "Zr91",
            abundance = 1
        ),
    ],
    density = 5.7
)
    
CinSiC = TSLMaterial(
    name = 'CinSiC',
    tsl_identifier = 'c_C_in_SiC',
    otf_identifier = 'tsl-CinSiC_OTF',
    composition = [
        Isotope(
            identifier = "C12",
            abundance = 1
        ),
        Isotope(
            identifier = "Si28",
            abundance = 1
        ),
    ],
    density = 3.16
)
    
SiinSiC = TSLMaterial(
    name = 'SiinSiC',
    tsl_identifier = 'c_Si_in_SiC',
    otf_identifier = 'tsl-SiinSiC_OTF',
    composition = [
        Isotope(
            identifier = "C12",
            abundance = 1
        ),
        Isotope(
            identifier = "Si28",
            abundance = 1
        ),
    ],
    density = 3.16
)
    
BeMetal = TSLMaterial(
    name = 'BeMetal',
    tsl_identifier = 'c_Be',
    otf_identifier = 'tsl-Be-metal_OTF',
    composition = [
        Isotope(
            identifier = "Be9",
            abundance = 1
        ),
    ],
    density = 1.85
)
    
BeinBeO = TSLMaterial(
    name = 'BeinBeO',
    tsl_identifier = 'c_Be_in_BeO',
    otf_identifier = 'tsl-BeinBeO_OTF',
    composition = [
        Isotope(
            identifier = "Be9",
            abundance = 1
        ),
        Isotope(
            identifier = "O16",
            abundance = 1
        ),
    ],
    density = 3
)
    
OinBeO = TSLMaterial(
    name = 'OinBeO',
    tsl_identifier = 'c_O_in_BeO',
    otf_identifier = 'tsl-OinBeO_OTF',
    composition = [
        Isotope(
            identifier = "Be9",
            abundance = 1
        ),
        Isotope(
            identifier = "O16",
            abundance = 1
        ),
    ],
    density = 3
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