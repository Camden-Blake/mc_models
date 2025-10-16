# Description of the model : https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nea4440-crocus.pdf

import numpy as np
import operator
from dataclasses import dataclass, field
from itertools import count

import mcnpy as mp
from mcnpy.elements import *

__id_counter = count(1)
def next_id():
    return next(__id_counter)

def sum_material_fraction(mat: mp.Material) -> float:
    return sum([nuclide.fraction for nuclide in mat.nuclides])

def combine_region(base: mp.Region | None, 
                   addition: mp.Region, 
                   op: operator = operator.and_) -> mp.Region:
    """Safely combine two regions using a given operator (default is intersection)."""
    if base == None:
        return addition
    return op(base, addition)

@dataclass
class CruciformPrism:
    """Class for holding the surfaces and regions of a cruciform prism."""
    surfaces: list[mp.Surface] = field(default_factory=list)
    internal_region: mp.Region = None


def cruciform_prism(distances:list[float], 
                    center:tuple[float, float]=(0,0), 
                    axis:str='z',
                    name: str | None = None) -> CruciformPrism:
    """
    Constructs the surfaces and regions of a generalized cruciform prism (cross-shaped).

    Parameters
    ----------
    distances : list of float
        A monotonically increasing iterable of distances in [cm] that form the planes of the generalized cruciform.
    center : (float, float)
        Center of the cruciform (in the plane perpendicular to the prism axis).
    axis : str
        Axis along which the prism extends ('x', 'y', or 'z').
    name : str | None
        Name for the cruciform prism to add to the surface comments.

    Returns
    -------
    CruciformPrism
        Object holding all surfaces and the combined regions.
    """

    base_comment = f"{name} - Cruciform Prism" if name else "Cruciform Prism"
    prism = CruciformPrism()
    c1, c2 = center

    axis_map = {
        "z": (mp.XPlane, "x0", mp.YPlane, "y0"),
        "y": (mp.XPlane, "x0", mp.ZPlane, "z0"),
        "x": (mp.YPlane, "y0", mp.ZPlane, "z0")
    }
    if axis not in axis_map:
        raise ValueError(f"Axis must be in {axis_map.keys()}.")
    surf_type_1, key_1, surf_type_2, key_2 = axis_map[axis]

    def make_surface(surf_type: mp.Surface, key: str, val, comment: str) -> mp.Surface:
        return surf_type(**{key: val}, comment=comment)

    def cruciform_arm(d1: float, d2: float, tag: str) -> tuple[list[mp.Surface],mp.Region]:
        s1 = make_surface(surf_type_1, key_1, c1 - d1, tag + " - axis 1 lower plane")
        s2 = make_surface(surf_type_1, key_1, c1 + d1, tag + " - axis 1 upper plane")
        s3 = make_surface(surf_type_2, key_2, c2 - d2, tag + " - axis 2 lower plane")
        s4 = make_surface(surf_type_2, key_2, c2 + d2, tag + " - axis 2 upper plane")
        return [s1, s2, s3, s4], (+s1 & -s2 & +s3 & -s4)

    iters = len(distances) // 2
    for i in range(iters):
        dl = distances[i]
        du = distances[-i - 1]

        arm1_surfs, arm1_region = cruciform_arm(dl, du, base_comment + f" - layer {i} - arm 1")
        arm2_surfs, arm2_region = cruciform_arm(du, dl, base_comment + f" - layer {i} - arm 2")

        prism.surfaces += arm1_surfs + arm2_surfs
        prism.internal_region = combine_region(prism.internal_region, arm1_region, operator.or_)
        prism.internal_region = combine_region(prism.internal_region, arm2_region, operator.or_)

    if len(distances) % 2:  # odd number of distances → center square
        d = distances[iters]
        square_surfs, square_region = cruciform_arm(d, d, base_comment + f" - layer {i} square")
        prism.surfaces += square_surfs
        prism.internal_region = combine_region(prism.internal_region, square_region, operator.or_)

    return prism

deck = mp.Deck()

# Materials
#===============================================================================

UO2 = mp.Material(
    U[235]@4.30565e-04 +
     O[16]@4.70902e-02 +
    U[238]@2.31145e-02 
)
UO2_density = sum_material_fraction(UO2)
deck += UO2

Umetal = mp.Material(
    U[235]@4.53160e-04 +
    U[238]@4.68003e-02
)
Umetal_density = sum_material_fraction(Umetal)
deck += Umetal

Clad_non_fuel = mp.Material(
    AL[27]@6.02611e-02
)
Clad_non_fuel_density = sum_material_fraction(Clad_non_fuel)
deck += Clad_non_fuel

Clad_uo2 = mp.Material(
    AL[27]@5.00614e-02
)
Clad_uo2_density = sum_material_fraction(Clad_uo2)
deck += Clad_uo2

Clad_umetal = mp.Material(
    AL[27]@5.17799e-2
)
Clad_umetal_density = sum_material_fraction(Clad_umetal)
deck += Clad_umetal

Moderator = mp.Material(
     H[1]@6.67578e-02 +
    O[16]@3.33789e-02
)
Moderator_density = sum_material_fraction(Moderator)
Moderator.s_alpha_beta = 'lwtr'
deck += Moderator

Support_Material = mp.Material(
    AL[27]@6.02611e-02
)
Support_Material_density = sum_material_fraction(Support_Material)
deck += Support_Material

Cadmium = mp.Material(
    CD[106]@0.0005768508299999999 +
    CD[108]@0.000411440592 +
    CD[110]@0.00577777498 +
    CD[111]@0.00592835853 +
    CD[112]@0.011170519406 +
    CD[113]@0.0056651848179999995 +
    CD[114]@0.013322705835999999 +
    CD[116]@0.003480565008
)
Cadmium_density = sum_material_fraction(Cadmium)
deck += Cadmium

Filler_Gas = mp.Material(
    HE[4]@1.64e-04
)
Filler_Gas_density = sum_material_fraction(Filler_Gas)
deck += Filler_Gas

Air = mp.Material(
    N[14]@4.14e-05 +
    O[16]@9.07e-06 
)
Air_density = sum_material_fraction(Air)
deck += Air

#===============================================================================
# Geometry
#===============================================================================

## Measured from the bottom of the tank
## ta = tank
## bp = base plate
## lg = lower grid
## ug = upper grid
## ca = cadmium
## pi = pin
## fu = fuel
## fg = fill gas
## wa = water

## top = top
## bot = bottom
## inr = inner
## out = outer

# Set the water hieght
wa_height = 96.51  # This measures from the bottom of the active fuel region with a maximum of 100

# Axial distances from tank bottom
ta_bot_dist     = 0
bp_top_dist     = ta_bot_dist + 3.0
lg_bot_bot_dist = bp_top_dist + 2.15
lg_bot_top_dist = lg_bot_bot_dist + .5
lg_top_bot_dist = lg_bot_top_dist + 0.05
lg_top_top_dist = lg_top_bot_dist + .5
ug_bot_bot_dist = lg_top_top_dist + 100
ug_bot_top_dist = ug_bot_bot_dist + .5
ug_top_bot_dist = ug_bot_top_dist + 0.05
ug_top_top_dist = ug_top_bot_dist + 1.5

wa_H0 = lg_top_bot_dist

pi_bot_dist = bp_top_dist
fu_bot_dist = pi_bot_dist + 2.7
fu_top_dist = fu_bot_dist + 100
uo2_fg_top_dist = fu_top_dist + 0.5
umet_fg_top_dist = fu_top_dist + 1.47
pi_top_dist = fu_top_dist + 17.3

ta_top_dist = pi_top_dist

# Shared axial Surfaces
ta_bot_surf      = mp.ZPlane(z0=ta_bot_dist,       comment = 'Tank Bottom')
bp_top_surf      = mp.ZPlane(z0=bp_top_dist,       comment = 'Base Plate Top')
lg_bot_bot_surf  = mp.ZPlane(z0=lg_bot_bot_dist,   comment = 'Lower Grid Bottom Bottom')
lg_bot_top_surf  = mp.ZPlane(z0=lg_bot_top_dist,   comment = 'Lower Grid Bottom Top')
lg_top_bot_surf  = mp.ZPlane(z0=lg_top_bot_dist,   comment = 'Lower Grid Top Bottom')
lg_top_top_surf  = mp.ZPlane(z0=lg_top_top_dist,   comment = 'Lower Grid Top Top')
ug_bot_bot_surf  = mp.ZPlane(z0=ug_bot_bot_dist,   comment = 'Upper Grid Bottom Bottom')
ug_bot_top_surf  = mp.ZPlane(z0=ug_bot_top_dist,   comment = 'Upper Grid Bottom Top')
ug_top_bot_surf  = mp.ZPlane(z0=ug_top_bot_dist,   comment = 'Upper Grid Top Bottom')
ug_top_top_surf  = mp.ZPlane(z0=ug_top_top_dist,   comment = 'Upper Grid Top Top')
wa_top_surf      = mp.ZPlane(z0=wa_H0 + wa_height, comment = 'Water Top')
# pi_bot_surf      = mp.ZPlane(z0=pi_bot_dist,       comment = 'Pin Bottom')
fu_bot_surf      = mp.ZPlane(z0=fu_bot_dist,       comment = 'Fuel Bottom')
fu_top_surf      = mp.ZPlane(z0=fu_top_dist,       comment = 'Fuel Top')
uo2_fg_top_surf  = mp.ZPlane(z0=uo2_fg_top_dist,   comment = 'UO2 Fill Gas Top')
umet_fg_top_surf = mp.ZPlane(z0=umet_fg_top_dist,  comment = 'UMET Fill Gas Top')
# pi_top_surf      = mp.ZPlane(z0=pi_top_dist,       comment = 'Pin Top')
ta_top_surf      = mp.ZPlane(z0=ta_top_dist,       comment = 'Tank Top')
deck += [
    ta_bot_surf,
    bp_top_surf,
    lg_bot_bot_surf,
    lg_bot_top_surf,
    lg_top_bot_surf,
    lg_top_top_surf,
    ug_bot_bot_surf,
    ug_bot_top_surf,
    ug_top_bot_surf,
    ug_top_top_surf,
    wa_top_surf,
    # pi_bot_surf,
    fu_bot_surf,
    fu_top_surf,
    uo2_fg_top_surf,
    umet_fg_top_surf,
    # pi_top_surf,
    ta_top_surf,
]

# Radial distances from tank center
uo2_pi_fu_rad = 0.526
uo2_pi_clad_rad = 0.63
uo2_pi_pitch = 1.8370
uo2_lat_ele_x = 22
uo2_lat_ele_y = 22

umet_pi_fu_rad = 0.85
umet_pi_clad_rad = 0.965
umet_pi_pitch = 2.9170
umet_lat_ele_x = 22
umet_lat_ele_y = 22

bp_incircle_rad = 36
bp_diagonal_intersect = float(np.sqrt(2)*bp_incircle_rad)
lg_incircle_rad = 38
lg_diagonal_intersect = float(np.sqrt(2)*lg_incircle_rad)
ug_incircle_rad = 42
ug_diagonal_intersect = float(np.sqrt(2)*ug_incircle_rad)

ta_rad = 65

# Radial surfaces
uo2_pi_fu_surf    = mp.ZCylinder(r = uo2_pi_fu_rad  ,  comment = "UO2 Pin Fuel")
uo2_pi_clad_surf  = mp.ZCylinder(r = uo2_pi_clad_rad,  comment = "UO2 Pin Clad")
umet_pi_fu_surf   = mp.ZCylinder(r = umet_pi_fu_rad  , comment = "UMET Pin Fuel")
umet_pi_clad_surf = mp.ZCylinder(r = umet_pi_clad_rad, comment = "UMET Pin Clad")
ta_surf           = mp.ZCylinder(r = ta_rad,           comment = "Tank")
deck += [
    uo2_pi_fu_surf, 
    uo2_pi_clad_surf, 
    umet_pi_fu_surf, 
    umet_pi_clad_surf, 
    ta_surf   
]

uo2_pin_xmin = mp.XPlane(x0=-0.5*uo2_pi_pitch, comment = 'UO2 Pin bound')
uo2_pin_xmax = mp.XPlane(x0= 0.5*uo2_pi_pitch, comment = 'UO2 Pin bound')
uo2_pin_ymin = mp.YPlane(y0=-0.5*uo2_pi_pitch, comment = 'UO2 Pin bound')
uo2_pin_ymax = mp.YPlane(y0= 0.5*uo2_pi_pitch, comment = 'UO2 Pin bound')
deck += [
    uo2_pin_xmin,
    uo2_pin_xmax,
    uo2_pin_ymin,
    uo2_pin_ymax
]
uo2_pin_region = +uo2_pin_xmin & -uo2_pin_xmax & +uo2_pin_ymin & -uo2_pin_ymax

umet_pin_xmin = mp.XPlane(x0=-0.5*umet_pi_pitch, comment = 'U_metal Pin bound')
umet_pin_xmax = mp.XPlane(x0= 0.5*umet_pi_pitch, comment = 'U_metal Pin bound')
umet_pin_ymin = mp.YPlane(y0=-0.5*umet_pi_pitch, comment = 'U_metal Pin bound')
umet_pin_ymax = mp.YPlane(y0= 0.5*umet_pi_pitch, comment = 'U_metal Pin bound')
deck += [
    umet_pin_xmin,
    umet_pin_xmax,
    umet_pin_ymin,
    umet_pin_ymax
]
umet_pin_region = +umet_pin_xmin & -umet_pin_xmax & +umet_pin_ymin & -umet_pin_ymax

# Assuming Regular Octagons
# Bottom Plate
bp_N = mp.YPlane(y0= bp_incircle_rad,                         comment = "Base Plate N")
bp_E = mp.XPlane(x0= bp_incircle_rad,                         comment = "Base Plate E")
bp_S = mp.YPlane(y0=-bp_incircle_rad,                         comment = "Base Plate S")
bp_W = mp.XPlane(x0=-bp_incircle_rad,                         comment = "Base Plate W")
bp_NE = mp.Plane(a= 1., b=1., c=0., d= bp_diagonal_intersect, comment = "Base Plate NE")
bp_SE = mp.Plane(a=-1., b=1., c=0., d=-bp_diagonal_intersect, comment = "Base Plate SE")
bp_SW = mp.Plane(a= 1., b=1., c=0., d=-bp_diagonal_intersect, comment = "Base Plate SW")
bp_NW = mp.Plane(a=-1., b=1., c=0., d= bp_diagonal_intersect, comment = "Base Plate NW")
bp_hex = -bp_N & +bp_S & +bp_W & -bp_E & -bp_NE & +bp_SW & -bp_NW & +bp_SE
bp_region = bp_hex & +ta_bot_surf & -bp_top_surf
deck += [
    bp_N,
    bp_E,
    bp_S,
    bp_W,
    bp_NE,
    bp_SE,
    bp_SW,
    bp_NW
]

# Lower Grid
lg_N = mp.YPlane(y0= lg_incircle_rad,                         comment = "Lower Grid N")
lg_E = mp.XPlane(x0= lg_incircle_rad,                         comment = "Lower Grid E")
lg_S = mp.YPlane(y0=-lg_incircle_rad,                         comment = "Lower Grid S")
lg_W = mp.XPlane(x0=-lg_incircle_rad,                         comment = "Lower Grid W")
lg_NE = mp.Plane(a= 1., b=1., c=0., d= lg_diagonal_intersect, comment = "Lower Grid NE")
lg_SE = mp.Plane(a=-1., b=1., c=0., d=-lg_diagonal_intersect, comment = "Lower Grid SE")
lg_SW = mp.Plane(a= 1., b=1., c=0., d=-lg_diagonal_intersect, comment = "Lower Grid SW")
lg_NW = mp.Plane(a=-1., b=1., c=0., d= lg_diagonal_intersect, comment = "Lower Grid NW")
lg_hex = -lg_N & +lg_S & +lg_W & -lg_E & -lg_NE & +lg_SW & -lg_NW & +lg_SE
lg_bot_region = lg_hex & +lg_bot_bot_surf & -lg_bot_top_surf
lg_cad_region = lg_hex & +lg_bot_top_surf & -lg_top_bot_surf
lg_top_region = lg_hex & +lg_top_bot_surf & -lg_top_top_surf
lg_region     = lg_hex & +lg_bot_bot_surf & -lg_top_top_surf
deck += [
    lg_N,
    lg_E,
    lg_S,
    lg_W,
    lg_NE,
    lg_SE,
    lg_SW,
    lg_NW
]

# Upper Grid
ug_N = mp.YPlane(y0= ug_incircle_rad,                         comment = "Upper Grid N")
ug_E = mp.XPlane(x0= ug_incircle_rad,                         comment = "Upper Grid E")
ug_S = mp.YPlane(y0=-ug_incircle_rad,                         comment = "Upper Grid S")
ug_W = mp.XPlane(x0=-ug_incircle_rad,                         comment = "Upper Grid W")
ug_NE = mp.Plane(a= 1., b=1., c=0., d= ug_diagonal_intersect, comment = "Upper Grid NE")
ug_SE = mp.Plane(a=-1., b=1., c=0., d=-ug_diagonal_intersect, comment = "Upper Grid SE")
ug_SW = mp.Plane(a= 1., b=1., c=0., d=-ug_diagonal_intersect, comment = "Upper Grid SW")
ug_NW = mp.Plane(a=-1., b=1., c=0., d= ug_diagonal_intersect, comment = "Upper Grid NW")
ug_hex = -ug_N & +ug_S & +ug_W & -ug_E & -ug_NE & +ug_SW & -ug_NW & +ug_SE
ug_bot_region = ug_hex & +ug_bot_bot_surf & -ug_bot_top_surf
ug_cad_region = ug_hex & +ug_bot_top_surf & -ug_top_bot_surf
ug_top_region = ug_hex & +ug_top_bot_surf & -ug_top_top_surf
ug_region     = ug_hex & +ug_bot_bot_surf & -ug_top_top_surf
deck += [
    ug_N,
    ug_E,
    ug_S,
    ug_W,
    ug_NE,
    ug_SE,
    ug_SW,
    ug_NW
]

grid_regions = lg_region | ug_region

# Cruciforms
inr_fu_crux      = cruciform_prism(distances=[3*uo2_pi_pitch, 6*uo2_pi_pitch, 9*uo2_pi_pitch, 11*uo2_pi_pitch], name="Inner Fuel")
out_inr_fu_crux  = cruciform_prism(distances=[umet_pi_pitch*2, umet_pi_pitch*4, umet_pi_pitch*6, umet_pi_pitch*7], name="Outer Inner Fuel")
out_out_fu_crux  = cruciform_prism(distances=[umet_pi_pitch*5,umet_pi_pitch*6,umet_pi_pitch*7,umet_pi_pitch*8,umet_pi_pitch*9,umet_pi_pitch*10,umet_pi_pitch*11], name="Outer Outer Fuel")
inr_fuel_region =  inr_fu_crux.internal_region
out_fuel_region = (~out_inr_fu_crux.internal_region & out_out_fu_crux.internal_region)
fuel_regions = inr_fuel_region | out_fuel_region
deck += inr_fu_crux.surfaces
deck += out_inr_fu_crux.surfaces
deck += out_out_fu_crux.surfaces

# Cells
# Create the cells for the uo2 pin cell
uo2_pin_cells = [
    mp.Cell(region=(-fu_bot_surf     & -uo2_pi_clad_surf),                                       material=Clad_non_fuel,    density=Clad_non_fuel_density,    comment='UO2 Pin Bottom'),
    mp.Cell(region=(+fu_bot_surf     & -fu_top_surf      & -uo2_pi_fu_surf),                     material=UO2,              density=UO2_density,              comment='UO2 Pin Fuel'),
    mp.Cell(region=(+fu_top_surf     & -uo2_fg_top_surf  & -uo2_pi_fu_surf),                     material=Filler_Gas,       density=Filler_Gas_density,       comment='UO2 Pin Filler Gas'),
    mp.Cell(region=(+fu_bot_surf     & -uo2_fg_top_surf  & +uo2_pi_fu_surf & -uo2_pi_clad_surf), material=Clad_uo2,         density=Clad_uo2_density,         comment='UO2 Pin Clad'),
    mp.Cell(region=(+uo2_fg_top_surf & -ta_top_surf      & -uo2_pi_clad_surf),                   material=Clad_non_fuel,    density=Clad_non_fuel_density,    comment='UO2 Pin Top'),
    mp.Cell(region=(+lg_bot_bot_surf & -lg_bot_top_surf  & +uo2_pi_clad_surf),                   material=Support_Material, density=Support_Material_density, comment='UO2 Pin Lower Grid Bottom'),
    mp.Cell(region=(+lg_bot_top_surf & -lg_top_bot_surf  & +uo2_pi_clad_surf),                   material=Cadmium,          density=Cadmium_density,          comment='UO2 Pin Lower Grid Cadmium'),
    mp.Cell(region=(+lg_top_bot_surf & -lg_top_top_surf  & +uo2_pi_clad_surf),                   material=Support_Material, density=Support_Material_density, comment='UO2 Pin Lower Grid Top'),
    mp.Cell(region=(+ug_bot_bot_surf & -ug_bot_top_surf  & +uo2_pi_clad_surf),                   material=Support_Material, density=Support_Material_density, comment='UO2 Pin Upper Grid Bottom'),
    mp.Cell(region=(+ug_bot_top_surf & -ug_top_bot_surf  & +uo2_pi_clad_surf),                   material=Cadmium,          density=Cadmium_density,          comment='UO2 Pin Upper Grid Cadmium'),
    mp.Cell(region=(+ug_top_bot_surf & -ug_top_top_surf  & +uo2_pi_clad_surf),                   material=Support_Material, density=Support_Material_density, comment='UO2 Pin Upper Grid Top'),
    mp.Cell(region=(-wa_top_surf     & ~grid_regions     & +uo2_pi_clad_surf),                   material=Moderator,        density=Moderator_density,        comment='UO2 Pin Moderator'),
    mp.Cell(region=(+wa_top_surf     & ~grid_regions     & +uo2_pi_clad_surf),                   material=Air,              density=Air_density,              comment='UO2 Pin Air')
]
uo2_u = mp.UniverseList(name=next_id(), cells=uo2_pin_cells)
deck += uo2_pin_cells

# Create the cells for the umet pin cell
Umet_pin_cells = [
    mp.Cell(region=(-fu_bot_surf      & -umet_pi_clad_surf),                                         material=Clad_non_fuel,    density=Clad_non_fuel_density,    comment="Umet Pin Bottom"),
    mp.Cell(region=(+fu_bot_surf      & -fu_top_surf       & -umet_pi_fu_surf),                      material=Umetal,           density=Umetal_density,           comment="Umet Pin Fuel"),
    mp.Cell(region=(+fu_top_surf      & -umet_fg_top_surf  & -umet_pi_fu_surf),                      material=Filler_Gas,       density=Filler_Gas_density,       comment="Umet Pin Filler Gas"),
    mp.Cell(region=(+fu_bot_surf      & -umet_fg_top_surf  & +umet_pi_fu_surf & -umet_pi_clad_surf), material=Clad_umetal,      density=Clad_umetal_density,      comment="Umet Pin Clad"),
    mp.Cell(region=(+umet_fg_top_surf & -ta_top_surf       & -umet_pi_clad_surf),                    material=Clad_non_fuel,    density=Clad_non_fuel_density,    comment="Umet Pin Top"),
    mp.Cell(region=(+lg_bot_bot_surf  & -lg_bot_top_surf   & +umet_pi_clad_surf),                    material=Support_Material, density=Support_Material_density, comment="Umet Pin Lower Grid Bottom"),
    mp.Cell(region=(+lg_bot_top_surf  & -lg_top_bot_surf   & +umet_pi_clad_surf),                    material=Cadmium,          density=Cadmium_density,          comment="Umet Pin Lower Grid Cadmium"),
    mp.Cell(region=(+lg_top_bot_surf  & -lg_top_top_surf   & +umet_pi_clad_surf),                    material=Support_Material, density=Support_Material_density, comment="Umet Pin Lower Grid Top"),
    mp.Cell(region=(+ug_bot_bot_surf  & -ug_bot_top_surf   & +umet_pi_clad_surf),                    material=Support_Material, density=Support_Material_density, comment="Umet Pin Upper Grid Bottom"),
    mp.Cell(region=(+ug_bot_top_surf  & -ug_top_bot_surf   & +umet_pi_clad_surf),                    material=Cadmium,          density=Cadmium_density,          comment="Umet Pin Upper Grid Cadmium"),
    mp.Cell(region=(+ug_top_bot_surf  & -ug_top_top_surf   & +umet_pi_clad_surf),                    material=Support_Material, density=Support_Material_density, comment="Umet Pin Upper Grid Top"),
    mp.Cell(region=(-wa_top_surf      & ~grid_regions      & +umet_pi_clad_surf),                    material=Moderator,        density=Moderator_density,        comment='Umet Pin Moderator'),
    mp.Cell(region=(+wa_top_surf      & ~grid_regions      & +umet_pi_clad_surf),                    material=Air,              density=Air_density,              comment='Umet Pin Air')
]
ume_u = mp.UniverseList(name=next_id(), cells=Umet_pin_cells)
deck += Umet_pin_cells

# Create the cells for the empty umet pin cell
Umet_empty_pin_cells = [
    mp.Cell(region=(+lg_bot_bot_surf & -lg_bot_top_surf), material=Support_Material, density=Support_Material_density, comment="Umet Empty Pin Lower Grid Bottom"),
    mp.Cell(region=(+lg_bot_top_surf & -lg_top_bot_surf), material=Cadmium,          density=Cadmium_density,          comment="Umet Empty Pin Lower Grid Cadmium"),
    mp.Cell(region=(+lg_top_bot_surf & -lg_top_top_surf), material=Support_Material, density=Support_Material_density, comment="Umet Empty Pin Lower Grid Top"),
    mp.Cell(region=(+ug_bot_bot_surf & -ug_bot_top_surf), material=Support_Material, density=Support_Material_density, comment="Umet Empty Pin Upper Grid Bottom"),
    mp.Cell(region=(+ug_bot_top_surf & -ug_top_bot_surf), material=Cadmium,          density=Cadmium_density,          comment="Umet Empty Pin Upper Grid Cadmium"),
    mp.Cell(region=(+ug_top_bot_surf & -ug_top_top_surf), material=Support_Material, density=Support_Material_density, comment="Umet Empty Pin Upper Grid Top"),
    mp.Cell(region=(-wa_top_surf & ~grid_regions),        material=Moderator,        density=Moderator_density,        comment='Umet Empty Pin Moderator'),
    mp.Cell(region=(+wa_top_surf & ~grid_regions),        material=Air,              density=Air_density,              comment='Umet Empty Pin Air')
]
emp_u = mp.UniverseList(name=next_id(), cells=Umet_empty_pin_cells)
deck += Umet_empty_pin_cells

# Create the cells that are outside of the lattices (includes the space between the lattices)
not_fuel_not_support_region = ~(fuel_regions | bp_region | lg_region | ug_region)
non_lattice_cells = [
    mp.Cell(region=(bp_region),                                                            material=Support_Material, density=Support_Material_density, comment="Non-Lattice Base Plate"),
    mp.Cell(region=(lg_bot_region & ~fuel_regions),                                        material=Support_Material, density=Support_Material_density, comment="Non-Lattice Lower Grid Bottom"),
    mp.Cell(region=(lg_cad_region & ~fuel_regions),                                        material=Cadmium,          density=Cadmium_density,          comment="Non-Lattice Lower Grid Cadmium"),
    mp.Cell(region=(lg_top_region & ~fuel_regions),                                        material=Support_Material, density=Support_Material_density, comment="Non-Lattice Lower Grid Top"),
    mp.Cell(region=(ug_bot_region & ~fuel_regions),                                        material=Support_Material, density=Support_Material_density, comment="Non-Lattice Upper Grid Bottom"),
    mp.Cell(region=(ug_cad_region & ~fuel_regions),                                        material=Cadmium,          density=Cadmium_density,          comment="Non-Lattice Upper Grid Cadmium"),
    mp.Cell(region=(ug_top_region & ~fuel_regions),                                        material=Support_Material, density=Support_Material_density, comment="Non-Lattice Upper Grid Top"),
    mp.Cell(region=(-ta_surf & +ta_bot_surf & -wa_top_surf & not_fuel_not_support_region), material=Moderator,        density=Moderator_density,        comment="Non-Lattice Moderator"),
    mp.Cell(region=(-ta_surf & +wa_top_surf & -ta_top_surf & not_fuel_not_support_region), material=Air,              density=Air_density,              comment="Non-Lattice Air")
]
deck += non_lattice_cells

# Create the UO2 fuel lattice
UO2_trans = mp.Transform([0.5*uo2_pi_pitch, 0.5*uo2_pi_pitch, 0])
UO2_numpy_lattice = np.tile(uo2_u.name, [1, 22, 22])
UO2_lattice = mp.Lattice(i=[-10,11], j=[-10,11], k=[0,0], lattice=UO2_numpy_lattice, type='REC', 
                         universes={uo2_u.name:uo2_u})
UO2_lattice_pin_cell = mp.Cell(region=uo2_pin_region, fill=UO2_lattice, transform=UO2_trans, comment='UO2 Lattice Pin Cell')
deck += UO2_lattice_pin_cell
UO2_lattice_universe = mp.UniverseList(name=next_id(), cells=UO2_lattice_pin_cell)
# UO2_lattice_cell = mp.Cell(region=(inr_fuel_region & +bp_top_surf & -ta_top_surf), fill=UO2_lattice_universe, transformation=UO2_trans, comment='UO2 Lattice')
UO2_lattice_cell = mp.Cell(region=(inr_fuel_region & +bp_top_surf & -ta_top_surf), fill=UO2_lattice_universe, comment='UO2 Lattice')
deck += UO2_lattice_cell

# Create the U Metal fuel lattice
Umet_trans = mp.Transform([0.5*umet_pi_pitch, 0.5*umet_pi_pitch, 0])
Umet_numpy_lattice = np.array([[
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name],
        [emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name],
        [emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name],
        [emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name],
        [emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name],
        [emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name],
        [emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name],
        [emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, ume_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
        [emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name, emp_u.name],
    ]])
Umet_lattice = mp.Lattice(i=[-10,11], j=[-10,11], k=[0,0], lattice=Umet_numpy_lattice, type='REC', 
                          universes={
                              ume_u.name:ume_u,
                              emp_u.name:emp_u
                            })
Umet_lattice_pin_cell = mp.Cell(region=umet_pin_region, fill=Umet_lattice, transform=Umet_trans, comment='U metal Lattice Pin Cell')
deck += Umet_lattice_pin_cell
Umet_lattice_universe = mp.UniverseList(name=next_id(), cells=Umet_lattice_pin_cell)
Umet_lattice_cell = mp.Cell(region=(~out_inr_fu_crux.internal_region & out_out_fu_crux.internal_region & +bp_top_surf & -ta_top_surf), fill=Umet_lattice_universe, comment='U Metal Lattice')
deck += Umet_lattice_cell

for cell in deck.cells.values():
    cell.importances = {'n' : 1.0}

deck += mp.Cell(name=999, region=(-ta_bot_surf | +ta_top_surf | +ta_surf), importances={'n':0}, comment='Outer Void')

#===============================================================================
# Settings
#===============================================================================

deck += mp.CriticalitySource(histories=1e6, keff_guess=1.0,
                             skip_cycles=50, cycles=650)

coords = [
    ( 0.5*uo2_pi_pitch,  0.5*uo2_pi_pitch, fu_bot_dist + 1),
    (-0.5*uo2_pi_pitch,  0.5*uo2_pi_pitch, fu_bot_dist + 1),
    ( 0.5*uo2_pi_pitch, -0.5*uo2_pi_pitch, fu_bot_dist + 1),
    (-0.5*uo2_pi_pitch, -0.5*uo2_pi_pitch, fu_bot_dist + 1),
]
deck += mp.CriticalitySourcePoints(coords)

deck.write('mcnpy_model.mcnp')

print("Don't directly use the deck from this script.")
print("Major changes need to be made to the input file.")
print("You need to manually change the SAB file to the correct ENDF8.0 version.")
print("You must also add the library number at the end of each ZAID.")