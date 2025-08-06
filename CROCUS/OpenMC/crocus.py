# Description of the model : https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/nea4440-crocus.pdf

import openmc
import numpy as np
import openmc.model
import openmc.source
import openmc.stats

model = openmc.Model()

#===============================================================================
# Materials
#===============================================================================

Reactor_Temperature = 20 + 273.15 # Kelvin

UO2 = openmc.Material(name="UO2")
UO2.add_nuclide('U235', percent=4.30565e-04, percent_type='ao')
UO2.add_nuclide('O16 ', percent=4.70902e-02, percent_type='ao')
UO2.add_nuclide('U238', percent=2.31145e-02, percent_type='ao')
UO2.set_density(units='sum')
UO2.temperature = Reactor_Temperature

Umetal = openmc.Material(name="UMetal")
Umetal.add_nuclide('U235', percent=4.53160e-04, percent_type='ao')
Umetal.add_nuclide('U238', percent=4.68003e-02, percent_type='ao')
Umetal.set_density(units="sum")
Umetal.temperature = Reactor_Temperature

Clad_non_fuel = openmc.Material(name="Cladding")
Clad_non_fuel.add_element("Al", percent=6.02611e-02, percent_type='ao')
Clad_non_fuel.set_density(units="sum")
Clad_non_fuel.temperature = Reactor_Temperature

Clad_uo2 = openmc.Material(name="Cladding")
Clad_uo2.add_element("Al", percent=5.00614e-02, percent_type='ao')
Clad_uo2.set_density(units="sum")
Clad_uo2.temperature = Reactor_Temperature

Clad_umetal = openmc.Material(name="Cladding")
Clad_umetal.add_element("Al", percent=5.17799e-2, percent_type='ao')
Clad_umetal.set_density(units="sum")
Clad_umetal.temperature = Reactor_Temperature

Moderator = openmc.Material(name="Moderator")
Moderator.add_nuclide("H1" , percent=6.67578e-02, percent_type='ao')
Moderator.add_nuclide("O16", percent=3.33789e-02, percent_type='ao')
Moderator.set_density(units="sum")
Moderator.add_s_alpha_beta('c_H_in_H2O')

Support_Material = openmc.Material(name="Support Material")
Support_Material.add_element("Al", percent=6.02611e-02, percent_type='ao')
Support_Material.set_density(units="sum")
Support_Material.temperature = Reactor_Temperature

Cadmium = openmc.Material(name="Cadmium")
Cadmium.add_element("Cd", percent=4.63334e-02, percent_type='ao')
Cadmium.set_density(units="sum")
Cadmium.temperature = Reactor_Temperature

Filler_Gas = openmc.Material(name="Filler Gas")
Filler_Gas.add_nuclide("He4", percent=1.64e-04, percent_type='ao')
Filler_Gas.set_density(units="sum")
Filler_Gas.temperature = Reactor_Temperature

Air = openmc.Material(name="Air")
Air.add_nuclide('N14', percent=4.14e-05, percent_type='ao')
Air.add_nuclide('O16', percent=9.07e-06, percent_type='ao')
Air.set_density(units="sum")
Air.temperature = Reactor_Temperature

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
ta_bot_surf      = openmc.ZPlane(ta_bot_dist, boundary_type = 'vacuum')
bp_top_surf      = openmc.ZPlane(bp_top_dist)
lg_bot_bot_surf  = openmc.ZPlane(lg_bot_bot_dist)
lg_bot_top_surf  = openmc.ZPlane(lg_bot_top_dist)
lg_top_bot_surf  = openmc.ZPlane(lg_top_bot_dist)
lg_top_top_surf  = openmc.ZPlane(lg_top_top_dist)
ug_bot_bot_surf  = openmc.ZPlane(ug_bot_bot_dist)
ug_bot_top_surf  = openmc.ZPlane(ug_bot_top_dist)
ug_top_bot_surf  = openmc.ZPlane(ug_top_bot_dist)
ug_top_top_surf  = openmc.ZPlane(ug_top_top_dist)
wa_top_surf      = openmc.ZPlane(wa_H0 + wa_height)
pi_bot_surf      = openmc.ZPlane(pi_bot_dist)
fu_bot_surf      = openmc.ZPlane(fu_bot_dist)
fu_top_surf      = openmc.ZPlane(fu_top_dist)
uo2_fg_top_surf  = openmc.ZPlane(uo2_fg_top_dist)
umet_fg_top_surf = openmc.ZPlane(umet_fg_top_dist)
pi_top_surf      = openmc.ZPlane(pi_top_dist, boundary_type = 'vacuum')
ta_top_surf      = openmc.ZPlane(ta_top_dist, boundary_type = 'vacuum')

# Radial distances from tank center
uo2_pi_fu_rad = 1.052/2.0
uo2_pi_clad_rad = uo2_pi_fu_rad + 0.085
uo2_pi_pitch = 1.8370

umet_pi_fu_rad = 1.7/2
umet_pi_clad_rad = umet_pi_fu_rad + 0.0975
umet_pi_pitch = 2.9170

ta_rad = 65

# Radial surfaces
uo2_pi_fu_surf = openmc.ZCylinder(r = uo2_pi_fu_rad)
uo2_pi_clad_surf = openmc.ZCylinder(r = uo2_pi_clad_rad)

umet_pi_fu_surf = openmc.ZCylinder(r = umet_pi_fu_rad)
umet_pi_clad_surf = openmc.ZCylinder(r = umet_pi_clad_rad)

ta_surf = openmc.ZCylinder(r = ta_rad, boundary_type = 'vacuum')

# Bottom Plate
bp_N  = openmc.YPlane(36.)
bp_NE = openmc.Plane(1., 1., 0., 54.)
bp_E  = openmc.XPlane(36.)
bp_SE = openmc.Plane(-1., 1., 0., -54.)
bp_S  = openmc.YPlane(-36.)
bp_SW = openmc.Plane(1., 1., 0., -54.)
bp_W  = openmc.XPlane(-36.)
bp_NW = openmc.Plane(-1., 1., 0., 54.)
bp_hex = -bp_N & +bp_S & +bp_W & -bp_E & -bp_NE & +bp_SW & -bp_NW & +bp_SE
bp_region = bp_hex & +ta_bot_surf & -bp_top_surf

# Lower Grid
lg_N  = openmc.YPlane(38.)
lg_NE = openmc.Plane(1., 1., 0., 57.)
lg_E  = openmc.XPlane(38.)
lg_SE = openmc.Plane(-1., 1., 0., -57.)
lg_S  = openmc.YPlane(-38.)
lg_SW = openmc.Plane(1., 1., 0., -57.)
lg_W  = openmc.XPlane(-38.)
lg_NW = openmc.Plane(-1., 1., 0., 57.)
lg_hex = -lg_N & +lg_S & +lg_W & -lg_E & -lg_NE & +lg_SW & -lg_NW & +lg_SE
lg_region = lg_hex & +lg_bot_bot_surf & -lg_top_top_surf

# Upper Grid
ug_N  = openmc.YPlane(42.)
ug_NE = openmc.Plane(1., 1., 0., 63.)
ug_E  = openmc.XPlane(42.)
ug_SE = openmc.Plane(-1., 1., 0., -63.)
ug_S  = openmc.YPlane(-42.)
ug_SW = openmc.Plane(1., 1., 0., -63.)
ug_W  = openmc.XPlane(-42.)
ug_NW = openmc.Plane(-1., 1., 0., 63.)
ug_hex = -ug_N & +ug_S & +ug_W & -ug_E & -ug_NE & +ug_SW & -ug_NW & +ug_SE
ug_region = ug_hex & +ug_bot_bot_surf & -ug_top_top_surf

# Cruciforms
inr_fu_crux_surf      = openmc.model.CruciformPrism(distances=[3*uo2_pi_pitch, 6*uo2_pi_pitch, 9*uo2_pi_pitch, 11*uo2_pi_pitch])
out_inr_fu_crux_surf  = openmc.model.CruciformPrism(distances=[umet_pi_pitch*2, umet_pi_pitch*4, umet_pi_pitch*6, umet_pi_pitch*7])
out_out_crux_surf     = openmc.model.CruciformPrism(distances=[umet_pi_pitch*5,umet_pi_pitch*6,umet_pi_pitch*7,umet_pi_pitch*8,umet_pi_pitch*9,umet_pi_pitch*10,umet_pi_pitch*11])

## Cells
# Create the cells for the uo2 pin cell
uo2_pi_bot_cell     = openmc.Cell(region=(+pi_bot_surf     & -fu_bot_surf     & -uo2_pi_clad_surf),                   fill = Clad_non_fuel)
uo2_pi_fu_cell      = openmc.Cell(region=(+fu_bot_surf     & -fu_top_surf     & -uo2_pi_fu_surf),                     fill = UO2)
uo2_pi_fg_cell      = openmc.Cell(region=(+fu_top_surf     & -uo2_fg_top_surf & -uo2_pi_fu_surf),                     fill = Filler_Gas)
uo2_pi_clad_cell    = openmc.Cell(region=(+fu_bot_surf     & -uo2_fg_top_surf & +uo2_pi_fu_surf & -uo2_pi_clad_surf), fill = Clad_uo2)
uo2_pi_top_cell     = openmc.Cell(region=(+uo2_fg_top_surf & -pi_top_surf     & -uo2_pi_clad_surf),                   fill = Clad_non_fuel)
uo2_pi_lg_bot_cell  = openmc.Cell(region=(+lg_bot_bot_surf & -lg_bot_top_surf & +uo2_pi_clad_surf),                   fill = Support_Material)
uo2_pi_lg_clad_cell = openmc.Cell(region=(+lg_bot_top_surf & -lg_top_bot_surf & +uo2_pi_clad_surf),                   fill = Cadmium)
uo2_pi_lg_top_cell  = openmc.Cell(region=(+lg_top_bot_surf & -lg_top_top_surf & +uo2_pi_clad_surf),                   fill = Support_Material)
uo2_pi_ug_bot_cell  = openmc.Cell(region=(+ug_bot_bot_surf & -ug_bot_top_surf & +uo2_pi_clad_surf),                   fill = Support_Material)
uo2_pi_ug_clad_cell = openmc.Cell(region=(+ug_bot_top_surf & -ug_top_bot_surf & +uo2_pi_clad_surf),                   fill = Cadmium)
uo2_pi_ug_top_cell  = openmc.Cell(region=(+ug_top_bot_surf & -ug_top_top_surf & +uo2_pi_clad_surf),                   fill = Support_Material)
uo2_pi_bot_mod_cell = openmc.Cell(region=(+bp_top_surf     & -lg_bot_bot_surf & +uo2_pi_clad_surf),                   fill = Moderator)
uo2_pi_top_mod_cell = openmc.Cell(region=(+lg_top_top_surf & -wa_top_surf     & +uo2_pi_clad_surf),                   fill = Moderator)
uo2_pi_bot_air_cell = openmc.Cell(region=(+wa_top_surf     & -ug_bot_bot_surf & +uo2_pi_clad_surf),                   fill = Air)
uo2_pi_top_air_cell = openmc.Cell(region=(+ug_top_top_surf & -pi_top_surf     & +uo2_pi_clad_surf),                   fill = Air)

# Create the up2 pin cell universe
uo2_pi_universe = openmc.Universe(cells=[
                                            uo2_pi_bot_cell,
                                            uo2_pi_fu_cell,
                                            uo2_pi_fg_cell,
                                            uo2_pi_clad_cell,
                                            uo2_pi_top_cell,
                                            uo2_pi_lg_bot_cell,
                                            uo2_pi_lg_clad_cell,
                                            uo2_pi_lg_top_cell,
                                            uo2_pi_ug_bot_cell,
                                            uo2_pi_ug_clad_cell,
                                            uo2_pi_ug_top_cell,
                                            uo2_pi_bot_mod_cell,
                                            uo2_pi_top_mod_cell,
                                            uo2_pi_bot_air_cell,
                                            uo2_pi_top_air_cell,
                                        ])

# Create the cells for the umet pin cell
umet_pi_bot_cell     = openmc.Cell(region=(+pi_bot_surf      & -fu_bot_surf      & -umet_pi_clad_surf),                    fill = Clad_non_fuel)
umet_pi_fu_cell      = openmc.Cell(region=(+fu_bot_surf      & -fu_top_surf      & -umet_pi_fu_surf),                      fill = Umetal)
umet_pi_fg_cell      = openmc.Cell(region=(+fu_top_surf      & -umet_fg_top_surf & -umet_pi_fu_surf),                      fill = Filler_Gas)
umet_pi_clad_cell    = openmc.Cell(region=(+fu_bot_surf      & -umet_fg_top_surf & +umet_pi_fu_surf & -umet_pi_clad_surf), fill = Clad_umetal)
umet_pi_top_cell     = openmc.Cell(region=(+umet_fg_top_surf & -pi_top_surf      & -umet_pi_clad_surf),                    fill = Clad_non_fuel)
umet_pi_lg_bot_cell  = openmc.Cell(region=(+lg_bot_bot_surf  & -lg_bot_top_surf  & +umet_pi_clad_surf),                    fill = Support_Material)
umet_pi_lg_clad_cell = openmc.Cell(region=(+lg_bot_top_surf  & -lg_top_bot_surf  & +umet_pi_clad_surf),                    fill = Cadmium)
umet_pi_lg_top_cell  = openmc.Cell(region=(+lg_top_bot_surf  & -lg_top_top_surf  & +umet_pi_clad_surf),                    fill = Support_Material)
umet_pi_ug_bot_cell  = openmc.Cell(region=(+ug_bot_bot_surf  & -ug_bot_top_surf  & +umet_pi_clad_surf),                    fill = Support_Material)
umet_pi_ug_clad_cell = openmc.Cell(region=(+ug_bot_top_surf  & -ug_top_bot_surf  & +umet_pi_clad_surf),                    fill = Cadmium)
umet_pi_ug_top_cell  = openmc.Cell(region=(+ug_top_bot_surf  & -ug_top_top_surf  & +umet_pi_clad_surf),                    fill = Support_Material)
umet_pi_bot_mod_cell = openmc.Cell(region=(+bp_top_surf      & -lg_bot_bot_surf  & +umet_pi_clad_surf),                    fill = Moderator)
umet_pi_top_mod_cell = openmc.Cell(region=(+lg_top_top_surf  & -wa_top_surf      & +umet_pi_clad_surf),                    fill = Moderator)
umet_pi_bot_air_cell = openmc.Cell(region=(+wa_top_surf      & -ug_bot_bot_surf  & +umet_pi_clad_surf),                    fill = Air)
umet_pi_top_air_cell = openmc.Cell(region=(+ug_top_top_surf  & -pi_top_surf      & +umet_pi_clad_surf),                    fill = Air)

# Create the umet pin cell universe
umet_pi_universe = openmc.Universe(cells=[
                                            umet_pi_bot_cell,
                                            umet_pi_fu_cell,
                                            umet_pi_fg_cell,
                                            umet_pi_clad_cell,
                                            umet_pi_top_cell,
                                            umet_pi_lg_bot_cell,
                                            umet_pi_lg_clad_cell,
                                            umet_pi_lg_top_cell,
                                            umet_pi_ug_bot_cell,
                                            umet_pi_ug_clad_cell,
                                            umet_pi_ug_top_cell,
                                            umet_pi_bot_mod_cell,
                                            umet_pi_top_mod_cell,
                                            umet_pi_bot_air_cell,
                                            umet_pi_top_air_cell,
                                        ])

# Create the cells for the empty umet pin cell
umet_empty_pi_mod_cell     = openmc.Cell(region=(+pi_bot_surf     & -wa_top_surf     & -umet_pi_clad_surf),                    fill = Moderator)
umet_empty_pi_air_cell     = openmc.Cell(region=(+wa_top_surf     & -ta_top_surf     & -umet_pi_clad_surf),                    fill = Air)
umet_empty_pi_lg_bot_cell  = openmc.Cell(region=(+lg_bot_bot_surf & -lg_bot_top_surf & +umet_pi_clad_surf),                    fill = Support_Material)
umet_empty_pi_lg_clad_cell = openmc.Cell(region=(+lg_bot_top_surf & -lg_top_bot_surf & +umet_pi_clad_surf),                    fill = Cadmium)
umet_empty_pi_lg_top_cell  = openmc.Cell(region=(+lg_top_bot_surf & -lg_top_top_surf & +umet_pi_clad_surf),                    fill = Support_Material)
umet_empty_pi_ug_bot_cell  = openmc.Cell(region=(+ug_bot_bot_surf & -ug_bot_top_surf & +umet_pi_clad_surf),                    fill = Support_Material)
umet_empty_pi_ug_clad_cell = openmc.Cell(region=(+ug_bot_top_surf & -ug_top_bot_surf & +umet_pi_clad_surf),                    fill = Cadmium)
umet_empty_pi_ug_top_cell  = openmc.Cell(region=(+ug_top_bot_surf & -ug_top_top_surf & +umet_pi_clad_surf),                    fill = Support_Material)
umet_empty_pi_bot_mod_cell = openmc.Cell(region=(+bp_top_surf     & -lg_bot_bot_surf & +umet_pi_clad_surf),                    fill = Moderator)
umet_empty_pi_top_mod_cell = openmc.Cell(region=(+lg_top_top_surf & -wa_top_surf     & +umet_pi_clad_surf),                    fill = Moderator)
umet_empty_pi_bot_air_cell = openmc.Cell(region=(+wa_top_surf     & -ug_bot_bot_surf & +umet_pi_clad_surf),                    fill = Air)
umet_empty_pi_top_air_cell = openmc.Cell(region=(+ug_top_top_surf & -pi_top_surf     & +umet_pi_clad_surf),                    fill = Air)

# Create the empty umet pin cell universe
umet_empty_pi_universe = openmc.Universe(cells=[
                                            umet_empty_pi_mod_cell,
                                            umet_empty_pi_air_cell,
                                            umet_empty_pi_lg_bot_cell,
                                            umet_empty_pi_lg_clad_cell,
                                            umet_empty_pi_lg_top_cell,
                                            umet_empty_pi_ug_bot_cell,
                                            umet_empty_pi_ug_clad_cell,
                                            umet_empty_pi_ug_top_cell,
                                            umet_empty_pi_bot_mod_cell,
                                            umet_empty_pi_top_mod_cell,
                                            umet_empty_pi_bot_air_cell,
                                            umet_empty_pi_top_air_cell,
                                        ])

# Create the cells that are outisde of the lattices (includes that space between the lattices)
bp_cell = openmc.Cell(region=(+ta_bot_surf & - bp_top_surf & bp_hex), fill=Support_Material)
out_lg_bot_cell  = openmc.Cell(region=(+lg_bot_bot_surf & -lg_bot_top_surf & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf)) & lg_hex), fill = Support_Material)
out_lg_clad_cell = openmc.Cell(region=(+lg_bot_top_surf & -lg_top_bot_surf & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf)) & lg_hex), fill = Cadmium)
out_lg_top_cell  = openmc.Cell(region=(+lg_top_bot_surf & -lg_top_top_surf & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf)) & lg_hex), fill = Support_Material)
out_ug_bot_cell  = openmc.Cell(region=(+ug_bot_bot_surf & -ug_bot_top_surf & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf)) & ug_hex), fill = Support_Material)
out_ug_clad_cell = openmc.Cell(region=(+ug_bot_top_surf & -ug_top_bot_surf & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf)) & ug_hex), fill = Cadmium)
out_ug_top_cell  = openmc.Cell(region=(+ug_top_bot_surf & -ug_top_top_surf & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf)) & ug_hex), fill = Support_Material)
out_mod_cell = openmc.Cell(region=(-ta_surf & +ta_bot_surf & -wa_top_surf & ~bp_region & ~lg_region & ~ug_region & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf))), fill = Moderator)
out_air_cell = openmc.Cell(region=(-ta_surf & +wa_top_surf & -ta_top_surf & ~bp_region & ~lg_region & ~ug_region & (+out_out_crux_surf | (-out_inr_fu_crux_surf & +inr_fu_crux_surf))), fill = Air)

# Create dummy universe
dummy_universe = openmc.Universe(name = "Dummy")
dummy_cell = openmc.Cell(name="Dummy")
dummy_universe.add_cell(dummy_cell)

# Create aliases for lattice definitions
uo2 = uo2_pi_universe
met = umet_pi_universe
emp = umet_empty_pi_universe
dum = dummy_universe

# Create the uo2 fuel lattice
inr_fuel_lattice = openmc.RectLattice()
inr_fuel_lattice.lower_left = [-11*uo2_pi_pitch, -11*uo2_pi_pitch, -ta_top_dist]
inr_fuel_lattice.pitch = [uo2_pi_pitch, uo2_pi_pitch, 2*ta_top_dist]
inr_fuel_lattice.universes = np.tile(uo2_pi_universe, [1, 22, 22])
inr_fuel_cell = openmc.Cell(region = (-inr_fu_crux_surf & +pi_bot_surf & -ta_top_surf), fill = inr_fuel_lattice)

# Create the umet fuel lattice
out_fuel_lattice = openmc.RectLattice()
out_fuel_lattice.lower_left = [-11*umet_pi_pitch, -11*umet_pi_pitch, -ta_top_dist]
out_fuel_lattice.pitch = [umet_pi_pitch, umet_pi_pitch, 2*ta_top_dist]
out_fuel_lattice.universes = \
    [[
        [dum, dum, dum, dum, dum, dum, emp, emp, emp, emp, emp, emp, emp, emp, emp, emp, dum, dum, dum, dum, dum, dum],
        [dum, dum, dum, dum, dum, emp, emp, emp, met, met, met, met, met, met, emp, emp, emp, dum, dum, dum, dum, dum],
        [dum, dum, dum, dum, emp, emp, met, met, met, met, met, met, met, met, met, met, emp, emp, dum, dum, dum, dum],
        [dum, dum, dum, emp, met, met, met, met, met, met, met, met, met, met, met, met, met, met, emp, dum, dum, dum],
        [dum, dum, emp, met, met, met, met, met, met, dum, dum, dum, dum, met, met, met, met, met, met, emp, dum, dum],
        [dum, emp, emp, met, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, met, emp, emp, dum],
        [emp, emp, met, met, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, met, met, emp, emp],
        [emp, emp, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, emp, emp],
        [emp, met, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, met, emp],
        [emp, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, emp],
        [emp, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, emp],
        [emp, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, emp],
        [emp, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, emp],
        [emp, met, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, met, emp],
        [emp, emp, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, emp, emp],
        [emp, emp, met, met, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, met, met, emp, emp],
        [dum, emp, emp, met, met, met, met, dum, dum, dum, dum, dum, dum, dum, dum, met, met, met, met, emp, emp, dum],
        [dum, dum, emp, met, met, met, met, met, met, dum, dum, dum, dum, met, met, met, met, met, met, emp, dum, dum],
        [dum, dum, dum, emp, met, met, met, met, met, met, met, met, met, met, met, met, met, met, emp, dum, dum, dum],
        [dum, dum, dum, dum, emp, emp, met, met, met, met, met, met, met, met, met, met, emp, emp, dum, dum, dum, dum],
        [dum, dum, dum, dum, dum, emp, emp, emp, met, met, met, met, met, met, emp, emp, emp, dum, dum, dum, dum, dum],
        [dum, dum, dum, dum, dum, dum, emp, emp, emp, emp, emp, emp, emp, emp, emp, emp, dum, dum, dum, dum, dum, dum],
    ]]
out_fuel_cell = openmc.Cell(region = (-out_out_crux_surf& +out_inr_fu_crux_surf & +pi_bot_surf & -ta_top_surf), fill = out_fuel_lattice)

# Add all cells to the root universe
Universe = openmc.Universe(cells = [
    bp_cell,
    out_lg_bot_cell,
    out_lg_clad_cell,
    out_lg_top_cell,
    out_ug_bot_cell,
    out_ug_clad_cell,
    out_ug_top_cell,
    out_mod_cell,
    out_air_cell,
    inr_fuel_cell,
    out_fuel_cell
])

# Fill in Geometry
model.geometry = openmc.Geometry(root = Universe)

#===============================================================================
# Settings
#===============================================================================

model.settings = openmc.Settings()

# Sources
source = openmc.IndependentSource()
source.space = openmc.stats.Box([-29.17, -29.17, 0.05], [29.17, 29.17, 96])
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt()
model.settings.source = [source]

# Tallies
flux_profile_tally_mesh = openmc.RegularMesh()
flux_profile_tally_mesh.lower_left = [-ta_rad, -ta_rad, fu_bot_dist]
flux_profile_tally_mesh.upper_right = [ta_rad, ta_rad, fu_top_dist]
flux_profile_tally_mesh.dimension = [500, 500, 100]
flux_profile_mesh_filter = openmc.MeshFilter(flux_profile_tally_mesh)

flux_profile_energy_filter = openmc.EnergyFilter(values=[1.E-5, 0.625, 20.E6])

flux_profile_tally = openmc.Tally(name='flux')
flux_profile_tally.scores = ['flux']
flux_profile_tally.filters = [flux_profile_energy_filter, flux_profile_mesh_filter]

flux_spectrum_tally = openmc.Tally(name='flux')
flux_spectrum_tally.scores = ['flux']

flux_spectrum_energy_filter = openmc.EnergyFilter(values=np.logspace(start=np.log10(1e-5), 
                                                       stop=np.log10(20e6), 
                                                       num=501))

flux_spectrum_tally.filters = [flux_spectrum_energy_filter]

model.tallies = openmc.Tallies([flux_profile_tally, flux_spectrum_tally])

# Entropy
model.settings.entropy_mesh = openmc.RegularMesh()
model.settings.entropy_mesh.lower_left = [-10*umet_pi_pitch, -10*umet_pi_pitch, fu_bot_dist]
model.settings.entropy_mesh.upper_right = [10*umet_pi_pitch, 10*umet_pi_pitch, fu_top_dist]
model.settings.entropy_mesh.dimension = [5, 5, 5]

# K-Eigenvalue
model.settings.run_mode = 'eigenvalue'
model.settings.particles = int(1e5)
model.settings.batches = 300
model.settings.inactive = 50

# Export/run
#=====================
model.export_to_xml()