import numpy as np

import mcnpy as mp
from mcnpy.elements import *

deck = mp.Deck()

# Materials
#=====================
fuel = mp.Material(
    U[235]*0.02644492 + 
    U[238]*0.85505247 + 
    O[16]*0.11850261
)
fuel_density = 10.3070
deck += [fuel]

burnable_poison = mp.Material(
     U[235]*0.00202753 + 
     U[238]*0.80898345 + 
      O[16]*0.11957908 +
    GD[152]*0.00013411 +
    GD[154]*0.00148108 +
    GD[155]*0.01012049 +
    GD[156]*0.01408805 +
    GD[157]*0.01083999 +
    GD[158]*0.01731511 +
    GD[160]*0.01543111 
)
burnable_poison_density = 10.3070
deck += [burnable_poison]

moderator = mp.Material(
     H[1]*6.663259e-1 +
    O[16]*3.330861e-1 +
    B[10]*7.186970e-5 +
    B[11]*2.892846e-4
)
moderator.s_alpha_beta = 'lwtr'
moderator_density = 0.79602
deck += [moderator]

helium = mp.Material(HE[4]@1)
helium_density = 0.0015981
deck += [helium]

fuel_cladding = mp.Material(
     O[16]*0.00125 +
    CR*0.00100 +
    FE*0.00210 +
    ZR*0.98115 +
    SN*0.01450 
)
fuel_cladding_density = 6.55
deck += [fuel_cladding]

control_rod_absorber = mp.Material(
    AG*0.80 +
    IN*0.15 +
    CD*0.05
)
control_rod_absorber_density = 10.16
deck += [control_rod_absorber]

control_rod_cladding = mp.Material(
    SI*0.0060 +
    CR*0.1900 +
    MN[55]*0.0200 +
    FE*0.6840 +
    NI*0.1000
)
control_rod_cladding_density = 8.03
deck += [control_rod_cladding]

# Geometry
#=====================

pitch = 1.26
fuel_rad = 0.3975
fuel_helium_gap_rad = 0.4125
fuel_clad_rad = 0.4750
control_rod_rad = 0.38227
control_rod_helium_rad = 0.38608
control_rod_clad_rad = 0.48387
guide_tube_inner_radius = 0.5725
guide_tube_outer_radius = 0.6125
assembly_length = 420
lattice_elements_x = 17
lattice_elements_y = 17
lattice_elements_z = 1

rod_insertion = 0
reactor_temperature = 565.0

if rod_insertion < 0 or rod_insertion > assembly_length:
    raise Exception(f"Rod insertion must be [0, {assembly_length}]")

# Surfaces
assembly_top_plane =    mp.ZPlane(z0 =  assembly_length/2, boundary_type = 'vacuum')
assembly_bottom_plane = mp.ZPlane(z0 = -assembly_length/2, boundary_type = 'vacuum')
rod_plane =             mp.ZPlane(z0 =  assembly_length/2 - rod_insertion)
deck += [
    assembly_top_plane,
    assembly_bottom_plane,
    rod_plane
]

lattice_xmin = mp.XPlane(x0=-lattice_elements_x*0.5*pitch, boundary_type='reflective')
lattice_xmax = mp.XPlane(x0= lattice_elements_x*0.5*pitch, boundary_type='reflective')
lattice_ymin = mp.YPlane(y0=-lattice_elements_y*0.5*pitch, boundary_type='reflective')
lattice_ymax = mp.YPlane(y0= lattice_elements_y*0.5*pitch, boundary_type='reflective')
deck += [
    lattice_xmin,
    lattice_xmax,
    lattice_ymin,
    lattice_ymax
]

pin_xmin = mp.XPlane(x0=-0.5*pitch)
pin_xmax = mp.XPlane(x0= 0.5*pitch)
pin_ymin = mp.YPlane(y0=-0.5*pitch)
pin_ymax = mp.YPlane(y0= 0.5*pitch)
deck += [
    pin_xmin,
    pin_xmax,
    pin_ymin,
    pin_ymax
]

fuel_cyl =                   mp.ZCylinder(x0 = 0, y0 = 0, r = fuel_rad)
fuel_helium_gap_cyl =        mp.ZCylinder(x0 = 0, y0 = 0, r = fuel_helium_gap_rad)
fuel_clad_cyl =              mp.ZCylinder(x0 = 0, y0 = 0, r = fuel_clad_rad)
control_rod_cyl =            mp.ZCylinder(x0 = 0, y0 = 0, r = control_rod_rad)
control_rod_helium_gap_cyl = mp.ZCylinder(x0 = 0, y0 = 0, r = control_rod_helium_rad)
control_rod_clad_cyl =       mp.ZCylinder(x0 = 0, y0 = 0, r = control_rod_clad_rad)
guide_tube_inner_cyl =       mp.ZCylinder(x0 = 0, y0 = 0, r = guide_tube_inner_radius)
guide_tube_outer_cyl =       mp.ZCylinder(x0 = 0, y0 = 0, r = guide_tube_outer_radius)
deck += [
    fuel_cyl,
    fuel_helium_gap_cyl,
    fuel_clad_cyl,
    control_rod_cyl,
    control_rod_helium_gap_cyl,
    control_rod_clad_cyl,
    guide_tube_inner_cyl,
    guide_tube_outer_cyl
]

# Regions and Cells
pin_region = +assembly_bottom_plane & -assembly_top_plane  & -pin_xmax & +pin_xmin & -pin_ymax & +pin_ymin
assembly_region = +assembly_bottom_plane & -assembly_top_plane & +lattice_xmin & -lattice_xmax & +lattice_ymin & -lattice_ymax
fuel_region = -fuel_cyl & +assembly_bottom_plane & -assembly_top_plane
burnable_poison_region = -fuel_cyl & +assembly_bottom_plane & -assembly_top_plane
fuel_helium_gap_region = -fuel_helium_gap_cyl & +fuel_cyl & +assembly_bottom_plane & -assembly_top_plane
fuel_clad_region = -fuel_clad_cyl & +fuel_helium_gap_cyl & +assembly_bottom_plane & -assembly_top_plane
pin_moderator_region = +fuel_clad_cyl & +assembly_bottom_plane & -assembly_top_plane  & -pin_xmax & +pin_xmin & -pin_ymax & +pin_ymin
guide_tube_moderator_fill_region = -guide_tube_inner_cyl & +assembly_bottom_plane & -rod_plane
control_rod_region = -control_rod_cyl & +rod_plane & -assembly_top_plane
control_rod_clad_region = -control_rod_clad_cyl & +control_rod_helium_gap_cyl & +rod_plane & -assembly_top_plane
control_rod_moderator_fill_region = -guide_tube_inner_cyl & +control_rod_clad_cyl & +rod_plane & -assembly_top_plane
guide_tube_region = -guide_tube_outer_cyl & +guide_tube_inner_cyl & +assembly_bottom_plane & -assembly_top_plane
control_rod_helium_gap_region = -control_rod_helium_gap_cyl & +control_rod_cyl & +rod_plane & -assembly_top_plane
guide_tube_moderator_region = +guide_tube_outer_cyl & +assembly_bottom_plane & -assembly_top_plane  & -pin_xmax & +pin_xmin & -pin_ymax & +pin_ymin

fuel_cell =                          mp.Cell(name= 1, region=fuel_region,                       material=fuel                , density=-fuel_density                )
burnable_poison_cell =               mp.Cell(name= 2, region=burnable_poison_region,            material=burnable_poison     , density=-burnable_poison_density     )
fuel_helium_gap_cell =               mp.Cell(name= 3, region=fuel_helium_gap_region,            material=helium              , density=-helium_density              )
burnable_poison_helium_gap_cell =    mp.Cell(name= 4, region=fuel_helium_gap_region,            material=helium              , density=-helium_density              )
fuel_clad_cell =                     mp.Cell(name= 5, region=fuel_clad_region,                  material=fuel_cladding       , density=-fuel_cladding_density       )
burnable_poison_clad_cell =          mp.Cell(name= 6, region=fuel_clad_region,                  material=fuel_cladding       , density=-fuel_cladding_density       )
fuel_pin_moderator_cell =            mp.Cell(name= 7, region=pin_moderator_region,              material=moderator           , density=-moderator_density           )
burnable_poison_pin_moderator_cell = mp.Cell(name= 8, region=pin_moderator_region,              material=moderator           , density=-moderator_density           )
guide_tube_moderator_fill_cell =     mp.Cell(name= 9, region=guide_tube_moderator_fill_region,  material=moderator           , density=-moderator_density           )
control_rod_cell =                   mp.Cell(name=10, region=control_rod_region,                material=control_rod_absorber, density=-control_rod_absorber_density)
control_rod_helium_gap_cell =        mp.Cell(name=11, region=control_rod_helium_gap_region,     material=helium              , density=-helium_density              )
control_rod_clad_cell =              mp.Cell(name=12, region=control_rod_clad_region,           material=control_rod_cladding, density=-control_rod_cladding_density)
control_rod_moderator_fill_cell =    mp.Cell(name=13, region=control_rod_moderator_fill_region, material=moderator           , density=-moderator_density           )
guide_tube_cell =                    mp.Cell(name=14, region=guide_tube_region,                 material=fuel_cladding       , density=-fuel_cladding_density       )
guide_tube_moderator_cell =          mp.Cell(name=15, region=guide_tube_moderator_region,       material=moderator           , density=-moderator_density           )

deck += [
    fuel_cell,
    burnable_poison_cell,
    fuel_helium_gap_cell,
    burnable_poison_helium_gap_cell,
    fuel_clad_cell,
    burnable_poison_clad_cell,
    fuel_pin_moderator_cell,
    burnable_poison_pin_moderator_cell,
    guide_tube_moderator_fill_cell,
    control_rod_cell,
    control_rod_helium_gap_cell,
    control_rod_clad_cell,
    control_rod_moderator_fill_cell,
    guide_tube_cell,
    guide_tube_moderator_cell
]

# Universes
fpu = mp.UniverseList(name=1, cells=[fuel_cell, fuel_helium_gap_cell, fuel_clad_cell, fuel_pin_moderator_cell])
bpu = mp.UniverseList(name=2, cells=[burnable_poison_cell, burnable_poison_helium_gap_cell, burnable_poison_clad_cell, burnable_poison_pin_moderator_cell])
gtu = mp.UniverseList(name=3, cells=[guide_tube_moderator_fill_cell, guide_tube_cell, guide_tube_moderator_cell, control_rod_moderator_fill_cell, control_rod_clad_cell, control_rod_helium_gap_cell, control_rod_cell])

# Lattice
numpy_lattice = np.array([[
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, bpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, gtu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          [fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name, fpu.name],
                          ]])

pin_lattice = mp.Lattice(i=[-8,8], j=[-8,8], k=[0,0], lattice=numpy_lattice, type='REC', 
                         universes={
                             fpu.name:fpu,
                             bpu.name:bpu,
                             gtu.name:gtu   
                         })
pin_cell = mp.Cell(region=assembly_region, fill=pin_lattice)
deck += pin_cell
# #NOTE: This may be needed for some reason.  example.py does this.
# pin_lattice_universe = mp.UniverseList(name=4, cells=pin_cell)
# assembly_cell = mp.Cell(region=assembly_region, fill=pin_lattice_universe)
# deck += assembly_cell

for cell in deck.cells.values():
    cell.importances = {'n' : 1.0}
    cell.temperature = reactor_temperature

deck += mp.Cell(name=99, region=~assembly_region, importances={'n':0})

# # Settings
# #=====================

deck += mp.CriticalitySource(histories=1e5, keff_guess=1.0,
                             skip_cycles=20, cycles=70)
deck += mp.CriticalitySourcePoints(
    [
        ( pitch/2, -pitch/2, 0),
        (-pitch/2,  pitch/2, 0),
        ( pitch/2,  pitch/2, 0),
        (-pitch/2, -pitch/2, 0),
    ]
)

# # Flux tally
# flux_tally_mesh = openmc.RegularMesh()
# flux_tally_mesh.lower_left = [-0.5*lattice_elements_x*pitch, -0.5*lattice_elements_y*pitch, -assembly_length/2]
# flux_tally_mesh.upper_right = [0.5*lattice_elements_x*pitch,  0.5*lattice_elements_y*pitch,  assembly_length/2]
# flux_tally_mesh.dimension = [100, 100, 1]
# flux_mesh_filter = openmc.MeshFilter(flux_tally_mesh)

# flux_energy_filter = openmc.EnergyFilter(values=[1.E-5, 0.625, 20.E6])

# flux_tally = openmc.Tally(name='flux')
# flux_tally.scores = ['flux']
# flux_tally.filters = [flux_energy_filter, flux_mesh_filter]

# # Heating tally
# heating_tally_mesh = openmc.RegularMesh()
# heating_tally_mesh.lower_left = [-0.5*lattice_elements_x*pitch, -0.5*lattice_elements_y*pitch, -assembly_length/2]
# heating_tally_mesh.upper_right = [0.5*lattice_elements_x*pitch,  0.5*lattice_elements_y*pitch,  assembly_length/2]
# heating_tally_mesh.dimension = [17, 17, 1]
# heating_mesh_filter = openmc.MeshFilter(heating_tally_mesh)

# heating_tally = openmc.Tally(name = 'heating')
# heating_tally.scores = ['heating-local']
# heating_tally.filters = [heating_mesh_filter]

# # Add tallies
# model.tallies = openmc.Tallies([flux_tally, heating_tally])

# ## Source
# source = openmc.IndependentSource()
# source.space = openmc.stats.Point((0,0,0))
# source.angle = openmc.stats.Isotropic()
# source.energy = openmc.stats.Watt()
# model.settings.source = [source]

# # Export/run
# #=====================
deck.write('model.mcnp')