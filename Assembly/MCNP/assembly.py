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
     H[1]@6.663259e-1 +
    O[16]@3.330861e-1 +
    B[10]@7.186970e-5 +
    B[11]@2.892846e-4
)
moderator.s_alpha_beta = ['lwtr']
moderator_density = 0.79602
deck += [moderator]

helium = mp.Material(
    HE[3]@0.000002 +
    HE[4]@0.999998
    )
helium_density = 0.0015981
deck += [helium]

fuel_cladding = mp.Material(
     O[16]*0.0012466835910117907 +
     O[17]*5.033555467821137e-07 +
     O[18]*2.813053441427251e-06 +
    CR[50]*4.173686622718313e-05 +
    CR[52]*0.000836993623513459 +
    CR[53]*9.673586818675183e-05 +
    CR[54]*2.4533642072606143e-05 +
    FE[54]*0.00011855672267684395 +
    FE[56]*0.0019299321043708116 +
    FE[57]*4.5367740887167366e-05 +
    FE[58]*6.143432065177235e-06 +
    ZR[90]*0.4975030719418411 +
    ZR[91]*0.10970127724609834 +
    ZR[92]*0.1695240925103877 +
    ZR[94]*0.1755385702549024 +
    ZR[96]*0.02888298804677048 +
    SN[112]*0.00013258696500899892 +
    SN[114]*9.18244938630758e-05 +
    SN[115]*4.771905899254868e-05 +
    SN[116]*0.0020584231441515095 +
    SN[117]*0.0010966473380910862 +
    SN[118]*0.003487981278597598 +
    SN[119]*0.0012475771051862219 +
    SN[120]*0.004771539508846636 +
    SN[122]*0.000689409486231945 +
    SN[124]*0.0008762916210303822
)
fuel_cladding_density = 6.55
deck += [fuel_cladding]

control_rod_absorber = mp.Material(
    AG[107]*0.41100940717766427 +
    AG[109]*0.3889905928223358 +
    IN[113]*0.006314443215116888 +
    IN[115]*0.14368555678488312 +
    CD[106]*0.0005864650095688713 +
    CD[108]*0.00042618832894341903 +
    CD[110]*0.00609573858367464 +
    CD[111]*0.006311586283065275 +
    CD[112]*0.011999697935677817 +
    CD[113]*0.0061401801606555464 +
    CD[114]*0.0145675033416658 +
    CD[116]*0.003872640356748635
)
control_rod_absorber_density = 10.16
deck += [control_rod_absorber]

control_rod_cladding = mp.Material(
    SI[28]*0.005512411036970206 +
    SI[29]*0.00028990501756718373 +
    SI[30]*0.00019768394546261053 +
    CR[50]*0.007930004583164795 +
    CR[52]*0.1590287884675572 +
    CR[53]*0.01837981495548285 +
    CR[54]*0.004661391993795167 +
    MN[55]*0.02 +
    FE[54]*0.03861561824331489 +
    FE[56]*0.6286064568522073 +
    FE[57]*0.014776921317534514 +
    FE[58]*0.0020010035869434425 +
    NI[58]*0.0671977052965899 +
    NI[60]*0.026775962899229563 +
    NI[61]*0.0011833590852967102 +
    NI[62]*0.0038348223034871697 +
    NI[64]*0.00100815041539667
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
xmin = -lattice_elements_x*0.5*pitch
xmax =  lattice_elements_x*0.5*pitch
ymin = -lattice_elements_y*0.5*pitch
ymax =  lattice_elements_y*0.5*pitch
zmin = -assembly_length/2
zmax =  assembly_length/2

rod_insertion = 375
reactor_temperature = 565.0

if rod_insertion < 0 or rod_insertion > assembly_length:
    raise Exception(f"Rod insertion must be [0, {assembly_length}]")

# Surfaces
assembly_top_plane =    mp.ZPlane(z0 = zmax, comment='Assembly top plane')
assembly_bottom_plane = mp.ZPlane(z0 = zmin, comment='Assembly bottom plane')
rod_plane =             mp.ZPlane(z0 = zmax - rod_insertion, comment = 'CR insertion plane')
deck += [
    assembly_top_plane,
    assembly_bottom_plane,
    rod_plane
]

lattice_xmin = mp.XPlane(x0=xmin, boundary_type='reflective', comment = 'Lattice bound')
lattice_xmax = mp.XPlane(x0=xmax, boundary_type='reflective', comment = 'Lattice bound')
lattice_ymin = mp.YPlane(y0=ymin, boundary_type='reflective', comment = 'Lattice bound')
lattice_ymax = mp.YPlane(y0=ymax, boundary_type='reflective', comment = 'Lattice bound')
deck += [
    lattice_xmin,
    lattice_xmax,
    lattice_ymin,
    lattice_ymax
]

pin_xmin = mp.XPlane(x0=-0.6*pitch, comment = 'Pin bound')
pin_xmax = mp.XPlane(x0= 0.6*pitch, comment = 'Pin bound')
pin_ymin = mp.YPlane(y0=-0.6*pitch, comment = 'Pin bound')
pin_ymax = mp.YPlane(y0= 0.6*pitch, comment = 'Pin bound')
deck += [
    pin_xmin,
    pin_xmax,
    pin_ymin,
    pin_ymax
]

fuel_cyl =                   mp.ZCylinder(x0 = 0, y0 = 0, r = fuel_rad,                comment = 'Fuel')
fuel_helium_gap_cyl =        mp.ZCylinder(x0 = 0, y0 = 0, r = fuel_helium_gap_rad,     comment = 'Fuel Helium')
fuel_clad_cyl =              mp.ZCylinder(x0 = 0, y0 = 0, r = fuel_clad_rad,           comment = 'Fuel Clad')
control_rod_cyl =            mp.ZCylinder(x0 = 0, y0 = 0, r = control_rod_rad,         comment = 'Control Rod')
control_rod_helium_gap_cyl = mp.ZCylinder(x0 = 0, y0 = 0, r = control_rod_helium_rad,  comment = 'Control Rod Helium')
control_rod_clad_cyl =       mp.ZCylinder(x0 = 0, y0 = 0, r = control_rod_clad_rad,    comment = 'Control Rod Clad')
guide_tube_inner_cyl =       mp.ZCylinder(x0 = 0, y0 = 0, r = guide_tube_inner_radius, comment = 'Guide Tube Inner')
guide_tube_outer_cyl =       mp.ZCylinder(x0 = 0, y0 = 0, r = guide_tube_outer_radius, comment = 'Guide Tube Outer')
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
fuel_region = -fuel_cyl
burnable_poison_region = -fuel_cyl
fuel_helium_gap_region = -fuel_helium_gap_cyl & +fuel_cyl
fuel_clad_region = -fuel_clad_cyl & +fuel_helium_gap_cyl
pin_moderator_region = +fuel_clad_cyl
guide_tube_moderator_fill_region = -guide_tube_inner_cyl & -rod_plane
control_rod_region = -control_rod_cyl & +rod_plane
control_rod_clad_region = -control_rod_clad_cyl & +control_rod_helium_gap_cyl & +rod_plane
control_rod_moderator_fill_region = -guide_tube_inner_cyl & +control_rod_clad_cyl & +rod_plane
guide_tube_region = -guide_tube_outer_cyl & +guide_tube_inner_cyl
control_rod_helium_gap_region = -control_rod_helium_gap_cyl & +control_rod_cyl & +rod_plane
guide_tube_moderator_region = +guide_tube_outer_cyl

fuel_cell =                          mp.Cell(name= 1, region=fuel_region,                       material=fuel                , density=-fuel_density                , comment = 'Fuel')
burnable_poison_cell =               mp.Cell(name= 2, region=burnable_poison_region,            material=burnable_poison     , density=-burnable_poison_density     , comment = 'Poison')
fuel_helium_gap_cell =               mp.Cell(name= 3, region=fuel_helium_gap_region,            material=helium              , density=-helium_density              , comment = 'Fuel Helium')
burnable_poison_helium_gap_cell =    mp.Cell(name= 4, region=fuel_helium_gap_region,            material=helium              , density=-helium_density              , comment = 'Poison Helium')
fuel_clad_cell =                     mp.Cell(name= 5, region=fuel_clad_region,                  material=fuel_cladding       , density=-fuel_cladding_density       , comment = 'Fuel Clad')
burnable_poison_clad_cell =          mp.Cell(name= 6, region=fuel_clad_region,                  material=fuel_cladding       , density=-fuel_cladding_density       , comment = 'Poison Clad')
fuel_pin_moderator_cell =            mp.Cell(name= 7, region=pin_moderator_region,              material=moderator           , density=-moderator_density           , comment = 'Fuel Mod')
burnable_poison_pin_moderator_cell = mp.Cell(name= 8, region=pin_moderator_region,              material=moderator           , density=-moderator_density           , comment = 'Poison Mod')
guide_tube_moderator_fill_cell =     mp.Cell(name= 9, region=guide_tube_moderator_fill_region,  material=moderator           , density=-moderator_density           , comment = 'Guide Tube Mod')
control_rod_cell =                   mp.Cell(name=10, region=control_rod_region,                material=control_rod_absorber, density=-control_rod_absorber_density, comment = 'Control Rod')
control_rod_helium_gap_cell =        mp.Cell(name=11, region=control_rod_helium_gap_region,     material=helium              , density=-helium_density              , comment = 'Control Rod Helium')
control_rod_clad_cell =              mp.Cell(name=12, region=control_rod_clad_region,           material=control_rod_cladding, density=-control_rod_cladding_density, comment = 'Control Rod Clad')
control_rod_moderator_fill_cell =    mp.Cell(name=13, region=control_rod_moderator_fill_region, material=moderator           , density=-moderator_density           , comment = 'Control Rod Mod')
guide_tube_cell =                    mp.Cell(name=14, region=guide_tube_region,                 material=fuel_cladding       , density=-fuel_cladding_density       , comment = 'Guide Tube')
guide_tube_moderator_cell =          mp.Cell(name=15, region=guide_tube_moderator_region,       material=moderator           , density=-moderator_density           , comment = 'Guide Tube Mod')

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

# NOTE: The i,j,k indexing makes no sense to me but whatever, it works.
pin_lattice = mp.Lattice(i=[0,0], j=[-8,8], k=[-8,8], lattice=numpy_lattice, type='REC', 
                         universes={
                             fpu.name:fpu,
                             bpu.name:bpu,
                             gtu.name:gtu   
                         })
pin_cell = mp.Cell(region=pin_region, fill=pin_lattice, comment = 'Pin Cell')
deck += pin_cell
pin_lattice_universe = mp.UniverseList(name=4, cells=pin_cell)
assembly_cell = mp.Cell(region=assembly_region, fill=pin_lattice_universe, comment = 'Assembly Cell')
deck += assembly_cell

for cell in deck.cells.values():
    cell.importances = {'n' : 1.0}
    # cell.temperature = reactor_temperature

deck += mp.Cell(name=99, region=~assembly_region, importances={'n':0}, comment = 'Outer Void')

# # Settings
# #=====================

deck += mp.CriticalitySource(histories=1e5, keff_guess=1.0,
                             skip_cycles=20, cycles=70)
# Place Source points in each pin cell
coords = [
    (pitch*(i - lattice_elements_x//2),
     pitch*(j - lattice_elements_y//2),
     assembly_length/2 - rod_insertion)
    for i in range(lattice_elements_x)
    for j in range(lattice_elements_y)
]
deck += mp.CriticalitySourcePoints(coords)
with open(".source_points.txt", "w") as f:
    for x, y, z in coords:
        f.write(f"{x:10.4f} {y:10.4f} {z:10.4f}\n")

flux_tally = mp.Tally.FMESH(name=14, geometry='xyz', origin=(xmin, ymin, zmin),
                            i_nodes=xmax, i_subdivisions=100,
                            j_nodes=ymax, j_subdivisions=100,
                            k_nodes=zmax,
                            energy_nodes=(0, 0.625e-6, 20),
                            format='xdmf'
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
deck.write('mcnpy_model.mcnp')

print("Don't use this script.  It was too frustrating to get it to do what I need.")
print("Major changes need to be made to the input file.")
print("You need to manually change the SAB file to the correct ENDF8.0 version.")
print("You must also add the library number at the end of each ZAID because the temp card causes the simulation to not run.")