import openmc
import openmc.stats

import numpy as np

# Materials
#=====================

reactor_temperature = 600

fuel = openmc.Material(name='Fuel')
fuel.add_nuclide(nuclide = 'U235', percent = 0.02644492, percent_type = 'wo')
fuel.add_nuclide(nuclide = 'U238', percent = 0.85505247, percent_type = 'wo')
fuel.add_nuclide(nuclide = 'O16' , percent = 0.11850261, percent_type = 'wo')
fuel.set_density('g/cc', 10.3070)
fuel.temperature = reactor_temperature

burnable_poison = openmc.Material(name='Burnable Poison')
burnable_poison.add_nuclide(nuclide = 'U235' , percent = 0.00202753, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'U238' , percent = 0.80898345, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'O16'  , percent = 0.11957908, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd152', percent = 0.00013411, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd154', percent = 0.00148108, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd155', percent = 0.01012049, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd156', percent = 0.01408805, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd157', percent = 0.01083999, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd158', percent = 0.01731511, percent_type = 'wo')
burnable_poison.add_nuclide(nuclide = 'Gd160', percent = 0.01543111, percent_type = 'wo')
burnable_poison.set_density('g/cc', 10.3070)
burnable_poison.temperature = reactor_temperature

moderator = openmc.Material(name='Borated Water Moderator')
moderator.add_nuclide(nuclide = 'H1' , percent = 6.663259e-1, percent_type = 'ao')
moderator.add_nuclide(nuclide = 'O16', percent = 3.330861e-1, percent_type = 'ao')
moderator.add_nuclide(nuclide = 'B10', percent = 7.186970e-5, percent_type = 'ao')
moderator.add_nuclide(nuclide = 'B11', percent = 2.892846e-4, percent_type = 'ao')
moderator.set_density('g/cc', 0.70602)
moderator.add_s_alpha_beta('c_H_in_H2O')
moderator.temperature = reactor_temperature

helium = openmc.Material(name='Helium')
helium.add_element(element = 'He', percent = 1, percent_type = 'ao')
helium.set_density('g/cc', 0.0015981)
helium.temperature = reactor_temperature

fuel_cladding = openmc.Material(name='Zircaloy-4 Cladding')
fuel_cladding.add_element(element = 'O' , percent = 0.00125, percent_type = 'wo')
fuel_cladding.add_element(element = 'Cr', percent = 0.00100, percent_type = 'wo')
fuel_cladding.add_element(element = 'Fe', percent = 0.00210, percent_type = 'wo')
fuel_cladding.add_element(element = 'Zr', percent = 0.98115, percent_type = 'wo')
fuel_cladding.add_element(element = 'Sn', percent = 0.01450, percent_type = 'wo')
fuel_cladding.set_density('g/cc', 6.55)
fuel_cladding.temperature = reactor_temperature

control_rod_absorber = openmc.Material(name='Control Rod Absorber')
control_rod_absorber.add_element(element = 'Ag', percent = 0.80, percent_type = 'wo')
control_rod_absorber.add_element(element = 'In', percent = 0.15, percent_type = 'wo')
control_rod_absorber.add_element(element = 'Cd', percent = 0.05, percent_type = 'wo')
control_rod_absorber.set_density('g/cc', 10.16)
control_rod_absorber.temperature = reactor_temperature

control_rod_cladding = openmc.Material(name='Control Rod Cladding')
control_rod_cladding.add_element(element = 'Si', percent = 0.0060, percent_type = 'wo')
control_rod_cladding.add_element(element = 'Cr', percent = 0.1900, percent_type = 'wo')
control_rod_cladding.add_element(element = 'Mn', percent = 0.0200, percent_type = 'wo')
control_rod_cladding.add_element(element = 'Fe', percent = 0.6840, percent_type = 'wo')
control_rod_cladding.add_element(element = 'Ni', percent = 0.1000, percent_type = 'wo')
control_rod_cladding.set_density('g/cc', 8.03)
control_rod_cladding.temperature = reactor_temperature

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

rod_insertion = 375

if rod_insertion < 0 or rod_insertion > assembly_length:
    raise Exception(f"Rod insertion must be [0, {assembly_length}]")

model = openmc.Model()

# Surfaces
assembly_top_plane =    openmc.ZPlane(z0 =  assembly_length/2, boundary_type = 'vacuum')
assembly_bottom_plane = openmc.ZPlane(z0 = -assembly_length/2, boundary_type = 'vacuum')
rod_plane = openmc.ZPlane(z0 = assembly_length/2 - rod_insertion)

xmin = openmc.XPlane(x0=-lattice_elements_x*0.5*pitch, boundary_type='reflective')
xmax = openmc.XPlane(x0= lattice_elements_x*0.5*pitch, boundary_type='reflective')
ymin = openmc.YPlane(y0=-lattice_elements_y*0.5*pitch, boundary_type='reflective')
ymax = openmc.YPlane(y0= lattice_elements_y*0.5*pitch, boundary_type='reflective')

fuel_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = fuel_rad)
fuel_helium_gap_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = fuel_helium_gap_rad)
fuel_clad_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = fuel_clad_rad)
control_rod_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = control_rod_rad)
control_rod_helium_gap_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = control_rod_helium_rad)
control_rod_clad_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = control_rod_clad_rad)
guide_tube_inner_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = guide_tube_inner_radius)
guide_tube_outer_cyl = openmc.ZCylinder(x0 = 0, y0 = 0, r = guide_tube_outer_radius)

# Regions and Cells
fuel_region = -fuel_cyl & +assembly_bottom_plane & -assembly_top_plane
fuel_cell = openmc.Cell(region=fuel_region, fill=fuel)

burnable_poison_region = -fuel_cyl & +assembly_bottom_plane & -assembly_top_plane
burnable_poison_cell = openmc.Cell(region=burnable_poison_region, fill=burnable_poison)

fuel_helium_gap_region = -fuel_helium_gap_cyl & +fuel_cyl & +assembly_bottom_plane & -assembly_top_plane
fuel_helium_gap_cell = openmc.Cell(region=fuel_helium_gap_region, fill = helium)
burnable_poison_helium_gap_cell = openmc.Cell(region=fuel_helium_gap_region, fill = helium)

fuel_clad_region = -fuel_clad_cyl & +fuel_helium_gap_cyl & +assembly_bottom_plane & -assembly_top_plane
fuel_clad_cell = openmc.Cell(region=fuel_clad_region, fill = fuel_cladding)
burnable_poison_clad_cell = openmc.Cell(region=fuel_clad_region, fill = fuel_cladding)

pin_moderator_region = +fuel_clad_cyl
fuel_pin_moderator_cell = openmc.Cell(region=pin_moderator_region, fill=moderator)
burnable_poison_pin_moderator_cell = openmc.Cell(region=pin_moderator_region, fill=moderator)

guide_tube_moderator_fill_region = -guide_tube_inner_cyl & +assembly_bottom_plane & -rod_plane
guide_tube_moderator_fill_cell = openmc.Cell(region=guide_tube_moderator_fill_region, fill = moderator)

control_rod_region = -control_rod_cyl & +rod_plane & -assembly_top_plane
control_rod_cell = openmc.Cell(region=control_rod_region, fill=control_rod_absorber)

control_rod_helium_gap_region = -control_rod_helium_gap_cyl & +control_rod_cyl & +rod_plane & -assembly_top_plane
control_rod_helium_gap_cell = openmc.Cell(region=control_rod_helium_gap_region, fill=helium)

control_rod_clad_region = -control_rod_clad_cyl & +control_rod_helium_gap_cyl & +rod_plane & -assembly_top_plane
control_rod_clad_cell = openmc.Cell(region=control_rod_clad_region, fill=control_rod_cladding)

control_rod_moderator_fill_region = -guide_tube_inner_cyl & +control_rod_clad_cyl & +rod_plane & -assembly_top_plane
control_rod_moderator_fill_cell = openmc.Cell(region=control_rod_moderator_fill_region, fill=moderator)

guide_tube_region = -guide_tube_outer_cyl & +guide_tube_inner_cyl & +assembly_bottom_plane & -assembly_top_plane
guide_tube_cell = openmc.Cell(region = guide_tube_region, fill=fuel_cladding)

guide_tube_moderator_region = +guide_tube_outer_cyl
guide_tube_moderator_cell = openmc.Cell(region=guide_tube_moderator_region, fill=moderator)

# Universes
fpu = openmc.Universe(name='Fuel Pin Universe', cells=[fuel_cell, fuel_helium_gap_cell, fuel_clad_cell, fuel_pin_moderator_cell])
bpu = openmc.Universe(name='Burnable Poison Pin Universe', cells=[burnable_poison_cell, burnable_poison_helium_gap_cell, burnable_poison_clad_cell, burnable_poison_pin_moderator_cell])
gtu = openmc.Universe(name='Guide Tube Universe', cells=[guide_tube_moderator_fill_cell, guide_tube_cell, guide_tube_moderator_cell, control_rod_moderator_fill_cell, control_rod_clad_cell, control_rod_helium_gap_cell, control_rod_cell])

# Lattice
pin_lattice = openmc.RectLattice()
pin_lattice.pitch = (pitch, pitch)
pin_lattice.lower_left = (-lattice_elements_x*0.5*pitch, -lattice_elements_y*0.5*pitch)
pin_lattice.universes = [
                        [fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, gtu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, gtu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, bpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, gtu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, gtu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, gtu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        [fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu, fpu],
                        ]

assembly_region = +xmin & -xmax & +ymin & -ymax & +assembly_bottom_plane & -assembly_top_plane
assembly_cell = openmc.Cell(region=assembly_region, fill=pin_lattice)

model.geometry = openmc.Geometry(root=[assembly_cell])

# Settings
#=====================

# K-eigenvalue
model.settings = openmc.Settings()
model.settings.run_mode = 'eigenvalue'
model.settings.particles = int(1e6)
model.settings.batches = 650
model.settings.inactive = 50

model.settings.temperature = {'method': 'interpolation'}

# Entropy mesh
model.settings.entropy_mesh = openmc.RegularMesh()
model.settings.entropy_mesh.lower_left = [-lattice_elements_x*0.5*pitch, -lattice_elements_y*0.5*pitch, -assembly_length/2]
model.settings.entropy_mesh.upper_right = [lattice_elements_x*0.5*pitch,  lattice_elements_y*0.5*pitch,  assembly_length/2]
model.settings.entropy_mesh.dimension = [17, 17, 17]

# Flux tally
flux_tally_mesh = openmc.RegularMesh()
flux_tally_mesh.lower_left = [-0.5*lattice_elements_x*pitch, -0.5*lattice_elements_y*pitch, -assembly_length/2]
flux_tally_mesh.upper_right = [0.5*lattice_elements_x*pitch,  0.5*lattice_elements_y*pitch,  assembly_length/2]
flux_tally_mesh.dimension = [100, 100, 1]
flux_mesh_filter = openmc.MeshFilter(flux_tally_mesh)

thermal_fast_energy_filter = openmc.EnergyFilter(values=[1.E-5, 0.625, 20.E6])
flux_spectrum_energy_filter = openmc.EnergyFilter(values=np.logspace(start=np.log10(1e-5), 
                                                       stop=np.log10(20e6), 
                                                       num=129))

flux_tally = openmc.Tally(name='Flux vs Position')
flux_tally.scores = ['flux']
flux_tally.filters = [thermal_fast_energy_filter, flux_mesh_filter]

flux_spectrum_tally = openmc.Tally(name='Flux Spectrum')
flux_spectrum_tally.filters = [flux_spectrum_energy_filter]
flux_spectrum_tally.scores = ['flux']

# Heating tally
heating_tally_mesh = openmc.RegularMesh()
heating_tally_mesh.lower_left = [-0.5*lattice_elements_x*pitch, -0.5*lattice_elements_y*pitch, -assembly_length/2]
heating_tally_mesh.upper_right = [0.5*lattice_elements_x*pitch,  0.5*lattice_elements_y*pitch,  assembly_length/2]
heating_tally_mesh.dimension = [17, 17, 1]
heating_mesh_filter = openmc.MeshFilter(heating_tally_mesh)

heating_tally = openmc.Tally(name = 'Heating')
heating_tally.scores = ['heating']
heating_tally.filters = [heating_mesh_filter]

# Add tallies
model.tallies = openmc.Tallies([flux_tally, flux_spectrum_tally, heating_tally])

## Source
source = openmc.IndependentSource()
source.space = openmc.stats.Box(
    (-lattice_elements_x*0.5*pitch, -lattice_elements_y*0.5*pitch, -assembly_length/2),
    ( lattice_elements_x*0.5*pitch,  lattice_elements_y*0.5*pitch,  assembly_length/2 - rod_insertion)
)
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Watt()
model.settings.source = [source]

# Export/run
#=====================
model.export_to_xml()
model.export_to_model_xml()