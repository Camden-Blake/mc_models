import openmc
from beavrs.builder import BEAVRS
import numpy as np

batches = 800
inactive = 300
particles = 5_000_000

entropy_mesh_x_elements = 15
entropy_mesh_y_elements = 15
entropy_mesh_z_elements = 10

lower_left = [-161.2773, -161.2773, 36.748]
upper_right = [161.2773, 161.2773, 402.508]

simulation_temp = 560

def get_beavrs_model() -> openmc.Model:
    model = BEAVRS(temperature=simulation_temp)
    return model

def set_settings(model: openmc.Model):
    set_keigen_settings(model)
    set_entropy_mesh(model)
    set_tallies(model)
    return

def set_keigen_settings(model: openmc.Model):
    model.settings.batches = batches
    model.settings.inactive = inactive
    model.settings.particles = particles
    return

def set_entropy_mesh(model: openmc.Model):
    mesh = openmc.RegularMesh()
    mesh.lower_left, mesh.upper_right = lower_left, upper_right
    mesh.dimension = (entropy_mesh_x_elements, entropy_mesh_y_elements, entropy_mesh_z_elements)
    model.settings.entropy_mesh = mesh
    return

def set_tallies(model: openmc.Model):
    flux_assembly_mesh = openmc.RegularMesh()
    flux_assembly_mesh.lower_left = lower_left
    flux_assembly_mesh.upper_right = upper_right
    flux_assembly_mesh.dimension = [15, 15, 1]
    flux_assembly_mesh_filter = openmc.MeshFilter(flux_assembly_mesh)

    thermal_fast_energy_filter = openmc.EnergyFilter(values=[1.E-5, 0.625, 20.E6])
    flux_spectrum_energy_filter = openmc.EnergyFilter(values=np.logspace(start=np.log10(1e-5), 
                                                           stop=np.log10(20e6), 
                                                           num=129))

    flux_assembly_tally = openmc.Tally(name='Assembly Flux')
    flux_assembly_tally.scores = ['flux']
    flux_assembly_tally.filters = [thermal_fast_energy_filter, flux_assembly_mesh_filter]

    flux_spectrum_tally = openmc.Tally(name='Flux Spectrum')
    flux_spectrum_tally.scores = ['flux']
    flux_spectrum_tally.filters = [flux_spectrum_energy_filter]

    model.tallies = openmc.Tallies([flux_assembly_tally, flux_spectrum_tally])
    return