import openmc
from beavrs.builder import BEAVRS

batches = 800
inactive = 300
particles = 5_000_000

entropy_mesh_x_elements = 15
entropy_mesh_y_elements = 15
entropy_mesh_z_elements = 10

simulation_temp = 560

def get_beavrs_model():
    model = BEAVRS(temperature=simulation_temp)
    return model

def set_settings(model):
    set_keigen_settings(model)
    set_entropy_mesh(model)
    return

def set_keigen_settings(model):
    model.settings.batches = batches
    model.settings.inactive = inactive
    model.settings.particles = particles
    return

def set_entropy_mesh(model):
    mesh = openmc.RegularMesh()
    mesh.lower_left, mesh.upper_right = [-161.2773, -161.2773, 36.748], [161.2773, 161.2773, 402.508]
    mesh.dimension = (entropy_mesh_x_elements, entropy_mesh_y_elements, entropy_mesh_z_elements)
    model.settings.entropy_mesh = mesh
    return