import openmc
import openmc.stats
import numpy as np

from dataclasses import dataclass
from typing import Protocol

import materials

DEFAULT_RUN_SETTINGS = openmc.Settings()
DEFAULT_RUN_SETTINGS.run_mode = 'eigenvalue'
DEFAULT_RUN_SETTINGS.particles = int(1e4)
DEFAULT_RUN_SETTINGS.batches = 105
DEFAULT_RUN_SETTINGS.inactive = 5

class ProblemGenerator(Protocol):
    def __call__(
        self,
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05,
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
        use_otf: bool = False,
    ) -> openmc.Model:
        ...

class CritSearchGenerator(Protocol):
    def __call__(
        self,
        U235_ao: float,
        tsl_material: materials.TSLMaterial = materials.LWTR,
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
    ) -> openmc.Model:
        ...

@dataclass
class GeneralizedProblem:
    name: str
    generator: ProblemGenerator
    crit_search_wrapper: CritSearchGenerator

    def __call__(
        self,
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05,
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
    ) -> openmc.Model:
        return self.generator(tsl_material, U235_ao, temperature, run_settings)

def generate_infinite_medium_model(
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05, 
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
        use_otf: bool = False,
        ) -> openmc.Model:

    model = openmc.Model()
    density = 10.0

    material = openmc.Material(name="Mixture")
    material.temperature = temperature
    material.set_density(units='g/cc', density=density)
    material.add_nuclide('B10' , percent=   0.01, percent_type='ao')
    material.add_nuclide('U235', percent=U235_ao, percent_type='ao')
    for isotope in tsl_material.composition:
        material.add_nuclide(nuclide=isotope.identifier, 
                             percent=isotope.abundance, 
                             percent_type='ao')
    material.add_s_alpha_beta(
        tsl_material.tsl_identifier if not use_otf else tsl_material.otf_identifier
    )

    dist = 5
    xl = openmc.XPlane(-dist, boundary_type='reflective')
    xu = openmc.XPlane( dist, boundary_type='reflective')
    yl = openmc.YPlane(-dist, boundary_type='reflective')
    yu = openmc.YPlane( dist, boundary_type='reflective')
    zl = openmc.ZPlane(-dist, boundary_type='reflective')
    zu = openmc.ZPlane( dist, boundary_type='reflective')
    box_region = +xl & -xu & +yl & -yu & +zl & -zu
    box_cell = openmc.Cell(region=box_region, fill = material)

    universe = openmc.Universe(cells=[box_cell])
    model.geometry = openmc.Geometry(root=universe)

    model.settings = run_settings

    source = openmc.IndependentSource()
    source.space = openmc.stats.Box([-dist, -dist, -dist], [dist, dist, dist])
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Watt()
    model.settings.source = [source]

    energy_filter = openmc.EnergyFilter(values=np.logspace(start=np.log10(1e-5), 
                                                       stop=np.log10(20e6), 
                                                       num=129))
    tally = openmc.Tally(name='Flux Spectrum')
    tally.filters = [energy_filter]
    tally.scores = ['flux']
    model.tallies = openmc.Tallies([tally])

    return model

def infinite_medium_wrapper(
        U235_ao: float, 
        tsl_material: materials.TSLMaterial = materials.LWTR,
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
        ) -> openmc.Model:
    return generate_infinite_medium_model(
        tsl_material=tsl_material,
        U235_ao=U235_ao,
        temperature=temperature,
        run_settings=run_settings
        )

INFINITE_MEDIUM = GeneralizedProblem(
    name="InfiniteMedium",
    generator=generate_infinite_medium_model,
    crit_search_wrapper=infinite_medium_wrapper,
)

def generate_slab_model(
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05, 
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
        use_otf: bool = False,
        ) -> openmc.Model:

    model = openmc.Model()

    fuel = openmc.Material(name="Fuel")
    fuel.temperature = temperature
    fuel.set_density(units='g/cc', density=10.97)
    fuel.add_nuclide('O16' , percent=2        , percent_type='ao')
    fuel.add_nuclide('U238', percent=1-U235_ao, percent_type='ao')
    fuel.add_nuclide('U235', percent=U235_ao  , percent_type='ao')


    moderator = openmc.Material(name="Moderator")
    moderator.temperature = temperature
    moderator.set_density(units='g/cc', density=tsl_material.density)
    for isotope in tsl_material.composition:
        moderator.add_nuclide(nuclide=isotope.identifier, 
                             percent=isotope.abundance, 
                             percent_type='ao')
    moderator.add_s_alpha_beta(
        tsl_material.tsl_identifier if not use_otf else tsl_material.otf_identifier
        )

    slab_thickness = 7.5
    xll = openmc.XPlane(-1.5*slab_thickness, boundary_type='vacuum')
    xlu = openmc.XPlane(-0.5*slab_thickness)
    xul = openmc.XPlane( 0.5*slab_thickness)
    xuu = openmc.XPlane( 1.5*slab_thickness, boundary_type='vacuum')

    dist = 5
    yl = openmc.YPlane(-dist, boundary_type='reflective')
    yu = openmc.YPlane( dist, boundary_type='reflective')
    zl = openmc.ZPlane(-dist, boundary_type='reflective')
    zu = openmc.ZPlane( dist, boundary_type='reflective')

    yz_bound_region = +yl & -yu & +zl & -zu
    mod_l_region = +xll & -xlu & yz_bound_region
    fuel_region  = +xlu & -xul & yz_bound_region
    mod_u_region = +xul & -xuu & yz_bound_region

    mod_l_cell = openmc.Cell(region = mod_l_region, fill = moderator)
    fuel_cell  = openmc.Cell(region = fuel_region,  fill = fuel     )
    mod_u_cell = openmc.Cell(region = mod_u_region, fill = moderator)

    universe = openmc.Universe(cells=[mod_l_cell, fuel_cell, mod_u_cell])
    model.geometry = openmc.Geometry(root=universe)

    model.settings = run_settings

    source = openmc.IndependentSource()
    source.space = openmc.stats.Box([-0.5*slab_thickness, -dist, -dist], [0.5*slab_thickness, dist, dist])
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Watt()
    model.settings.source = [source]

    flux_position_mesh = openmc.RegularMesh()
    flux_position_mesh.lower_left  = [-1.5*slab_thickness, -dist, -dist]
    flux_position_mesh.upper_right = [ 1.5*slab_thickness,  dist,  dist]
    flux_position_mesh.dimension = [150, 50, 50]
    flux_position_mesh_filter = openmc.MeshFilter(flux_position_mesh)

    thermal_fast_energy_filter = openmc.EnergyFilter(values=[1.E-5, 0.625, 20.E6])
    spectrum_energy_filter = openmc.EnergyFilter(values=np.logspace(start=np.log10(1e-5), 
                                                       stop=np.log10(20e6), 
                                                       num=129))

    flux_position_tally = openmc.Tally(name='Flux vs Position')
    flux_position_tally.scores = ['flux']
    flux_position_tally.filters = [thermal_fast_energy_filter, flux_position_mesh_filter]

    spectrum_tally = openmc.Tally(name='Flux Spectrum')
    spectrum_tally.scores = ['flux']
    spectrum_tally.filters = [spectrum_energy_filter]

    model.tallies = openmc.Tallies([flux_position_tally, spectrum_tally])

    return model
    
def slab_wrapper(
        U235_ao: float, 
        tsl_material: materials.TSLMaterial = materials.LWTR,
        temperature: float = 293.15,
        run_settings: openmc.Settings = DEFAULT_RUN_SETTINGS,
        ) -> openmc.Model:
    return generate_slab_model(
        tsl_material=tsl_material,
        U235_ao=U235_ao,
        temperature=temperature,
        run_settings=run_settings
        )

SLAB = GeneralizedProblem(
    name="Slab",
    generator=generate_slab_model,
    crit_search_wrapper=slab_wrapper,
)


PROBLEMS = [
    INFINITE_MEDIUM,
    SLAB,
]

if __name__ == "__main__":
    # model = generate_slab_model(materials.LWTR, U235_ao=0.5)
    # model.export_to_model_xml()
    pass