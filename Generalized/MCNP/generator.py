import numpy as np
import argparse

import mcnpy as mp
from mcnpy.elements import *

from dataclasses import dataclass
from typing import Protocol, Optional
from pathlib import Path

import materials

DEFAULT_RUN_SETTINGS = mp.CriticalitySource(
    histories=1e6,
    keff_guess=1.0,
    skip_cycles=5,
    cycles=105,
)

class ProblemGenerator(Protocol):
    def __call__(
        self,
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05,
        run_settings: mp.CriticalitySource = DEFAULT_RUN_SETTINGS,
        use_otf: bool = False,
    ) -> mp.Deck:
        ...

@dataclass
class GeneralizedProblem:
    name: str
    generator: ProblemGenerator

    def __call__(
        self,
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05,
        run_settings: mp.CriticalitySource = DEFAULT_RUN_SETTINGS,
    ) -> mp.Deck:
        return self.generator(tsl_material, U235_ao, run_settings)

def generate_infinite_medium_model(
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05, 
        run_settings: mp.CriticalitySource = DEFAULT_RUN_SETTINGS,
        ) -> mp.Deck:

    deck = mp.Deck()
    density = 10.0

    U235 = mp.Nuclide(name = 'U235', fraction = U235_ao, unit = materials.UNIT, library = materials.LIBRARY)
    B10  = mp.Nuclide(name = 'B10',  fraction = 0.01,    unit = materials.UNIT, library = materials.LIBRARY)
    material = mp.Material(name = 1, nuclides = [U235, B10])
    for nuclide in tsl_material.nuclides:
        material.nuclides.append(nuclide)
    material.s_alpha_beta = tsl_material.mcnp_placeholder_tsl
    deck += [material]

    dist = 5
    xl = mp.XPlane(x0=-dist, boundary_type='reflective')
    xu = mp.XPlane(x0= dist, boundary_type='reflective')
    yl = mp.YPlane(y0=-dist, boundary_type='reflective')
    yu = mp.YPlane(y0= dist, boundary_type='reflective')
    zl = mp.ZPlane(z0=-dist, boundary_type='reflective')
    zu = mp.ZPlane(z0= dist, boundary_type='reflective')
    deck += [xl, xu, yl, yu , zl, zu]
    box_region = +xl & -xu & +yl & -yu & +zl & -zu
    box_cell = mp.Cell(name=1, region=box_region, material=material, density=density)
    deck += box_cell

    for cell in deck.cells.values():
        cell.importances = {'n' : 1.0}

    deck += mp.Cell(name=99, region=~box_region, importances={'n':0}, comment = 'Outer Void')

    deck += run_settings

    deck += mp.CriticalitySourcePoints([[0,0,0]])
    deck += run_settings

    energies = np.logspace(start=np.log10(1e-5), stop=np.log10(20e6), num=129)
    energies = (float(i) for i in energies)
    origin = mp.Point(-dist, -dist, -dist)
    particle = mp.Particle('N')
    spectrum_tally = mp.Tally.FMESH(name=14, origin=origin, geometry='xyz',
                                    i_nodes=[dist], i_subdivisions=[1], 
                                    j_nodes=[dist], j_subdivisions=[1], 
                                    k_nodes=[dist], k_subdivisions=[1], 
                                    energy_nodes=energies,
                                    particles = [particle])
    deck += spectrum_tally

    return deck

INFINITE_MEDIUM = GeneralizedProblem(
    name="InfiniteMedium",
    generator=generate_infinite_medium_model,
)

def generate_slab_model(
        tsl_material: materials.TSLMaterial = materials.LWTR,
        U235_ao: float = 0.05, 
        run_settings: mp.CriticalitySource = DEFAULT_RUN_SETTINGS,
        ) -> mp.Deck:

    deck = mp.Deck()

    O16 =  mp.Nuclide(name = 'O16' , fraction = 2        , unit = materials.UNIT, library = materials.LIBRARY)
    U238 = mp.Nuclide(name = 'U238', fraction = 1-U235_ao, unit = materials.UNIT, library = materials.LIBRARY)
    U235 = mp.Nuclide(name = 'U235', fraction = U235_ao  , unit = materials.UNIT, library = materials.LIBRARY)
    fuel = mp.Material(name = 1, nuclides = [O16, U238, U235])

    moderator = mp.Material(name = 2, nuclides = tsl_material.nuclides)
    moderator.s_alpha_beta = tsl_material.mcnp_placeholder_tsl

    deck += [fuel, moderator]

    slab_thickness = 7.5
    xll = mp.XPlane(x0=-1.5*slab_thickness, boundary_type='vacuum')
    xlu = mp.XPlane(x0=-0.5*slab_thickness)
    xul = mp.XPlane(x0= 0.5*slab_thickness)
    xuu = mp.XPlane(x0= 1.5*slab_thickness, boundary_type='vacuum')
    deck += [xll, xlu, xul, xuu]

    dist = 5
    yl = mp.YPlane(y0=-dist, boundary_type='reflective')
    yu = mp.YPlane(y0= dist, boundary_type='reflective')
    zl = mp.ZPlane(z0=-dist, boundary_type='reflective')
    zu = mp.ZPlane(z0= dist, boundary_type='reflective')
    deck += [yl, yu, zl, zu]

    mod_l_region = +xll & -xlu & +yl & -yu & +zl & -zu
    fuel_region  = +xlu & -xul & +yl & -yu & +zl & -zu
    mod_u_region = +xul & -xuu & +yl & -yu & +zl & -zu
    outside_region = -xll | +xuu | -yl | +yu | -zl | +zu

    mod_l_cell = mp.Cell(name = 1, region = mod_l_region, material = moderator, density = tsl_material.density)
    fuel_cell  = mp.Cell(name = 2, region = fuel_region,  material = fuel     , density = 10.97)
    mod_u_cell = mp.Cell(name = 3, region = mod_u_region, material = moderator, density = tsl_material.density)

    deck += [mod_l_cell, fuel_cell, mod_u_cell]

    for cell in deck.cells.values():
        cell.importances = {'n' : 1.0}

    deck += mp.Cell(name=99, region=outside_region, importances={'n':0}, comment = 'Outer Void')

    deck += run_settings

    deck += mp.CriticalitySourcePoints([[0,0,0]])
    deck += run_settings

    energies = np.logspace(start=np.log10(1e-5), stop=np.log10(20e6), num=129)
    energies = (float(i) for i in energies)
    origin = mp.Point(-1.5*slab_thickness, -dist, -dist)
    particle = mp.Particle('N')
    spectrum_tally = mp.Tally.FMESH(name=14, origin=origin, geometry='xyz',
                                    i_nodes=[1.5*slab_thickness], i_subdivisions=[1], 
                                    j_nodes=[dist]              , j_subdivisions=[1], 
                                    k_nodes=[dist]              , k_subdivisions=[1], 
                                    energy_nodes=energies,
                                    particles = [particle])
    deck += spectrum_tally

    flux_position_tally = mp.Tally.FMESH(name=24, origin=origin, geometry='xyz',
                                    i_nodes=[1.5*slab_thickness], i_subdivisions=[150], 
                                    j_nodes=[dist]              , j_subdivisions=[50], 
                                    k_nodes=[dist]              , k_subdivisions=[50], 
                                    energy_nodes=(0.625, 20.E6),
                                    particles = [particle])
    deck += flux_position_tally

    return deck
    
SLAB = GeneralizedProblem(
    name="Slab",
    generator=generate_slab_model,
)

PROBLEMS = [
    INFINITE_MEDIUM,
    SLAB,
]

PROBLEM_NAMES = [problem.name for problem in PROBLEMS]

def get_problem(name:str = 'InfiniteMedium') -> GeneralizedProblem:
    for problem in PROBLEMS:
        if name == problem.name:
            return problem
    raise KeyError(f"{name}: Unknown problem name. Allowable names: {PROBLEM_NAMES}")

def get_material(name: str = 'LWTR') -> materials.TSLMaterial:
    for material in materials.MATERIALS:
        if name == material.name:
            return material
    raise KeyError(f"{name}: Unknown material name. Allowable names: {materials.MATERIAL_NAMES}")

def fix_sab_name(loc: Path, old:str, new:str):
    text = loc.read_text()
    text = text.replace(old, new)
    loc.write_text(text)
    return

def write_slurm_script(path: Path, otf_name: Optional[str] = None):
    lines = [
        "#!/bin/bash",
        "#",
        "# ===== SLURM Job Settings =====",
        "#SBATCH --job-name=mcnp             # Job name",
        "#SBATCH --partition=debug           # Partition",
        "#SBATCH --ntasks=1                  # Single process",
        "#SBATCH --cpus-per-task=11          # Number of OpenMP threads",
        "#SBATCH --time=0                    # Max runtime",
        "#SBATCH --output=slurm-%j.out       # Stdout",
        "#SBATCH --error=slurm-%j.err        # Stderr",
        "",
        "# ===== Environment Setup =====",
        "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK",
        "",
        "# Export MCNP input files or other variables",
        'export DATAPATH="/home/camden/Data/MCNP/MCNP_DATA"',
        "",
        'EXE="/home/camden/MyPrograms/MCNP-thermal-otf/MCNP_6.3_SOURCE/mcnp-src/mcnp-6.3.0-Source/mcnp6/build/mcnp6"',
        'INP="model.mcnp"',
    ]
    if otf_name is not None:
        lines.append(f'OTF="/home/camden/Data/OTF/Test_ENDF8_0_OTF_Files/{otf_name}.h5"')
    lines.extend([
        "",
        "# ===== Run MCNP =====",
    ])
    if otf_name is not None:
        lines.append('$EXE inp=$INP otf $OTF tasks $OMP_NUM_THREADS')
    else:
        lines.append('$EXE inp=$INP tasks $OMP_NUM_THREADS')
    path.write_text("\n".join(lines) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Creates a generalized MCNP input file along with its slurm file.")
    parser.add_argument("--problem_name",
                        type=str,
                        default='InfiniteMedium',
                        help="Sets the name of the problem to run.")
    parser.add_argument("--material_name",
                        type=str,
                        default='LWTR',
                        help="Sets the name of the material to use.")
    parser.add_argument("--u235_ao",
                        type=float,
                        default=0.05,
                        help="Sets the concentration of U235.")
    parser.add_argument("--use-otf",
                        action="store_true",
                        help="Sets the command line options to run MCNP with the OTF file.")
    parser.add_argument("--save_loc",
                        type=Path,
                        default=Path('.'),
                        help="Sets the location to write the input and slurm file.")
    parser.add_argument("--input_name",
                        type=Path,
                        default=Path('model.mcnp'),
                        help="Sets the name of the input file.")
    parser.add_argument("--slurm_name",
                        type=Path,
                        default=Path('mcnp.slurm'),
                        help="Sets the name of the slurm file.")

    args = parser.parse_args()

    problem_name = args.problem_name
    material_name = args.material_name
    u235_ao = args.u235_ao
    use_otf = args.use_otf
    save_loc = args.save_loc
    input_name = args.input_name
    slurm_name = args.slurm_name

    input_file = Path(save_loc, input_name)
    slurm_file = Path(save_loc, slurm_name)

    problem = get_problem(problem_name)
    material = get_material(material_name)

    deck = problem.generator(material, u235_ao)
    deck.write(input_file)
    fix_sab_name(input_file, material.mcnp_placeholder_tsl, material.tsl_identifier)
    if use_otf:
        write_slurm_script(slurm_file,  material.otf_identifier)
    else:
        write_slurm_script(slurm_file)


if __name__ == "__main__":
    main()