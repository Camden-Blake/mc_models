import openmc

from pathlib import Path

import materials
import generator
import optimizer

run_settings = openmc.Settings()
run_settings.run_mode = 'eigenvalue'
run_settings.particles = int(1e6)
run_settings.batches = 105
run_settings.inactive = 5

def read_optimized_results() -> dict[str, optimizer.OptimizedResults]:
    files = [file for file in optimizer.OPTIMIZER_SAVE_LOCATION.glob('*') 
             if file.is_file()]
    results: dict[str, optimizer.OptimizedResults] = {}
    for file in files:
        problem = file.stem.removeprefix(optimizer.OPTIMIZER_PREFIX)
        results[problem] = optimizer.load_optimized_results(file)
    return results

def create_input(
        problem:generator.GeneralizedProblem,
        material:materials.TSLMaterial,
        value:float,
        location:Path,
        use_otf:bool = False,
        ):
    model = problem.generator(
        tsl_material=material,
        U235_ao=value,
        run_settings=run_settings,
        use_otf=use_otf
        )
    model.export_to_model_xml(location)
    pass

def create_inputs(
        optimized_results:dict[str, optimizer.OptimizedResults],
        materials:list[materials.TSLMaterial],
        problems:list[generator.GeneralizedProblem],
        ):
    for problem in problems:
        prob_loc = Path(problem.name)
        for material in materials:
            value = optimized_results[problem.name][material.name].value
            
            standard_loc = prob_loc / f"{problem.name}_{material.name}" / "Standard"
            standard_loc.mkdir(parents=True, exist_ok=True)

            otf_loc = prob_loc / f"{problem.name}_{material.name}" / "OTF"
            otf_loc.mkdir(parents=True, exist_ok=True)

            create_input(problem, material, value, standard_loc)
            create_input(problem, material, value, otf_loc, use_otf=True)

if __name__ == "__main__":
    results = read_optimized_results()
    create_inputs(
        results,
        materials.MATERIALS,
        generator.PROBLEMS,
    )