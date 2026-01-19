import os
import time

import openmc

from pathlib import Path

import materials
import generator
import optimizer

def read_optimized_results() -> dict[str, optimizer.OptimizedResults]:
    files = [file for file in optimizer.OPTIMIZER_SAVE_LOCATION.glob('*') 
             if file.is_file()]
    results: dict[str, optimizer.OptimizedResults] = {}
    for file in files:
        problem = file.stem.removeprefix(optimizer.OPTIMIZER_PREFIX)
        results[problem] = optimizer.load_optimized_results(file)
    return results

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

            os.system(f"python generator.py --problem_name {problem.name} --material_name {material.name} --u235_ao {value} --save_loc {standard_loc}")
            os.system(f"python generator.py --problem_name {problem.name} --material_name {material.name} --u235_ao {value} --save_loc {otf_loc} --use-otf")

if __name__ == "__main__":
    results = read_optimized_results()
    create_inputs(
        results,
        materials.MATERIALS,
        generator.PROBLEMS,
    )