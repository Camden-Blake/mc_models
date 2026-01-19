import json

from uncertainties import UFloat
from dataclasses import dataclass
from dataclasses import asdict
from pathlib import Path

import openmc

import materials
import generator

OPTIMIZER_SETTINGS = openmc.Settings()
OPTIMIZER_SETTINGS.run_mode = 'eigenvalue'
OPTIMIZER_SETTINGS.particles = int(1e4)
OPTIMIZER_SETTINGS.batches = 105
OPTIMIZER_SETTINGS.inactive = 5

OPTIMIZER_PREFIX = "optimized_"
OPTIMIZER_SAVE_LOCATION = Path("Optimized")
OPTIMIZER_SAVE_LOCATION.mkdir(parents=True, exist_ok=True)

def serialize_ufloats(vals:list[UFloat]):
    return [(float(k.n), float(k.s)) for k in vals]

@dataclass
class OptimizedMaterial:
    name: str
    value: float
    guesses: list[float]
    keffs: list[tuple[float, float]]
    temperature: float = 293.15

OptimizedResults = dict[str, OptimizedMaterial]

def load_optimized_results(filename: Path) -> OptimizedResults:
    with open(filename, "r") as f:
        raw = json.load(f)
    results: OptimizedResults = {}
    for key, data in raw.items():
        mat = OptimizedMaterial(
            name=data["name"],
            value=data["value"],
            guesses=data["guesses"],
            keffs=data["keffs"],
            temperature=data["temperature"],
        )
        results[key] = mat
    return results

def optimize_problem(problem: generator.GeneralizedProblem, out_filename: Path):
    results = {}
    for mat in materials.MATERIALS:
        print(f"Working on {mat.name} ({problem.name})")

        value, guesses, keffs = openmc.search_for_keff(
            problem.crit_search_wrapper,
            bracket=[0.01, 0.9],
            tol=0.01,
            model_args={
                "tsl_material": mat,
                "run_settings": OPTIMIZER_SETTINGS,
                },
            print_iterations=True,
            run_args={"output": False},
        )

        results[mat.name] = OptimizedMaterial(
            name=mat.name,
            value=value,
            guesses=guesses,
            keffs=serialize_ufloats(keffs),
        )

    with open(out_filename, "w") as f:
        json.dump({k: asdict(v) for k, v in results.items()}, f, indent=2)


def main():
    for problem in generator.PROBLEMS:
        outfile = OPTIMIZER_SAVE_LOCATION / f"{OPTIMIZER_PREFIX}{problem.name}.json"
        optimize_problem(problem, outfile)


if __name__ == "__main__":
    main()