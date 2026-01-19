import json

from uncertainties import UFloat
from dataclasses import dataclass
from dataclasses import asdict
from pathlib import Path

import materials
import generator

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
    raise RuntimeError("Not Implemented")


def main():
    for problem in generator.PROBLEMS:
        outfile = OPTIMIZER_SAVE_LOCATION / f"{OPTIMIZER_PREFIX}{problem.name}.json"
        optimize_problem(problem, outfile)


if __name__ == "__main__":
    main()