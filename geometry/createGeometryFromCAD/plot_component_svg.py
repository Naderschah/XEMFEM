#!/usr/bin/env python3
import importlib.util
from pathlib import Path


def _load_root_plotter():
    root_script = Path(__file__).resolve().parents[1] / "plot_component_svg.py"
    spec = importlib.util.spec_from_file_location("geometry_plot_component_svg", root_script)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise RuntimeError(f"failed to load plotter from {root_script}")
    spec.loader.exec_module(module)
    return module


def main():
    module = _load_root_plotter()
    module.main()


if __name__ == "__main__":
    main()
