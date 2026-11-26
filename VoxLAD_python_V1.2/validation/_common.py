import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import voxlad                                                     # noqa: E402

EXAMPLE = os.path.join(HERE, "..", "examples", "shea_tree", "shea_tree.toml")


def example_scans():
    """The two example scans as Scan objects, and the settings of the example."""
    settings, specs, _ = voxlad.load_config(EXAMPLE)
    return voxlad.load_scans(specs), settings
