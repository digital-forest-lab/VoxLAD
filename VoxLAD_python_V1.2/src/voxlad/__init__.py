"""VoxLAD: voxel-level leaf area density from ground-based (terrestrial) lidar.

Python implementation of the model of Beland, Widlowski & Fournier (2014),
Environmental Modelling & Software 51, 184-189, doi:10.1016/j.envsoft.2013.09.034
"""
__version__ = "1.2.0"

import os                                                          # noqa: E402

from .io import InputError, LEAF, WOOD, NOISE                      # noqa: E402,F401
from .model import (Scan, Settings, Result, compute, FLAG_NAMES, FLAG_OK,   # noqa: E402,F401
                    FLAG_LOW_SAMPLING, FLAG_UNRESOLVED, FLAG_OCCLUDED, FLAG_NOT_SEEN)
from .config import load_config, load_scans, load_scan, TEMPLATE  # noqa: E402,F401
from .report import summary_text, write_results, vertical_profile  # noqa: E402,F401
from .viewer import view, write_viewer, write_viewer_from_csv      # noqa: E402,F401
from . import orientation                                          # noqa: E402,F401
from .orientation import estimate_axis                             # noqa: E402,F401


def run(config_path, progress=None, write=True, view=False, **overrides):
    """Run a configuration file end to end.  Returns the Result.

    Keyword overrides replace settings from the file, e.g. run("a.toml",
    voxel_size=0.2, compat="matlab").  ``view=True`` opens the 3-D viewer of
    the result in the web browser when the run is done.
    """
    from . import viewer as _viewer
    settings, specs, out = load_config(config_path)
    for k, v in overrides.items():
        if not hasattr(settings, k):
            raise InputError(f"unknown setting '{k}'")
        setattr(settings, k, v)
    scans = load_scans(specs, progress)
    res = compute(scans, settings, progress)
    if write:
        res.files = write_results(res, out["prefix"], out["per_scan"], out["npz"], out["matlab_text"], out["viewer"])
    if view:
        if write and out["viewer"]:
            _viewer.open_in_browser(res.files[-1])
        else:
            _viewer.view(res, title=os.path.basename(out["prefix"]))
    return res
