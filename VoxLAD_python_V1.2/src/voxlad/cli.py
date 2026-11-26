"""Command line: voxlad run | view | inspect | template."""
from __future__ import annotations

import argparse
import sys

from . import __version__
from .io import InputError


def _cmd_run(args):
    from . import config, model, report
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    settings, specs, out = config.load_config(args.config)
    if args.compat:
        settings.compat = args.compat
    if args.voxel_size:
        settings.voxel_size = args.voxel_size
    if args.threads:
        settings.threads = args.threads
    if args.output:
        out["prefix"] = args.output
    scans = config.load_scans(specs, say)
    res = model.compute(scans, settings, say)
    files = report.write_results(res, out["prefix"], out["per_scan"], out["npz"], out["matlab_text"], out["viewer"])
    if not args.quiet:
        print()
        print(report.summary_text(res))
    for f in files:
        say("wrote " + f)
    if out["viewer"] and out["open_viewer"] and not args.no_view:
        from . import viewer
        if viewer.open_in_browser(files[-1]):
            say("opened the 3-D viewer in the web browser")
        else:
            say("no web browser could be started here; open " + files[-1] + " in a browser to see the result in 3-D")
    return 0


def _cmd_view(args):
    """Open the 3-D viewer of a finished run, writing it from the voxel table if needed."""
    import os
    from . import config, viewer
    path = args.path
    low = path.lower()
    if low.endswith(".toml"):
        _, _, out = config.load_config(path)
        page, table = out["prefix"] + "_viewer.html", out["prefix"] + "_voxels.csv"
        if not os.path.isfile(page) and not os.path.isfile(table):
            raise InputError(f"no result found for {path}: run 'voxlad run {path}' first")
    elif not os.path.isfile(path):
        raise InputError(f"file not found: {path}")
    elif low.endswith(".csv"):
        table = path
        page = (path[:-len("_voxels.csv")] if low.endswith("_voxels.csv") else os.path.splitext(path)[0]) + "_viewer.html"
    elif low.endswith((".html", ".htm")):
        page = path
        table = path[:-len("_viewer.html")] + "_voxels.csv" if low.endswith("_viewer.html") else ""
    else:
        raise InputError("give a configuration file (.toml), a voxel table (_voxels.csv) or a viewer file (.html)")
    have_table = bool(table) and os.path.isfile(table)
    stale = have_table and os.path.isfile(page) and os.path.getmtime(table) > os.path.getmtime(page)
    if not os.path.isfile(page) or stale or args.rebuild:
        if not have_table:
            raise InputError(f"cannot write {page}: the voxel table {table or '(_voxels.csv)'} was not found")
        viewer.write_viewer_from_csv(table, page)
        print("wrote " + page)
    if args.no_open:
        return 0
    if viewer.open_in_browser(page):
        print("opened " + page)
    else:
        print("no web browser could be started here; open " + page + " in a browser")
    return 0


def _cmd_inspect(args):
    from . import config, scaninfo
    settings, specs, out = config.load_config(args.config)
    for spec in specs:
        print(scaninfo.describe(config.load_scan(spec)))
        print()
    return 0


def _cmd_template(args):
    from .config import TEMPLATE
    if args.file:
        import os
        if os.path.exists(args.file) and not args.force:
            raise InputError(f"{args.file} exists; use --force to overwrite it")
        with open(args.file, "w") as fh:
            fh.write(TEMPLATE)
        print("wrote " + args.file)
    else:
        sys.stdout.write(TEMPLATE)
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="voxlad",
        description="Leaf area density per voxel from classified terrestrial lidar scans "
                    "(Beland, Widlowski & Fournier 2014).")
    p.add_argument("--version", action="version", version="voxlad " + __version__)
    sub = p.add_subparsers(dest="command")
    r = sub.add_parser("run", help="run the model described by a configuration file")
    r.add_argument("config")
    r.add_argument("-o", "--output", help="output prefix (overrides [output] prefix)")
    r.add_argument("--compat", choices=["none", "matlab"], help="override [model] compat")
    r.add_argument("--voxel-size", type=float, help="override [grid] voxel_size (m)")
    r.add_argument("--threads", type=int, help="number of threads (default: all)")
    r.add_argument("--no-view", action="store_true", help="do not open the 3-D viewer at the end of the run")
    r.add_argument("-q", "--quiet", action="store_true", help="print nothing but errors")
    r.set_defaults(func=_cmd_run)
    v = sub.add_parser("view", help="open the 3-D viewer of a finished run")
    v.add_argument("path", help="configuration file, <prefix>_voxels.csv or <prefix>_viewer.html")
    v.add_argument("--rebuild", action="store_true", help="write the page again from the voxel table even if it exists")
    v.add_argument("--no-open", action="store_true", help="write the viewer file if needed, without opening it")
    v.set_defaults(func=_cmd_view)
    i = sub.add_parser("inspect", help="describe the scans of a configuration file and check their settings")
    i.add_argument("config")
    i.set_defaults(func=_cmd_inspect)
    t = sub.add_parser("template", help="print a commented configuration file, or write it to FILE")
    t.add_argument("file", nargs="?")
    t.add_argument("--force", action="store_true")
    t.set_defaults(func=_cmd_template)
    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 1
    try:
        return args.func(args)
    except InputError as err:
        print("voxlad: error: " + str(err), file=sys.stderr)
        return 2
    except BrokenPipeError:            # output piped into a program that stopped reading
        return 0


if __name__ == "__main__":
    sys.exit(main())
