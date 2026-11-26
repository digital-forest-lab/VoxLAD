"""Build runnable GNU Octave scripts from the ORIGINAL voxlad_v7.m.

    python make_reference.py /path/to/voxlad_v7.m /path/to/scans OUTDIR

writes into OUTDIR
    a1s?_t10_thin<k>.asc      the two example scans thinned 1:k (angular step x sqrt(k))
    ref_thin<k>_d.m, _s.m     v7 on the thinned scans, double / single precision storage
    ref_full_<nn>.m           v7 on the full scans for 24 voxels, three per script
    run_all.sh                runs everything (several hours in Octave)

Only the edits listed in build() are applied to the script: file paths, the
number of scans, one array size that otherwise stops v7 with an index error,
optionally 'single' -> double storage, an extra debug file, and an optional
restriction to some voxels.  No line that computes a result is changed.

Each run writes <tag>_dbg.txt with one line per voxel and scan:
    x_min y_min z_min scan Aproj Re(G) Im(G) d_bar Nt Nb Ni_l Ni_w Ni_n theta
and the result files of the original script (<tag>_out_1.txt, <tag>_out_2.txt).
"""
import math
import os
import sys

import numpy as np

SCANS = ("a1s1_t10", "a1s2_t10")


def build(src, tag, f1, f2, subset=None, double=False, xrange_=None, fastload=False):
    s = src

    def sub(old, new):
        nonlocal s
        assert old in s, old
        s = s.replace(old, new, 1)

    # 1. inputs and outputs
    sub("nscan=6;", "nscan=2;")
    sub("file1='F:\\TLS data processing\\LAD\\Tram\\pos8_tram_area_voxlad_input.txt' ;", f"file1='{f1}';")
    sub("file2='F:\\TLS data processing\\LAD\\Tram\\pos9_tram_area_voxlad_input.txt' ;", f"file2='{f2}';")
    sub("result_file='C:\\MBeland\\Postdoc\\Data Analysis\\LAD\\tram_area_voxlad_output.txt' ;", f"result_file='{tag}_out.txt';")
    sub("result_file='C:\\MBeland\\Postdoc\\Data Analysis\\LAD\\tram_area_voxlad_output_1.txt';", f"result_file='{tag}_out_1.txt';")
    sub("result_file='C:\\MBeland\\Postdoc\\Data Analysis\\LAD\\tram_area_voxlad_output_2.txt';", f"result_file='{tag}_out_2.txt';")
    # 2. debug file (extra output only)
    sub("parfor split_mat=1:2", f"fiddbg=fopen('{tag}_dbg.txt','w'); tic;\nparfor split_mat=1:2")
    sub("                    entering_pts(file_number)=Nt(file_number)-Nb(file_number);",
        "                    entering_pts(file_number)=Nt(file_number)-Nb(file_number);\n"
        "                    fprintf(fiddbg,'%.6f %.6f %.6f %d %.17g %.17g %.17g %.17g %d %d %d %d %d %.17g\\n',"
        "x_min,y_min,z_min,file_number,Aproj,real(G(file_number)),imag(G(file_number)),d_bar(file_number),"
        "Nt(file_number),Nb(file_number),Ni_l(file_number),Ni_w(file_number),Ni_n(file_number),theta);")
    sub("LA=total_LAD*vox_size^3", "elapsed=toc\nfclose(fiddbg);\nLA=total_LAD*vox_size^3")
    s = s.replace("beep\n", "")
    # 3. one more cell per axis in the occupancy map: without it v7 stops with an index error
    sub("vox_map=zeros(num_x_iteration,num_y_iteration,num_z_iteration, 'uint8');",
        "vox_map=zeros(num_x_iteration+1,num_y_iteration+1,num_z_iteration+1, 'uint8');")
    # 4. optional: double precision storage, as in voxlad_v3.m
    if double:
        for v in "xyz":
            sub(f"PC(file_number).{v}=zeros(1,1, 'single');", f"PC(file_number).{v}=zeros(1,1);")
    # 5. optional: only some voxels
    if subset is not None:
        rows = ";".join("%.6f %.6f %.6f" % tuple(r) for r in subset)
        sub("parfor split_mat=1:2", f"sel=[{rows}];\nparfor split_mat=1:2")
        sub("                z_max=z_min+vox_size;",
            "                if ~any(all(abs(sel-[x_min y_min z_min])<1e-6,2)), continue, end\n                z_max=z_min+vox_size;")
    if xrange_ is not None:
        old = "            y_max=y_min+vox_size;"
        sub(old, f"            if x_min < {xrange_[0]!r} || x_min >= {xrange_[1]!r}, continue, end\n" + old)
    # 6. optional: fill the workspace with a fast loader and use the script's own skip_read switch
    if fastload:
        loader = f"""
files={{'{f1}','{f2}'}};
for file_number=1:2
  fid=fopen(files{{file_number}});
  x_tls(file_number)=fscanf(fid,'%g',[1]); y_tls(file_number)=fscanf(fid,'%g',[1]); z_tls(file_number)=fscanf(fid,'%g',[1]);
  x_1(file_number)=fscanf(fid,'%g',[1]); y_1(file_number)=fscanf(fid,'%g',[1]);
  x_2(file_number)=fscanf(fid,'%g',[1]); y_2(file_number)=fscanf(fid,'%g',[1]);
  angle_units=fscanf(fid,'%s',[1]);
  switch angle_units
    case 'd'
      alpha(file_number)=deg2rad(fscanf(fid,'%g',[1]));
    case 'r'
      alpha(file_number)=fscanf(fid,'%g',[1]);
  end
  npts(file_number)=fscanf(fid,'%u',[1]);
  M=fscanf(fid,'%g',[4 Inf])';
  fclose(fid);
  keep=(M(:,3)<12) & (M(:,3)>-1);       % the script's height filter (max/min_canopy_height)
  PC(file_number).x=M(keep,1)'; PC(file_number).y=M(keep,2)'; PC(file_number).z=M(keep,3)';
  PC(file_number).type=uint8(M(keep,4))';
  npts(file_number)=sum(keep);
end
clear M keep;
max_canopy_height=12; min_canopy_height=-1;
"""
        sub("skip_read=0;", "skip_read=1;" + loader)
    with open(tag + ".m", "w", encoding="latin-1") as fh:
        fh.write(s)


def thin(scan_dir, out_dir, k):
    """Every k-th return, angular step multiplied by sqrt(k) (same pulse density per return)."""
    out = []
    for f in SCANS:
        lines = open(os.path.join(scan_dir, f + ".asc")).read().split("\n")
        hdr, pts = lines[:5], [l for l in lines[5:] if l.strip()]
        pts = pts[::k]
        hdr[3] = "r %.8f" % (float(hdr[3].split()[1]) * math.sqrt(k))
        hdr[4] = str(len(pts))
        fn = os.path.join(out_dir, f"{f}_thin{k}.asc")
        with open(fn, "w") as fh:
            fh.write("\n".join(hdr + pts) + "\n")
        out.append(fn)
    return out


def voxel_subset(scan_dir, n=24, seed=1):
    """How reference/subset_voxels.txt was drawn: n leaf-occupied voxels of v7's
    grid on the full scans, the six with most leaf returns and n-6 at random."""
    rng = np.random.default_rng(seed)
    allp = np.vstack([np.loadtxt(os.path.join(scan_dir, f + ".asc"), skiprows=5) for f in SCANS])
    allp = allp[(allp[:, 2] < 12) & (allp[:, 2] > -1)]
    lo = np.floor(allp[:, :3].min(0))
    leaf = allp[allp[:, 3] == 1]
    ijk = np.floor((leaf[:, :3] - lo) / 0.3).astype(int)
    u, cnt = np.unique(ijk, axis=0, return_counts=True)
    order = np.argsort(cnt)
    pick = np.concatenate([order[-6:], rng.choice(order[:-6], n - 6, replace=False)])
    return lo + u[pick] * 0.3


def main(v7_path, scan_dir, out_dir):
    scan_dir, out_dir = os.path.abspath(scan_dir), os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    src = open(v7_path, encoding="latin-1").read().replace("\r\n", "\n")
    os.chdir(out_dir)
    jobs = []
    for k in (400, 100):
        a, b = thin(scan_dir, out_dir, k)
        build(src, f"ref_thin{k}_d", a, b, double=True)
        build(src, f"ref_thin{k}_s", a, b, double=False)
        jobs += [f"ref_thin{k}_d", f"ref_thin{k}_s"]
    here = os.path.dirname(os.path.abspath(__file__))
    sel = np.loadtxt(os.path.join(here, "reference", "subset_voxels.txt"))
    f1, f2 = (os.path.join(scan_dir, f + ".asc") for f in SCANS)
    for c in range(0, len(sel), 3):
        tag = f"ref_full_{c // 3:02d}"
        build(src, tag, f1, f2, subset=sel[c:c + 3], double=True, fastload=True)
        jobs.append(tag)
    with open("run_all.sh", "w") as fh:
        fh.write("#!/bin/bash\n# each line can be restarted on its own\n")
        for t in jobs:
            fh.write(f"[ -f {t}.done ] || ( octave --no-gui --quiet {t}.m > {t}.log 2>&1 && touch {t}.done )\n")
        fh.write("cat ref_full_??_dbg.txt > full_dbg.txt; cat ref_full_??_out_1.txt > full_out_1.txt; cat ref_full_??_out_2.txt > full_out_2.txt\n")
    print(f"{len(jobs)} Octave scripts written to {out_dir}; start them with: bash run_all.sh")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
