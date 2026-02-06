from pathlib import Path
import numpy as np
import pytest
from mqc_runner import MQCArgs, run_case

# Usage: 
# Run pytest --markers to see registered markers
# e.g. Test all mqc runs in the Shin-Metiu model: pytest -m mqc -v
# e.g. Test SHXF runs in the Shin-Metiu model: pytest -m shxf -v

# For a specific mqc run by invoking the exact case ID: pytest -m mqc -k TEST-SHXF-FG-a-
# The ID format is TEST-{MQC name}-{width scheme}_{momentum jump scheme}
# MQC name - BOMD, Eh, SH, SHXF, EhXF, CT
# width scheme - FG, TD
# momentum jump scheme - e+, e-, v+, v-, p+, p-, a+, a-

REF_ROOT = Path("reference")

KEY_RESCALE = ["e", "v", "p", "a"]
KEY_REJECT = ["+", "-"]
KEY_WIDTH = ["FG", "TD"]

# Define the test matrix you care about
ALL_CASES = [
    # BOMD, Eh
    pytest.param(
        MQCArgs(md=0), "TEST-BOMD", 
        marks=(pytest.mark.mqc, pytest.mark.bomd)
    ),
    
    pytest.param(
        MQCArgs(md=1), "TEST-Eh", 
        marks=(pytest.mark.mqc, pytest.mark.eh)
    ),

    # SH
    *[ 
        pytest.param(
            MQCArgs(md=2, rescale=r, reject=j), f"TEST-SH-{KEY_RESCALE[r]}{KEY_REJECT[j]}", 
            marks=(pytest.mark.mqc, pytest.mark.sh)
        ) 
        for r in range(4) for j in range(2)
     ],

    # SHXF
    *[
        pytest.param(
            MQCArgs(md=3, width=w, rescale=r, reject=j), f"TEST-SHXF-{KEY_WIDTH[w]}-{KEY_RESCALE[r]}{KEY_REJECT[j]}",
            marks=(pytest.mark.mqc, pytest.mark.shxf)
        )
        for w in (0, 1) for r in range(4) for j in range(2)
     ],

    # EhXF
    *[
        pytest.param(
            MQCArgs(md=4, width=w, rescale=r, reject=j), f"TEST-EhXF-{KEY_WIDTH[w]}-{KEY_RESCALE[r]}{KEY_REJECT[j]}",
            marks=(pytest.mark.mqc, pytest.mark.ehxf)
        )
        for w in (0, 1) for r in range(4) for j in range(2)
     ],

    # CT
    pytest.param(
        MQCArgs(md=5), "TEST-CT", marks=(pytest.mark.mqc, pytest.mark.ct)
    )
]

def _load_numeric(path: Path, is_xyz: bool):
    if is_xyz:
        return np.loadtxt(path, skiprows=2, usecols=(1, 2))
    return np.loadtxt(path, skiprows=1)

def _compare_file(out_file: Path, ref_file: Path, tg: str):

    is_xyz = tg.endswith(".xyz")
    out_data = _load_numeric(out_file, is_xyz)
    ref_data = _load_numeric(ref_file, is_xyz)

    assert out_data.shape == ref_data.shape, f"Shape mismatch: {tg}"

    # Compare absolute values for phase-dependent quantities (eigenvector phase convention)
    if tg in ["NACME", "BOCOH"]:
        out_data = np.abs(out_data)
        ref_data = np.abs(ref_data)

    np.testing.assert_allclose(out_data, ref_data, rtol=1e-7, atol=1e-9)

@pytest.mark.parametrize("args,case_id", ALL_CASES)
def test_mqc_case(args, case_id):

    run_case(args)
    #return
    # Determine what to compare
    targets = ["MDENERGY", "FINAL.xyz"]

    if args.md != 0:    # not BOMD
        targets += ["BOPOP", "NACME"]

    if args.md in (2, 3, 4):   # have the SH feature
        targets += ["SHSTATE", "SHPROB"]

    # Compare test results and the reference
    case_id = Path(case_id)
    for tg in targets:
        if args.md != 5:    # not CT
            _compare_file(case_id / "md" / tg, REF_ROOT / case_id / "md" / tg, tg)
        else:    # CT
            _compare_file(case_id / "TRAJ_1" / "md" / tg, REF_ROOT / case_id / "TRAJ_1" / "md" / tg, tg)
            _compare_file(case_id / "TRAJ_2" / "md" / tg, REF_ROOT / case_id / "TRAJ_2" / "md" / tg, tg)

