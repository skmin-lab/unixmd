# mqc_runner.py (recommended new module)
from dataclasses import dataclass
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path
import os, copy
import random

from molecule import Molecule
import qm, mqc
from misc import data

@dataclass(frozen=True)
class MQCArgs:
    md: int = 3      # the MQC method (0: BOMD, 1: Eh, 2: SH, 3: SHXF, 4: EhXF, 5: CT)
    rescale: int = 3 # the hop rescale option (0: energy, 1: velocity, 2: momentum, 3: augment)
    reject: int = 1  # the hop reject option (0: keep, 1: reverse)
    width: int = 0   # the width scheme in XF (0: frozen Gaussian as 0.1 Bohr, 1: TD)

MOMENTUM_JUMP_SCHEME = {
    (0, 0): ("-e+", "energy",   "keep"),
    (1, 0): ("-v+", "velocity", "keep"),
    (2, 0): ("-p+", "momentum", "keep"),
    (3, 0): ("-a+", "augment",  "keep"),
    (0, 1): ("-e-", "energy",   "reverse"),
    (1, 1): ("-v-", "velocity", "reverse"),
    (2, 1): ("-p-", "momentum", "reverse"),
    (3, 1): ("-a-", "augment",  "reverse"),
}

def build_mqc(args: MQCArgs, nsteps: int = 10):

    # Define initial conditions for the Shin-Metiu model
    data["X1"] = 1836  # au
    geom = """
1
Shin-Metiu model
X1       -2.0     0.02
"""
    mol = Molecule(geometry=geom, ndim=1, nstates=2, ndof=1, unit_pos="au", l_model=True)
    qm_model = qm.model.Shin_Metiu(molecule=mol)
    
    # Decide which MQC is used
    out_dir = "TEST"
    if args.md == 0:    # BOMD
        out_dir += "-BOMD"
        md = mqc.BOMD(molecule=mol, nsteps=nsteps, dt=5.0, unit_dt="au", istate=1)

    elif args.md == 1:  # Eh
        out_dir += "-Eh"
        md = mqc.Eh(molecule=mol, nsteps=nsteps, dt=5.0, unit_dt="au", istate=1)

    elif args.md == 2:  # SH
        out_dir += "-SH"
        md = mqc.SH(molecule=mol, nsteps=nsteps, nesteps=1, dt=5.0, unit_dt="au", istate=1)
        suffix, md.hop_rescale, md.hop_reject = MOMENTUM_JUMP_SCHEME[(args.rescale, args.reject)]
        out_dir += suffix

    elif args.md == 3:  # SHXF
        out_dir += "-SHXF"
        if args.width == 0:
            md = mqc.SHXF(molecule=mol, nsteps=nsteps, nesteps=1, dt=5.0, unit_dt="au", sigma=0.1, istate=1)
            out_dir += "-FG"
        elif args.width == 1:
            md = mqc.SHXF(molecule=mol, nsteps=nsteps, nesteps=1, dt=5.0, unit_dt="au", l_td_sigma=True, istate=1)
            out_dir += "-TD"
        else:
            raise ValueError(f"Invalid width={args.width} for SHXF")
        suffix, md.hop_rescale, md.hop_reject = MOMENTUM_JUMP_SCHEME[(args.rescale, args.reject)]
        out_dir += suffix

    elif args.md == 4:  # EhXF
        out_dir += "-EhXF"
        if args.width == 0:
            md = mqc.EhXF(molecule=mol, nsteps=nsteps, nesteps=1, dt=5.0, unit_dt="au", sigma=0.1, istate=1)
            out_dir += "-FG"
        elif args.width == 1:
            md = mqc.EhXF(molecule=mol, nsteps=nsteps, nesteps=1, dt=5.0, unit_dt="au", l_td_sigma=True, istate=1)
            out_dir += "-TD"
        else:
            raise ValueError(f"Invalid width={args.width} for EhXF")
        suffix, md.hop_rescale, md.hop_reject = MOMENTUM_JUMP_SCHEME[(args.rescale, args.reject)]
        out_dir += suffix

    elif args.md == 5:  # CT
        out_dir += "-CT"
        mol1 = copy.deepcopy(mol)
        mol1.pos[0, 0] = -1.9
        md = mqc.CT(molecules=[mol, mol1], nsteps=nsteps, nesteps=1, dt=5.0, unit_dt="au", istates=[1, 1])

    else:
        raise ValueError(f"Invalid md={args.md}")

    return qm_model, md, out_dir

def run_case(args: MQCArgs, nsteps: int = 10):
    qm, md, out_dir = build_mqc(args, nsteps=nsteps)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "run.log"

    cwd = os.getcwd()
    try:
        random.seed(1000)  # keep random numbers for SH tests
        
        with open(log_path, "w") as f, redirect_stdout(f), redirect_stderr(f):
            md.run(qm=qm, output_dir=str(out_dir))
        #md.run(qm=qm, output_dir=out_dir)
    finally:
        os.chdir(cwd)

