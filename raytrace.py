import tempfile
import subprocess
import os

# Goal: Wrap vm_trace so that I can raytrace a velocity model directly from Python.


def vm_trace(
    vm_file,
    vtop,
    drp,
    isec,
    itop,
    ibot,
    ins,
    xin,
    zin,
    pfile,
    npx,
    sfile,
    nsh,
    new,
    rayfile,
):
    """
    vm_trace parameters:
        vm_file - path to vm file
        vtop - max slowness
        drp - graph grid spacing
        isec - seed for RNG
        itop,ibot - shallowest, deepest phases
        ins - instrument number
        xin,zin - x, z position of instrument in model space
        pfile - path to picks file
        npx - number of picks
        sfile - shot file
        nsh - number of shots
        new - 1 to start new rayfile, 0 to append to existing rayfile
        rayfile - path to rayfile
    """
    with tempfile.TemporaryFile("w") as input_f:
        for arg in [
            vm_file,
            vtop,
            drp,
            isec,
            itop,
            ibot,
            ins,
            xin,
            zin,
            pfile,
            npx,
            sfile,
            nsh,
            new,
            rayfile,
        ]:
            input_f.write(f"{arg}\n")
        input_f.seek(0)
        subprocess.Popen("./fortran_codes/vm_trace", stdin=input_f)
