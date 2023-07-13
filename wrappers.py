import tempfile
import subprocess
import os


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
        subprocess.Popen(
            "./fortran_codes/vm_trace", stdin=input_f, stdout=subprocess.DEVNULL
        ).wait()


def vm_tomo(
    vm_file,
    rayfile,
    cmax,
    nxc,
    nzc,
    dz1,
    itop,
    ibot,
    vscal,
    zscal,
    xreg,
    zreg,
    asr,
    reg0,
    reg1,
    reg2,
    crf,
    ch2n,
    vpmin,
    dwsv,
    dwsz,
    mask_file,
    vmnew,
):
    """
    vm_tomo parameters:
        vm_file - path to vm file
        rayfile - path to rayfile
        cmax - max slowness
        nxc - number of horizontal grid points
        nzc - number of vertical grid points
        dz1 - vertical grid spacing near top of model
        itop - shallowest phase
        ibot - deepest phase
        vscal - velocity dimension scaling
        zscal - depth dimension scaling
        xreg - horizontal reach of regularization (km)
        zreg - vertical reach of regularization (km)
        asr - aspect ratio horizontal vs. vertical regularization
        reg0 - relative strength dampening
        reg1 - relative strength flattening
        reg2 - relative strength smoothing
        crf - strength regularization boundaries vs velocities
        ch2n - target chi-squared
        vpmin - (plotting) illuminate all velocities lower than vpmin
        dwsv - derivative weight sum (illumination)
        dwsz -
        mask_file - (plotting) mask part of model without ray coverage
        vmnew - path to write output model
    """
    with tempfile.TemporaryFile("w") as input_f:
        for arg in [
            vm_file,
            rayfile,
            cmax,
            nxc,
            nzc,
            dz1,
            itop,
            ibot,
            vscal,
            zscal,
            xreg,
            zreg,
            asr,
            reg0,
            reg1,
            reg2,
            crf,
            ch2n,
            vpmin,
            dwsv,
            dwsz,
            mask_file,
            vmnew,
        ]:
            input_f.write(f"{arg}\n")
        input_f.seek(0)
        subprocess.Popen("./fortran_codes/vm_tomo", stdin=input_f).wait()
