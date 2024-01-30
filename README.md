# USP-JEWEL
Jet propagation (currently done using an unofficial modification of JEWEL) in realistic viscous hydrodynamics modeled using vUSPhydro

---

Disclaimers: 
- Send an email to leonardo.barreto.campos@usp.br for compatible medium profiles. Files are large and we are developing an easier way to distribute them.

---

## Installation 
### Requirements
- [LHAPDF 6](https://lhapdf.hepforge.org/install.html) and desired PDF sets.
- [JEWEL 2.4.0 source code](https://jewel.hepforge.org/downloads/).

### Set up
1. Clone this repository and add the following original JEWEL `meix.f` and `pythia6425mod-lhapdf6.f` files to the repo directory, as we do not distribute those,

```bash
git clone git@github.com:leo-barreto/USP-JEWEL.git
cp -t USP-JEWEL path/jewel-x.y.z/meix.f path/jewel-x.y.z/pythia6425mod-lhapdf6.f
```

2. Edit the `LHAPDF_PATH := ` line in `Makefile` to set your LHAPDF library path (e.g. `path/lhapdf/lib`).

3. Compile USP-JEWEL

```bash
cd USP-JEWEL
make
```

4. Make sure to setup LHAPDF environment paths (just like JEWEL installation)

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/lhapdf/lib
export LHAPATH=/path/lhapdf/share/lhapdf
```

### Running the example
Let's see if everything is setup correctly by running an example. The folder `example/` has parameters files for the JEWEL (`params_example.dat` and the medium `medium_params_example.dat` (similar to original JEWEL), and a *TReNTov-USPhydro* PbPb 0-10% 5.02 TeV medium profile. For this test, the nuclear PDF *EPPS21nlo_CT18Anlo_Pb208* (referenced by the parameter `PDFSET 904400`) must be provided by LHAPDF (get it from [here](https://lhapdf.hepforge.org/pdfsets.html) and follow LHAPDF set up instructions). The file is rather large so be sure to decompress it (should be .dat file!) before usage. 

```bash
cd example/
unxz 0-10_example.dat.xz
.././usp-jewel params_example.dat
```

After the execution is completed, the expected outputs are `out.log`, logging multiple events, and the final-particle distribution `out.hepmc`. And it is done! You now can simulate **JEWEL parton shower evolution in realistic 2+1D hydrodynamic profiles**!
