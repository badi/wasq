{pkgs}:

with pkgs; with python27Packages;
buildPythonPackage rec {
  version = "0.9.0";
  baseName = "mdtraj";
  name = "${baseName}-${version}";

  src = fetchurl {
    url = "https://github.com/rmcgibbo/${baseName}/archive/${version}.tar.gz";
    sha256 = "e687c31e32dc7e4486ee9277d7df51ed8c303ce6112ae03f9e6e4514d90690d8";
  };

  buildInputs = [ cython ];
  propagatedBuildInputs = [ numpy scipy pandas nose ];

  meta = {
    description = "Read, write and analyze Molecular Dynamics trajectories";
    homepage = "http://mdtraj.org/${version}";
    longDescription = ''
      MDTraj is a python library that allows users to manipulate
      molecular dynamics (MD) trajectories and perform a variety of
      analyses, including fast RMSD, solvent accessible surface
      area, hydrogen bonding, etc. A highlight of MDTraj is the wide
      variety of molecular dynamics trajectory file formats which
      are supported, including RCSB pdb, GROMACS xtc and trr, CHARMM
      / NAMD dcd, AMBER binpos, AMBER NetCDF, AMBER mdcrd, TINKER
      arc and MDTraj HDF5.

      MDTraj is based on numpy and is both easy to use and
      fast. Core routines like RMSD are written in C with explicit
      SSE vectorization and multicore parallelization. The RMSD
      code, in particular, is based on the Theobald QCP method and
      is 4x the speed of the original Theobald code and over 3x as
      fast as Theobald code modified to use GotoBLAS.

      The library also ships with a flexible command-line
      application for converting trajectories between formats. When
      you install MDTraj, the script will be installed under the
      name mdconvert.
      '';

    license = pkgs.lib.licenses.lgpl2;
    platforms = pkgs.lib.platforms.all;
    maintainers = [ "Badi' Abdul-Wahid <abdulwahidc@gmail.com>" ];
  };
}