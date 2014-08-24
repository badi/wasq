
let
  pkgs = import <nixpkgs> {};

  licenses = pkgs.lib.licenses;
  platforms = pkgs.lib.platforms;

in

with pkgs;
with python27Packages;

let
  pythonEnv = {
      inherit stdenv buildPythonPackage fetchurl;
      inherit licenses platforms;
      inherit python;
    };

  buildEnv = {
    inherit stdenv fetchurl;
    inherit licenses platforms;
  };

  # depencencies
  mdtraj    = callPackage ./nix/mdtraj.nix  (pythonEnv      // {inherit cython numpy scipy pandas nose;});
  cctools   = callPackage ./nix/cctools.nix (buildEnv       // {inherit which swig libzip zlib perl fuse openssl python;});
  prody     = callPackage ./nix/prody.nix   (pythonEnv      // {inherit numpy;});
  pxul      = callPackage ./nix/pxul.nix     pythonEnv;
  mdprep    = callPackage ./nix/mdprep.nix  (pythonEnv      // { inherit pyyaml mdtraj prody; });
  pwq       = callPackage ./nix/pwq.nix     (pythonEnv      // {inherit cctools pyyaml pyzmq;});
  mdq       = callPackage ./nix/mdq.nix     (pythonEnv      // {inherit cctools mdprep pwq guamps prody mdtraj pyyaml pyzmq;});
  gromacs   = callPackage ./nix/gromacs.nix {inherit stdenv fetchurl cmake; fftw = fftwSinglePrec; };
  guamps    = callPackage ./nix/guamps.nix  (buildEnv       // {inherit cmake gromacs;});

  # profiling
  guppy     = callPackage ./nix/guppy.nix    pythonEnv;

in
{
  devEnv = buildPythonPackage rec {
    version = "0.1.0";
    basename = "wasq";
    name = "${basename}-${version}";
    src = ./.;
    buildInputs = [
      # provided by nixpkgs
      ipython numpy scipy matplotlib pyyaml
      tables pandas
      tkinter

      # externals
      imagemagick
      mdtraj prody cctools
      gromacs
      guppy

      # my libraries
      guamps
      pxul mdprep pwq mdq
    ];
  };
}
