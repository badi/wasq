
let
  pkgs = import <nixpkgs> {};

  licenses = pkgs.lib.licenses;
  platforms = pkgs.lib.platforms;

  pythonEnv = {
      inherit (pkgs) stdenv buildPythonPackage fetchurl;
      inherit (pkgs.lib) licenses platforms;
      inherit python;
    };

  buildEnv = {
    inherit (pkgs) stdenv fetchurl;
    inherit (pkgs.lib) licenses platforms;
  };

  pythonPackages = pkgs.python27Packages;
  python = pythonPackages.python;
  pyyaml = pythonPackages.pyyaml;
  pyzmq = pythonPackages.pyzmq;
  cython = pythonPackages.cython;
  numpy = pythonPackages.numpy;
  scipy = pythonPackages.scipy;
  pandas = pythonPackages.pandas;
  nose = pythonPackages.nose;


  mdtraj = pkgs.callPackage ./nix/mdtraj.nix (pythonEnv // {inherit cython numpy scipy pandas nose;});
  cctools = pkgs.callPackage ./nix/cctools.nix (buildEnv // {
    inherit (pkgs) which swig libzip zlib perl fuse openssl;
    inherit python;
  });
  prody = pkgs.callPackage ./nix/prody.nix (pythonEnv // {inherit numpy;});
  pxul = pkgs.callPackage ./nix/pxul.nix pythonEnv;
  mdprep = pkgs.callPackage ./nix/mdprep.nix (pythonEnv // { inherit pyyaml mdtraj prody; });
  pwq = pkgs.callPackage ./nix/pwq.nix (pythonEnv // {inherit cctools pyyaml pyzmq;});
  mdq = pkgs.callPackage ./nix/mdq.nix (pythonEnv // {inherit cctools mdtraj mdprep pwq;});
  gromacs = with pkgs;
            pkgs.callPackage ./nix/gromacs.nix {inherit stdenv fetchurl cmake;
	                                        fftw = fftwSinglePrec; };
  guamps = with pkgs; pkgs.callPackage ./nix/guamps.nix (buildEnv // {inherit cmake gromacs;});

in
with pkgs;
with python27Packages;
{
  devEnv = buildPythonPackage rec {
    version = "0.1.0";
    basename = "wasq";
    name = "${basename}-${version}";
    src = ./.;
    buildInputs = [
      # provided by nixpkgs
      ipython numpy scipy matplotlib pyyaml

      # externals
      mdtraj prody cctools
      gromacs

      # my libraries
      guamps
      pxul mdprep pwq
    ];
  };
}
