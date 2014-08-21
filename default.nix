
let
  pkgs = import <nixpkgs> {};
  mdtraj = pkgs.callPackage ./nix/mdtraj.nix { inherit pkgs; };
in
with pkgs;
with python27Packages;
{
  devEnv = buildPythonPackage rec {
    version = "0.1.0";
    basename = "wasq";
    name = "${basename}-${version}";
    src = ./.;
    buildInputs = [ ipython numpy scipy matplotlib mdtraj ];
  };
}
