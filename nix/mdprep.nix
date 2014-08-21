{ buildPythonPackage, fetchgit,
python, pyyaml, mdtraj, prody, ... }:

buildPythonPackage rec {
  version = "2014.07.25.42f4f58";
  basename = "mdprep";
  name = "${basename}-${version}";

  src = fetchgit {
    url = "git://github.com/badi/{basename}.git";
    rev = "42f4f583c849281362d971a28fc576a8060040a8";
    sha256 = "0hkhak8hqdl65xp03qpclbvjddlfldj47qxqzxnzpgcplgryn1dc";
  };

  buildInputs = [ python pyyaml mdtraj prody ];

  meta = {
    description = "Setup for Molecular Dynamics projects";
    homepage = "http://github.com/badi/${basename}";
  };
}