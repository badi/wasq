{ buildPythonPackage, fetchurl,
python, pyyaml, mdtraj, prody, ... }:

buildPythonPackage rec {
  basename = "mdprep";
  revision = "42f4f58";
  version = "2014.07.25.${revision}";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.zip";
    sha256 = "014zqxldka23k5ma0gabddz8zm9d2xa82p8xd0916nci2cwvj4kg";
  };

  buildInputs = [ python pyyaml mdtraj prody ];

  meta = {
    description = "Setup for Molecular Dynamics projects";
    homepage = "http://github.com/badi/${basename}";
  };
}