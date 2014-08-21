{ buildPythonPackage, fetchurl,
python, cctools, mdprep, pwq,
mdtraj, prody, pyyaml, pyzmq,
...}:

buildPythonPackage rec {
  basename = "mdq";
  revision = "3d73933";
  version = "2014.04.05.${revision}";
  name = "${basename}-${version}";
  doCheck = true;

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.tar.gz";
    sha256 = "11ws2l7c9xdizkw61paxsf2gp14bagndhpr5d8kjvkn0q30xags8";
  };

  buildInputs = [ python cctools mdprep pwq mdtraj prody pyyaml pyzmq ];

  meta = {
    description = "WorkQueue for running MD Tasks";
    homepage = "http://github.com/badi/${basename}";
  };
}