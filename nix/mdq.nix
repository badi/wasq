{ buildPythonPackage, fetchurl,
python, cctools, mdtraj, prody, mdprep, pwq,
...}:

buildPythonPackage rec {
  basename = "mdq";
  revision = "3d73933";
  version = "2014.04.05.${revision}";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.zip";
    sha256 = "05cdyp7rf4gz5640y7ryyfnyj1696fkmpj4amdlkjnpdchq30vsj";
  };

  buildInputs = [ python cctools mdtraj prody mdprep pwq ];

  meta = {
    description = "WorkQueue for running MD Tasks";
    homepage = "http://github.com/badi/${basename}";
  };
}