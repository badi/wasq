{buildPythonPackage, fetchurl,
python,  numpy,
...}:

buildPythonPackage rec {
  version = "1.5.1";
  basename = "ProDy";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "http://pypi.python.org/packages/source/P/${basename}/${basename}-${version}.tar.gz";
    sha256 = "9fccfcff0df016c69e525a706a2b017db9cc54544cb81685047478c5aae53bb7";
  };

  buildInputs = [ python numpy ];

  meta = {
    description = "A Python Package for Protein Dynamics Analysis";
    homepage = "http://prody.csb.pitt.edu";
  };
}