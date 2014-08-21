{buildPythonPackage, fetchurl,
python, cctools, pyyaml, pyzmq,
...}:

buildPythonPackage rec {
  basename = "pwq";
  revision = "1ccfc79";
  version = "2014.04.04.${revision}";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.zip";
    sha256 = "3534ea88d86601e8ecf11fce09dfcc9e0e24fcd25cc8a05bc47f33b277c4fba4";
  };

  buildInputs = [ python cctools pyyaml pyzmq ];

  meta = {
    description = "Process-safe Python interface to the CCL WorkQueue library.";
    homepage = "http://github.com/badi/${basename}";
  };
}