{buildPythonPackage, fetchgit,
python,
...}:

buildPythonPackage rec {
  version = "2014.08.19.30b4caa";
  baseName = "mdtraj";
  name = "${baseName}-${version}";

  src = fetchgit {
    url = "git://github.com/badi/pxul.git";
    rev = "30b4caa2a35c441c01dab7daf16abc3638049f23";
    sha256 = "0hkhak8hqdl65xp03qpclbvjddlfldj47qxqzxnzpgcplgryn1dc";
  };

  buildInputs = [ python ];

  meta = {
    description = "Python eXtras and UtiLities";
    homepage = "http://github.com/badi/pxul";
  };
}