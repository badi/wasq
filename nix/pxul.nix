{buildPythonPackage, fetchurl,
python,
...}:

buildPythonPackage rec {
  basename = "pxul";
  revision = "30b4caa";
  version = "2014.08.19.${revision}";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.zip";
    sha256 = "09zwws9vfaakc90cbkin2xrbyng071413yfl76psmk9gilkbr6hb";
  };

  buildInputs = [ python ];

  meta = {
    description = "Python eXtras and UtiLities";
    homepage = "http://github.com/badi/${basename}";
  };
}