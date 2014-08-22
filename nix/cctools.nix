{stdenv, fetchurl, licenses,
which, swig, python, libzip, zlib, perl, fuse, openssl,
...}:

stdenv.mkDerivation rec {
  version = "3.7.5.badi-${revision}";
  basename = "cctools";
  name = "${basename}-${version}";
  revision = "1512c37";

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.tar.gz";
    sha256 = "1gbwx8r094ga143bkam9qfpg7jddpzjvnzz3226cpnygcfxrry8g";
  };

  buildInputs = [ stdenv which swig python libzip zlib perl fuse openssl ];

  configurePhase = ''
    ./configure --prefix $out \
       --with-python-path ${python} --with-perl-path ${perl} --with-zlib-path ${zlib}
  '';

  meta = {
    description = "Collection of software tools from the Cooperative Computing Lab at the University of Notre Dame";
    homepage = "http://ccl.cse.nd.edu";
    license = licenses.gpl2;
    maintainers = [ "Badi' Abdul-Wahid <abdulwahidc@gmail.com>" ];
  };
}