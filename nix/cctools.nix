{pkgs}:

with pkgs;
stdenv.mkDerivation rec {
  version = "4.2.2";
  basename = "cctools";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "http://ccl.cse.nd.edu/software/files/cctools-4.2.2-source.tar.gz";
    sha256 = "89e14601e258885f39fd9e3cb85d8314e37b4869535c2416b6319039a98f14f8";
  };

  buildInputs = [ stdenv which swig python27Packages.python libzip zlib perl fuse openssl ];

  configurePhase = ''
    ./configure --prefix $out \
       --with-python-path ${python27Packages.python} --with-perl-path ${perl} --with-zlib-path ${zlib}
  '';

  meta = {
    description = "Collection of software tools from the Cooperative Computing Lab at the University of Notre Dame";
    homepage = "http://ccl.cse.nd.edu";
    license = pkgs.lib.licenses.gpl2;
    maintainers = [ "Badi' Abdul-Wahid <abdulwahidc@gmail.com>" ];
  };
}