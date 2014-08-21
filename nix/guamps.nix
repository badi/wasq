{stdenv, fetchurl, licenses,
cmake, gromacs,
...}:

stdenv.mkDerivation rec {
  basename = "guamps";
  version = "2014.07.29.${revision}";
  revision = "e132aee";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "https://github.com/badi/${basename}/archive/${revision}.tar.gz";
    sha256 = "1rb4zyviqxp7c23gh4gm5ccl7py7cqadj5f0bnhqqqmps3hxwhkq";
  };

  buildInputs = [ cmake gromacs ];

  cmakeFlags = "-DGROMACS_LIBRARY=${gromacs}/lib/libgmx.a -DGROMACS_INCLUDE_DIR=${gromacs}/include";

  meta = {
    description = "Gromacs Utilities Are a Major Pain in the Shins";
    homepage = "http://github.com/badi/${basename}";
    maintainers = [ "Badi' Abdul-Wahid <abdulwahidc@gmail.com>" ];
  };
}