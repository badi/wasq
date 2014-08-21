{ stdenv, fetchurl, cmake,
  singlePrec ? true,
  staticLib ? true,
  fftw
}:


stdenv.mkDerivation rec {
  basename = "gromacs";
  version = "4.5.7";
  name = "${basename}-${version}";

  src = fetchurl {
    url = "ftp://ftp.gromacs.org/pub/gromacs/${name}.tar.gz";
    md5 = "24febafaf51be785b1c755ef679bea08";
  };

  buildInputs = [cmake fftw];

  flagPrecision = if singlePrec then "-DGMX_DOUBLE=OFF" else "-DGMX_DOUBLE=ON -DGMX_DEFAULT_SUFFIX=OFF";
  flagStatic = if staticLib then "-DBUILD_SHARED_LIBS=OFF -DGMX_PREFER_STATIC_LIBS=ON" else "";

  cmakeFlags = "${flagPrecision} ${flagStatic}";

  meta = {
    homepage    = "http://www.gromacs.org";
    licence     = "GPLv2";
    description = "The GROMACS molecular dynamics software package";
    longDescription = ''
      GROMACS is a versatile package to perform molecular dynamics,
      i.e. simulate the Newtonian equations of motion for systems
      with hundreds to millions of particles.

      It is primarily designed for biochemical molecules like
      proteins, lipids and nucleic acids that have a lot of
      complicated bonded interactions, but since GROMACS is
      extremely fast at calculating the nonbonded interactions (that
      usually dominate simulations) many groups are also using it
      for research on non-biological systems, e.g. polymers.

      GROMACS supports all the usual algorithms you expect from a
      modern molecular dynamics implementation, (check the online
      reference or manual for details), but there are also quite a
      few features that make it stand out from the competition.

      See: http://www.gromacs.org/About_Gromacs for details.
    '';
  };
}
