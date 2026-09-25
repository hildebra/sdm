#!/usr/bin/env bash
# Builds a fully static sdm with HTSlib (SAM/BAM/CRAM input) and libdeflate from a micromamba env.
#
#   ./build_static.sh -j8           # build ./sdm
#   ./build_static.sh -j8 test      # build and run all test suites with the same toolchain
#
# Arguments go to make. The first run creates the env (conda-forge: GCC with a static glibc
# sysroot, make, cmake, pkg-config, curl, and zlib/bzip2/xz with static archives) and builds
# libdeflate and HTSlib from source into it: conda-forge's libdeflate ships no libdeflate.a, and
# bioconda's htslib needs libcurl/OpenSSL, which conda does not provide as static archives.
# HTSlib is built without remote-URL (libcurl, S3, GCS) and plugin support. Later runs only
# rebuild sdm. Start over with: micromamba env remove -y -n sdm-static
#
# Environment:
#   SDM_ENV             env name (default sdm-static); an existing env gets missing packages
#   SDM_MICROMAMBA      micromamba executable (default: $MAMBA_EXE, else micromamba on PATH)
#   HTSLIB_VERSION      default 1.24; set HTSLIB_SHA256 to verify the download of other versions
#   LIBDEFLATE_VERSION  default 1.26; LIBDEFLATE_SHA256 likewise
set -euo pipefail

die() { echo "build_static.sh: $*" >&2; exit 1; }

env_name=${SDM_ENV:-sdm-static}
htslib_version=${HTSLIB_VERSION:-1.24}
libdeflate_version=${LIBDEFLATE_VERSION:-1.26}
# Checksums of the default versions, as in the bioconda htslib and conda-forge libdeflate recipes.
if [[ $htslib_version == 1.24 ]]; then
    : "${HTSLIB_SHA256:=28a8de191381c7a97a35675ceac76fa1ea95e7b678d6a2e9d600a7874e4077de}"
fi
if [[ $libdeflate_version == 1.26 ]]; then
    : "${LIBDEFLATE_SHA256:=bba03fffc5538576213675ce6968fcff6ce2e67d82e4d5febea2d05f9f13cf85}"
fi
packages=(gcc_linux-64 gxx_linux-64 sysroot_linux-64 make cmake pkg-config curl
          zlib bzip2 liblzma-devel liblzma-static)
target=x86_64-conda-linux-gnu
required=(bin/make bin/cmake bin/pkg-config bin/curl bin/bzip2
          "bin/$target-gcc" "bin/$target-g++" "bin/$target-readelf"
          lib/libz.a lib/libbz2.a lib/liblzma.a include/zlib.h include/bzlib.h include/lzma.h)

[[ $(uname -s)-$(uname -m) == Linux-x86_64 ]] || die "needs x86-64 Linux (conda-forge *_linux-64 compilers)"
cd "$(dirname "${BASH_SOURCE[0]}")"

# The micromamba shell function and aliases (mima) are not visible in a script, but
# "micromamba shell init" exports the executable's path as MAMBA_EXE.
micromamba=${SDM_MICROMAMBA:-${MAMBA_EXE:-micromamba}}
command -v "$micromamba" > /dev/null || die "micromamba not found (put it on PATH or set SDM_MICROMAMBA)"

# Rows of "micromamba env list" are "name [*] path".
env_prefix() { "$micromamba" env list 2> /dev/null | awk -v n="$env_name" '$1 == n { print $NF; exit }'; }
missing_files() {
    local f
    for f in "${required[@]}"; do [[ -e $prefix/$f ]] || printf '%s ' "$f"; done
}

prefix=$(env_prefix)
if [[ -z $prefix || ! -d $prefix/conda-meta ]]; then
    echo "== creating micromamba env $env_name"
    "$micromamba" create -y -n "$env_name" --override-channels -c conda-forge "${packages[@]}"
    prefix=$(env_prefix)
    [[ -n $prefix && -d $prefix/conda-meta ]] || die "cannot locate micromamba env $env_name"
fi
missing=$(missing_files)
if [[ -n $missing ]]; then
    echo "== adding packages to micromamba env $env_name (missing: $missing)"
    "$micromamba" install -y -n "$env_name" --override-channels -c conda-forge "${packages[@]}"
    missing=$(missing_files)
    [[ -z $missing ]] || die "still missing in $prefix: $missing"
fi

# Use the env directly rather than activating it: this replaces the compiler flags an activated
# env exports (-O2, -DNDEBUG, ...) with the Makefile's defaults.
export CONDA_PREFIX=$prefix PATH=$prefix/bin:$PATH PKG_CONFIG_PATH=$prefix/lib/pkgconfig
export CC=$prefix/bin/$target-gcc CXX=$prefix/bin/$target-g++
export AR=$prefix/bin/$target-ar RANLIB=$prefix/bin/$target-ranlib
export CPPFLAGS="-I$prefix/include" CFLAGS=-O2 LDFLAGS="-L$prefix/lib"
unset CXXFLAGS LDLIBS

fetch() {  # url file sha256
    curl -fL --retry 3 -o "$2" "$1"
    if [[ -n $3 ]]; then
        echo "$3  $2" | sha256sum -c - || die "checksum mismatch for $2"
    else
        echo "warning: $2 is not verified (no checksum given)" >&2
    fi
}

deps="libdeflate-$libdeflate_version htslib-$htslib_version"
stamp=$prefix/share/sdm-static-deps
if [[ $(cat "$stamp" 2> /dev/null) != "$deps" ]]; then
    ncpu=$(nproc 2> /dev/null || echo 4)
    work=$(mktemp -d)
    trap 'rm -rf "$work"' EXIT
    (
        cd "$work"
        echo "== building static libdeflate $libdeflate_version"
        fetch "https://github.com/ebiggers/libdeflate/archive/v$libdeflate_version.tar.gz" \
            libdeflate.tar.gz "${LIBDEFLATE_SHA256:-}"
        tar -xzf libdeflate.tar.gz
        # PIC because HTSlib also builds a shared library, which links libdeflate.a in.
        cmake -S "libdeflate-$libdeflate_version" -B libdeflate-build -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_INSTALL_LIBDIR=lib \
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
            -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF
        cmake --build libdeflate-build -j "$ncpu"
        cmake --install libdeflate-build

        echo "== building static HTSlib $htslib_version"
        fetch "https://github.com/samtools/htslib/releases/download/$htslib_version/htslib-$htslib_version.tar.bz2" \
            htslib.tar.bz2 "${HTSLIB_SHA256:-}"
        tar -xjf htslib.tar.bz2
        cd "htslib-$htslib_version"
        ./configure --prefix="$prefix" --with-libdeflate --disable-libcurl --disable-gcs \
            --disable-s3 --disable-plugins --disable-ref-cache
        make -j "$ncpu" AR="$AR" RANLIB="$RANLIB"
        make install AR="$AR" RANLIB="$RANLIB"
    )
    mkdir -p "$(dirname "$stamp")"
    echo "$deps" > "$stamp"
fi

hts_cflags=$(pkg-config --cflags htslib) || die "pkg-config cannot find htslib in $prefix"
hts_libs=$(pkg-config --static --libs htslib) || die "pkg-config cannot resolve htslib's static libraries"
echo "== building sdm (HTS_LIBS: $hts_libs)"
# Objects go to .build/conda-<the Makefile's configuration name>, apart from builds with other
# compilers in this checkout. The single quotes leave the $(...) to make.
make STATIC=1 WITH_HTS=1 WITH_LIBDEFLATE=1 CC="$CC" CXX="$CXX" CXXFLAGS=-O3 \
    'BUILD_DIR=.build/conda-$(if $(filter 1,$(WITH_HTS)),hts,portable)-$(STATIC)' \
    HTS_CFLAGS="$hts_cflags" HTS_LIBS="$hts_libs" "$@"

if [[ -e sdm ]]; then
    if [[ $("$prefix/bin/$target-readelf" -d sdm) == *NEEDED* ]]; then
        die "sdm is dynamically linked (make STATIC=0 given?)"
    fi
    echo "== static sdm: $(./sdm -v 2>&1 | head -n 1)"
fi
