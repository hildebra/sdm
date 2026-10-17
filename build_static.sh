#!/usr/bin/env bash
# Builds a fully static sdm with HTSlib (SAM/BAM/CRAM input), ISA-L and zstd from a micromamba env.
#
#   ./build_static.sh -j8           # build ./sdm
#   ./build_static.sh -j8 test      # build and run all test suites with the same toolchain
#
# Arguments go to make. The first run creates the env (conda-forge: GCC with a static glibc
# sysroot, make, cmake, pkg-config, curl, nasm, and zlib/bzip2/xz with static archives, which
# HTSlib needs) and builds HTSlib from source into it: bioconda's htslib needs libcurl/OpenSSL,
# which conda does not provide as static archives. HTSlib is built without libdeflate (optional for
# it); sdm itself uses neither zlib nor libdeflate.
# ISA-L (all gzip input, BGZF output) and zstd are built as static libraries by build_deps.sh into
# <env>/sdm-deps, apart from the env's own zstd (a dependency of its cmake).
# HTSlib is built without remote-URL (libcurl, S3, GCS) and plugin support. Later runs only
# rebuild sdm. Start over with: micromamba env remove -y -n sdm-static
#
# Environment:
#   SDM_ENV             env name (default sdm-static); an existing env gets missing packages
#   SDM_MICROMAMBA      micromamba executable (default: $MAMBA_EXE, else micromamba on PATH)
#   HTSLIB_VERSION      default 1.24; set HTSLIB_SHA256 to verify the download of other versions
#   ISAL_VERSION, ZSTD_VERSION and their *_SHA256 and *_URL: see build_deps.sh
set -euo pipefail

die() { echo "build_static.sh: $*" >&2; exit 1; }

env_name=${SDM_ENV:-sdm-static}
htslib_version=${HTSLIB_VERSION:-1.24}
# The checksum of the default version, as in the bioconda htslib recipe.
if [[ $htslib_version == 1.24 ]]; then
    : "${HTSLIB_SHA256:=28a8de191381c7a97a35675ceac76fa1ea95e7b678d6a2e9d600a7874e4077de}"
fi
packages=(gcc_linux-64 gxx_linux-64 sysroot_linux-64 make cmake pkg-config curl nasm
          zlib bzip2 liblzma-devel liblzma-static)
target=x86_64-conda-linux-gnu
required=(bin/make bin/cmake bin/pkg-config bin/curl bin/bzip2 bin/nasm
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

deps="htslib-$htslib_version without-libdeflate"
stamp=$prefix/share/sdm-static-deps
if [[ $(cat "$stamp" 2> /dev/null) != "$deps" ]]; then
    ncpu=$(nproc 2> /dev/null || echo 4)
    work=$(mktemp -d)
    trap 'rm -rf "$work"' EXIT
    (
        cd "$work"
        # The libdeflate earlier versions of this script built into the env.
        rm -rf "$prefix"/lib/libdeflate.* "$prefix"/include/libdeflate.h "$prefix"/lib/pkgconfig/libdeflate.pc \
            "$prefix"/lib/cmake/libdeflate

        echo "== building static HTSlib $htslib_version"
        fetch "https://github.com/samtools/htslib/releases/download/$htslib_version/htslib-$htslib_version.tar.bz2" \
            htslib.tar.bz2 "${HTSLIB_SHA256:-}"
        tar -xjf htslib.tar.bz2
        cd "htslib-$htslib_version"
        ./configure --prefix="$prefix" --without-libdeflate --disable-libcurl --disable-gcs \
            --disable-s3 --disable-plugins --disable-ref-cache
        make -j "$ncpu" AR="$AR" RANLIB="$RANLIB"
        make install AR="$AR" RANLIB="$RANLIB"
    )
    mkdir -p "$(dirname "$stamp")"
    echo "$deps" > "$stamp"
fi

# Static ISA-L and zstd in a prefix of their own (build_deps.sh does nothing once they are there).
deps_prefix=$prefix/sdm-deps
bash ./build_deps.sh --prefix "$deps_prefix"

hts_cflags=$(pkg-config --cflags htslib) || die "pkg-config cannot find htslib in $prefix"
hts_libs=$(pkg-config --static --libs htslib) || die "pkg-config cannot resolve htslib's static libraries"
echo "== building sdm (HTS_LIBS: $hts_libs)"
# Objects go to .build/conda-<the Makefile's configuration name>, apart from builds with other
# compilers in this checkout. The single quotes leave the $(...) to make. ISA-L's and zstd's
# headers come before the env's, which may hold another zstd.h (CPPFLAGS from the environment, so
# that the Makefile still adds its own).
CPPFLAGS="-I$deps_prefix/include $CPPFLAGS" \
make STATIC=1 WITH_HTS=1 CC="$CC" CXX="$CXX" CXXFLAGS=-O3 \
    'BUILD_DIR=.build/conda-$(if $(filter 1,$(WITH_HTS)),hts,portable)-$(STATIC)' \
    HTS_CFLAGS="$hts_cflags" HTS_LIBS="$hts_libs" \
    ISAL_CFLAGS= ISAL_LIBS="$deps_prefix/lib/libisal.a" \
    ZSTD_CFLAGS= ZSTD_LIBS="$deps_prefix/lib/libzstd.a -pthread" "$@"

if [[ -e sdm ]]; then
    if [[ $("$prefix/bin/$target-readelf" -d sdm) == *NEEDED* ]]; then
        die "sdm is dynamically linked (make STATIC=0 given?)"
    fi
    echo "== static sdm: $(./sdm -v 2>&1 | head -n 1)"
fi
