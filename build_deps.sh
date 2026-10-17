#!/usr/bin/env bash
# Builds static ISA-L and zstd from their release tarballs into a prefix, for a static sdm:
# Debian/Ubuntu's libisal-dev ships no libisal.a, and a conda env has neither archive. The
# Makefile uses the default prefix (.deps, `make deps`) whenever it holds them;
# build_static.sh builds them into its micromamba env.
#
#   ./build_deps.sh                  # into .deps
#   ./build_deps.sh --prefix DIR     # into DIR (lib/libisal.a, lib/libzstd.a, include/)
#
# Needs a C compiler, make, cmake (3.16 or later) and curl or wget. ISA-L's x86-64 code is nasm
# assembly: a nasm of 2.14.01 or later on the PATH is used, else nasm is built from source first.
# Both libraries are built position-independent (they also link into STATIC=0 builds), without
# -march: they choose their SSE/AVX2/AVX-512 code at run time, so sdm still runs on any x86-64.
# Later runs do nothing unless a version changes; remove the prefix to start over.
#
# Environment:
#   CC                    C compiler (default: cc)
#   ISAL_VERSION          default 2.32.1; set ISAL_SHA256 to verify the download of others
#   ZSTD_VERSION          default 1.5.7; ZSTD_SHA256 likewise
#   NASM_VERSION          default 2.16.03; NASM_SHA256 likewise
#   ISAL_URL, ZSTD_URL, NASM_URL   a URL or local path of the tarball instead of the download
#                         (without network access, download them elsewhere)
set -euo pipefail

die() { echo "build_deps.sh: $*" >&2; exit 1; }

cd "$(dirname "${BASH_SOURCE[0]}")"
prefix=.deps
while [[ $# -gt 0 ]]; do
    case $1 in
        --prefix) [[ $# -ge 2 ]] || die "--prefix needs a directory"; prefix=$2; shift 2 ;;
        --prefix=*) prefix=${1#--prefix=}; shift ;;
        -h|--help) sed -n '2,23p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) die "unknown argument: $1 (see --help)" ;;
    esac
done
mkdir -p "$prefix"
prefix=$(cd "$prefix" && pwd)

isal_version=${ISAL_VERSION:-2.32.1}
zstd_version=${ZSTD_VERSION:-1.5.7}
nasm_version=${NASM_VERSION:-2.16.03}
# Checksums of the default versions (as protal pins them, lib/static-deps.cmake).
[[ $isal_version == 2.32.1 ]] && : "${ISAL_SHA256:=d9f7179ab0e14a3db9b610fac22793854a1435e8423ec9ce07f4cbedc5f92f5e}"
[[ $zstd_version == 1.5.7 ]] && : "${ZSTD_SHA256:=eb33e51f49a15e023950cd7825ca74a4a2b43db8354825ac24fc1b7ee09e6fa3}"
[[ $nasm_version == 2.16.03 ]] && : "${NASM_SHA256:=1412a1c760bbd05db026b6c0d1657affd6631cd0a63cddb6f73cc6d4aa616148}"
# ISA-L publishes no release tarballs; GitHub's tag archive is the release's source.
isal_url=${ISAL_URL:-https://github.com/intel/isa-l/archive/refs/tags/v$isal_version.tar.gz}
zstd_url=${ZSTD_URL:-https://github.com/facebook/zstd/releases/download/v$zstd_version/zstd-$zstd_version.tar.gz}
nasm_url=${NASM_URL:-https://www.nasm.us/pub/nasm/releasebuilds/$nasm_version/nasm-$nasm_version.tar.xz}

cc=${CC:-cc}
command -v "$cc" > /dev/null || die "no C compiler ($cc; set CC)"
command -v cmake > /dev/null || die "cmake is needed (Debian/Ubuntu: apt install cmake; conda: cmake)"
command -v make > /dev/null || die "make is needed"
ncpu=$(nproc 2> /dev/null || sysctl -n hw.ncpu 2> /dev/null || echo 4)

deps="isa-l-$isal_version zstd-$zstd_version"
stamp=$prefix/share/sdm-deps
if [[ $(cat "$stamp" 2> /dev/null) == "$deps" && -e $prefix/lib/libisal.a && -e $prefix/lib/libzstd.a ]]; then
    echo "== static ISA-L $isal_version and zstd $zstd_version are in $prefix"
    exit 0
fi

sha256() {
    if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

fetch() {  # url-or-path file sha256
    if [[ -f $1 ]]; then
        cp "$1" "$2"
    elif command -v curl > /dev/null; then
        curl -fL --retry 3 -o "$2" "$1" || die "cannot download $1"
    elif command -v wget > /dev/null; then
        wget -O "$2" "$1" || die "cannot download $1"
    else
        die "neither curl nor wget: download $1 and pass its path (see --help)"
    fi
    if [[ -n $3 ]]; then
        [[ $(sha256 "$2") == "$3" ]] || die "checksum mismatch for $1"
    else
        echo "warning: $1 is not verified (no checksum given)" >&2
    fi
}

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cmake_args=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="$cc" -DCMAKE_INSTALL_PREFIX="$prefix"
            -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_INSTALL_INCLUDEDIR=include
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON)

echo "== building static zstd $zstd_version"
fetch "$zstd_url" "$work/zstd.tar.gz" "${ZSTD_SHA256:-}"
mkdir "$work/zstd" && tar -xzf "$work/zstd.tar.gz" -C "$work/zstd" --strip-components=1
cmake -S "$work/zstd/build/cmake" -B "$work/zstd-build" "${cmake_args[@]}" \
    -DZSTD_BUILD_STATIC=ON -DZSTD_BUILD_SHARED=OFF -DZSTD_BUILD_PROGRAMS=OFF \
    -DZSTD_BUILD_TESTS=OFF -DZSTD_BUILD_CONTRIB=OFF -DZSTD_MULTITHREAD_SUPPORT=ON
cmake --build "$work/zstd-build" -j "$ncpu"
cmake --install "$work/zstd-build"

# nasm for ISA-L on x86-64: the PATH's if it is 2.14.01 or later (ISA-L's minimum), else built here.
isal_nasm=()
if [[ $(uname -m) == x86_64 || $(uname -m) == amd64 ]]; then
    nasm=$(command -v nasm || true)
    have=""
    [[ -n $nasm ]] && have=$("$nasm" -v 2> /dev/null | sed -n 's/^NASM version \([0-9][0-9.]*\).*/\1/p')
    if [[ -n $have ]] && [[ $(printf '%s\n%s\n' 2.14.01 "$have" | sort -V | head -n 1) == 2.14.01 ]]; then
        echo "== nasm for ISA-L: $nasm ($have)"
    else
        echo "== building nasm $nasm_version for ISA-L (none of 2.14.01 or later on the PATH)"
        fetch "$nasm_url" "$work/nasm.tar.xz" "${NASM_SHA256:-}"
        mkdir "$work/nasm" && tar -xJf "$work/nasm.tar.xz" -C "$work/nasm" --strip-components=1
        (cd "$work/nasm" && ./configure CC="$cc" > configure.log && make -j "$ncpu" nasm > make.log) ||
            die "building nasm failed (logs in $work/nasm)"
        nasm=$work/nasm/nasm
    fi
    isal_nasm=(-DCMAKE_ASM_NASM_COMPILER="$nasm")
fi

echo "== building static ISA-L $isal_version"
fetch "$isal_url" "$work/isa-l.tar.gz" "${ISAL_SHA256:-}"
mkdir "$work/isa-l" && tar -xzf "$work/isa-l.tar.gz" -C "$work/isa-l" --strip-components=1
# Without its tests (which need zlib), programs and zlib shim.
cmake -S "$work/isa-l" -B "$work/isa-l-build" "${cmake_args[@]}" "${isal_nasm[@]}" \
    -DBUILD_SHARED_LIBS=OFF -DISAL_BUILD_TESTS=OFF -DISAL_BUILD_PERF_TESTS=OFF \
    -DISAL_BUILD_ISAL_SHIM=OFF -DISAL_BUILD_IGZIP_CLI=OFF
cmake --build "$work/isa-l-build" -j "$ncpu"
cmake --install "$work/isa-l-build"

[[ -e $prefix/lib/libisal.a && -e $prefix/include/isa-l/igzip_lib.h ]] || die "ISA-L did not install into $prefix"
[[ -e $prefix/lib/libzstd.a && -e $prefix/include/zstd.h ]] || die "zstd did not install into $prefix"
mkdir -p "$(dirname "$stamp")"
echo "$deps" > "$stamp"
echo "== static ISA-L $isal_version and zstd $zstd_version installed in $prefix"
