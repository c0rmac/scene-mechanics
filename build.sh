#!/bin/sh
# Build the demo against MuJoCo. MUJOCO_DIR is a directory holding include/mujoco/mujoco.h and the libmujoco shared library;
# by default the one shipped with the Python package (pip install mujoco).
set -e
cd "$(dirname "$0")"
: "${MUJOCO_DIR:=$(python3 -c 'import mujoco, os; print(os.path.dirname(mujoco.__file__))')}"
LIB=$(ls "$MUJOCO_DIR"/libmujoco*.dylib "$MUJOCO_DIR"/libmujoco.so* 2>/dev/null | head -1)
[ -n "$LIB" ] || { echo "no libmujoco in $MUJOCO_DIR"; exit 1; }
mkdir -p build
c++ -O2 -std=c++20 -Wno-dangling-else -Iinclude -I"$MUJOCO_DIR/include" examples/demo.cpp -o build/demo "$LIB" -Wl,-rpath,"$MUJOCO_DIR" -lpthread
if [ "$(uname)" = "Darwin" ]; then   # the pip package's library names itself as a framework; point the binary at the file itself
  NAME=$(otool -L build/demo | awk '/libmujoco/ {print $1; exit}')
  case "$NAME" in *mujoco.framework*) install_name_tool -change "$NAME" "@rpath/$(basename "$LIB")" build/demo; codesign -f -s - build/demo >/dev/null 2>&1 || true;; esac
fi
echo "built build/demo (MuJoCo from $MUJOCO_DIR)"
