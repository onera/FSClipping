# FSClipping

FSClipping is a C++20 library that clips two meshes against each other at a shared boundary. It is built on the FSDM data-management framework.

Given two meshes and one boundary marker on each, it computes the intersection of the two boundary face sets. It then builds a new clipped mesh for each side, so both meshes conform to a common interface. It is used for non-conforming mesh coupling in CFD, such as sliding or chimera interfaces.

## Features

- Sutherland–Hodgman polygon clipping of boundary faces, with BVH-accelerated candidate search.
- Volume reconstruction: cells next to the interface become polyhedral cells that conform to the clipped faces.
- MPI-parallel execution on distributed meshes.
- Integrates with FSDM as a `FSDataManager` operation (`ClippedMesh`), with Python bindings (SWIG).

## Requirements

- CMake ≥ 3.22
- A C++20 compiler
- MPI (OpenMPI ≥ 4.0)
- FSDM ≥ 25.08 (found through `pkg-config`)
- Python ≥ 3.10 with development headers
- SWIG (for the Python bindings)

## Build

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build -j
cmake --install build
```

## Tests

The tests use GoogleTest (bundled in `thirdparty/`) and are on by default (`BUILD_TESTING=ON`).

```bash
cd build
ctest --output-on-failure
```

The sequential tests run on 1 MPI process. The parallel tests run on 3, 5 or 6 processes and are prefixed `np3.`, `np5.` and `np6.`.

## Usage (Python)

The clipping runs on `1 + n1 + n2` MPI processes: one process does the matching, `n1` processes hold mesh 1 and `n2` processes hold mesh 2.

```python
dataManagerOps = (("ClippedMesh", {"MeshKeyOrig1"    : "original1",
                                   "MeshKeyOrig2"    : "original2",
                                   "MeshKeyClipped1" : "clippedMesh1",
                                   "MeshKeyClipped2" : "clippedMesh2",
                                   "ClippedMarker1"  : 1,
                                   "ClippedMarker2"  : 1,
                                   "Tolerance"       : 1e-8}),)
dm.DoOps(dataManagerOps) or FSError.PrintAndExit()
```

| Parameter | Description |
|---|---|
| `MeshKeyOrig1`, `MeshKeyOrig2` | Keys of the input meshes in the data manager |
| `MeshKeyClipped1`, `MeshKeyClipped2` | Keys of the output clipped meshes |
| `ClippedMarker1`, `ClippedMarker2` | Boundary marker that selects the faces to clip on each mesh |
| `Tolerance` | Geometric tolerance of the clipping (a length) |

For a complete example, including mesh loading and the MPI setup, see [example/meshOp/howToComputeClippedMesh.py](example/meshOp/howToComputeClippedMesh.py).

## Project layout

```
include/FSClipping/   Public headers
src/                  Library sources
py/                   Python bindings (SWIG)
test/                 Unit tests and test meshes
example/              Python examples
tools/                pkg-config template
```

## License

FSClipping is distributed under the GNU Lesser General Public License v3.0. See [LICENSE](LICENSE).
