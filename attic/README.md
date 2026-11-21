# attic

Folders that are not part of the build. The root `CMakeLists.txt` does not add anything from
here, and the contents are not maintained or tested.

| Folder | What it was |
|---|---|
| `extractandpreprocessgeobodies/` | extraction tool, not wired into the root build |
| `sparsesteepest/` | sparse variant of the steepest-descent gradient tool, not wired into the root build |
| `test_convolution/` | standalone convolution experiments with their own CMake project |

A folder leaves `attic/` only by being wired into the build: move it under `tools/` or `tests/`
and add its option and `add_subdirectory` line to the root `CMakeLists.txt`.
