# Third-party submodule patches

Vendored submodules under `third-party/` are pinned to upstream commits, so a
fresh `git clone --recurse-submodules` (or `git submodule update`) checks out
**pristine upstream** and drops any local edits. We keep those edits here as
patch files and replay them at build time.

- **Layout:** `patches/<submodule>/*.patch` → applied inside `third-party/<submodule>`
  (patch paths are relative to that submodule's root).
- **Apply:** `scripts/apply_patches.sh` (idempotent — already-applied patches are
  skipped). Called automatically at the top of `scripts/build_all.sh`, after
  submodule checkout and before configure.
- **Add a new patch:** make the edit in the submodule working tree, then
  `git -C third-party/<sub> diff -- <file> > patches/<sub>/<short-name>.patch`,
  and add a row below.
- **On upstream bump:** if a patch stops applying, `apply_patches.sh` fails loudly
  (`FAIL ...`) — refresh or drop the patch, since upstream may have fixed it.

## Patch list

| Submodule | Patch | Reason | Upstream status |
|-----------|-------|--------|-----------------|
| vtk | `scn-libcxx-wrap_iter-base.patch` | VTK's vendored `vtkscn` calls `std::__wrap_iter<char*>::base()`, which LLVM ≥18 libc++ removed; switched to `operator->()` (returns `std::__to_address(__i_)`). Only triggers on the libc++ build (`-DTABERNA_LIBCXX=ON`). | Not reported upstream; libstdc++ builds are unaffected. |
