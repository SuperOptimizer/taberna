#!/usr/bin/env bash
# Re-apply taberna's local patches to vendored git submodules.
#
# Submodules (third-party/*) are pinned to upstream commits, so a fresh
# `git clone --recurse-submodules` or `git submodule update` checks out
# pristine upstream and DROPS any working-tree edits. We keep those edits as
# patch files under patches/<submodule>/*.patch and replay them here.
#
# Layout:  patches/<name>/*.patch  ->  applied inside third-party/<name>
#          (patch paths are relative to that submodule's root)
#
# Idempotent: a patch that is already applied (reverse-applies cleanly) is
# skipped, so this is safe to run on every build. Run after submodule checkout
# and before configuring. Called automatically by scripts/build_all.sh.
set -euo pipefail
ROOT=/home/forrest/taberna
cd "$ROOT"

shopt -s nullglob
any=0
for dir in patches/*/; do
  name="$(basename "$dir")"
  sub="third-party/$name"
  if [ ! -d "$sub" ]; then
    echo "skip $name: $sub not checked out"; continue
  fi
  for p in "$dir"*.patch; do
    any=1
    rel="${p#"$ROOT/"}"
    if git -C "$sub" apply --reverse --check "$ROOT/$p" 2>/dev/null; then
      echo "ok   $name <- $(basename "$p") (already applied)"
    elif git -C "$sub" apply --check "$ROOT/$p" 2>/dev/null; then
      git -C "$sub" apply "$ROOT/$p"
      echo "APPL $name <- $(basename "$p")"
    else
      echo "FAIL $name <- $(basename "$p"): does not apply cleanly (upstream moved?)" >&2
      exit 1
    fi
  done
done
[ "$any" = 1 ] || echo "no patches found under patches/"
