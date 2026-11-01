#!/usr/bin/env bash
set -euo pipefail

# Bump the shared workspace version and publish the salmon-* crates to crates.io
# in dependency order. Modeled on piscem's bump_and_publish.sh, adapted for a
# multi-crate workspace whose members all inherit `version.workspace = true`.
#
# This script publishes ONLY the salmon-* workspace crates. The external
# COMBINE-lab deps (cf1-rs, piscem-rs, ksw2rs) are published interactively from
# their own repositories first; salmon depends on their published crates.io
# versions (see the release plan, R2/R3).

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./scripts/bump_and_publish.sh <version> [--publish] [--dry-run]

Bumps [workspace.package] version (all salmon-* crates inherit it), commits,
tags, pushes, then optionally publishes the salmon-* crates to crates.io in
dependency order.

Options:
  --publish  Publish to crates.io after bumping, committing, tagging, pushing
  --dry-run  Do not modify this checkout, commit, tag, push, or publish. Applies
             the bump in a scratch git worktree and runs
             `cargo publish --workspace --dry-run` there, so every crate is
             packaged and verified at the new version without uploading.
  -h, --help Show this help message

Publish order (dependency-topological):
  salmon-core -> salmon-eqclass -> salmon-model -> salmon-index -> salmon-infer
  -> salmon-map -> salmon-rad -> salmon-align -> salmon-quant -> salmon-cli
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish) PUBLISH=true ;;
        --dry-run) DRY_RUN=true ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1" ;;
        *)
            [[ -z "$VERSION" ]] || die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || { usage; exit 1; }

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$WORKSPACE_ROOT"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
TAG="v${VERSION}"

# Dependency-topological publish order. crates.io will not accept a crate until
# the versions of its path+version deps are already in the index; cargo publish
# blocks (wait-for-publish) until each appears, so a sequential loop is enough.
CRATES=(
    salmon-core
    salmon-eqclass
    salmon-model
    salmon-index
    salmon-infer
    salmon-map
    salmon-rad
    salmon-align
    salmon-quant
    salmon-cli
)

MANIFEST_BACKUP=""
LOCKFILE_BACKUP=""
MANIFEST_UPDATED=false
COMMIT_CREATED=false

cleanup() {
    local status=$?
    if [[ "$status" -ne 0 && "$DRY_RUN" == false && "$MANIFEST_UPDATED" == true && "$COMMIT_CREATED" == false ]]; then
        [[ -n "$MANIFEST_BACKUP" && -f "$MANIFEST_BACKUP" ]] && cp "$MANIFEST_BACKUP" "$ROOT_CARGO"
        [[ -n "$LOCKFILE_BACKUP" && -f "$LOCKFILE_BACKUP" ]] && cp "$LOCKFILE_BACKUP" "$LOCKFILE"
        echo "restored $ROOT_CARGO and $LOCKFILE after failure" >&2
    fi
    [[ -n "$MANIFEST_BACKUP" && -f "$MANIFEST_BACKUP" ]] && rm -f "$MANIFEST_BACKUP"
    [[ -n "$LOCKFILE_BACKUP" && -f "$LOCKFILE_BACKUP" ]] && rm -f "$LOCKFILE_BACKUP"
    return "$status"
}
trap cleanup EXIT

[[ -f "$ROOT_CARGO" ]] || die "not found: $ROOT_CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE (run 'cargo build' once to generate it)"

# Current workspace version lives under [workspace.package].
CURRENT_VERSION="$(sed -n '/^\[workspace.package\]/,/^\[/{s/^version = "\(.*\)"/\1/p}' "$ROOT_CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine [workspace.package] version from $ROOT_CARGO"

[[ "$CURRENT_VERSION" != "$VERSION" ]] || die "workspace version is already $VERSION"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi
if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi
git remote get-url origin >/dev/null 2>&1 || die "git remote 'origin' is not configured"

echo "Current workspace version : $CURRENT_VERSION"
echo "New workspace version     : $VERSION"
echo "Tag                       : $TAG"
echo "Publish                   : $([[ "$PUBLISH" == true ]] && echo yes || echo no)"
echo "Dry-run                   : $([[ "$DRY_RUN" == true ]] && echo yes || echo no)"
echo "Crates (in order)         : ${CRATES[*]}"
echo

echo "Preflight: cargo check"
cargo check -q

echo "Updating [workspace.package] version: $CURRENT_VERSION -> $VERSION"

bump_manifest() {
    # Bump the version line within the [workspace.package] block, and the
    # version="X" requirements on the internal salmon-* path deps in
    # [workspace.dependencies] (so the published deps resolve to the new version).
    local manifest="$1"
    sed -i.bak "/^\[workspace.package\]/,/^\[/{s/^version = \".*\"/version = \"${VERSION}\"/}" "$manifest"
    sed -i.bak "/^\[workspace.dependencies\]/,/^\[profile/{s/^\(salmon-[a-z]* = { path = \"crates\/salmon-[a-z]*\", version = \"\).*\(\" }\)/\1${VERSION}\2/}" "$manifest"
    rm -f "${manifest}.bak"
}

if [[ "$DRY_RUN" == true ]]; then
    # Validate exactly what --publish would upload: apply the same bump in a
    # throwaway worktree of HEAD and run a workspace publish dry-run there.
    # `cargo publish --workspace` (cargo >= 1.90) verifies every crate in
    # dependency order, supplying the not-yet-published salmon-* versions from
    # a temporary local registry, so the dependents of a bumped crate are
    # checked too. Nothing in this checkout is modified.
    command -v git >/dev/null || die "git is required"
    SIM_DIR="$(mktemp -d "${TMPDIR:-/tmp}/salmon-release-sim.XXXXXX")"
    trap 'git worktree remove --force "$SIM_DIR" >/dev/null 2>&1 || true; rm -rf "$SIM_DIR"; cleanup' EXIT
    git worktree add --quiet --detach "$SIM_DIR" HEAD
    echo
    echo "Dry-run: validating the v${VERSION} packages in a scratch worktree ($SIM_DIR)"
    (
        cd "$SIM_DIR"
        bump_manifest "$ROOT_CARGO"
        cargo check -q
        cargo publish --workspace --dry-run --allow-dirty
    )
    echo
    echo "Dry-run: would rewrite [workspace.package] version and salmon-* dep versions in $ROOT_CARGO"
    print_cmd git add "$ROOT_CARGO" "$LOCKFILE"
    print_cmd git commit -m "chore(release): bump salmon workspace to v${VERSION}"
    print_cmd git tag -a "$TAG" -m "Release ${VERSION}"
    print_cmd git push origin HEAD
    print_cmd git push origin "$TAG"
    if [[ "$PUBLISH" == true ]]; then
        for crate in "${CRATES[@]}"; do
            print_cmd cargo publish -p "$crate"
        done
    fi
    echo
    echo "Dry-run complete: all crates packaged and verified at v${VERSION}; no commit, tag, push, or publish performed"
    exit 0
fi

MANIFEST_BACKUP="$(mktemp "${TMPDIR:-/tmp}/salmon-Cargo.toml.XXXXXX")"
LOCKFILE_BACKUP="$(mktemp "${TMPDIR:-/tmp}/salmon-Cargo.lock.XXXXXX")"
cp "$ROOT_CARGO" "$MANIFEST_BACKUP"
cp "$LOCKFILE" "$LOCKFILE_BACKUP"
bump_manifest "$ROOT_CARGO"
MANIFEST_UPDATED=true

# Refresh Cargo.lock for the new versions.
cargo check -q

UPDATED_VERSION="$(sed -n '/^\[workspace.package\]/,/^\[/{s/^version = "\(.*\)"/\1/p}' "$ROOT_CARGO" | head -1)"
[[ "$UPDATED_VERSION" == "$VERSION" ]] || die "workspace version update failed"

run git add "$ROOT_CARGO" "$LOCKFILE"
run git commit -m "chore(release): bump salmon workspace to v${VERSION}"
COMMIT_CREATED=true

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

if [[ "$PUBLISH" == true ]]; then
    for crate in "${CRATES[@]}"; do
        echo "Publishing $crate ..."
        run cargo publish -p "$crate"
    done
    echo "All salmon-* crates published for v${VERSION}"
else
    echo "Skipping crates.io publish; re-run with --publish to publish v${VERSION}"
fi

echo
echo "Release bump complete for v${VERSION}"
