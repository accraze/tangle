# Releasing

Releases are built and published by CI (`.github/workflows/ci.yml`) when a `v*` tag is pushed.
Nothing is built or uploaded from a local machine.

## Cutting a release

1. Bump `version` in `Cargo.toml` (semver: patch for fixes/docs/licensing, minor for new features
   or sound changes, major once presets/automation could break).
2. In `CHANGELOG.md`, move the `[Unreleased]` entries into a new `## [X.Y.Z] - YYYY-MM-DD` section
   (user-facing changes only), and update the compare links at the bottom.
3. Refresh the lockfile and commit:
   ```sh
   cargo metadata --format-version 1 -q > /dev/null   # updates Cargo.lock's version entry
   git commit -am "Release X.Y.Z"
   git push origin main
   ```
4. Wait for CI on `main` to go green. A tag on a broken commit still runs the checks, but you'd
   rather find out before the tag exists.
5. Tag with an **annotated** tag. Its message becomes the release notes verbatim:
   ```sh
   git tag -a vX.Y.Z -F notes.md     # or -m "..." for short notes
   git push origin vX.Y.Z
   ```
6. Watch the run (`gh run list`, `gh run watch <id>`), then check the release:
   ```sh
   gh release view vX.Y.Z
   gh release download vX.Y.Z --dir /tmp/check && unzip -l /tmp/check/*.zip
   ```
   The zip should contain `Tangle.clap`, `Tangle.vst3/` and `LICENSE`.

### Release notes

First line is the title (`Tangle vX.Y.Z`), then a one-line summary, then the version's
`CHANGELOG.md` entries, then an Install section. See `git tag -n99 v0.1.0` for the shape.

## What CI does on a tag

Same job as every push, plus one step at the end:

1. Format check, clippy, `cargo test --workspace --lib` (plugin and `core` tests)
2. `cargo xtask bundle <crate> --release`
3. clap-validator 0.4.1 against the CLAP bundle (four state tests excluded, see the comment in
   the workflow for why)
4. Upload the bundles as a workflow artifact
5. **Tag only:** zip the bundles plus `LICENSE` as `<crate>-vX.Y.Z-linux.zip` and run
   `gh release create --notes-from-tag`

If any step before 5 fails, no release is created.

Gotchas baked into the workflow:
- `actions/checkout` fetches the tag as a lightweight ref, so the release step re-fetches the
  annotated tag object before `--notes-from-tag` can read its message.
- The job needs `permissions: contents: write` to create releases.
- GPL: the VST3 build links nih-plug's GPLv3 bindings, so release zips must ship `LICENSE`.

## Fixing a bad release

Prefer a new patch release over rewriting history. Don't move or re-push an existing tag.

- Cut `vX.Y.Z+1` with the fix.
- Prepend a note to the old release pointing at the new one:
  ```sh
  { echo '> **Superseded by [vX.Y.Z+1](...)**, ...'; echo; gh release view vX.Y.Z --json body -q .body; } > notes.md
  gh release edit vX.Y.Z --notes-file notes.md
  ```
- Only delete a release if it must not be downloaded at all (e.g. a broken binary). Deleting
  breaks any links already shared.

## Reusing this for a new plugin

The workflow is generic apart from names. When copying it to another nih-plug repo:

| What | Where | Change |
|---|---|---|
| Crate name | `cargo xtask bundle tangle` in `ci.yml` | your crate's `[package] name` |
| Bundle names | `Tangle.clap` / `Tangle.vst3` in `ci.yml` (validate, upload, zip steps) | the `name` in your `bundler.toml` |
| Artifact and zip names | `tangle-linux`, `tangle-$GITHUB_REF_NAME-linux.zip` | your crate name |
| Release title | `--title "Tangle ..."` | your plugin name |
| Plugin IDs | `CLAP_ID`, `VST3_CLASS_ID` in `src/lib.rs` | must be unique per plugin, and never change after release |
| License | `LICENSE`, `license` in `Cargo.toml` | GPL-3.0 if you export VST3; anything if CLAP only |
| Validator excludes | `-x ...` flags in `ci.yml` | re-check: drop excludes once nih-plug/clap-validator agree again |

Also copy `rust-toolchain.toml` and the apt package list (needed for any nih_plug_egui GUI).
