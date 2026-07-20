#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  printf 'usage: %s JOB_ID CPU BINARY ORIGINAL_INPUT RESULT_DIR CHECKPOINT PREPARER\n' "$0" >&2
  exit 64
fi

job_id=$1
cpu=$2
binary=$3
original_input=$4
result_dir=$5
checkpoint=$6
preparer=$7
root=$(cd "$(dirname "$0")/.." && pwd)
launcher="$root/launch_pinned_cp2k.sh"

[[ -x $binary && -r $original_input && -d $result_dir ]] || exit 66
[[ -s $checkpoint && -x $preparer && -x $launcher ]] || exit 66
if [[ -s $result_dir/cp2k.out ]] && grep -q 'PROGRAM ENDED AT' "$result_dir/cp2k.out"; then
  printf 'refusing to restart an already completed result: %s\n' "$result_dir" >&2
  exit 65
fi

attempt_root="$result_dir/attempts"
stamp=$(date --utc +%Y%m%dT%H%M%SZ)
attempt_dir="$attempt_root/$stamp"
mkdir -p "$attempt_dir"

# The new run may rotate PROJECT-RESTART.kp before the initial guess is read.
# Preserve an immutable source copy under the attempt directory and point the
# restart input there, so input and output checkpoint paths cannot alias.
checkpoint_source="$result_dir/restart-source-${stamp}.kp"
source_checkpoint_sha=$(sha256sum "$checkpoint" | awk '{print $1}')
cp -p --reflink=auto "$checkpoint" "$checkpoint_source"
archived_checkpoint_sha=$(sha256sum "$checkpoint_source" | awk '{print $1}')
if [[ $source_checkpoint_sha != "$archived_checkpoint_sha" ]]; then
  printf 'checkpoint copy failed integrity verification\n' >&2
  exit 74
fi
printf '%s  %s\n' "$archived_checkpoint_sha" "$(basename "$checkpoint_source")" \
  > "$attempt_dir/source-checkpoint.sha256"
printf '%s\n' "$checkpoint_source" > "$attempt_dir/source-checkpoint.path"

for name in cp2k.out exit_status affinity_preexec.txt binary.sha256 input.sha256; do
  if [[ -e $result_dir/$name ]]; then
    mv "$result_dir/$name" "$attempt_dir/$name"
  fi
done

restart_input="$result_dir/restart-${stamp}.inp"
"$preparer" "$original_input" "$checkpoint_source" "$restart_input" \
  > "$result_dir/restart-${stamp}.provenance.stdout.json"

bash "$launcher" "${job_id}-resume-${stamp}" "$cpu" "$binary" \
  "$restart_input" "$result_dir"

restart_mode=
if grep -q 'KPOINT_RESTART| Strict same-mesh restart accepted' \
  "$result_dir/cp2k.out"; then
  restart_mode=strict_same_mesh
elif grep -q 'KPOINT_RESTART| Validated BvK mesh transfer accepted' \
  "$result_dir/cp2k.out"; then
  restart_mode=validated_bvk_transfer
else
  printf 'FAIL\n' > "$result_dir/restart_acceptance"
  printf 'CP2K completed without accepting a validated k-point restart\n' >&2
  exit 86
fi
grep -q 'PROGRAM ENDED AT' "$result_dir/cp2k.out"
printf 'PASS\n' > "$result_dir/restart_acceptance"
printf '%s\n' "$restart_mode" > "$result_dir/restart_acceptance_mode"
