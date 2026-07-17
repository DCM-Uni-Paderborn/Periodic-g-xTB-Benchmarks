#!/usr/bin/env bash
set -u

if [[ $# -ne 1 ]]; then
  echo "usage: $0 ORDERED_PE_RESERVATION" >&2
  echo "example: $0 96,97,... (exactly 49 distinct available CPUs)" >&2
  exit 64
fi
reservation=$1
IFS=, read -r -a reserved_cpus <<< "$reservation"
required_cpus=49
if [[ ${#reserved_cpus[@]} -ne $required_cpus ]]; then
  echo "ordered PE reservation needs exactly $required_cpus CPUs" >&2
  exit 65
fi
declare -A seen_cpus=()
for cpu in "${reserved_cpus[@]}"; do
  if [[ ! $cpu =~ ^[0-9]+$ ]] || [[ -n ${seen_cpus[$cpu]+x} ]]; then
    echo "ordered PE reservation contains an invalid or duplicate CPU: $cpu" >&2
    exit 66
  fi
  seen_cpus[$cpu]=1
done

root=/home/kuehne88/work/cp2k-gxtb-streamed-reverse-consumer-20260716
envdir=/home/kuehne88/work/codex-gxtb-pbc-20260714T1038Z-18d37c-1449feb/env
runner="$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/run_with_rss.py"
matrix="$root/evidence/linux_matrix"

mkdir -p "$matrix"
exec 9> "$matrix/.run_linux_matrix.lock"
if ! flock -n 9; then
  echo "another linux-matrix writer already owns $matrix" >&2
  exit 75
fi
cpu_lock_root="/tmp/periodic-gxtb-cpu-reservations-$(id -u)"
mkdir -p "$cpu_lock_root"
reservation_lock_fds=()
for cpu in "${reserved_cpus[@]}"; do
  exec {lock_fd}>> "$cpu_lock_root/cpu-${cpu}.lock"
  if ! flock -n "$lock_fd"; then
    echo "logical CPU $cpu is already reserved by another production launcher" >&2
    exit 76
  fi
  printf '{"cpu":%s,"pid":%s,"source":"run_linux_matrix.sh"}\n' \
    "$cpu" "$$" >&"$lock_fd"
  reservation_lock_fds+=("$lock_fd")
done

export PATH="$envdir/bin:$PATH"
export LD_LIBRARY_PATH="$root/build/src:$envdir/lib:${LD_LIBRARY_PATH:-}"
export CP2K_DATA_DIR="$root/source/data"
export CP2K_GXTB_EXCHANGE_GRADIENT_MODE=QUALIFY
export CP2K_GXTB_EXCHANGE_IMAGE_BATCH_SIZE=2
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export LSAN_OPTIONS=detect_leaks=0
for inherited_mca_key in ${!OMPI_MCA_@} ${!PRTE_MCA_@}; do
  unset "$inherited_mca_key"
done
ulimit -c 0

names=(
  k290_rks_3d_fd
  shifted_spglib_rks_3d
  tr_rks_3d_fd
  fullmesh_rks_3d_fd
  fullmesh_uks_3d_fd
  rks_1d_fd
  spglib_rks_2d_fd
)
inputs=(
  "$root/source/tests/xTB/regtest-tblite-gxtb/CH4_gxtb_kp_k290_force_stress.inp"
  "$root/source/tests/xTB/regtest-tblite-gxtb-spglib/Si_prim_gxtb_kp_shifted_spglib.inp"
  "$root/source/tests/xTB/regtest-tblite-gxtb/H2_gxtb_kp_311_tr_force_stress.inp"
  "$root/source/tests/xTB/regtest-tblite-gxtb/H2_gxtb_kp_311_force_stress.inp"
  "$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/inputs/O2_gxtb_uks_kp_311_force_stress.inp"
  "$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/inputs/Ar2_chain_gxtb_x_kp_debug.inp"
  "$root/source/GXTB_STREAMED_REVERSE_EVIDENCE/inputs/Ar_layer_gxtb_yz_kp_debug.inp"
)

run_one() {
  local name=$1 input=$2 ranks=$3 pe_list=$4
  local run="$matrix/${name}_p${ranks}"
  mkdir -p "$run"
  cd "$run" || exit 90
  rm -f PASS FAIL RUNNING run.out launcher.log rusage.json returncode.txt \
    selected_evidence.txt SHA256SUMS
  printf '%s\n' "$input" > input.path
  printf '%s\n' "$pe_list" > ordered_pe_list.txt
  printf '%s\n' "$ranks" > mpi_ranks.txt
  sha256sum "$input" > input.sha256
  date --iso-8601=seconds > started.txt
  : > RUNNING

  set +e
  python3 "$runner" --mpi-ranks "$ranks" --ordered-pe-list "$pe_list" \
    --cp2k "$root/build/bin/cp2k.pdbg" --launcher-log launcher.log \
    rusage.json timeout 1200 \
    "$envdir/bin/mpiexec" --map-by "pe-list=${pe_list}:ordered" --bind-to core \
    --report-bindings -n "$ranks" \
    "$root/build/bin/cp2k.pdbg" -i "$input" -o run.out > launcher.log 2>&1
  local rc=$?
  set -e
  printf '%s\n' "$rc" > returncode.txt
  date --iso-8601=seconds > finished.txt
  if [[ -f run.out ]]; then
    grep -E "GXTB-QUALIFICATION_ONLY STREAMED-REVERSE|GXTB STREAMED-REVERSE|DEBUG.*Sum of differences|PROGRAM ENDED" \
      run.out > selected_evidence.txt || true
  fi
  rm -f RUNNING
  if [[ $rc -eq 0 ]] && grep -q "PROGRAM ENDED" run.out 2> /dev/null &&
    grep -q "GXTB-QUALIFICATION_ONLY STREAMED-REVERSE" run.out 2> /dev/null; then
    : > PASS
  else
    : > FAIL
    if [[ $rc -eq 0 ]]; then
      rc=98
      printf '%s\n' "$rc" > returncode.txt
    fi
  fi
  sha256sum input.path input.sha256 ordered_pe_list.txt mpi_ranks.txt launcher.log rusage.json \
    returncode.txt started.txt finished.txt > SHA256SUMS 2> /dev/null || true
  [[ -f run.out ]] && sha256sum run.out >> SHA256SUMS
  [[ -f PASS ]] && sha256sum PASS >> SHA256SUMS
  [[ -f FAIL ]] && sha256sum FAIL >> SHA256SUMS
  return "$rc"
}

next_cpu=0
pids=()
for index in "${!names[@]}"; do
  for ranks in 1 2 4; do
    pe_list=""
    for ((offset = 0; offset < ranks; offset++)); do
      if [[ -n $pe_list ]]; then
        pe_list+=,
      fi
      pe_list+=${reserved_cpus[$((next_cpu + offset))]}
    done
    next_cpu=$((next_cpu + ranks))
    run_one "${names[$index]}" "${inputs[$index]}" "$ranks" "$pe_list" &
    pids+=("$!")
  done
done

overall=0
for pid in "${pids[@]}"; do
  wait "$pid" || overall=1
done

find "$matrix" -mindepth 2 -maxdepth 2 -type f \( -name PASS -o -name FAIL \) -print | sort \
  > "$matrix/status.txt"
date --iso-8601=seconds > "$matrix/finished.txt"
exit "$overall"
