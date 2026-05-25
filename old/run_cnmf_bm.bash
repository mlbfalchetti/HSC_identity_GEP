
NAME="sk_2022_adult_bm_01_hpc_old"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_01_hpc_old.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_01_hpc_old"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done



NAME="sk_2022_adult_bm_01_hpc_young"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_01_hpc_young.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_01_hpc_young"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done



NAME="sk_2022_adult_bm_01_hsc_old"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_01_hsc_old.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_01_hsc_old"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done



NAME="sk_2022_adult_bm_01_hsc_young"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_01_hsc_young.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_01_hsc_young"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done



NAME="sk_2022_adult_bm_02_hpc_old"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_02_hpc_old.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_02_hpc_old"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done



NAME="sk_2022_adult_bm_02_hpc_young"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_02_hpc_young.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_02_hpc_young"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done



NAME="sk_2022_adult_bm_02_hsc_young"

COUNTS="/home/mlbfalchetti/cnmf_inputs_1.5/sk_2022_adult_bm_02_hsc_young.txt"
OUTDIR="/home/mlbfalchetti/cnmf_inputs_1.5/cnmf_runs/sk_2022_adult_bm_02_hsc_young"

mkdir -p "$OUTDIR"

for k in 5 6 7 8 9 10
do
  cnmf prepare \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    --counts "$COUNTS" \
    -k "$k" \
    --n-iter 100 \
    --numgenes 5000

  cnmf factorize \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf combine \
    --name "$NAME" \
    --output-dir "$OUTDIR"

  cnmf consensus \
    --name "$NAME" \
    --output-dir "$OUTDIR" \
    -k "$k"
done
