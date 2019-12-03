function write_base_options {
  ILU=$1
cat > tenstream.options.base <<EOF
-solar_dir_ksp_converged_reason
#-solar_dir_ksp_view
-solar_diff_ksp_converged_reason
#-solar_diff_ksp_view
-thermal_diff_ksp_converged_reason
#-thermal_diff_ksp_view

#-twostr_ratio 1.5
#-rrtm_bands 1,5

-force_theta
-theta 60
-force_phi
-phi 180

-solar_diff_ksp_type bcgs
-solar_diff_pc_type bjacobi
-solar_diff_sub_pc_type ilu
-solar_diff_sub_pc_factor_levels $ILU

-thermal_diff_ksp_type bcgs
-thermal_diff_pc_type bjacobi
-thermal_diff_sub_pc_type ilu
-thermal_diff_sub_pc_factor_levels $ILU

-solar_dir_ksp_type bcgs
-solar_dir_pc_type bjacobi
-solar_dir_sub_pc_type ilu
-solar_dir_sub_pc_factor_levels 1
-solar_dir_pc_factor_nonzeros_along_diagonal
-solar_dir_pc_factor_reuse_ordering
EOF
}

function write_mg_options {
  write_base_options 1
  cp tenstream.options.base tenstream.options.mg
cat >> tenstream.options.mg <<EOF
-max_solution_time $3

-da_refine_x $2
-da_refine_y 2
-da_refine_z 2

-solar_diff_pc_type mg
-solar_diff_pc_mg_galerkin
-solar_diff_pc_mg_levels $1
-solar_diff_mg_levels_ksp_type richardson
-solar_diff_mg_levels_pc_type sor
-solar_diff_mg_levels_ksp_max_it 1
-solar_diff_mg_coarse_ksp_type preonly
-solar_diff_mg_coarse_ksp_max_it 1
-solar_diff_mg_coarse_pc_type sor
-solar_diff_ksp_reuse_preconditioner true

-thermal_diff_pc_type mg
-thermal_diff_pc_mg_galerkin
-thermal_diff_pc_mg_levels $1
-thermal_diff_mg_levels_ksp_type richardson
-thermal_diff_mg_levels_pc_type sor
-thermal_diff_mg_levels_ksp_max_it 1
-thermal_diff_mg_coarse_ksp_type preonly
-thermal_diff_mg_coarse_ksp_max_it 1
-thermal_diff_mg_coarse_pc_type sor
-thermal_diff_ksp_reuse_preconditioner true
EOF
}

function write_job_file {
  BINARY=$1
  WDIR=$2
  OPT=$3
  SBATCH=$4
if [[ $PETSC_ARCH = *"amd" ]]; then
  CONSTRAINT="-C amd"
fi
  cat > slurm.job <<EOF
#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem=250G
#SBATCH -n 64
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH $SBATCH
#SBATCH $CONSTRAINT

. $HOME/.profile
lpl $PETSC_ARCH

cd $WDIR
srun --mpi=pmix $BINARY namoptions.000 -log_view $OPT
EOF
}

function run_ex {
  BINARY=$1
  IDENT=$2
  IRAD=$3
  OPTFILE=$4
  OPTFLAG=$5
  SBATCH=$6

  RUNDIR=run.fix_angles.$IDENT
  if [ ! -d $RUNDIR ]; then
    echo "Generating Job for $IDENT in [$RUNDIR]"
    mkdir $RUNDIR
    cd $RUNDIR
    cp ../afglus_100m.dat ../*.000 ../backrad.inp.000.nc ../rrtmg_*.nc .

    [ -e ../$OPTFILE ] && (cp ../$OPTFILE tenstream.options) || (echo "Warning: tenstream options file: [$OPTFILE] does not exist?")

    sed -i "s/iradiation.*/iradiation = $IRAD/" namoptions.000

    write_job_file $BINARY $(pwd) "$OPTFLAG" "$SBATCH"
    sbatch slurm.job
    cd ..
  fi
}

write_base_options 1

. $WHOME/.profile
lpl prod_single_gcc_amd

BIN=$WORK/lib/dales_menno/build_$PETSC_ARCH/src/dales4

run_ex $BIN dales_rrtmg 4 "no_opt"
run_ex $BIN rrtmg 5 tenstream.options.base "-rrtmg_only"
run_ex $BIN twostr 5 tenstream.options.base "-twostr_only"
run_ex $BIN twostr_schwarz 5 tenstream.options.base "-twostr_only -schwarzschild -schwarzschild_Nmu 3"


for SOLTIME in 0 300; do
  for ILU in 0 1; do
    write_base_options $ILU
    run_ex $BIN ilu${ILU}_$SOLTIME 5 tenstream.options.base "-max_solution_time $SOLTIME"
  done
done

for SOLTIME in 0 300; do
for LVL in 1 2 3 4; do
for ZREFINE in 1 2; do
  write_mg_options $LVL $ZREFINE $SOLTIME
  run_ex $BIN mg_${LVL}_${ZREFINE}_${SOLTIME} 5 tenstream.options.mg ""
done; done; done
