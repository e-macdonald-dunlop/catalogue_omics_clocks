#!/bin/bash

mkdir -p data
mkdir -p logs

mkdir -p p01_basic_qc

cd p01_basic_qc

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a01_basic_qc.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a01_log.log

cd ../

mkdir -p p02_correct_for_covariates

cd p02_correct_for_covariates

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a02_correct_for_covariates.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a02_log.log

cd ../

mkdir -p p03_make_clock

cd p03_make_clock

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a03_make_clock.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a03_log.log

cd ../


mkdir -p p04_clock_500_iterations

cd p04_clock_500_iterations

qsub /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a04_clock_500_iterations.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" #&> ../logs/a04_clock_500_iterations_log.log

cd ../

mkdir -p p05_predictor_inclusion

cd p05_predictor_inclusion

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a05_predictor_inclusion.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a05_predictor_inclusion_log.log

cd ../


mkdir -p p06_core_model_prediction

cd p06_core_model_prediction

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a06_core_model_prediction.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a06_core_model_prediction_log.log

cd ../


mkdir -p p18_make_pcs

cd p18_make_pcs

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a18_make_pcs.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" 3 5 10 20 &> ../logs/a18_make_pcs_log.log

cd ../

mkdir -p p19_testing_ukb_range

cd p19_testing_ukb_range

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/rep_ukb_age_range.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a19_testing_ukb_range.log

cd ../


mkdir -p p05b_predictor_inclusion_negative_effects_only

cd p05b_predictor_inclusion_negative_effects_only

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a05b_predictor_inclusion_negative_effects_only.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a05b_predictor_inclusion_negative_effects_only_log.log

cd ../


mkdir -p p05c_predictor_inclusion_positive_effects_only

cd p05c_predictor_inclusion_positive_effects_only

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a05c_predictor_inclusion_positive_effects_only.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a05c_predictor_inclusion_positive_effects_only_log.log

cd ../
