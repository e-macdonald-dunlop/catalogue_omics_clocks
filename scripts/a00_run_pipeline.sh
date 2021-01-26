#!/bin/bash

mkdir -p data
mkdir -p logs

mkdir -p a_minus_b
mkdir -p b_only
mkdir -p a_only


cd a_minus_b/

mkdir data
mkdir logs
cp -r ../control_files/ .


cd b_only
mkdir data
mkdir logs
cp -r ../control_files/ .

mkdir 3_pcs
cd 3_pcs
mkdir data
mkdir logs
cp -r ../control_files/ .


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


mkdir -p p03b_pca_clustering

cd p03b_pca_clustering

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a03b_pca_clustering.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ../logs/a03b_pca_clusteringlog.log

cd ../

module load igmm/apps/R/3.6.0

mkdir -p p04_clock_500_iterations

cd p04_clock_500_iterations

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a04_clock_500_iterations.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a04_clock_500_iterations_log.log

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


mkdir -p p07_model_error_figures

cd p07_model_error_figures

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a07_model_error_figures.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a07_model_error_figures_log.log

cd ../

mkdir -p p08_smoking_status_by_cohort

cd p08_smoking_status_by_cohort

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a08_smoking_status_by_cohort.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ../logs/a08_smoking_status_by_cohort_log.log

cd ../


mkdir -p p10_correlation_distribution

cd p10_correlation_distribution

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a10_correlation_distribution.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a10_correlation_distribution_log.log

cd ../


mkdir -p p13_smoking_figures

cd p13_smoking_figures

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a13_smoking_figures.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ../logs/a13_smoking_figures_log.log

cd ../



mkdir -p p14_make_b_only

cd p14_make_b_only

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a14_make_b_only.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a14_make_b_only_log.log

cd ../


mkdir -p p14b_make_a_only

cd p14b_make_a_only

Rscript /exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/a14b_make_a_only.R "/exports/igmm/eddie/wilson-lab/projects/prj_107_omics_clock_pipeline/scripts/eddie/" &> ../logs/a14b_make_a_only_log.log

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



mkdir p05_predictor_inclusion_negative

cd p05_predictor_inclusion_negative

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a05b_predictor_inclusion_negative_effects_only.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/"

cd ../

cd /opt/working/wilson/projects/prj_086_omics_clocks/pipeline/age_at_vene/figures/

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a11_compare_omics_within.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ./logs/a11_compare_omics_within_log.log



cd /opt/working/wilson/projects/prj_086_omics_clocks/pipeline/age_at_vene/figures/

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a12_r_suared_figure_within.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ./logs/a12_r_suared_figure_within_log.log



cd /opt/working/wilson/projects/prj_086_omics_clocks/pipeline/age_at_vene/compare_across_iterations/figures/

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a11_compare_omics_within.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ./logs/a11_compare_omics_within_log.log



cd /opt/working/wilson/projects/prj_086_omics_clocks/pipeline/age_at_vene/compare_across_iterations/figures/

Rscript /opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/a12_r_suared_figure_within.R "/opt/working/wilson/projects/prj_107_omics_clock_pipeline/scripts/" > ./logs/a12_r_suared_figure_within_log.log

