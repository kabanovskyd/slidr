#!/bin/bash


#========================== Define Conda Paths =======================
#PROFILE_CONDA_PATH=/home/tools/miniconda3/envs/trekker/
#source /home/tools/miniconda3/etc/profile.d/conda.sh

#conda activate ${PROFILE_CONDA_PATH}


TREKKEROUT_ROOT="${OUT_DIR}/${ANALYSIS_DATE}_${SAMPLE_ID}/"
TREKKEROUT_SAMP="${TREKKEROUT_ROOT}/trekker_${SAMPLE_ID}/"
TREKKEROUT_MAIN="${TREKKEROUT_SAMP}/output/"
TREKKEROUT_MISC="${TREKKEROUT_SAMP}/misc/"
TREKKEROUT_TILE="${TREKKEROUT_MISC}/${TILE_ID}/"


mkdir -p "${TREKKEROUT_ROOT}"
mkdir -p "${TREKKEROUT_SAMP}"
mkdir -p "${TREKKEROUT_MAIN}"
mkdir -p "${TREKKEROUT_MISC}"
mkdir -p "${TREKKEROUT_TILE}"


if [[ "$SC_PLATFORM" == "TrekkerQ_P" ]]; then
   if [[ ! -f "$SC_OUTDIR/$barcodes" || ! -f "$SC_OUTDIR/$features" || ! -f "$SC_OUTDIR/$matrix" ]]; then
      if [[ -f "$SC_OUTDIR/count_matrix.mtx" && -f "$SC_OUTDIR/cell_metadata.csv" && -f "$SC_OUTDIR/all_genes.csv" ]] || [[ -f "$SC_OUTDIR/count_matrix.mtx.gz" && -f "$SC_OUTDIR/cell_metadata.csv.gz" && -f "$SC_OUTDIR/all_genes.csv.gz" ]]; then
         echo "Start parse_preprocessing"
         python "${SCRIPT_DIR}/parse_preprocessing.py" \
            "${SC_OUTDIR}" \
            "${SCRIPT_DIR}" >& \
            "${LOG_DIR}/parse_preprocessing.log"
         echo "parse_preprocessing Done"
      fi
   fi
fi


echo "Start fastq_parser"
python "${SCRIPT_DIR}/fastq_parser.py" \
   "${SAMPLE_ID}" \
   "${FASTQ_CB}" \
   "${FASTQ_TAGS}" \
   "${SC_OUTDIR}" \
   "${SCRIPT_DIR}/cellbarcode_whitelists/" \
   "${SC_PLATFORM}" \
   "${TREKKEROUT_MISC}" >& \
   "${LOG_DIR}/fastq_parser.log"
echo "fastq_parser Done"


echo "Start splitspatialbarcodes"
python "${SCRIPT_DIR}/splitspatialbarcodes.py" \
   "${TILE_ID_PATH}" \
   "${TREKKEROUT_TILE}" >& \
   "${LOG_DIR}/beadbarcode_splitter.log"
echo "splitspatialbarcodes Done"


echo "Start bead_matching"
python "${SCRIPT_DIR}/bead_matching.py" \
   -b "${TREKKEROUT_MISC}/reads_per_SB_${SAMPLE_ID}.txt" \
   -i "${TILE_ID}" \
   -d "${TREKKEROUT_MISC}" \
   -o "${TREKKEROUT_MISC}/matching_result_${SAMPLE_ID}.csv" >& \
   "${LOG_DIR}/bead_matching.log"
echo "bead_matching Done"


echo "Start spatial"
Rscript --vanilla "${SCRIPT_DIR}/spatial.R" \
   "${SAMPLE_ID}" \
   "${TREKKEROUT_MISC}" \
   "${TREKKEROUT_MAIN}" \
   "${SC_PLATFORM}" \
   "${SUBSAMPLE_UPDATE}" \
   "${CORES}" >& \
   "${LOG_DIR}/spatial.log"
echo "spatial Done"

echo "Start analysis"
Rscript --vanilla "${SCRIPT_DIR}/analysis.R" \
   "${SAMPLE_ID}" \
   "${SC_PLATFORM}" \
   "${SCRIPT_DIR}" \
   "${SC_OUTDIR}" \
   "${TREKKEROUT_MISC}" \
   "${TREKKEROUT_MAIN}" \
   "${PROFILE_CONDA_PATH}" \
   "${SCMULTI_OUTDIR}" \
   "${SCMULTI_PREFIX}" >& \
   "${LOG_DIR}/analysis.log"
echo "analysis Done"


echo "Start mismatch_analysis"
python "${SCRIPT_DIR}/mismatch_analysis.py" \
   "${SAMPLE_ID}" \
   "${TREKKEROUT_MISC}" >& \
   "${LOG_DIR}/mismatch_analysis.log"
echo "mismatch_analysis Done"


echo "Start genmetrics"
python "${SCRIPT_DIR}/genmetrics.py" \
   "${SAMPLE_ID}" \
   "${TILE_ID}" \
   "${SC_PLATFORM}" \
   "${LOG_DIR}" \
   "${TREKKEROUT_MISC}" \
   "${TREKKEROUT_MAIN}" >& \
   "${LOG_DIR}/genmetrics.log"
echo "genmetrics Done"


echo "Start genreport"
Rscript --vanilla "${SCRIPT_DIR}/genreport.R" \
   "${SAMPLE_ID}" \
   "${TILE_ID}" \
   "${SCRIPT_DIR}/htmlrender.Rmd" \
   "${TREKKEROUT_MAIN}" \
   "${SC_PLATFORM}" \
   "${TREKKEROUT_MAIN}/${SAMPLE_ID}_summary_metrics.csv" \
   "${TREKKEROUT_MAIN}/${SAMPLE_ID}_variable_features_clusters.csv" \
   "${TREKKEROUT_MAIN}/${SAMPLE_ID}_variable_features_spatial_moransi.txt" \
   "${TREKKEROUT_MAIN}/intermediates/${SAMPLE_ID}_seurat_spatial.rds" \
   "${TREKKEROUT_MAIN}/${SAMPLE_ID}_ConfPositioned_seurat_spatial.rds" \
   "${TREKKEROUT_MISC}/${SAMPLE_ID}_mismatchfreq.csv" >& \
   "${LOG_DIR}/genreport.log"
echo "genreport Done"
