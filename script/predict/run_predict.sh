SCRIPT=$(dirname $(readlink -f "$0"))
contig='chr1-22'
PARALLEL='parallel'
WHATSHAP='whatshap'
SAMTOOLS='samtools'
THREADS=24
RETRIES=3
phase=1
GEN_PY='generate_hap_both_onehot_refalt.py'

help_message="
Usage: $(basename "$0") [-o <out_path>] [-b <bam_file>] [-v <vcf>] [-m <model>] [-r <ref_path>] [-p <phase>] [-h <help>]
Options:
  -o, -out_path    Output path.
  -b, -bam_file    BAM file path.
  -v, -vcf         VCF file path.
  -m, -model       Model path.
  -r, -ref_path    Reference file path.
  -t, -threads     Number of threads.The default is 32.
  -c, -contig      Contig to process.The default is chr1-22.
  -p, -phase       Enable phase.
  -h, -help        Display this help message.
"

while getopts ":o:b:v:m:r:t:c:d:ph" opt; do
  case $opt in
    o|out_path)
      out_path="$OPTARG"
      ;;
    b|bam_file)
      BAM_FILE_PATH="$OPTARG"
      ;;
    v|vcf)
      aim_vcf="$OPTARG"
      ;;
    m|model)
      model="$OPTARG"
      ;;
    r|ref_path)
      REFERENCE_FILE_PATH="$OPTARG"
      ;;
    t|threads?)
      THREADS="$OPTARG"
      ;;
    c|contig?)
      contig="$OPTARG"
      ;;
    d|model_path)
      model_path="$OPTARG"
      ;;
    p|phase?)
      phase=2
      ;;
    h|help)
      echo $help_message
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo $help_message
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo $help_message
      exit 1
      ;;
    *)
      echo "Invalid option: -$OPTARG" >&2
      echo $help_message
      exit 1
      ;;
  esac
done

if [ $# -eq 0 ]; then
  echo "No options provided. Use -h or -help for usage information." >&2
  exit 1
fi
if [[ $# -lt 5 ]]; then
  echo "Insufficient parameters" >&2
  echo $help_message
  exit 1
fi

echo "out_path: $out_path"
echo "BAM_FILE_PATH: $BAM_FILE_PATH"
LOG_PATH=${out_path}'/tmp/logs'
CHR=${out_path}'/tmp/CONTIGS'
PHASE_VCF_PATH=${out_path}'/tmp/phase_vcf'
PHASE_BAM_PATH=${out_path}'/tmp/phase_bam'
OUTPUT=${out_path}'/data'
OUTPUT_DATA=${OUTPUT}'/onehot'

mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}
mkdir -p ${LOG_PATH}
mkdir -p ${OUTPUT}
mkdir -p ${out_path}
mkdir -p ${OUTPUT_DATA}

cp ${aim_vcf} ${OUTPUT}
if [[ "$aim_vcf" == *.gz ]]; then
    gzip -dc ${aim_vcf} > ${OUTPUT}/original.vcf
    aim_vcf=${OUTPUT}/original.vcf
fi

python3 ${SCRIPT}'/select_snp_for_phase.py' --vcf_fn  ${aim_vcf} --split_folder ${PHASE_VCF_PATH}
python3 ${SCRIPT}'/generate_chr.py' ${contig} ${CHR}
time ${PARALLEL}  --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_1_phase.log -j${THREADS} \
        "${WHATSHAP} phase \
            --output ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
            --reference ${REFERENCE_FILE_PATH} \
            --chromosome {1} \
            --distrust-genotypes \
            --ignore-read-groups \
            ${PHASE_VCF_PATH}/{1}.vcf \
            ${BAM_FILE_PATH}" :::: ${CHR} |& tee ${LOG_PATH}/phase.log
time ${PARALLEL} -j${THREADS} tabix -f -p vcf ${PHASE_VCF_PATH}/phased_{1}.vcf.gz :::: ${CHR}

time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_2_haplotag.log -j${THREADS} \
        "${WHATSHAP} haplotag \
            --output ${PHASE_BAM_PATH}/{1}.bam \
            --reference ${REFERENCE_FILE_PATH} \
            --ignore-read-groups \
            --regions {1} \
            ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
            ${BAM_FILE_PATH}" :::: ${CHR} |& tee ${LOG_PATH}/haplotag.log
${PARALLEL} -j${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/{1}.bam :::: ${CHR}

python3 ${SCRIPT}'/split_chr.py' ${aim_vcf} ${OUTPUT} ${contig}

ls ${OUTPUT}/chr*_* > ${OUTPUT}/batch.txt
awk -F '/' '{print $NF}' ${OUTPUT}/batch.txt > ${OUTPUT}/cut_batch.txt
awk -F '_' '{print $1}' ${OUTPUT}/cut_batch.txt > ${OUTPUT}/batch_chr.txt
time ${PARALLEL} --bar --link --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_gen_data_f.log -j${THREADS} \
      "python3 ${SCRIPT}'/'${GEN_PY} \
          --bam_path ${PHASE_BAM_PATH}/{1}.bam \
          --ref_path ${REFERENCE_FILE_PATH} \
          --chromosome {1} \
          -o ${OUTPUT_DATA} \
          --pos_path ${OUTPUT}/{2}" :::: ${OUTPUT}/batch_chr.txt :::: ${OUTPUT}/cut_batch.txt |& tee ${LOG_PATH}/gen_data_f.log
ls ${OUTPUT_DATA}/*.pileup > ${OUTPUT_DATA}/all_files.txt
PTHREADS=$(( (${THREADS} / 8) * 3 ))
time ${PARALLEL} --bar -j ${PTHREADS} python3 ${SCRIPT}/pred_batch.py {1} {2} :::: ${OUTPUT_DATA}/all_files.txt ::: ${model_path}/${model}
python3 ${SCRIPT}/batch_con.py -i ${OUTPUT_DATA}/all_files.txt -o ${out_path} -f ${aim_vcf} -a ${OUTPUT_DATA}/miss.pos -c ${contig}

rm -rf ${PHASE_BAM_PATH}
