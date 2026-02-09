SCRIPT=$(dirname $(readlink -f "$0"))

object='train'
mode='both'
THREADS=32
RETRIES=3
help_message="
Usage: $(basename "$0") [-w <work_path>] [-a <aim_vcf>] [-b <bam_file>] [-r <ref>] [-q <vcf>] [-m <mode>] [-j <object>] [-p <phase>]
Options:
  -w, --work_path   Working directory path.
  -a, --aim_vcf     Aim VCF file path.
  -b, --bam_file    BAM file path.
  -r, --ref         Reference file path.
  -q, --vcf         Benchmark VCF file path.
  -m, --mode        Mode of operation.You can choose snp or indel. The default is both.
  -j, --object      Object to process.
  -p, --phase       Enable phase.
"

while getopts ":w:a:b:r:q:m:j:ph" opt; do
  case $opt in
    w|--work_path)
      work_path="$OPTARG"
      ;;
    a|--aim_vcf)
      aim_vcf="$OPTARG"
      ;;
    b|--bam_file)
      BAM_FILE_PATH="$OPTARG"
      ;;
    r|--ref)
      REFERENCE_FILE_PATH="$OPTARG"
      ;;
    q|--vcf)
      B_VCF="$OPTARG"
      ;;
    m|--mode)
      mode="$OPTARG"
      ;;
    j|--object)
      object="$OPTARG"
      ;;
    p|--phase)
      phase=1
      ;;
    h|help)
      echo $help_message
      exit 0
      ;;
    train)
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


PARALLEL='parallel'
LOG_PATH=${work_path}'/tmp/logs'
WHATSHAP='whatshap'
SAMTOOLS='samtools'
CHR=${work_path}'/tmp/CONTIGS'
PHASE_VCF_PATH=${work_path}'/tmp/phase_vcf'
PHASE_BAM_PATH=${work_path}'/tmp/phase_bam'
pileup_output=${work_path}'/tmp/pileup_output'
OUTPUT=${work_path}'/data'
REFERENCE_FILE_PATH='/home/gyc/ont/seq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
GEN_PY='generate_hap_both_onehot_refalt.py'
OUTPUT_DATA=${OUTPUT}'/onehot'


rm -rf ${OUTPUT}

mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}
mkdir -p ${LOG_PATH}
mkdir -p ${OUTPUT}
mkdir -p ${OUTPUT_DATA}
(time (

echo "work_path: $work_path"
echo "BAM_FILE_PATH: $BAM_FILE_PATH"
echo "B_VCF: $B_VCF"
echo "mode: $mode"
echo "object: $object"

if [[ "$aim_vcf" == *.gz ]]; then
    gzip -dc ${aim_vcf} > ${OUTPUT}/original.vcf
    aim_vcf=${OUTPUT}/original.vcf
fi
if [ "$phase" -eq 1 ]; then
  python3 ${SCRIPT}'/select_snp_for_phase.py' --vcf_fn  ${OUTPUT}/original.vcf --split_folder ${PHASE_VCF_PATH}
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
  echo "[INFO] 4/7 Haplotag input BAM file using Whatshap"
      time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_2_haplotag.log -j${THREADS} \
      "${WHATSHAP} haplotag \
          --output ${PHASE_BAM_PATH}/{1}.bam \
          --reference ${REFERENCE_FILE_PATH} \
          --ignore-read-groups \
          --regions {1} \
          ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
          ${BAM_FILE_PATH}" :::: ${CHR} |& tee ${LOG_PATH}/haplotag.log
      ${PARALLEL} -j${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/{1}.bam :::: ${CHR}
fi
if [ ${mode} = 'snp' ] || [ ${mode} = 'both' ]
then
  snp_mode='snp'
  python3 ${SCRIPT}'/select_pf_snp.py' ${OUTPUT}/original.vcf ${B_VCF} ${snp_mode}
  echo "[INFO] 位点分类完毕"
  if test ${object} = 'train'
  then
    num_f=`wc -l ${OUTPUT}'/f_'${snp_mode}'.txt' | awk '{print $1}'`
    echo $num_f
    shuf -n $num_f ${OUTPUT}'/p_'${snp_mode}'.txt' > ${OUTPUT}'/p_s.txt'
    rm ${OUTPUT}'/p_'${snp_mode}'.txt'
    mv ${OUTPUT}'/p_s.txt' ${OUTPUT}'/p_'${snp_mode}'.txt'
  fi
#  echo "[INFO] 文件清除完毕"
  python3 ${SCRIPT}'/split_chr.py' ${OUTPUT}'/f_'${snp_mode}'.txt' ${snp_mode}
  python3 ${SCRIPT}'/split_chr.py' ${OUTPUT}'/p_'${snp_mode}'.txt' ${snp_mode}

  ls ${OUTPUT}/f_chr*_${snp_mode}* > ${OUTPUT}/f_file_${snp_mode}.txt
  ls ${OUTPUT}/p_chr*_${snp_mode}* > ${OUTPUT}/p_file_${snp_mode}.txt
  awk -F '/' '{print $NF}' ${OUTPUT}/f_file_${snp_mode}.txt > ${OUTPUT}/f_cut_${snp_mode}.txt
  awk -F '/' '{print $NF}' ${OUTPUT}/p_file_${snp_mode}.txt > ${OUTPUT}/p_cut_${snp_mode}.txt
  awk -F '_' '{print $2}' ${OUTPUT}/f_cut_${snp_mode}.txt > ${OUTPUT}/f_file_num_${snp_mode}.txt
  awk -F '_' '{print $2}' ${OUTPUT}/p_cut_${snp_mode}.txt > ${OUTPUT}/p_file_num_${snp_mode}.txt
  time ${PARALLEL} --link --retries ${RETRIES} --bar --joblog ${LOG_PATH}/parallel_gen_data_f.log -j${THREADS} \
        "python3 ${SCRIPT}'/'${GEN_PY} \
            --bam_path ${PHASE_BAM_PATH}/{1}.bam \
            --ref_path ${REFERENCE_FILE_PATH} \
            --chromosome {1} \
            -m {3} \
            -o ${OUTPUT_DATA} \
            -l f \
            --pos_path ${OUTPUT}/{2}" :::: ${OUTPUT}/f_file_num_${snp_mode}.txt :::: ${OUTPUT}/f_cut_${snp_mode}.txt ::: ${snp_mode} |& tee ${LOG_PATH}/gen_data_f.log
  time ${PARALLEL} --link --retries ${RETRIES} --bar --joblog ${LOG_PATH}/parallel_gen_data_p.log -j${THREADS} \
        "python3 ${SCRIPT}'/'${GEN_PY} \
            --bam_path ${PHASE_BAM_PATH}/{1}.bam \
            --ref_path ${REFERENCE_FILE_PATH} \
            --chromosome {1} \
            -m {3} \
            -o ${OUTPUT_DATA} \
            -l p \
            --pos_path ${OUTPUT}/{2}" :::: ${OUTPUT}/p_file_num_${snp_mode}.txt :::: ${OUTPUT}/p_cut_${snp_mode}.txt ::: ${snp_mode} |& tee ${LOG_PATH}/gen_data_p.log

fi
if [ ${mode} = 'indel' ] || [ ${mode} = 'both' ]
then
  indel_mode='indel'
  python3 ${SCRIPT}'/select_pf_indel.py' ${OUTPUT}/original.vcf ${B_VCF} ${indel_mode}

#  del='rm '${OUTPUT}'/*chr*'
  if test ${object} = 'train'
  then
    num_f=`wc -l ${OUTPUT}'/f_'${indel_mode}'.txt' | awk '{print $1}'`
    echo $num_f

    shuf -n $num_f ${OUTPUT}'/p_'${indel_mode}'.txt' > ${OUTPUT}'/p_s.txt'
    rm ${OUTPUT}'/p_'${indel_mode}'.txt'
    mv ${OUTPUT}'/p_s.txt' ${OUTPUT}'/p_'${indel_mode}'.txt'
  fi

  python3 ${SCRIPT}'/split_chr.py' ${OUTPUT}'/f_'${indel_mode}'.txt' ${indel_mode}
  python3 ${SCRIPT}'/split_chr.py' ${OUTPUT}'/p_'${indel_mode}'.txt' ${indel_mode}

  ls ${OUTPUT}/f_chr*_${indel_mode}* > ${OUTPUT}/f_file_${indel_mode}.txt
  ls ${OUTPUT}/p_chr*_${indel_mode}* > ${OUTPUT}/p_file_${indel_mode}.txt
  awk -F '/' '{print $NF}' ${OUTPUT}/f_file_${indel_mode}.txt > ${OUTPUT}/f_cut_${indel_mode}.txt
  awk -F '/' '{print $NF}' ${OUTPUT}/p_file_${indel_mode}.txt > ${OUTPUT}/p_cut_${indel_mode}.txt
  awk -F '_' '{print $2}' ${OUTPUT}/f_cut_${indel_mode}.txt > ${OUTPUT}/f_file_num_${indel_mode}.txt
  awk -F '_' '{print $2}' ${OUTPUT}/p_cut_${indel_mode}.txt > ${OUTPUT}/p_file_num_${indel_mode}.txt
  time ${PARALLEL} --link --retries ${RETRIES} --bar --joblog ${LOG_PATH}/parallel_gen_data_f.log -j${THREADS} \
        "python3 ${SCRIPT}'/'${GEN_PY} \
            --bam_path ${PHASE_BAM_PATH}/{1}.bam \
            --ref_path ${REFERENCE_FILE_PATH} \
            --chromosome {1} \
            -m {3} \
            -o ${OUTPUT_DATA} \
            -l f \
            --pos_path ${OUTPUT}/{2}" :::: ${OUTPUT}/f_file_num_${indel_mode}.txt :::: ${OUTPUT}/f_cut_${indel_mode}.txt ::: ${indel_mode} |& tee ${LOG_PATH}/gen_data_f.log
  time ${PARALLEL} --link --retries ${RETRIES} --bar --joblog ${LOG_PATH}/parallel_gen_data_p.log -j${THREADS} \
        "python3 ${SCRIPT}'/'${GEN_PY} \
            --bam_path ${PHASE_BAM_PATH}/{1}.bam \
            --ref_path ${REFERENCE_FILE_PATH} \
            --chromosome {1} \
            -m {3} \
            -o ${OUTPUT_DATA} \
            -l p \
            --pos_path ${OUTPUT}/{2}" :::: ${OUTPUT}/p_file_num_${indel_mode}.txt :::: ${OUTPUT}/p_cut_${indel_mode}.txt ::: ${indel_mode} |& tee ${LOG_PATH}/gen_data_p.log

fi

rm -rf ${PHASE_BAM_PATH}
rm -rf ${pileup_output}

cat "${OUTPUT_DATA}"/*.pileup > "${OUTPUT_DATA}/pileup.data"
shuf "${OUTPUT_DATA}/pileup.data" > "${OUTPUT_DATA}/shufed.data"
python3 ${SCRIPT}/model_train.py "${OUTPUT_DATA}/shufed.data" ${work_path}

)) |& tee ${LOG_PATH}/gen_data.log