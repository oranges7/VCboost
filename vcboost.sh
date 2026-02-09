SCRIPT=$(dirname $(readlink -f "$0"))
if [[ "$@" == *"-train"* ]]; then
    # 调用训练脚本，并传递所有参数
    sh ${SCRIPT}/script/train/run_train.sh "$@"
else
    # 调用预测脚本，并传递所有参数
    sh ${SCRIPT}/script/predict/run_predict.sh "$@"
fi