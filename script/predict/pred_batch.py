import tensorflow as tf
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Layer, Flatten, Dense, Bidirectional, Dropout, LSTM, \
    LayerNormalization, MultiHeadAttention, Concatenate
from tensorflow.keras.models import Model
from tensorflow.keras import backend as K
import sys
import inspect
from tensorflow.keras.layers import Input, Conv2D, BatchNormalization, Activation, Add
import pandas as pd
import numpy as np
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'


def get_available_path_log(path):

    if not os.path.exists(path):
        return path

    index = 1
    while True:

        path1 = path.split('.')[0]
        path2 = path.split('.')[1]
        new_path = f"{path1}_{index}.{path2}"
        if not os.path.exists(new_path):
            return new_path
        index += 1


def get_available_path(path):
    # 如果文件夹不存在，直接返回原始路径
    if not os.path.exists(path):
        return path

    # 文件夹存在，尝试在原始路径后面加 1
    index = 1
    while True:
        new_path = f"{path}_{index}"
        if not os.path.exists(new_path):
            return new_path
        index += 1

def calculate_metrics(confusion_matrix):
    num_classes = confusion_matrix.shape[0]
    precision = np.zeros(num_classes)
    recall = np.zeros(num_classes)
    f1 = np.zeros(num_classes)

    for i in range(num_classes):
        tp = confusion_matrix[i, i]
        fp = np.sum(confusion_matrix[:, i]) - tp
        fn = np.sum(confusion_matrix[i, :]) - tp

        precision[i] = tp / (tp + fp)
        recall[i] = tp / (tp + fn)
        f1[i] = 2 * (precision[i] * recall[i]) / (precision[i] + recall[i])

    return precision, recall, f1


def parse_data_line(line):
    # 这里假设数据行是以空格分隔的数字
    # 首先使用 tf.strings.split() 函数将字符串分割成列表
    data = tf.strings.split(line, sep=' ')
    data = data[2:]
    data = tf.strings.to_number(data, out_type=tf.int32)
    # 重新调整数据的维度形状，这里假设你希望数据的形状为 (72, 5, 4)
    data = tf.reshape(data, (72, 5, 4))

    return data


# 划分数据集为训练集和测试集


e = 40
batch = 350

def residual_block(input_tensor, filters):
    x = Conv2D(filters, kernel_size=(3, 3), padding='same')(input_tensor)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)
    x = Conv2D(filters, kernel_size=(3, 3), padding='same')(x)
    x = BatchNormalization()(x)
    x = Add()([x, input_tensor])
    x = Activation('relu')(x)
    return x


class DifferenceLayer(Layer):
    def __init__(self, **kwargs):
        super(DifferenceLayer, self).__init__(**kwargs)

    def call(self, inputs):
        if inputs.shape[-1] != 4:
            raise ValueError("DifferenceLayer expects the last dimension to be 4.")

        diff1 = K.abs(inputs[:, :, :, 0] - inputs[:, :, :, 2])
        diff2 = K.abs(inputs[:, :, :, 1] - inputs[:, :, :, 3])
        diff3 = K.abs(inputs[:, :, :, 0] - inputs[:, :, :, 1])
        diff4 = K.abs(inputs[:, :, :, 2] - inputs[:, :, :, 3])
        merged_diff = K.stack([diff1, diff2, diff3, diff4], axis=-1)

        return merged_diff


def residual_convnet_model(input_shape, num_classes):
    inputs = Input(shape=input_shape)
    # 添加差异层（放在卷积层之前）
    # x_diff = DifferenceLayer()(inputs)
    # merged = Concatenate(axis=-1)([inputs, x_diff])
    x = Conv2D(128, kernel_size=(7, 5), padding='same')(inputs)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)

    x = make_basic_block_layer(filter_num=128, blocks=1, stride=1)(x)

    x = Conv2D(128, kernel_size=(7, 5), strides=(1, 1), padding="same")(x)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)

    x = make_basic_block_layer(filter_num=128, blocks=1, stride=1)(x)

    # x = MaxPooling2D(pool_size=(3, 3), padding='same')(x)
    # 添加SPP层，指定池化尺度

    # 将SPP层应用于x

    x = Flatten()(x)
    x = tf.keras.layers.Reshape((72, -1))(x)

    x = Bidirectional(LSTM(128, return_sequences=True))(x)

    # 添加LSTM层

    x = MultiHeadAttention(key_dim=128, num_heads=8)(x, x)
    x = Bidirectional(LSTM(128, return_sequences=True))(x)
    x = Bidirectional(LSTM(128, return_sequences=True))(x)
    x = LSTM(128)(x)
    # x = Dense(units=128, activation='relu')(x)
    # x = Dropout(0.5)(x)
    x = Dense(units=256, activation='relu')(x)
    x = Dropout(0.5)(x)

    outputs = Dense(units=num_classes, activation='softmax')(x)

    model = Model(inputs=inputs, outputs=outputs)

    return model


L2_regularizers = tf.keras.regularizers.l2(1e-7)


class BasicBlock(tf.keras.layers.Layer):

    def __init__(self, filter_num, stride=1):
        super(BasicBlock, self).__init__()
        # conv = tf.keras.layers.SeparableConv2D if SeparableConv else tf.keras.layers.Conv2D
        self.filter_num = filter_num
        self.stride = stride
        self.conv1 = tf.keras.layers.Conv2D(filters=filter_num,
                                            kernel_size=(7, 5),
                                            strides=stride,
                                            padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn1 = tf.keras.layers.BatchNormalization()
        self.conv2 = tf.keras.layers.Conv2D(filters=filter_num,
                                            kernel_size=(7, 5),
                                            strides=1,
                                            padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn2 = tf.keras.layers.BatchNormalization()
        if stride != 1:
            self.downsample = tf.keras.Sequential()
            self.downsample.add(tf.keras.layers.Conv2D(filters=filter_num,
                                                       kernel_size=(7, 5),
                                                       strides=stride,
                                                       kernel_regularizer=L2_regularizers))
            self.downsample.add(tf.keras.layers.BatchNormalization())
        else:
            self.downsample = lambda x: x

    def call(self, inputs):
        residual = self.downsample(inputs)

        x = self.conv1(inputs)
        x = self.bn1(x, )
        x = tf.nn.relu(x)
        x = self.conv2(x)
        x = self.bn2(x, )

        output = tf.nn.relu(tf.keras.layers.add([residual, x]))

        return output

    def get_config(self):
        config = super().get_config()
        config.update({
            'filter_num': self.filter_num,
            'stride': self.stride,
        })
        return config


def make_basic_block_layer(filter_num, blocks, stride=1):
    res_block = tf.keras.Sequential()

    res_block.add(BasicBlock(filter_num, stride=stride))

    for _ in range(1, blocks):
        res_block.add(BasicBlock(filter_num, stride=1))

    return res_block


def create_model(input_shape, num_classes):
    inputs = tf.keras.Input(shape=input_shape)

    conv1 = Conv2D(64, kernel_size=(5, 5), strides=(1, 1), padding="same")(inputs)
    conv1 = BatchNormalization()(conv1)
    conv1 = Activation('relu')(conv1)

    # Residual block 1
    res_block1 = make_basic_block_layer(filter_num=64, blocks=1, stride=1, SeparableConv=False)(conv1)

    conv2 = Conv2D(128, kernel_size=(5, 5), strides=(1, 1), padding="same")(res_block1)
    conv2 = BatchNormalization()(conv2)
    conv2 = Activation('relu')(conv2)

    # Residual block 2
    res_block2 = make_basic_block_layer(filter_num=128, blocks=1, stride=1, SeparableConv=False)(conv2)

    conv3 = Conv2D(256, kernel_size=(5, 5), strides=(1, 1), padding="same")(res_block2)
    conv3 = BatchNormalization()(conv3)
    conv3 = Activation('relu')(conv3)

    # Residual block 3
    res_block3 = make_basic_block_layer(filter_num=256, blocks=1, stride=1)(conv3)

    flatten = Flatten()(res_block3)

    dropout = Dropout(0.4)(flatten)

    dense1 = Dense(256, activation='relu')(dropout)
    dropout1 = Dropout(0.4)(dense1)

    dense2 = Dense(num_classes, activation='softmax')(dropout1)

    model = Model(inputs=inputs, outputs=dense2)

    return model


pred = sys.argv[1]
# Model parameters
input_shape = (72, 5, 4)
# print(f"特征维度: {input_shape}")
num_classes = 2

source_code = inspect.getsource(residual_convnet_model)
# print(f"模型代码:\n{source_code}")

model = residual_convnet_model(input_shape, num_classes)

path = sys.argv[2]
# output = sys.argv[3]
model.load_weights(path)
data = pd.read_csv(pred, sep=' ', header=None)
memory_usage = data.memory_usage().sum()
memory_usage_gib = memory_usage / (1024 ** 3)
# print(f"Memory usage: {memory_usage_gib}GiB")
# print("HG004预测读取数据完成")
# data = np.array(data)
# labels=data.iloc[:,-1].values.astype(np.int32)
chr = data.iloc[:, 0]
pos = data.iloc[:, 1]
data = data.iloc[:, 2:].values.astype(np.int32)
# labels_one = np.eye(2)[labels.astype(int)]
batch = 256
# 调整数据形状为(样本数, 5,128,2)
data = data.reshape(-1, 72, 5, 4)

pred_labels_one = model.predict(data)
chr = chr.to_numpy().reshape(-1, 1)
pos = pos.to_numpy().reshape(-1, 1)
output_file_path=pred+'.txt'

out_array=np.concatenate([chr,pos,pred_labels_one], axis=1)

df = pd.DataFrame(out_array)

df.to_csv(output_file_path, sep=' ', index=False, header=False)


