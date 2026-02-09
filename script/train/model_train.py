import sys
import logging
import os
import tensorflow as tf
import tensorflow_addons as tfa
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Layer, Flatten, Dense, Bidirectional, Dropout, LSTM, \
    LayerNormalization, MultiHeadAttention, Concatenate
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv2D, BatchNormalization, Activation, Add
from datetime import datetime
import numpy as np
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
gpus = tf.config.experimental.list_physical_devices('GPU')

if gpus:
    try:
        for gpu in gpus:
            tf.config.experimental.set_memory_growth(gpu, True)
    except RuntimeError as e:
        print(e)


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
    if not os.path.exists(path):
        return path

    index = 1
    while True:
        new_path = f"{path}_{index}"
        if not os.path.exists(new_path):
            return new_path
        index += 1


original_stdout = sys.stdout

logger = logging.getLogger('print_logger')
logger.setLevel(logging.INFO)
log_path = "train.log"
log_path = get_available_path_log(log_path)
file_handler = logging.FileHandler(log_path)
formatter = logging.Formatter('%(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


class PrintToLogger:
    def __init__(self):
        self.buffer = ''

    def write(self, message):
        self.buffer += message
        if '\n' in self.buffer:
            lines = self.buffer.split('\n')
            for line in lines[:-1]:
                logger.info(line.strip())
            self.buffer = lines[-1]

    def flush(self):
        if self.buffer:
            logger.info(self.buffer.strip())
            self.buffer = ''


sys.stdout = PrintToLogger()


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
    data = tf.strings.split(line, sep=' ')
    data = data[2:]

    data = tf.strings.to_number(data, out_type=tf.int32)

    last_column = data[-1]
    data = data[:-1]

    data = tf.reshape(data, (72, 5, 4))
    one_hot_encoded = tf.one_hot(last_column, depth=2)

    return data, one_hot_encoded


current_time = datetime.now()
formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
print(os.path.basename(__file__))
print(f"Time: {current_time}")
current_directory = os.getcwd()
print(f"Current Directory: {current_directory}")

ft = sys.argv[1]

dataset = tf.data.TextLineDataset(ft)
your_total_dataset_size = dataset.reduce(tf.constant(0, dtype=tf.int32), lambda x, _: x + 1).numpy()

data = dataset.map(parse_data_line)

num_samples = your_total_dataset_size

train_size = int(0.8 * your_total_dataset_size)
test_size = your_total_dataset_size - train_size

train_data = data.take(train_size)
test_data = data.skip(train_size)

repeat_count = 40
train_data = train_data.repeat(repeat_count)


batch_size = 256
train_data = train_data.batch(batch_size)
test_data = test_data.batch(batch_size)


e = 40
batch = 256
print(f"epoch: {e}")
print(f"batch size: {batch}")

num_trains = train_size


def residual_convnet_model(input_shape, num_classes):
    inputs = Input(shape=input_shape)

    x = Conv2D(128, kernel_size=(7, 5), padding='same')(inputs)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)

    x = make_basic_block_layer(filter_num=128, blocks=1, stride=1)(x)

    x = Conv2D(128, kernel_size=(7, 5), strides=(1, 1), padding="same")(x)
    x = BatchNormalization()(x)
    x = Activation('relu')(x)

    x = make_basic_block_layer(filter_num=128, blocks=1, stride=1)(x)

    x = Flatten()(x)
    x = tf.keras.layers.Reshape((128, -1))(x)

    x = Bidirectional(LSTM(128, return_sequences=True))(x)

    # 添加LSTM层

    x = MultiHeadAttention(key_dim=128, num_heads=8)(x,x)
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


# Model parameters
input_shape = (72, 5, 4)
num_classes = 2
learning_rate = 0.001
total_steps = num_trains * e // batch

strategy = tf.distribute.MirroredStrategy()
model = residual_convnet_model(input_shape, num_classes)
optimizer = tfa.optimizers.Lookahead(tfa.optimizers.RectifiedAdam(
    lr=learning_rate,
    total_steps=total_steps,
    warmup_proportion=0.1,
    min_lr=learning_rate * 0.8,
))

model.compile(optimizer=optimizer,
              loss='categorical_crossentropy',
              metrics=['acc'])


path = sys.argv[2]
path = get_available_path(path)
early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5, mode="min")
model_save_callback = tf.keras.callbacks.ModelCheckpoint("./" + path + "/Select" + ".{epoch:02d}", period=1,
                                                         save_weights_only=False)
model_best_callback = tf.keras.callbacks.ModelCheckpoint("./" + path + "/best_val_loss", monitor='val_loss',
                                                         save_best_only=True, mode="min")
train_log_callback = tf.keras.callbacks.CSVLogger("./" + path + "/training.log", separator='\t')
tf_callback = tf.keras.callbacks.TensorBoard(log_dir="./" + path + "/logs")
csv_logger = tf.keras.callbacks.CSVLogger("./" + path + '/training_log.csv')

sys.stdout = original_stdout
model.summary()

steps_per_epoch = train_size // batch
train_history = model.fit(train_data,
                          epochs=e,
                          validation_data=test_data,
                          callbacks=[early_stop_callback,
                                     model_save_callback,
                                     model_best_callback,
                                     tf_callback,
                                     train_log_callback,
                                     csv_logger],
                          steps_per_epoch=steps_per_epoch,
                          verbose=1,
                          shuffle=True)

sys.stdout = PrintToLogger()

if 'val_loss' in train_history.history:
    best_validation_epoch = np.argmin(np.array(train_history.history["val_loss"])) + 1
    print("[INFO] Best validation loss at epoch: %d" % best_validation_epoch)
    info_str = str(best_validation_epoch)
else:
    best_train_epoch = np.argmin(np.array(train_history.history["loss"])) + 1
    print("[INFO] Best train loss at epoch: %d" % best_train_epoch)
    info_str = str(best_train_epoch)

current_time = datetime.now()
formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
print(f"Time: {current_time}")

