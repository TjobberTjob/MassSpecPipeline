from keras.layers import Dense, Input, Flatten, Conv2D, MaxPooling2D, Concatenate, GlobalAveragePooling1D, \
    GlobalAveragePooling2D
from keras.models import Model
from classes.utilities import r2
from keras.callbacks import ModelCheckpoint, EarlyStopping

class Network_Model():
    def __init__(self, whichMS, classification, n_classes, ms1size, ms2size, n_channels, lenMS2, metapath, imageclass, patience):
        self.patience = patience
        self.imageclass = imageclass
        self.metapath = metapath
        self.lenMS2 = lenMS2
        self.n_channels = n_channels
        self.ms2size = ms2size
        self.ms1size = ms1size
        self.whichMS = whichMS
        self.classification = classification
        self.n_classes = n_classes

    def get_phase1_ms1_net(self):
        input = Input(shape=(self.ms1size[0], self.ms1size[1], self.n_channels,))
        x = Conv2D(32, kernel_size=(5, 5), activation='relu', padding='same')(input)
        x = MaxPooling2D(pool_size=(5, 5), strides=(1, 1), padding='same')(x)
        x1 = Conv2D(32, kernel_size=(3, 3), activation='relu', padding='same')(input)
        x1 = MaxPooling2D(pool_size=(3, 3), strides=(1, 1), padding='same')(x1)
        x2 = Conv2D(32, kernel_size=(2, 2), activation='relu', padding='same')(input)
        x2 = MaxPooling2D(pool_size=(2, 2), strides=(1, 1), padding='same')(x2)
        x = Concatenate()([x, x1, x2])
        x = GlobalAveragePooling2D()(x)

        return input, x

    def get_phase1_ms2_net(self):
        input = Input(shape=(self.lenMS2, self.ms2size[1]))
        x = Dense(128, activation='relu')(input)
        x = Dense(64, activation='relu')(input)
        x = Flatten()(x)

        return input, x

    def get_phase2_net(self):
        network_input = Dense(128, activation='relu')(self.p1_output_layer)
        network_output = Dense(64, activation='relu')(network_input)
        return network_input, network_output

    def get_phase3_net(self):
        if self.classification:
            return Dense(self.n_classes, activation='softmax')
        else:
            return Dense(self.n_classes, activation='linear')


    def get_network(self):
        if self.whichMS == 'ms1':
            (self.p1_input_layer, self.p1_output_layer) = self.get_phase1_ms1_net()

        elif self.whichMS == 'ms2':
            (self.p1_input_layer, self.p1_output_layer) = self.get_phase1_ms2_net()

        else:
            (ms1_input, ms1_output) = self.get_phase1_ms1_net()
            (ms2_input, ms2_output) = self.get_phase1_ms2_net()
            self.p1_input_layer = [ms1_input, ms2_input]

            self.p1_output_layer = Concatenate()([ms1_output, ms2_output])

        (p2_input_layer, p2_output_layer) = self.get_phase2_net()

        self.p3_output_layer = self.get_phase3_net()(p2_output_layer)

        self.model = Model(inputs=self.p1_input_layer, outputs=self.p3_output_layer)

        if self.classification:
            self.model.compile(loss='categorical_crossentropy', metrics=['accuracy'], optimizer='adam')

        else:
            self.model.compile(loss='mse', metrics=[r2], optimizer='rmsprop')

        return self.model


    def get_callbacks(self):
        if self.classification:
            checkpoint = ModelCheckpoint(f'{self.metapath}Best-{self.imageclass}.h5', monitor='val_accuracy', save_best_only=True)
        else:
            checkpoint = ModelCheckpoint(f'{self.metapath}Best-{self.imageclass}.h5', monitor='val_r2', save_best_only=True)
        early_stopping = EarlyStopping(monitor='val_loss', patience=self.patience)
        return [checkpoint, early_stopping]