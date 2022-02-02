import torch.nn as nn
import torch.nn.functional as F
from torch import randn
import numpy as np

'''
ResNet 50 implementation for Armadillo. 
Adapted from Aladdin Persson's code (https://github.com/aladdinpersson/Machine-Learning-Collection) and The Bioinformatics Repository (https://github.com/bioinform)
'''
class ArmBlock(nn.Module):
    def __init__(self, in_channels, intermediate_channels, identity_downsample=None, stride=1):
        super(ArmBlock, self).__init__()
        self.expansion = 4
        self.conv1 = nn.Conv2d(in_channels, intermediate_channels, kernel_size=1, stride=1, padding=0)
        self.bn1 = nn.BatchNorm2d(intermediate_channels)
        self.conv2 = nn.Conv2d(intermediate_channels, intermediate_channels, kernel_size=3, stride=stride, padding=1)
        self.bn2 = nn.BatchNorm2d(intermediate_channels)
        self.conv3 = nn.Conv2d(intermediate_channels, intermediate_channels * self.expansion, kernel_size=1, stride=1, padding=0)
        self.bn3 = nn.BatchNorm2d(intermediate_channels * self.expansion)
        self.relu = nn.ReLU()
        self.identity_downsample = identity_downsample
        self.stride = stride

    def forward(self, x):
        identity = x.clone()

        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.conv2(x)
        x = self.bn2(x)
        x = self.relu(x)
        x = self.conv3(x)
        x = self.bn3(x)

        if self.identity_downsample is not None:
            identity = self.identity_downsample(identity)

        x += identity
        x = self.relu(x)
        return x


class ArmNet(nn.Module):
    def __init__(self, ArmBlock, layers, channels = 7, num_classes = 3):
        super(ArmNet, self).__init__()
        self.in_channels = 64
        self.conv1 = nn.Conv2d(channels, self.in_channels, kernel_size=(1,3), stride=1, padding=(0,1))
        self.bn1 = nn.BatchNorm2d(self.in_channels)
        self.relu = nn.ReLU()
        self.maxpool = nn.MaxPool2d((1, 3), padding=(0, 1), stride=(1, 1))

        self.layer1 = self._make_layer(ArmBlock, layers[0], intermediate_channels=64, stride=1)
        self.layer2 = self._make_layer(ArmBlock, layers[1], intermediate_channels=128, stride=2)
        self.layer3 = self._make_layer(ArmBlock, layers[2], intermediate_channels=256, stride=2)
        self.layer4 = self._make_layer(ArmBlock, layers[3], intermediate_channels=512, stride=2)

        self.avgpool = nn.AdaptiveAvgPool2d((1, 1))
        self.fc = nn.Linear(512 * 4, num_classes)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x = self.layer4(x)
        x = self.avgpool(x)
        x = x.reshape(x.shape[0], -1)

        return self.fc(x)

    def _make_layer(self, ArmBlock, num_residual_blocks, intermediate_channels, stride):
        identity_downsample = None
        layers = []

        if stride != 1 or self.in_channels != intermediate_channels * 4:
            identity_downsample = nn.Sequential(nn.Conv2d(self.in_channels, intermediate_channels * 4, kernel_size=1, stride=stride), 
            nn.BatchNorm2d(intermediate_channels * 4))

        layers.append(ArmBlock(self.in_channels, intermediate_channels, identity_downsample, stride))

        self.in_channels = intermediate_channels * 4

        for i in range(num_residual_blocks - 1):
            layers.append(ArmBlock(self.in_channels, intermediate_channels))

        return nn.Sequential(*layers)

