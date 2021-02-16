#!/bin/python3
from features import model_features
import CNN_model
from argparse import ArgumentParser
from torch import load, stack, tensor, no_grad, set_num_threads, max as torch_max, device as torch_device
from torch.utils.data import DataLoader, TensorDataset

from sys import stderr

def arguments():
    parser = ArgumentParser(description='Train the model for Armadillo')
    parser.add_argument(
    '-i', '--input', required = True, metavar = 'input', type = str,
    help = 'TSV file with: CHROM, POS, ALT, GROUP, TUMOR_BAM, CONTROL_BAM columns')
    parser.add_argument(
    '-f', '--fasta', required = True, metavar = 'fasta', type = str,
    help = 'Path to fasta file')
    parser.add_argument(
    '-m', '--model', required = True, metavar = 'model', type = str,
    help = 'CNN trained model for Armadillo results')
    parser.add_argument(
    '-t', '--threads', required = False, default = 10, metavar = 'threads', type = int,
    help = 'Number of threads [10]')
    return parser.parse_args()

def fitArmNet(chrom, pos, alt, tumor, control, fasta, model):
    device = torch_device("cpu")
    features = model_features(tumor, control, chrom, int(pos), alt, fasta)

    net = CNN_model.ArmNet(CNN_model.ArmBlock, layers = [3, 8, 36, 3], channels = 5, num_classes = 2).to(device)   
    net.load_state_dict(load(model))
    net.eval()

    outputs = net(features).to(device)
    _, predicted = torch_max(outputs, 1)

    keep = True if predicted.item() == 0 else False

    return(keep) 

def main():
    args = arguments()
    set_num_threads(args.threads)
    device = torch_device("cpu")

    print("Initializing ArmNet!", file = stderr)
    net = CNN_model.ArmNet(CNN_model.ArmBlock, layers = [3, 8, 36, 3], channels = 5, num_classes = 2).to(device)
    net.load_state_dict(load(args.model))
    net.eval()

    print("Extracting features!", file = stderr)
    for line in open(args.input):
        chrom, pos, alt, tumor_bam, control_bam = line.strip().split()
        features = model_features(tumor_bam, control_bam, chrom, int(pos), alt, args.fasta)
        features = stack([features]).float()
        with no_grad():
            outputs = net(features).to(device)
            _, predicted = torch_max(outputs, 1)
        print(chrom, pos, alt, predicted.item(), sep = "\t")

if __name__ == "__main__":
    main()