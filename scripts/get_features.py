#!/bin/python3
from features import model_features
from torch import save, stack, tensor, no_grad, nonzero, set_num_threads, squeeze, max as torch_max, device as torch_device
from torch.utils.data import DataLoader, TensorDataset
from torch.utils.data.sampler import SubsetRandomSampler
from torch.nn import CrossEntropyLoss
from torch.optim import SGD
from torch.optim.lr_scheduler import CyclicLR
from torch.multiprocessing import set_sharing_strategy, Pool
from numpy.random import seed, shuffle
from numpy import floor
import CNN_model
from argparse import ArgumentParser
from sys import argv, stderr, exit
from tqdm import tqdm
from itertools import repeat
from mmap import mmap
import random

def arguments():
    parser = ArgumentParser(description='Train the model for Armadillo')
    parser.add_argument(
    '-i', '--input', required = True, metavar = 'file', type = str,
    help = 'TSV file with: CHROM, POS, ALT, GROUP, TUMOR_BAM, CONTROL_BAM columns')
    parser.add_argument(
    '-f', '--fasta', required = True, metavar = 'file', type = str,
    help = 'Path to fasta file')
    parser.add_argument(
    '-t', '--threads', required = False, default = 10, metavar = 'int', type = int,
    help = 'Number of threads [10]')

    if len(argv) < 2:
        parser.print_help()
        exit(0)

    return parser.parse_args()

def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

def run(input):
    line, fasta = input[:2]
    chrom, pos, alt, group, tumor_bam, control_bam = line.strip().split()
    r = model_features(tumor_bam, control_bam, chrom, int(pos), alt, fasta)
    return(r, int(group))
    
def main():
    args = arguments()
    samples, groups = [], []
    validation_fraction = 0.20
    set_sharing_strategy('file_system')
    set_num_threads(args.threads)

    print("Loading data!", file = stderr)
    lines = list(open(args.input))     
    with Pool(processes=args.threads) as pool:
        with tqdm(total=len(lines)) as pbar:
            for _, res in enumerate(pool.imap_unordered(run, zip(lines, repeat(args.fasta)))):
                pbar.update()
                if res[1] == 1:
                    if random.random() < 0.25:
                        continue
                samples.append(res[0].double())
                groups.append(res[1])
    pool.join()
    
    print("Preparing datasets!", len(samples), file = stderr)
    device = torch_device("cpu")
    groups = tensor(groups).to(device)
    samples = stack(samples)
    
    #Split into validation / test sets
    indices = list(range(len(samples)))
    seed(123)
    shuffle(indices)
    split_range = int(floor(validation_fraction * len(samples)))
    train_indices, val_indices = indices[split_range:], indices[:split_range]
    trainingsamp = SubsetRandomSampler(train_indices)
    validationsamp = SubsetRandomSampler(val_indices)
    train = TensorDataset(samples.float(), groups)

    batch_size = 30
    num_classes = 2
    train_loader = DataLoader(train, batch_size=batch_size, sampler=trainingsamp, num_workers=args.threads)
    val_loader = DataLoader(train, batch_size=batch_size, sampler=validationsamp, num_workers=args.threads)
    save(train_loader, "/home/pbousquets/models/training_set.pt")
    save(val_loader, "/home/pbousquets/models/val_set.pt")

    print("Initializing ArmNet!", file = stderr)
    net = CNN_model.ArmNet(CNN_model.ArmBlock, layers = [3, 8, 36, 3], channels = 5, num_classes = num_classes).to(device)
    criterion = CrossEntropyLoss()
    optimizer = SGD(net.parameters(), lr=0.001, momentum=0.9) 
    scheduler = CyclicLR(optimizer, base_lr=0.0000001, max_lr=0.001, step_size_up=3, step_size_down = 7, mode="triangular2")

    epochs = 250
    last_file = f'/home/pbousquets/models/armNet_epoch1.pth'

    print("Training!", file = stderr)
    for epoch in range(epochs):  # loop over the dataset multiple times
        validation_loss = 0.00
        validation_total = 0
        training_loss = 0.00
        training_total, training_correct = 0, 0
        print ('Epoch [{}/{}]...\t'.format(epoch+1, epochs), end = "")

        ## TRAIN ##
        net.train() # Set training mode 
        for _, data in enumerate(train_loader, 0):
            # get the inputs; data is a list of [inputs, labels]
            inputs, labels = data

            # zero the parameter gradients
            optimizer.zero_grad()

            # forward + backward + optimize
            outputs = net(inputs).to(device)

            loss = criterion(outputs, labels).to(device)
            loss.backward()
            optimizer.step()
            # print statistics
            training_loss += loss.item() * len(data)
            _, predicted = torch_max(outputs, 1)
            training_total += labels.size(0)
            training_correct += (predicted == labels).sum().item()

        ## VALIDATE ##
        class_correct = [0.] * num_classes
        class_total = [0.] * num_classes
        class_failed = [0.] * num_classes
        class_missed = [0.] * num_classes
        validation_total, validation_correct = 0, 0
        net.eval() # Set validation mode

        with no_grad():
            for _, data in enumerate(val_loader, 0):
                # get the inputs; data is a list of [inputs, labels]
                inputs, labels = data

                optimizer.zero_grad()
                outputs = net(inputs).to(device)
                loss = criterion(outputs, labels).to(device)

                # print statistics
                validation_loss += loss.item() * len(data)
                _, predicted = torch_max(outputs, 1)
                validation_total += labels.size(0)
                validation_correct += (predicted == labels).sum().item()

                for each in range(num_classes):
                    idx = (labels == each).nonzero()
                    T = len(labels[idx])
                    t = (predicted == each).sum()
                    c = (predicted[idx] == each).sum().item()
                    class_total[each] += T
                    class_correct[each] += c
                    class_failed[each] += t - c
                    class_missed[each] += T - c

        print ('Tr. loss: {:.4f} \t Tr. accuracy: {:.2f}% \t Val. loss: {:.4f} \t Val. accuracy: {:.2f}%'.format(training_loss/len(train_loader), 100*training_correct/training_total, validation_loss/len(val_loader), 100*validation_correct/validation_total))
        scheduler.step()
        for each in range(num_classes):
            class_sens = 100 * class_correct[each] / class_total[each]
            print(f"Class: {each:.0f} \t Sensitivity: {class_sens:.2f} \t True: {class_correct[each]:.0f} \t False: \t {class_failed[each]:.0f} \t FN: \t {class_missed[each]:.0f}")
        
        last_file = f'/home/pbousquets/models/armNet_epoch{epoch}.pth'
        save(net.state_dict(), last_file)

main()

