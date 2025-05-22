import matplotlib.pyplot as plt
import numpy as np


def plot_performance(results, algorithms, file_sizes):
    file_sizes = np.array(file_sizes)/(1024*1024)
    sorted_indices = np.argsort(file_sizes)
    results = results[:,sorted_indices,:]
    file_sizes = file_sizes[sorted_indices]

    fig, (ax0, ax1) = plt.subplots(1,2, figsize=(10, 4))
    for i in range(results.shape[0]):
        ax0.plot(file_sizes, results[i,:,0], '-x', label=algorithms[i])
        ax0.set_title("Runtime")
        ax0.set_ylabel("Time [s]")
        ax0.set_xlabel("Filesize [Mb]")
    ax0.legend()
    for i in range(results.shape[0]):
        ax1.plot(file_sizes, results[i,:,1], '-x', label=algorithms[i])
        ax1.set_title("Memory")
        ax1.set_ylabel("Size [Mb]")
        ax1.set_xlabel("Filesize [Mb]")
    ax1.legend()
    plt.show()


def plot_accuracy_FM(results, algorithms):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.errorbar(np.arange(results.shape[1]), results[1,:,2], np.abs(results[0,:,2] - results[1,:,2]), fmt='o', linewidth=2, capsize=6, label="Flajolet Martin")
    ax.set_title("Accuracy")
    ax.set_ylabel("Error")
    ax.set_xlabel("Dataset")
    ax.set_xticks(np.arange(results.shape[1]))
    # ax.legend()
    plt.show()


def plot_accuracies(results, algorithms):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.errorbar(np.arange(results.shape[1])-0.1, results[1,:,2], np.abs(results[0,:,2] - results[1,:,2]), fmt='o', linewidth=2, capsize=6, label="FM")
    ax.errorbar(np.arange(results.shape[1])+0.1, results[2,:,2], np.abs(results[0,:,2] - results[2,:,2]), fmt='o', linewidth=2, capsize=6, label="HLL")
    ax.set_title("Accuracy")
    ax.set_ylabel("Error")
    ax.set_xlabel("Dataset")
    ax.set_xticks(np.arange(results.shape[1]))
    ax.legend()
    plt.show()


def plot_accuracy_runtime(errors, runtimes, probabilities):
    fig, (ax0, ax1) = plt.subplots(1,2, figsize=(10, 4))
    ax0.plot(probabilities, errors)
    ax0.set_title("Accuracy")
    ax0.set_ylabel("Error")
    ax0.set_xlabel("Sample Ratio")
    ax0.set_ylim(0, 1)
    ax1.plot(probabilities, runtimes)
    ax1.set_title("Runtime")
    ax1.set_ylabel("Time [s]")
    ax1.set_xlabel("Sample Ratio")

def plot_minhash(errors, runtimes, permutations):
    fig, (ax0, ax1) = plt.subplots(1,2, figsize=(8, 3.5))
    ax0.plot(permutations, errors)
    ax0.set_title("Accuracy")
    ax0.set_ylabel("Error")
    ax0.set_xlabel("Number permutations")
    ax0.set_ylim(0, 1)
    ax1.plot(permutations, runtimes)
    ax1.set_title("Runtime")
    ax1.set_ylabel("Time [s]")
    ax1.set_xlabel("Number permutations")


def plot_minhashs(hashfunctions, errors, runtimes, permutations):
    fig, (ax0, ax1) = plt.subplots(1,2, figsize=(8, 3.5))
    for hashfunction in hashfunctions:
        ax0.plot(permutations, errors[hashfunction], label=hashfunction)
    ax0.set_title("Accuracy")
    ax0.set_ylabel("Error")
    ax0.set_xlabel("Number permutations")
    ax0.set_ylim(0, 1)
    ax0.legend()
    for hashfunction in hashfunctions:
        ax1.plot(permutations, runtimes[hashfunction], label=hashfunction)
    ax1.set_title("Runtime")
    ax1.set_ylabel("Time [s]")
    ax1.set_xlabel("Number permutations")
    ax1.legend()


def plot_fracminhash(errors, runtimes, scaling_factors):
    fig, (ax0, ax1) = plt.subplots(1,2, figsize=(8, 3.5))
    ax0.plot(scaling_factors, errors)
    ax0.set_title("Accuracy")
    ax0.set_ylabel("Error [%]")
    ax0.set_xlabel("scaling factor")
    # ax0.set_ylim(0, 1)
    ax1.plot(scaling_factors, runtimes)
    ax1.set_title("Runtime")
    ax1.set_ylabel("Time [s]")
    ax1.set_xlabel("scaling factor")


    