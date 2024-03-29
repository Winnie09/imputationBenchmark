n_epochs_all = None
import dataclasses
import os
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from scvi.models import *
from scvi.dataset import CsvDataset
from scvi.inference import UnsupervisedTrainer
import sys
import torch

gene_dataset = CsvDataset(sys.argv[1], save_path='',new_n_genes=False)

n_epochs=400 if n_epochs_all is None else n_epochs_all
lr=1e-3
use_batches=False
use_cuda=True

vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches)
trainer = UnsupervisedTrainer(vae,
                              gene_dataset,
                              train_size=0.75,
                              use_cuda=use_cuda)

trainer.train(n_epochs=n_epochs, lr=lr)

full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))

imputed_values = full.sequential().imputation()

latent, batch_indices, labels = full.sequential().get_latent()

np.savetxt(sys.argv[2], imputed_values,delimiter=",")
np.savetxt(sys.argv[3], latent, delimiter=",")
