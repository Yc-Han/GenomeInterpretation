# GenomeInterpretation

This is an ongoing project in collaboration with LMU Munich, BIFO \@
HZI, and Harvard T.H. Chan School of Public Health. The goal is to
develop a set of tools or modelling standards that can enhance the interpretability
of phenotype-genome models. Furthermore, we want to automate the process of biomarker discovery
from genome sequences.

## Background

See `reports/`. In the [`deepG`](https://deepg.de/index.html) package,
we have implemented *Intergrated Gradients (IG)*. However, the results
do not correspond to biological significant features.

We typically train CNN models on one-hot encoded sequences. While IG is
known to work well on images, our results clearly showcase the need to
adjust our approach to the biological context.

## Project Plan

The project consists of three main components:

1.  **Synthetic Data Generation**: we are developing an algorithm that
    generates sequences based on Hidden Markov Models (HMMs) and real bacterial genomes.
    The sequences then
    undergo a mutator to ensure heterogeneity and variability in the data. Finally, a random set of
    homologous genes are inserted randomly to represent the presence of
    known biomarkers.

2.  **GENTOMA**: a permutation-based explanation strategy that utlizes
    the advantages of genetic algorithms to generate masked sequences
    that can efficiently flip the prediction of a model. The randomness
    of the permutation process is controlled by a follow-up neural
    network mask predictor.

3.  **Robust Training**: inspired by [Wang et al.
    2021](https://arxiv.org/abs/2103.11257), we are developing a robust
    model training process for one-hot encoded sequences that can
    enhance interpretability and thus allow the use of gradient-based
    interpretations, such as BIG.

The project is expected to have more progress by May 2025.

If you have any questions or are interested in more details, please
contact me at `yichen.han AT campus.lmu.de`
