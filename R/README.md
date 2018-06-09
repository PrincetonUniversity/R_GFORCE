## GFORCE: An R package for high-dimensional clustering and inference in cluster-based graphical models

Author: [Carson Eisenach](http://princeton.edu/~eisenach)

Please send all correspondence to eisenach [AT] princeton.edu.

### Summary
This package provides implementations of state-of-the-art clustering algorithms and inference procedures introduced in
 - Eisenach, C. and Liu, H. (2017). Efficient, Certifiably Optimal High-Dimensional Clustering. [arXiv:1806.00530](https://arxiv.org/abs/1806.00530).
 - Eisenach, C., Bunea, F., Ning, Y. and Dinicu, C. (2018). Efficient, High-Dimensional Inference for Cluster-Based Graphical Models. *Manuscript submitted for publication*.

The new methods implemented include:
 1. FORCE - a fast solver for a semi-definite programming (SDP) relaxation of the K-means problem. For certain data generating distributions it produces a certificate of optimality with high probability, and
 2. Inferential procedures and FDR control for cluster based graphical models.

It also includes high quality implementations of traditional clustering algorithms:
 - Lloyd's algorithm,
 - kmeans++ initializations,
 - hierarchical clustering