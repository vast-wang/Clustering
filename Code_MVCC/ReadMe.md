Title: Multi-View Clustering via Concept Factorization with Local Manifold Regularization (MVCC)
--------------------------------------------------------------

The project deals with multi-view clustering problem. The work has been published in IEEE ICDM 2016. If you have any questions, please do not hesitate to contact us. (Email: hwang@my.swjtu.edu.cn)

To run the MVCC algorithm you can simply execute the script: run_MVCC.m

The details of the main folders are given below,
- Datasets: Contains various datasets used for testing the models.
- run_MVCC: Startup for some examples. Here, you can choose the dataset and set the number of iterations.
- MVCC    : Conatins implementation of MVCC.
- bestmap : permute labels of L2 match L1 as good as possible.
- CalcMeasures: get the clustering performance results, including ACC, NMI and F-measures.

Citation:
@inproceedings{Wang:7837980,
author={H. Wang and Y. Yang and T. Li},
booktitle={2016 IEEE 16th International Conference on Data Mining (ICDM)},
title={Multi-view Clustering via Concept Factorization with Local Manifold Regularization},
year={2016},
pages={1245-1250},
doi={10.1109/ICDM.2016.0167},
month={December},
address={Bacerlona, Spain}
}
