"""
Two data imputation functions that combine MAGIC imputation with BBKNN batch correction:

BBKNN creates a neighborhood network and accounts for batch effects by identifying neighbors
in each batch separately. This is currently used for combined visualization of different batches
in a UMAP embedding
MAGIC uses a diffusion-like process to denoise single-cell RNA seq data based on a kNN graph.

These functions calculate the BBKNN neighborhood graph and supply it to MAGIC data imputation
in order to batch-integrate the given spatiomolecular matrix.

Author: Marius Klein, February 2023

"""

from typing import Optional, Union
import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import bbknn
import magic
import graphtools


def bb_magic(
    adata: ad.AnnData,
    batch_key: str = 'batch',
    neighbors_within_batch: int = 3,
    n_pcs: int = 100,
    copy = False,
    layer_from: str = None,
    layer_to: str = None,
    t: int = 3,
    n_jobs: int = None,
    bbknn_kws: dict = None,
    **magic_kws
) -> Optional[ad.AnnData]:
    """Perform data imputation and batch correction on annotated data matrix

    Args:
        adata (ad.AnnData): 
        batch_key (str, optional): Columns of adata.obs that distinguishes separate batches. 
        Defaults to 'batch'.
        neighbors_within_batch (int, optional): How many top neighbours to report for each batch; 
        total number of neighbours in the initial k-nearest-neighbours computation will be this 
        number times the number of batches. This then serves as the basis for the construction of a 
        symmetrical matrix of connectivities. Defaults to 3.
        n_pcs (int, optional): How many dimensions (in case of PCA, principal components) to use in 
        the analysis. Defaults to 100.
        copy (bool, optional): If input adata should be copied and returned. Otherwise it is updated
        inplace. Defaults to False.
        layer_from (str, optional): Layer of adata to take input data from. Takes adata.X if None given. Defaults to None.
        layer_to (str, optional): Layer of adata to save inputed data matrix to. Overwrites adata.X if None given. Defaults to None.
        t (int, optional): power to which the diffusion operator is powered. This sets the level of 
        diffusion. If 'auto', t is selected according to the Procrustes disparity of the diffused 
        data. Defaults to 3.
        n_jobs (int, optional): The number of jobs to use for the computation. If -1 all CPUs are used. 
        If 1 is given, no parallel computing code is used at all, which is useful for debugging. For n_jobs 
        below -1, (n_cpus + 1 + n_jobs) are used. Thus for n_jobs = -2, all CPUs but one are used
        bbknn_kws (dict, optional): Dictionary of keyword arguments passed to bbknn.bbknn. 
        Defaults to None.
        **magic_kws: All other keyword arguments are passed to magic.MAGIC
    Returns:
        Optional[ad.AnnData]: Returns a copy of the supplied adata if copy=True. Otherwise alters
        supplied adata inplace and returns None
    """

    # running pca, bbknn etc. should not influence the original adata
    adata_analysis = adata.copy()
    del adata_analysis.obsm, adata_analysis.obsp
    if layer_from is not None:
        adata_analysis.X = adata.layers[layer_from]

    # initialize kws to hand over to bbknn
    bbknn_kws = {} if bbknn_kws is None else bbknn_kws.copy()
    n_jobs = 1 if n_jobs is None else n_jobs
    
    # run pca and bbknn on input data
    sc.tl.pca(
        adata_analysis,
        n_comps=n_pcs
    )
    bbknn.bbknn(
        adata_analysis,
        batch_key=batch_key,
        neighbors_within_batch = neighbors_within_batch,
        n_pcs = n_pcs,
        **bbknn_kws
    )

    # use precomputed connectivities from bbknn as adjacency-matrix in magic
    graph = graphtools.graphs.TraditionalGraph(
        adata_analysis.obsp['connectivities'].todense(),
        precomputed='adjacency',
        decay=None # TODO: check if removing changes anything
    )
    
    # calculate diffusion operator (magic does this internally when creating the graph)
    graph.build_kernel()
    graph._kernel = np.asarray(graph._kernel)

    # Initialize magic and fit input data with precomputed graph
    magic_op = magic.MAGIC(t=t, n_jobs=n_jobs, **magic_kws)
    magic_op.fit(adata_analysis.X, graph = graph)

    # transform input data 
    imputed_matrix = magic_op.transform(adata_analysis.X)
    
    if copy:
        adata = adata.copy()

    if layer_to is None:
        adata.X = imputed_matrix
    else:
        adata.layers[layer_to] = imputed_matrix

    if copy:
        return adata



def bb_magic_array(
    array: Union[pd.DataFrame, np.array],
    batch_array: Union[np.array, list, pd.Series],
    neighbors_within_batch: int = 3,
    n_pcs: int = 100,
    t: int = 3,
    bbknn_kws: dict = None,
    **magic_kws
) -> Union[pd.DataFrame, np.array]:
    """Perform data imputation and batch correction on dataframe or array
    Args:
        array (Union[pd.DataFrame, np.array]): Pandas dataframe or numpy array with cells
        representing rows and ions/genes representing columns.
        batch_array (Union[np.array, list]): 1-dimensional array or list of the corresponding batch 
        keys of the cells given in array. Length of batch_array has to equal the number of rows of 
        array.
        neighbors_within_batch (int, optional): How many top neighbours to report for each batch; 
        total number of neighbours in the initial k-nearest-neighbours computation will be this 
        number times the number of batches. This then serves as the basis for the construction of a 
        symmetrical matrix of connectivities. Defaults to 3.
        n_pcs (int, optional): How many dimensions (in case of PCA, principal components) to use in 
        the analysis. Defaults to 100.
        t (int, optional): power to which the diffusion operator is powered. This sets the level of 
        diffusion. If 'auto', t is selected according to the Procrustes disparity of the diffused 
        data. Defaults to 3.
        bbknn_kws (dict, optional): Dictionary of keyword arguments passed to bbknn.bbknn. 
        Defaults to None.
        **magic_kws: All other keyword arguments are passed to magic.MAGIC

    Returns:
        Union[pd.DataFrame, np.array]: Array or Dataframe with the inputed and batch-corrected data
        matrix. Type matches the type of the imput array.
    """

    # create dummy adata to run actual function
    adata = ad.AnnData(X=array)
    adata.obs['batch'] = np.array(batch_array)
    
    # perform imputation and batch correction 
    impute_MAGIC_with_BBKNN(
        adata=adata,
        batch_key='batch',
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        t=t,
        bbknn_kws=bbknn_kws,
        **magic_kws
    )

    # return result in type that was supplied
    if isinstance(array, pd.DataFrame):
        return adata.to_df()
    else:
        return adata.X