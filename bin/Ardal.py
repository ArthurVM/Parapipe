import pandas as pd
import numpy as np
import re
from sys import stdout, stderr
from scipy.spatial.distance import pdist
from collections import defaultdict

class Ardal(object):
    """ Class for handling Ardal-neighbourhood objects. """


    def __init__( self, csv : str, ref : bool=False ):
        """ Ardal constructor
        """
        super(Ardal, self).__init__()

        ## if the csv argument is a dataframe the just assign it as the matrix
        if isinstance(csv, pd.DataFrame):
            self.matrix = csv
        ## otherwise, assume it is the path to a csv and read it in
        else:
            self.matrix = self._parseMatrix(csv, ref)

        self.ref = ref


    def _parseMatrix( self, csv_path : str, ref : bool=False ) -> pd.DataFrame:
        """ parse a sparse presence absence matrix in csv format.
        """
        wg_matrix = pd.read_csv(csv_path, index_col=0)

        if ref:
            wg_matrix.loc["REF"] = 0

        ## check all GUIDs are unique
        if wg_matrix.index.duplicated().any():
            raise KeyError(f"Multiple instances of GUIDs found in {csv_path}.")
        
        return wg_matrix
    

    def subsetbyGUID( self, guid_list : list ):
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        Returns an Ardal object with the subsetted matrix.
        """
        ## Filter rows based on the provided list of indexes
        subset_matrix = self.matrix[self.matrix.index.isin(guid_list)]
    
        ## Remove alleles which do not exist within the subset
        subset_matrix = subset_matrix.loc[:, (subset_matrix != 0).any(axis=0)]

        ## return a new Ardal object
        return Ardal(subset_matrix)

    
    def neighbourhood( self, guid : str, n : int ) -> set:
        """ get the SNP neighbourhood of a GUID
        """
        if guid not in self.matrix.index:
            return f"Variant ID '{guid}' not found in the matrix."

        query_alleles = self.matrix.loc[guid]

        close_variants = set()
        for variant, alleles in self.matrix.iterrows():
            dist = Ardal._absD(query_alleles, alleles)

            if dist <= n:
                close_variants.add(variant)

        return set(close_variants)
    

    def toDict( self ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_mapping = defaultdict(list)

        for variant, alleles in self.matrix.iterrows():
            allele_signature = Ardal._getAlleleSignature(alleles)
            allele_mapping[variant] = allele_signature

        return allele_mapping


    def unique( self, guids : list ) -> set:
        """ Take a set of guids and return alleles unique to this subset.
        """
        intersection, symmetric = self._getIntersectionAndSymmetric(guids)

        return intersection.difference(symmetric)
    

    def common( self, guids : list ) -> set:
        """ Take a set of guids and return alleles common to this subset.
        """
        intersection, symmetric = self._getIntersectionAndSymmetric(guids)

        return intersection
    

    def allele( self, alleles : list ) -> set:
        """ Take a set of alleles and return all GUIDs containing all of those alleles.
        """

        ## check guids are present within the matrix and construct a list of present guids to proceed with
        present_alleles = self._checkAlleles(alleles)
        allele_count = len(present_alleles)

        ## initialise a dictionary for storing present alleles
        allele_dict = { guid : [] for guid in self.matrix.index.values }

        for allele_id in present_alleles:
            mask = self.matrix[allele_id].astype(bool)
            for guid in self.matrix.loc[mask].index.values:
                allele_dict[guid].append(allele_id)

        return set([guid for guid, alleles in allele_dict.items() if len(alleles) == allele_count])


    def old_pairwise( self, guids : list=[], metric : str="euclidean" ) -> pd.DataFrame:
        """ Takes a set of GUIDs and calculates the pairwise distance between them, returning a distance matrix.
        Pairwise distance can be calculated using Jaccard, Euclidean, or absolute distance functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """

        ## check the specified distance function is valid
        accepted_dist_functions = ["absolute", "euclidean", "jaccard"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")

        present_guids = self._checkGUIDs(guids)

        if len(guids) == 0:
            subset_matrix = self.matrix
            present_guids = self.matrix.index.values
        elif len(present_guids) == 0:
            raise KeyError(f"No input GUIDs found within the matrix.")
        elif len(present_guids) == 1:
            raise ValueError(f"Cannot construct a distance matrix using fewer than 2 valid GUIDs.")
        else:
            subset_matrix = self.matrix.loc[guids]

        ## initialise an empty distance matrix
        dist_matrix = np.zeros((len(present_guids), len(present_guids)))

        ## iterate through the subset matrix in an O(n**2) fashion
        ## TODO: optimise this big-O horror show
        for i, (vi, ai) in enumerate(subset_matrix.iterrows()):
            for j, (vj, aj) in enumerate(subset_matrix.iterrows()):
                
                ## calculate distance between two samples using the given distance function
                if metric == "absolute":
                    d = Ardal._absD(ai, aj)
                elif metric == "euclidean":
                    d = Ardal._eD(ai, aj)
                elif metric == "jaccard":
                    d = Ardal._jD(ai, aj)

                dist_matrix[i][j] = d
        
        return pd.DataFrame(dist_matrix, columns=present_guids, index=present_guids)


    def pairwise(self, guids: list = [], metric: str = "euclidean", chunk_size : int=100 ) -> pd.DataFrame:
        """ Takes a set of GUIDs and calculates the pairwise distance between them, returning a distance matrix.
        Pairwise distance can be calculated using Jaccard, Euclidean, or absolute distance functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """

        ## check the specified distance function is valid
        accepted_dist_functions = ["absolute", "euclidean", "jaccard"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")

        present_guids = self._checkGUIDs(guids)

        if len(guids) == 0:
            subset_matrix = self.matrix
            present_guids = self.matrix.index.values
        elif len(present_guids) == 0:
            raise KeyError(f"No input GUIDs found within the matrix.")
        elif len(present_guids) == 1:
            raise ValueError(f"Cannot construct a distance matrix using fewer than 2 valid GUIDs.")
        else:
            subset_matrix = self.matrix.loc[guids]

        ## initialize an empty distance matrix
        dist_matrix = np.zeros((len(present_guids), len(present_guids)))

        ## Calculate the number of chunks
        n_chunks = len(subset_matrix) // chunk_size + 1

        ## Calculate pairwise distances using vectorized operations in chunks
        for chunk_idx in range(n_chunks):

            ## track progress
            # stdout.write(f"Processing chunk {chunk_idx}/{n_chunks}\r")
            # stdout.flush()
            
            start_idx = chunk_idx * chunk_size
            end_idx = min((chunk_idx + 1) * chunk_size, len(subset_matrix))
            chunk_data = subset_matrix.iloc[start_idx:end_idx]

            if metric == "absolute":
                dist_chunk = np.sum(np.abs(chunk_data.values[:, None] - subset_matrix.values), axis=2)
            elif metric == "euclidean":
                dist_chunk = np.sqrt(((chunk_data.values[:, None, :] - subset_matrix.values) ** 2).sum(axis=2))
            elif metric == "jaccard":
                ## TODO: implement your Jaccard distance function using vectorized operations
                
                for guid_i, a_i in chunk_data.iterrows():
                    for guid_j, a_j in chunk_data.iterrows():
                        print(Ardal._absD(a_i, a_j))
                        
                # intersection = chunk_data.values[:, None] & chunk_data.values[:, None]  # Count the number of common 1s
                # union = chunk_data.values[:, None] | chunk_data.values[:, None]        # Count the number of total 1s
                # jaccard_similarity = intersection / union
                # dist_chunk = 1 - jaccard_similarity

            dist_matrix[start_idx:end_idx] = dist_chunk

        return pd.DataFrame(dist_matrix, columns=present_guids, index=present_guids)


    def removeGUID( self, guid : str ):
        """ Removes a sample from the matrix, defined by a unique GUID
        """
        None


    def matchAlleleNames(self, expression: str) -> list:
        """ Return all allele names that match the given expression with wildcards.
        """
        if not isinstance(expression, str):
            raise ValueError("Expression must be a string.")

        pattern = re.compile(expression.replace('*', '.*'))
        return set([allele for allele in self.matrix.columns if pattern.match(allele)])
    
    
    def _getIntersectionAndSymmetric( self, guids : list ) -> list:
        """ Take a set of guids and return both the allelic intersection and the symmetric set of the union.
        """

        present_guids = self._checkGUIDs(guids)

        if len(present_guids) == 0:
            present_guids = list(self.matrix.index)

        ## get allele signature for the first guid to use as a reference
        intersection = set(Ardal._getAlleleSignature(self.matrix.loc[present_guids[0]]))

        ## list to store (A U B) / A
        symmetric = []

        for variant, alleles in self.matrix.iterrows():
            if variant in present_guids:
                intersection.intersection_update(Ardal._getAlleleSignature(alleles))
            else:
                symmetric.append(Ardal._getAlleleSignature(alleles))
        
        symmetric = set([x for xs in symmetric for x in xs])
        
        return intersection, symmetric
    

    def _checkGUIDs( self, guids : list ) -> list:
        """ Check guids are present within the matrix and construct a list of present guids to proceed with
        """

        present_guids = []
        for id in guids:
            if self._checkIndex(id):
                present_guids.append(id)
            else:
                print(f"{id} not present in allele matrix.")

        return present_guids
    

    def _checkAlleles( self, alleles : list ) -> list:
        """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
        """
        present_alleles = []
        for allele_id in alleles:
            if self._checkCol(allele_id):
                present_alleles.append(allele_id)
            else:
                print(f"{allele_id} not present in allele matrix.")

        return present_alleles

    
    @staticmethod
    def _absD( ai : pd.DataFrame, aj : pd.DataFrame ) -> int:
        """ return the absolute SNP distance of two rows
        """
        diff_alleles = Ardal._getAlleleSignature(ai) ^ Ardal._getAlleleSignature(aj)
        return len(diff_alleles)


    @staticmethod
    def _eD( ai : pd.DataFrame, aj : pd.DataFrame ) -> pd.DataFrame:
        """ return the Euclidean distance of two rows
        """
        return pdist([ai, aj], metric='euclidean')


    @staticmethod
    def _jD( ai : pd.DataFrame, aj : pd.DataFrame ) -> pd.DataFrame:
        """ return the Jaccard distance of two rows
        """
        return pdist([ai, aj], metric='jaccard')
    

    def jaccard_distance(X, Y):
        """
        Calculate the Jaccard distance between two sets of binary vectors.

        Parameters:
        - X: 2D numpy array, shape (n_samples, n_features)
            Binary feature matrix for the first set of samples.
        - Y: 2D numpy array, shape (n_samples, n_features)
            Binary feature matrix for the second set of samples.

        Returns:
        - distances: 1D numpy array, shape (n_samples,)
            Jaccard distances between corresponding samples in X and Y.
        """
        intersection = np.sum(X & Y, axis=1)  # Count the number of common 1s
        union = np.sum(X | Y, axis=1)         # Count the number of total 1s
        jaccard_similarity = intersection / union
        distances = 1 - jaccard_similarity    # Convert similarity to distance
        return distances


    @staticmethod
    def _getAlleleSignature( df : pd.DataFrame ) -> set:
        """ Use boolean indexing to filter alleles by presence and return the allele IDs
        """
        return set(df[df.eq(1)].index)
    

    def _checkIndex( self, id : str ) -> bool:
        """ take an index ID and check it exists within the matrix.
        """
        return id in list(self.matrix.index)
    

    def _checkCol( self, id : str ) -> bool:
        """ take a column ID and check it exists within the matrix.
        """
        return id in list(self.matrix.columns.values)