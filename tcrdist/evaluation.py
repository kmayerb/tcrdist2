import pickle
import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.model_selection import LeavePGroupsOut


class EvalTCR():
    def __init__(self, tcrrep):
        self.tcrrep  = tcrrep
        self.roc_dict = {}
        self._validate_tcrrep_index_cols()
        self._validate_tcrrep_clone_df()

    def __repr__(self):
        return 'tcrdist.evaluation.EvalTCR WITH A TRACTOR BEAM FOCUSED ON:\n\t {}'.format(str(self.tcrrep))

    def _validate_tcrrep_index_cols(self):
        """
        Checks that index columns contain epitope and subject
        """
        assert("subject" in self.tcrrep.index_cols),\
         "EvalTCR needs a TCRrep with index columns 'subject' and 'epitope'"
        assert("epitope" in self.tcrrep.index_cols),\
        "EvalTCR needs a TCRrep with index columns 'subject' and 'epitope'"

    def _validate_tcrrep_clone_df(self):
        """
        Checks that clone_df is present in the TCRrep instance passed to EvalTCR()
        """
        assert(isinstance(self.tcrrep.__dict__['clone_df'], pd.DataFrame )),\
                "EvalTCR needs a TCRrep instance with a valid .clone_df DataFrame attribute"


    def hold_one_subject_knn(self,
                             tcrrep = None,
                             max_k = 100,
                             inverse_distance_weight = True,
                             account_for_prior = True,
                             power = -1.0):
        """
        Parameters
        ----------
        tcrrep  : tcrdist.repertoire.TCRrep
            with clone_df attribute with 'subject' and 'epitope' features
        max_k : int
            integer for the number of nearest neighbor
        inverse_distance_weight : boolean
            set to True to inverse weight contributions of neighbors
        account_for_prior : boolean
            set to True to weight contributions of neighbors based on their prior
            probabilities in tr.clone_df.epitope
        power : float
            the exponent used for the inverse distance weighting

        Returns
        -------
        roc_dict : dictionary of pandas DataFrames with roc_parameters

        """
        if tcrrep is None:
            tcrrep  = self.tcrrep
            tr = tcrrep


        # extract classes (epitopes) from TCRrep
        my_cl = np.array(tr.clone_df.epitope)

        # define unique epitopes
        unique_epitopes = pd.Series(my_cl).unique().tolist()

        # extract distances from TCRrep
        my_d = tr.paired_tcrdist
        # extract subjects from TCRrep
        my_sub = np.array(tr.clone_df.subject)

        # define cross-validation partions
        cv_partions = self._cv_split_hold_out_by_subject_using_sklearn(tr)
        # apply over cross-valdidation folds
        knn_pred = [ self.knn( class_labels = my_cl, d = my_d, cv = x, k = max_k) for x in cv_partions ]
        # extract cv-fold results two knn DataFrames with classes of neighbors and distance to neighbots
        knn_pred_class = pd.concat([x['knn_class'] for x in knn_pred])
        knn_pred_dist  = pd.concat([x['knn_dist'] for x in knn_pred])
        # supervised labels
        knn_obsv_class = pd.concat([x['knn_label'] for x in knn_pred])
        knn_pred_class_truth = pd.get_dummies(knn_obsv_class, prefix = '', prefix_sep='')




        # calc_inverse_distance
        idw = self._calc_inverse_distance(knn_pred_dist, k = max_k, power = -1)

        if inverse_distance_weight:
            knn_pred_class_prob =  self._calc_idw_weighted_score(knn_pred_dist,
                                                        knn_pred_class,
                                                        k = max_k,
                                                        power = power,
                                                        classes = unique_epitopes)
        else:
            # If distance is ignored probs are just frequency
            knn_pred_class_prob = self._calc_frequency(knn_pred_class, k = max_k)



        if account_for_prior:
            # prior probabilities
            priors = self._calc_prior_probabilities(knn_obsv_class[0])

            # map priors to a DF of same dimensions as knn_pred_class_prob
            priors = pd.DataFrame(priors, index = knn_pred_class_prob.index.copy() )

            # adjust probabilities by prior probability
            knn_pred_class_prob2 = knn_pred_class_prob / priors

            # normalize over the row to force between [0,1]
            knn_pred_class_posterier = knn_pred_class_prob2.div(knn_pred_class_prob2.sum(axis=1), axis=0)

        else:
            knn_pred_class_posterier = knn_pred_class_prob



        roc_dict = {k:{} for k in unique_epitopes}

        for epitope in unique_epitopes:
            roc_dict[epitope] = self.roc_pddf(epitope,
                                        obs = knn_pred_class_truth,
                                        pr = knn_pred_class_posterier)

        return roc_dict




    def _cv_split_hold_out_by_subject_using_sklearn(self, tcrrep = None):
        """
        returns a generator with train and test set indices based on hold on
        subject out cross-validation. This is based on the LeavePGroupsOut


        Parameters
        ----------
        tcrrep : TCRrep class instance
            TCRrep class instance, with TCRrep.clone_df.subject and TCRrep.clone_df.epitope fields

        Returns
        -------
        partitions : generator object BaseCrossValidator.split from sklearn

        """
        if tcrrep  is None:
            tcrrep  = self.tcrrep
        # unique epitope mapped to unique numbers
        encoder_epitope = preprocessing.LabelEncoder()
        encoder_epitope.fit(list(tcrrep.clone_df.epitope.unique()))

        # `y` target vector
        y = encoder_epitope.transform(tcrrep.clone_df.epitope)

        # `X` distance matrix (metric = 'precomputed')
        X = tcrrep.paired_tcrdist

        # Cross Validation Split
        # unique subjects mapped to unique numbers
        encoder_subjects = preprocessing.LabelEncoder()
        encoder_subjects = encoder_subjects.fit(list(tcrrep.clone_df.subject.unique()))

        # define groups based on subject
        groups = list(encoder_subjects.transform(tcrrep.clone_df.subject))

        # Leave P Groups Out
        lpgo = LeavePGroupsOut(n_groups=1)
        lpgo.get_n_splits(X, y, groups)
        partitions  = lpgo.split(X, y, groups)
        return partitions



    def _cv_split_hold_out_by_subject(self, tcrrep = None):
        """
        returns a generator with train and test set indices based on hold on
        subject out cross-validation.


        Parameters
        ----------
        tcrrep : TCRrep class instance
            TCRrep class instance, with TCRrep.clone_df.subject and TCRrep.clone_df.epitope fields

        Returns
        -------
        partitions : generator object BaseCrossValidator.split from sklearn

        """
        if tcrrep is None:
            tcrrep  = self.tcrrep

        subs = pd.Series(tcrrep.clone_df.subject)
        subs_index = pd.Categorical(subs)
        partitions = list()
        for i in range(len(subs_index.categories)):
            test_index = np.where(subs_index.codes == i)
            train_index = np.where(subs_index.codes != i)
            partitions.append([train_index[0], test_index[0]])
        return partitions


    def _train_test(self, chunks):
        """
        Get training and testing CV indicies [(train_index, test_index)] from a set of chunks

        Parameters
        ----------
        chunks : list
            list of list containing indcies (i.e., [[1, 2], [3, 4], [5, 6]])

        Returns
        -------
        r : list
            list of tuples

        Examples
        --------
        >>> _train_test([[1, 2], [3, 4], [5, 6])
        [([3, 4, 5, 6], [1, 2]), ([1, 2, 5, 6], [3, 4]), ([1, 2, 3, 4], [5, 6])]

        """
        r = list()
        for i in range(len(chunks)):
            train = _flatten(chunks[:i] + chunks[(i+1):])
            test = chunks[i]
            r.append((train, test))
        return(r)

    def _flatten(self, l):
        """
        return([item for sublist in l for item in sublist])

        """

        return([item for sublist in l for item in sublist])

    def _partition(self, l, n):
        """
        Function that takes a list and maximum number of elements,
        and break the list into sublists of length n

        Parameters
        ----------
        l : list

        n : int
            size of the sublist
        Returns
        -------
        r : list of lists

        Examples
        --------
        >>> _partition([1,2,3,4,5,6], 2)
        [[1, 2], [3, 4], [5, 6]]
        >>> _partition([1,2,3,4,5,6], 3)
        [[1, 2, 3], [4, 5, 6]]
            """
        n = int(n)
        r = [l[i:i + n] for i in range(0, len(l), n)]
        return r

    def _slice_partition(self, l, n):
        """
        Function that takes a list and maximum number of chunks
        and splits puting the 1st index in first partition and second
        index in the second partition. Achieves better randomization than
        _partition

        Parameters
        ----------
        l : list

        n : int
            number of sublists
        Returns
        -------
        r : list of lists

        Examples
        --------
        >>> _slice_partition([1,2,3,4,5,6], 2)
        [[1, 3, 5], [2, 4, 6]]

        >>> _slice_partition([1,2,3,4,5,6], 3)
        [[1, 4], [2, 5], [3, 6]]

        """
        n = int(n)
        r =  [l[i::n] for i in range(n)]
        return r


    def knn(self, class_labels, d, cv, k ):
        """
        k-nearest neighbor classifer enabled for cross validation

        Parameters
        ----------
        class_labels : list
            list of class labels of lenght (n) in the same order as rows and columns (n x n)
            square distance_matrix

        d : np.array
            square distance matrix (n x n)

        cv : list
            list of tuples containing (train_index, test_index) which can be produced
            by _cv_split_hold_out_by_subject_using_sklearn or other methods (see Notes)

        k : number of k nearest neighbors to consider

        Returns
        -------
        r : dictionary
            dictionary with pandas DataFrames,
            'knn_class' - (n x k) classes labels of k nearest neighobrs
            'knn_dist'  - (n x k) distances to k nearest neighobrs ,
            'knn_label' - (n x 1) supervised label

        Notes
        -----
        Multiple suggestsed methods for generating cv


        _cv_split_hold_out_by_subject_using_sklearn()
        # or #
        _cv_split_hold_out_by_subject()
        # or #
        parts = _partition(list(range(len(my_cl))), 50)
        cv    = _train_test(part
        # or #
        parts =  _slice_partition(list(range(len(my_cl))), 50)
        cv    = _train_test(part


        """

        train_index, test_index = cv
        d_cv = d[test_index,:][:, train_index]
        class_labels_train = class_labels[train_index]
        class_labels_test  = class_labels[test_index]
        knn_ind   = np.argsort(d_cv, axis = 1)[:,1:k+1]
        knn_dist  = np.sort(d_cv, axis = 1)[:,1:k+1]
        knn_class = class_labels_train[knn_ind]
        r = {'knn_class' : pd.DataFrame(knn_class,         index = test_index),
             'knn_dist'  : pd.DataFrame(knn_dist,          index = test_index),
             'knn_label' : pd.DataFrame(class_labels_test, index = test_index)}
        return r

    def _normalize_rows(self, df):
        """
        row normalize pd.DataFrame

        Parameters
        ----------
        df : DataFrame
            pandas DataFrame

        Returns
        -------
        df_row_norm : DataFrame
            pandas DataFrame that has been row normalized

        Examples
        --------
        >>> print(_normalize_rows(pd.DataFrame({"A":[1,2,3],"B":[1,1,1]})))
                  A         B
        0  0.500000  0.500000
        1  0.666667  0.333333
        2  0.750000  0.250000

        """
        df_row_norm = df.div(df.sum(axis=1), axis=0)
        return df_row_norm


    def _calc_counts(self, df, k):
        """
        calculate the number of times a catagorical value occurs per row


        Parameters
        ----------
        df : DataFrame
            pandas DataFrame

        k : nearest neigbors to consider

        Returns:
        r : DataFrame
            pandas DataFrame with counts of each values applied row wise

        Examples
        -------
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        ... print(_calc_counts(pddf , k = 1 ))
             A    B
        0  1.0  0.0
        1  0.0  1.0
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        ... print(_calc_counts(pddf , k = 2 ))
           A    B
        0  2.0  0.0
        1  0.0  2.0
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        ... print(_calc_counts(pddf , k = 3 ))
             A    B    C
        0  2.0  0.0  1.0
        1  1.0  2.0  0.0
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        print(_calc_counts(pddf , k = 4 ))
            A    B    C    D
        0  3.0  0.0  1.0  0.0
        1  1.0  2.0  0.0  1.0
        """
        dfk = df.iloc[:,0:k].copy()
        r = dfk.apply(pd.Series.value_counts, axis=1).fillna(0)
        return r

    def _calc_frequency(self, df, k):
        """
        Parameters
        ----------
        df : DataFrame
            pandas DataFrame

        Returns
        -------
        r : DataFrame
            pandas DataFrame with counts of each values reported as a frequency

        Examples
        --------
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        ... print(_calc_frequency(pddf , k = 1 ))
             A    B
        0  1.0  0.0
        1  0.0  1.0
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        ... print(_calc_frequency(pddf , k = 2 ))
             A    B
        0  1.0  0.0
        1  0.0  1.0
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        ... print(_calc_frequency(pddf , k = 3 ))
                  A         B         C
        0  0.666667  0.000000  0.333333
        1  0.333333  0.666667  0.000000
        >>> pddf = pd.DataFrame({0:["A","B"],1:["A","B"],2:["C","A"],3:["A","D"]})
        print(_calc_frequency(pddf , k = 4 ))
              A    B     C     D
        0  0.75  0.0  0.25  0.00
        1  0.25  0.5  0.00  0.25

        """
        dfk = df.iloc[:,0:k].copy()
        r = dfk.apply(pd.Series.value_counts, axis=1, normalize=True).fillna(0)
        return r

    def _calc_prior_probabilities(self, series):
        """
        Parameters
        ----------
        series : Series
            pandas Seies of categorical values

        Returns
        -------
        r : dictionary
            dictionary of prior probabilites for each class

        Examples
        --------
        >>> _calc_prior_probabilities(pd.Series(["A","A","B"]))
        {'A': 0.6666666666666666, 'B': 0.3333333333333333}

        """
        r = series.value_counts(normalize = True).to_dict()
        return r


    def _calc_inverse_distance(self, df , k, power = -1.):
        """
        Parameters
        ----------
        df : DataFrame
            pandas DataFrame
        k : int
            intereger for number of neareest neighbors
        power : float
            exponentiation factor for inverse distance weighting where (dist)^power

        Returns
        -------
        idw : DataFrame
            row wise normalized inverse distance weightings

        Example
        -------
        >>> pddf = pd.DataFrame({0:[1.],1:[2.],2:[2.]})
        ... print(_calc_inverse_distance(pddf, k=3))
             0     1     2
        0  0.5  0.25  0.25
        """
        df = df.iloc[:,0:k].copy()
        df = df.pow(power)
        # normalize
        idw = df.div(df.sum(axis=1), axis=0)
        return idw

    def _calc_idw_weighted_score(self, knn_pred_dist, knn_pred_class, k, power, classes):
        """
        Parameters
        ----------
        knn_pred_dist : DataFrame
            pandas DataFrame
        knn_pred_class : DataFrame
            pandas DataFrame
        k : int
            number of nearest neighbors
        power : float
            exponent for inverse distance
        classes : list
            list of unique class labels (e.g., epitopes)

        Returns
        -------
        wsn : DataFrame
            pandas data frame with inverse distance weighted scores (row normalized)


        Examples
        --------
        Given:

        distances
               1      2    3
        1    1.0    2.0  2.0
        2    1.0    1.0  1.0
        3  100.0  100.0  1.0

        classes
            1    2    3
        1   PA   PA   PA
        2  PB1  PB1  M45
        3  PB1  PB1  M45

        >>> k_p_d = pd.DataFrame({"1":{"1":1.0,"2":1.0,"3":100.0},"2":{"1":2.0,"2":1.0,"3":100.0},"3":{"1":2.0,"2":1.0,"3":1.0}})
        ... k_p_c = pd.DataFrame({"1":{"1":"PA","2":"PB1","3":"PB1"},"2":{"1":"PA","2":"PB1","3":"PB1"},"3":{"1":"PA","2":"M45","3":"M45"}})
        ... print(_calc_idw_weighted_score(k_p_d, k_p_c, k = 3, classes = ['M45', 'PA', 'PB1']))
                M45   PA       PB1
        1  0.000000  1.0  0.000000
        2  0.333333  0.0  0.666667
        3  0.961905  0.0  0.038095
        """
        l = []
        idw = self._calc_inverse_distance(knn_pred_dist + 1, k = k, power = power) # the +1 avoids infinite weight when distance = 0 to closest point
        for e in classes:
            weighted_score = (idw * (knn_pred_class.iloc[:,0:k] == e).astype(int))
            l.append(weighted_score.sum(axis = 1))
        ws = pd.concat(l, axis  = 1)
        ws.columns = classes
        wsn = self._normalize_rows(ws)
        return wsn


    def roc_pddf(self, epi, obs, pr):
        """
        Parameters
        ----------
        df : DataFrame
            pandas DataFrame

        Returns
        -------
        roc_df : DataFrame


        """
        fpr, tpr, threshold = roc_curve(obs[epi], pr[epi])
        roc_auc = auc(fpr, tpr)
        roc_df = pd.DataFrame({'fpr': fpr, 'tpr' : tpr, 'thr' : threshold, 'roc_auc': roc_auc})
        return roc_df
