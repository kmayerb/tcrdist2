import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import roc_curve
from sklearn.metrics import auc


class TCRproximity():
    '''
    Percentile Nearest Neighbor distance score:
    Focusing on the TCRs repertoire specific to a certain ("target") epitope,
    the sampling density of the target epitope-specific TCRs nearby each TCR
    is estimated by taking the weighted average distance to a defined percentile
    of nearest neighbours. A small nearest-neighbours distance (NN-distance) score
    indicates that there are many other nearby receptors and hence greater
    local sampling density.

    To calculate the NN-distance of a TCR (target epitope-specific or "other"),
    in respect to the epitope-specific repertoire, only the distances to the
    TCRs in the epitope-specific repertoire are considered. The distances to the
    nearest percentile of epitope-specific TCRs are summed with a weight that
    decreases from nearest to farthest neighbours by their rank.
    Importantly, NN-distance scores are calculated after removing all
    neighbors from the same subject as the receptor being scored
    (effectively a leave-one-subject-out control).
    Or, if a subjects_folds dictionary is provided, NN-distance scores are
    calculated for each subject's TCR only based on subjects not in its fold
    (i.e., based on subjects in its training folds).

    To compute AUROC scores for a NN-distance classifier,
    epitope-specific TCRs (positives) and other receptors (negatives)
    are sorted by NN-distance for the target epitope (normalized
    to values between 0-1); an ROC curve is constructed by plotting sensitivity
    (fractional recovery of epitope-specific receptors, or true positive rate) versus
    1-specificity (fractional recovery of background receptors, or false positive
    rate) as the NN-distance threshold increases; and the area under this ROC curve
    is measured (AUC).
    '''

    def __init__(self, tcrrep, target_epitope, other_epitopes,
                 nn_percentile=10, chain='alpha-beta', cdrs='all', subjects_folds=None):
        '''
        :param tcrrep: TCRrep object with the clone_df TCRs dataframe
                      and TCRdist distances matrix to be analyzed
        :param target_epitope: string (epitope name).
                      For each TCR, an NN-distance score is computed with respect
                      to the target epitope repertoire.
                      Also, TCRs specific to the target epitope will be considered
                      as "positive" examples in the prediction task.
        :param other_epitopes: list of strings (epitope names).
                      list of epitopes to be considered as "other"
                      or "negative" examples in the prediction task.
        :param nn_percentile: int. the percentile of target epitope samples
                      considered for NN-distance score calculation.
                      Negative percentile means take exactly -percentile top n
        :param chain: string. Chain/s distance matrix to consider. Possible values:
                      'alpha', 'beta', 'alpha-beta', 'gamma', 'delta' or 'gamma-delta'.
        :param cdrs:  'cdr3' or 'all'. CDRs distance matrix to consider. Possible values,
                      depending on chain being 'alpha' / 'beta' / 'alpha-beta' / ... :
                      if 'cdr3': consider tcrrep.cdr3_a_aa_pw / tcrrep.cdr3_b_aa_pw / tcrrep.cdr3_a_aa_pw + tcrrep.cdr3_b_aa_pw / ...
                      if 'all': consider tcrrep.dist_a / tcrrep.dist_b / tcrrep.paired_tcrdist / tcrrep.dist_g / tcrrep.dist_d / tcrrep.paired_tcrdist
        :param subjects_folds: a dictionary (keys are ints. Each value is a list of
                      strings - subjects, as described next). For each fold there is a key
                      numbered from 0 to the number_of_folds-1. A key holds a list with
                      all subjects in the fold. For each TCR of a subject in fold x,
                      only subjects in other folds will be used for calculating the TCR's
                      NN-score ("training"). Basically, this is a leave-some-subjects-out
                      cross-validation procedure.
                      Default is None. If none - NN-score is calculated after removing all
                      neighbors from the same subject (effectively a leave-one-subject-out
                      cross-validation procedure).
        '''

        self.tcrrep  = tcrrep
        self.target_epitope = target_epitope
        self.other_epitopes = other_epitopes
        self.nn_percentile = nn_percentile
        self.chain = chain
        self.cdrs = cdrs
        self.subjects_folds = subjects_folds
        self.dist_matrix = self._choose_dist_matrix()

        # Get only relevant indices - of target epitope and other epitopes
        indices = tcrrep.clone_df['epitope'].isin([self.target_epitope] + self.other_epitopes)

        # Isolate only relevant TCRs and 'group' other epitopes
        self.prox_df = tcrrep.clone_df.loc[indices, ['epitope', 'subject']].copy()
        self.prox_df['epitope_group'] = self.prox_df['epitope']
        self.prox_df.loc[self.prox_df['epitope'] != self.target_epitope, 'epitope_group'] = 'other'

        # calc NN-distance scores and predicted probabilities
        self._calc_nn_scores()
        self._calc_predictions()


    def _choose_dist_matrix(self):
        '''
        Choose the relevant distance matrix from tcrrep object,
        according to the 'chain' and 'cdrs' defined. (see __init__ documentation)
        '''
        if self.cdrs == 'cdr3':
            if self.chain == 'alpha':
                return self.tcrrep.cdr3_a_aa_pw
            elif self.chain == 'beta':
                return self.tcrrep.cdr3_b_aa_pw
            elif self.chain == 'alpha-beta':
                return self.tcrrep.cdr3_a_aa_pw + self.tcrrep.cdr3_b_aa_pw
            elif self.chain == 'gamma':
                return self.tcrrep.cdr3_g_aa_pw
            elif self.chain == 'delta':
                return self.tcrrep.cdr3_d_aa_pw
            elif self.chain == 'gamma-delta':
                return self.tcrrep.cdr3_g_aa_pw + self.tcrrep.cdr3_d_aa_pw
            else:
                raise ValueError('Dist2Rep: chain value must be alpha / beta / alpha-beta / gamma / delta / gamma-delta.')

        elif self.cdrs == 'all':
            if self.chain == 'alpha':
                return self.tcrrep.dist_a
            elif self.chain == 'beta':
                return self.tcrrep.dist_b
            elif self.chain == 'gamma':
                return self.tcrrep.dist_g
            elif self.chain == 'delta':
                return self.tcrrep.dist_d
            elif self.chain == 'alpha-beta' or self.chain == 'gamma-delta':
                return self.tcrrep.paired_tcrdist
            else:
                raise ValueError('Dist2Rep: chain value must be alpha / beta / alpha-beta / gamma / delta / gamma-delta.')

        else:
            raise ValueError('Dist2Rep: cdrs value must be cdrs or all.')


    def _calc_nn_scores(self):
        '''
        Calculate the NN-distance scores for all TCRs, in respect to the
        target epitope-specific TCRs only. For each TCR, only consider
        TCRs from other subjects (or subjects in other folds, if self.subjects_folds
        was provided) as neighbors. Calculation is made based on the nearest
        neighbors percentile, summing their weighted distances (by order rank).
        '''
        self.prox_df['nn_score'] = None

        # Get indices of target epitope-specific TCRs
        ind_only_epitope = self.prox_df.index[self.prox_df['epitope_group'] == self.target_epitope]

        # For each TCR (target epitope-specific and others), calc nn_score
        for ind_TCR, row in self.prox_df.iterrows():
            # Get distances to non-subject/training folds, target epitope-specific TCRs
            if self.subjects_folds is None:
                ind_non_subject = self.prox_df.index[self.prox_df['subject'] != row['subject']]  # Index hold out all data from that subject
            else:
                subject_fold = self._find_subjects_fold(row['subject'])
                subjects_in_other_folds = self._subjects_in_other_folds_to_list(subject_fold)
                ind_non_subject = self.prox_df.index[self.prox_df['subject'].isin(subjects_in_other_folds)]  # Index use only data from that subject's training set

            ind_epitope_non_subject = ind_only_epitope.intersection(ind_non_subject)
            distances_ep_non_subject = self.dist_matrix[ind_TCR, ind_epitope_non_subject]  # Get Distances from the ith row, holding out subject

            # Calculate weighted average nn-distance score to % of these TCRs
            nn_score = self._sort_and_compute_weighted_nbrdist_from_distances(distances_ep_non_subject)
            self.prox_df.loc[ind_TCR, 'nn_score'] = nn_score


    def _sort_and_compute_weighted_nbrdist_from_distances(self, dist_list):
        '''
        Compute the NN-distance score of a single TCR, to TCRs of the given distances
        (i.e., the epitope-specific TCRs), using the closest percentile of them.

        :param dist_list: list. TCRdist distances to the TCRs (in the epitope-specific
                          repertoire)
        :param percentile: int/float. Percent of TCRs to consider (from the epitope-
                          specific repertoire).
                          negative percentile means take exactly -percentile top n
        :return: float. The NN-distance score of the TCR
        '''
        dist_list.sort()
        assert dist_list[0] <= dist_list[-1]
        if self.nn_percentile < 0:
            n = int(max(1, min(len(dist_list), -1 * self.nn_percentile)))
        else:
            n = int(max(1, (self.nn_percentile * len(dist_list)) / 100))

        total_wt = 0.0
        nbrdist = 0.0
        for i, dist in enumerate(dist_list[:n]):
            wt = 1.0 - float(i) / n
            total_wt += wt
            nbrdist += wt * dist

        return nbrdist / total_wt


    def _calc_predictions(self):
        ''' in prox_df,
        1. create column 'predicted_p':  NN-distance scores normalized to values between 0-1.
                                      (predicted probabilities)
        2. create column 'ground_truth': target epitope specific TCRs are 1s,
                                      other TCRs are 0s. (observed values)

        '''
        min_score = self.prox_df['nn_score'].min()
        max_score = self.prox_df['nn_score'].max()

        self.prox_df['predicted_p'] = 1 - ((self.prox_df['nn_score'] - min_score) /
                                           (max_score - min_score))

        self.prox_df['ground_truth'] = 0
        self.prox_df.loc[self.prox_df['epitope_group'] == self.target_epitope, 'ground_truth'] = 1


    def plot_ROC(self, ax=None):
        '''
        Get receiver operating characteristic (ROC) curve plot for observed
        (prox_df['ground_truth']) and predicted probabilities (prox_df['predicted_p']).

        :param ax: matplotlib axes object. If given, plot will be created over it.

        :return: Dict holding ROC figure axes object, auc score,
                 fpr (false positive rate) and tpr (true positive rate)
        '''

        fpr, tpr, thresholds = roc_curve(self.prox_df['ground_truth'], self.prox_df['predicted_p'])
        roc_auc = auc(fpr, tpr)

        ax = TCRproximity._plot_roc(roc_auc, fpr, tpr,
                                    title='ROC - Epitope {}'.format(self.target_epitope), ax=ax)

        return {'axes': ax, 'auc': roc_auc, 'rates': {'fpr': fpr, 'tpr': tpr}}


    def _find_subjects_fold(self, subject_name):
        '''
        For a given subject, find its fold in the self.subjects_folds dict.
        :param subject_name: string.
        :return: int - fold number
        '''
        for fold in self.subjects_folds.keys():
            if subject_name in self.subjects_folds[fold]:
                return fold


    def _subjects_in_other_folds_to_list(self, fold):
        '''
        For a given fold number, return a list with all subjects
        in all other folds.

        :param fold: int.
        :return: A list of strings - subjects in all other folds.
        '''
        new_list = []
        for check_fold in self.subjects_folds.keys():
            if check_fold != fold:
                new_list += self.subjects_folds[check_fold]
        return new_list


    def plot_NN_score_distribution(self, ax=None):
        '''
        Plot NN-distance score distribution of target epitope TCRs and other TCRs.

        :param ax: matplotlib axes object. If given, plot will be created over it.

        :return: matplotlib axes object with the distribution plot.
        '''
        ax = sns.distplot(self.prox_df.loc[self.prox_df['epitope_group'] == self.target_epitope, 'nn_score'].to_list(),
                          color="blue", label=self.target_epitope, ax=ax)
        sns.distplot(self.prox_df.loc[self.prox_df['epitope_group'] == 'other', 'nn_score'].to_list(),
                     color="red", label="Other", ax=ax)
        ax.legend()
        ax.set_title('NN score distribution - {}'.format(self.target_epitope), color='#800000')
        ax.set_xlabel('NN_score')
        ax.set_ylabel('Frequency')
        return ax


    @staticmethod
    def _plot_roc(auc_val, fpr, tpr, title='ROC curve', ax=None):
        '''
        Plot receiver operating characteristic (ROC) curve.

        :param auc_val: float between 0-1. precalculated AUC value
        :param fpr: list/array of floats. vector of false positive rates
        :param tpr: list/array of floats. vector of true positive rates
        :param title: string. figure title string
        :param ax: matplotlib axes object. If given, plot will be created over it.
        :return: matplotlib axes object containing the ROC plot.
        '''

        if ax is None:
            ax = plt.gca()

        ax.plot([0, 1], [0, 1], 'k--')
        ax.plot(fpr, tpr, label='AUC = {:.3f}'.format(auc_val))
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(title, color='#800000')
        ax.legend(loc='best')

        return ax



