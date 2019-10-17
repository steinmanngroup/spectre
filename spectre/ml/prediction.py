import numpy
import qml
from qml import Compound
import qml.fchl
import qml.math


def representations_from_files(files, max_size, cut_distance):
    """ Constructs FCHL representations directly from files

        basically wraps :func:`compounds_from_files` and :func:`representations_from_compounds`

    :param files:
    :type files: list[str]
    :param max_size: the maximum size of the representation
    :type max_size: int
    :param cut_distance: cutoff distance
    :type cut_distance: float
    :return: The representation of all compounds given as input
    :rtype: numpy.ndarray
    """

    compounds = compounds_from_files(files)
    return representations_from_compounds(compounds, max_size, cut_distance)


def compounds_from_files(files):
    """ Generates a list of :class:`qml.Compound`

    :param files: list of input filenames
    :type: list[str]
    :return: list of qml compounds
    :rtype: list[Compound]
    """
    return [Compound(file) for file in files]


def representations_from_compounds(compounds, max_size, cut_distance):
    """ Generates representations from compounds

        :param compounds: compounds used to generate represenatation vector
        :type compounds: list[Compound]
        :param max_size: the maximum size of the representation
        :type max_size: int
        :param cut_distance: cutoff distance
        :type cut_distance: float
        :return: The representation of all compounds given as input
        :rtype: numpy.ndarray
    """
    for compound in compounds:
        compound.generate_fchl_representation(max_size=max_size, cut_distance=cut_distance)
    return restructure_representation(numpy.array([mol.representation for mol in compounds]))


def restructure_representation(X):
    """ Restructures the representations for when constructing the kernel

        :param X: the representation to restructure
        :type X: numpy.ndarray
        :return: the restructured representation
        :rtype: numpy.ndarray
    """
    number_of_molecules, nat, dim, cut = numpy.shape(X)
    return X.reshape(number_of_molecules * nat, dim, cut)


def atomic_properties(filenames, ml_data):
    alpha_p = ml_data['alpha_p']
    alpha_q = ml_data['alpha_q']
    repr_train = ml_data['representation']
    sigmas = ml_data['sigma']
    max_size = ml_data['max_size'][0]
    cut_distance = ml_data['cut_distance'][0]

    # construct representation from files
    repr_predict = representations_from_files(filenames, max_size, cut_distance)
    kernel_predict = qml.fchl.get_atomic_kernels(repr_train, repr_predict, sigmas, alchemy='off')[0]
    all_pols = alpha_p.dot(kernel_predict)
    all_charges = alpha_q.dot(kernel_predict)
    return all_charges, all_pols, max_size
