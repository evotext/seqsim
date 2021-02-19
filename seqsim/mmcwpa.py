# encoding: utf-8

"""
Implementaiton of the Modified Moving Contracting Window Pattern Algorithm (MMCWPA).

See details as explained in Tresoldi (2017).
"""

def _mmcwpa(f_x, f_y, ssnc):
    """
    An implementation of the Modified Moving Contracting
    Window Pattern Algorithm (MMCWPA) to calculate string
    similarity, returns a list of non-overlapping,
    non-contiguous fields Fx, a list of non-overlapping,
    non-contiguous fields Fy, and a number indicating the
    Sum of the Square of the Number of the same
    Characters. This function is intended to be "private",
    called from the "public" stringcomp() function below.

    @param f_x: A C{list} of C{strings}.
    @param f_y: A C{list} of C{strings}.
    @param ssnc: A C{float}.

    @return: A C{list} of C{strings} with non-overlapping,
             non-contiguous subfields for Fx, a C{list} of
             C{strings} with non-overlapping,
             non-contiguous subfields for Fy, and a C{float}
             with the value of the SSNC collected so far.
    @rtype: C{list} of C{strings}, C{list} of strings,
            C{float}
    """

    # The boolean variable indicating if a total or partial
    # match was found between subfields Fx and Fy; when
    # a match is found, the variable is used to cascade
    # out of the loops of the function (while such coding approach 
    # might be less than elegant, it makes the loop much clearer
    # to understand for other people).
    match = False

    # The variables where to store the new collections of
    # subfields, if any match is found; if these values
    # are not changed and the empty lists are returned,
    # the encapsulating function will break the loop of comparison,
    # calculate the similarity ratio and return its value.
    new_f_x, new_f_y = [], []

    # Search for patterns in all subfields of Fx; the index of
    # the subfield in the list is used for upgrading the
    # list in case a pattern is a found.
    for idx_x, sf_x in enumerate(f_x):
        # `length` stores the length of the sliding window,
        # from full length down to a single character
        for length in range(len(sf_x), 0, -1):
            # `i` stores the starting index of the sliding
            # window in Fx
            for i in range(len(sf_x)-length+1):
                # extract the pattern for matching
                pattern = sf_x[i:i+length]

                # look for the pattern in Fy
                for idx_y, sf_y in enumerate(f_y):
                    # 'j' stores the starting index in Fy; the
                    # Python find() function returns -1 if there
                    # is no match
                    j = sf_y.find(pattern)
                    if j > -1:
                        # the pattern was found; set 'new_fx' and
                        # 'new_fy' to version of 'fx' and 'fy' with
                        # the patterns removed, update the SSNC and
                        # set 'match' as True, in order to cascade
                        # out of the loops
                        new_f_x = \
                            f_x[:idx_x] + \
                            [sf_x[:i], sf_x[i+length:]] + \
                            f_x[idx_x+1:]
                        new_f_y = \
                            f_y[:idx_y] + \
                            [sf_y[:j], sf_y[j+length:]] + \
                            f_y[idx_y+1:]
                        
                        # Remove any empty subfields introduced by
                        # pattern removal and return
                        new_f_x = [sf for sf in new_f_x if sf]
                        new_f_y = [sf for sf in new_f_y if sf]                        

                        # Update the `ssnc` and return
                        ssnc += (2*length)**2
                        return new_f_x, new_f_y, ssnc

    # The loop will only get here if empty strings are provided;
    # this should never happen in such internal function, but we
    # provide the return nonetheless so that the library does not break
    return new_f_x, new_f_y, ssnc


# TODO: add function docstring
def mmcwpa(str_x, str_y):
    len_x, len_y = len(str_x), len(str_y)

    # The internal function operates on lists of strings, as
    # initialized here
    f_x, f_y = [str_x], [str_y]

    ssnc = 0.0
    while True:
        f_x, f_y, ssnc = _mmcwpa(f_x, f_y, ssnc)

        # If any one of the list of substrings was entirely consumed,
        # we can break out and return
        if not len(f_x) or not len(f_y):
            break

    # Normalize the `ssnc` and return
    return (ssnc / ((len_x+len_y)**2.))**0.5
