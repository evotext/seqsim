# encoding: utf-8

"""
Module with the implementation of the Modified Moving Contracting Window
Pattern Algorithm (MMCWPA)
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

    # the boolean value indicating if a total or partial
    # match was found between subfields Fx and Fy; when
    # a match is found, the variable is used to cascade
    # out of the loops of the function
    match = False

    # the variables where to store the new collections of
    # subfields, if any match is found; if these values
    # are not changed and the empty lists are returned,
    # stringcomp() will break the loop of comparison,
    # calculate the similarity ratio and return its value
    new_f_x, new_f_y = [], []

    # search patterns in all subfields of Fx; the index of
    # the subfield in the list is used for upgrading the
    # list, if a pattern is a found
    for idx_x, sf_x in enumerate(f_x):
        # 'length' stores the length of the sliding window,
        # from full length to a single character
        for length in range(len(sf_x), 0, -1):
            # 'i' stores the starting index of the sliding
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
                        tmp_x = [sf_x[:i], sf_x[i+length:]]
                        tmp_y = [sf_y[:j], sf_y[j+length:]]
                        new_f_x = f_x[:idx_x] + tmp_x + f_x[idx_x+1:]
                        new_f_y = f_y[:idx_y] + tmp_y + f_y[idx_y+1:]

                        ssnc += (2*length)**2

                        match = True

                        break

                    # if the current match was found, end search
                    if match:
                        break

                # if a match was found, end the sliding window
                if match:
                    break

            # if a match was found, end Fx subfield enumeration
            if match:
                break

        # remove any empty subfields due to pattern removal
        new_f_x = [sf for sf in new_f_x if sf]
        new_f_y = [sf for sf in new_f_y if sf]

        return new_f_x, new_f_y, ssnc


