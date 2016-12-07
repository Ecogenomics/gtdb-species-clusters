###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import itertools

def aai_thresholds(seq1, seq2, max_mismatches, min_matches):
    """Calculate AAI between sequences.

    Mismatches are only calculate across
    positions where both sequences have
    an amino acid. The max_mismatches
    threshold is used to quickly stop
    comparisons between divergent
    sequences.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : float
        Second sequence.
    max_mismatches : int
        Maximum allowed mismatches between sequences.
    min_matches : int
        Minimum required matches between sequences.

    Returns
    -------
    float
        AAI between sequences, or 0 if criteria is not meet.
    """

    mismatches = 0
    matches = 0
    for c1, c2 in itertools.izip(seq1, seq2):
        if c1 == '-' or c2 == '-':
            continue
        elif c1 != c2:
            mismatches += 1
            if mismatches > max_mismatches:
                return 0
        else:
            matches += 1
            
    if matches < min_matches:
        return 0

    aai = float(matches) / (matches + mismatches)
    return aai


def aai(seq1, seq2, threshold):
    """Determine AAI between sequences.

    The identify is only calculate across
    positions where both sequences have
    an amino acid. A threshold value is
    used to allow quick identification
    of highly similar sequences.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : float
        Second sequence.
    threshold : float
        Minimum AAI required for reporting.

    Returns
    -------
    float
        AAI between sequences if it is greater than the
        specified threshold, else None.
    """

    assert len(seq1) == len(seq2)

    max_mismatches = (1.0 - threshold) * len(seq1)

    mismatches = 0
    matches = 0
    for c1, c2 in itertools.izip(seq1, seq2):
        if c1 == '-' or c2 == '-':
            continue
        elif c1 != c2:
            mismatches += 1
            if mismatches >= max_mismatches:
                return None
        else:
            matches += 1

    aai = float(matches) / max(1, (matches + mismatches))
    if aai < threshold:
        return None

    return aai


def aai_test(seq1, seq2, threshold):
    """Test AAI between sequences.

    The identify is only calculate across
    positions where both sequences have
    an amino acid. A threshold value is
    used to allow quick identification
    of highly similar sequences.

    Parameters
    ----------
    seq1 : str
        First sequence.
    seq2 : float
        Second sequence.
    threshold : float
        Minimum AAI required for reporting.

    Returns
    -------
    bool
        True if AAI >= threshold, else False.
    """

    if aai(seq1, seq2, threshold):
        return True

    return False
