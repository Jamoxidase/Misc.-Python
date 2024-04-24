'''
Python implimentation of PBWT (Positional Burrows-Wheeler Transform) algorithm based
on on Durbin et al. paper. (doi:10.1093/bioinformatics/btu014). This script contains 
functions to construct the reverse prefix sort matrix, construct Y from X, construct 
X from Y, construct common suffix matrix, and efficiently get long matches from a 
set of binary strings from transformed representation.


Given M sequences with N bi-allelic variable sites, an O(NM) algorithm to derive a representation
of the data based on positional prefix arrays is given, which is termed the positional 
Burrows–Wheeler transform (PBWT). On large datasets this compresses with run-length encoding by 
more than a factor of a hundred smaller than using gzip on the raw data. Using this representation 
a method is given to find all maximal haplotype matches within the set in O(NM) time rather than 
O(NM2) as expected from naive pairwise comparison, and also a fast algorithm, empirically 
independent of M given sufficient memory for indexes, to find maximal matches between a new 
sequence and the set.
'''

import numpy as np

# Construct Reverse Prefix Sort Matrix
def constructReversePrefixSortMatrix(X):
    """
    Computes the reverse prefix sort matrix A for a given set of binary strings X. For 
    integer 0<=i<=N we define an ith prefix sort as a lexicographic sort (here 0 precedes 1)
    of the set of ith prefixes: { x[:i] | x in X }. Similarly an ith reverse prefix sort 
    is a lexicographic sort of the set of ith prefixes after each prefix is reversed.

    Parameters:
    - X: A list of binary strings.

    Returns:
    - A: An Mx(N+1) Numpy integer array, where M is the number of strings in X
         and N is the length of each string.


    Let X be a list of M binary strings (over the alphabet { 0, 1 }) each of length 
    N. 

    """
    n = len(X)  # length of set
    k = len(X[0])  # number of biallelic positions
    A = np.empty(shape=(len(X), 1 if len(X) == 0 else len(X[0]) + 1), dtype=int)

    # Initialize the first column of A with range(n)
    A[:, 0] = np.arange(n)

    # Process each position
    for position in range(k):
        a = []  # List for 0s
        b = []  # List for 1s
        # Process each haplotype sample
        for haplotypeIndex in A[:, position]:
            currentHaplotype = X[haplotypeIndex]
            allele = int(currentHaplotype[position])
            if allele == 0:
                a.append(haplotypeIndex)
            else:
                b.append(haplotypeIndex)
        # Update A with the values of a+b in the current position
        A[:, position + 1] = a + b
    return A


# Construct Y from X
def constructYFromX(X):
    """
    Constructs the matrix Y from the matrix X.

    Parameters:
    - X: A list of binary strings.

    Returns:
    - Y: An MxN Numpy integer array, where M is the number of strings in X
         and N is the length of each string.
    """
    A = constructReversePrefixSortMatrix(X)
    M, N = A.shape
    Y = np.empty((M, N-1), dtype=int)

    for i in range(M):
        for j in range(N-1):
            # Extract character from string at index A[i, j]
            Y[i, j] = int(X[A[i, j]][j])

    return Y


# Construct X from Y
def constructXFromY(Y):
    """
    Constructs the list of binary strings X from the matrix Y.

    Parameters:
    - Y: An MxN Numpy integer array.

    Returns:
    - X: A list of binary strings.
    """
    M, N = Y.shape
    X = np.empty(shape=(M, N), dtype=str)

    A_i = [i for i in range(M)]

    for i in range(N):
        for j in range(M):
            X[A_i[j], i] = str(Y[j, i])

        A_i = sorted(A_i, key=lambda k: X[k, i])

    X = ["".join([str(i) for i in j]) for j in X]
    return X


# Construct Common Suffix Matrix
def computeCSL(s1, s2):
    """
    Computes the length of the common suffix between two strings.

    Parameters:
    - s1, s2: Strings for comparison.

    Returns:
    - sl: Length of the common suffix.
    """
    sl = 0
    while sl < len(s1) and sl < len(s2) and s1[-1 - sl] == s2[-1 - sl]:
        sl += 1
    return sl


def constructCommonSuffixMatrix(A, X):
    """
    Constructs the common suffix matrix D for a given reverse prefix sort matrix A and
    set of binary strings X.

    Parameters:
    - A: An Mx(N+1) Numpy integer array.
    - X: A list of binary strings.

    Returns:
    - D: An MxN Numpy integer array representing the common suffix matrix.
    """
    M, N = A.shape
    D = np.zeros(shape=(M, N), dtype=int)

    for i in range(1, M):
        for j in range(N):
            s1 = X[A[i, j]][:j]
            s2 = X[A[i - 1, j]][:j]
            suffix_length = computeCSL(s1, s2)
            D[i, j] = suffix_length

    return D


# Get Long Matches
def getLongMatches(X, minLength):
    """
    Enumerates all long matches between all substrings of X.

    Parameters:
    - X: A list of binary strings.
    - minLength: The minimum length threshold for a long match.

    Yields:
    - Tuple: A tuple containing the indices of strings in X that form a long match
             and the position at which the match ends.
    """
    assert minLength > 0

    A = constructReversePrefixSortMatrix(X)
    D = constructCommonSuffixMatrix(A, X)

    for j in range(1, 0 if len(X) == 0 else len(X[0])):
        b, c = [], []

        for i in range(len(X)):
            if D[i, j] < minLength:
                for x in b:
                    for y in c:
                        yield (x, y, j) if x < y else (y, x, j)
                b, c = [], []

            if X[A[i, j]][j] == '0':
                b.append(A[i, j])
            else:
                c.append(A[i, j])

        for x in b:
            for y in c:
                yield (x, y, j) if x < y else (y, x, j)


# Example usage
X = np.array([
    "101010101010101010",
    "010101010101010101",
    "110011001100110011",
    "001100110011001100",
    "111100001111000011",
    "000011110000111100",
    "111111111111111111",
    "000000000000000000",
    "101101011010110101",
    "010010100101001010",
    "110000001111111100",
    "001111110000000011",
    "100011110000111100",
    "011100001111000011",
    "101111111110000000",
    "010000000001111111"
])
X = np.array(['010101', '110001', '111111', '011110', '000000', '100010', '110001', '010110'])
print('X')
print(X)

# Construct reverse prefix sort matrix
A = constructReversePrefixSortMatrix(X)
print("Reverse Prefix Sort Matrix (A):")
print(A)

# Construct Y from X
Y = constructYFromX(X)
print("\nY constructed from X:")
print(Y)

# Construct X from Y
X_reconstructed = constructXFromY(Y)
print("\nX reconstructed from Y:")
print(X_reconstructed)

# Construct common suffix matrix
D = constructCommonSuffixMatrix(A, X)
print("\nCommon Suffix Matrix (D):")
print(D)

# Get long matches
min_length = 1
long_matches = list(getLongMatches(X_reconstructed, min_length))
print("\nLong Matches (x, y, j):")
print(long_matches)

