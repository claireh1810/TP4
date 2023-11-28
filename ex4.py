# Remplissez les "..." ci-dessous pour préciser les complexités des opérations
# d'addition en fonction du type des matrices, du nombre de non-zéros de la première
# matrice, noté n1, et du nombre de non-zéros de la seconde matrice, noté n2.

# DOK + DOK: O(...)
# DOK + COO: O(...)
# COO + COO: O(...)
# COO + DOK: O(...)

from bisect import bisect_left


class SparseMatrix:
    """a generic class to represent sparse matrices"""

    def __init__(self, m, n):
        if not (type(m) == int and type(n) == int and m > 0 and n > 0):
            raise ValueError("init: wrong matrix dimension")
        self.m = m  # number of rows
        self.n = n  # number of columns
        self.nnz = 0  # number of non zero values

    def check_key(self, key):
        """Check that a key (= a pair of positions or of slices) is correct"""
        if not (len(key) == 2):
            raise ValueError("wrong matrix indexing")
        if type(key[0]) == int and type(key[1]) == int:
            if not (0 <= key[0] < self.m and 0 <= key[1] < self.n):
                raise ValueError("wrong matrix indexing")
            else:
                return int
        elif type(key[0]) == slice and type(key[1]) == slice:
            if key[0].step is not None or key[1].step is not None:
                raise ValueError("slices with steps are not supported")
            if not (
                0 <= key[0].start < self.m
                and 0 <= key[0].stop <= self.m
                and key[0].start < key[0].stop
                and 0 <= key[1].start < self.n
                and 0 <= key[1].stop <= self.n
                and key[1].start < key[1].stop
            ):
                raise ValueError("wrong matrix indexing")
            else:
                return slice
        else:
            raise ValueError("wrong matrix indexing")

    def getnbnz(self):
        """Return the number of non zero values in the matrix"""
        return self.nnz

    def __str__(self):
        """Return a string representation of the matrix"""
        if self.m <= 10 and self.n <= 10:
            repr = ""
            for i in range(self.m):
                for j in range(self.n):
                    repr += f" {self.__getitem__((i,j)):5.2f}"
                if i != self.m - 1:
                    repr += "\n"
        else:
            repr = f"SparseMatrix(m:{self.m},n:{self.n},nnz:{self.nnz})"
        return repr


class SparseDOK(SparseMatrix):
    """a class for sparse matrices using the dictionary of keys (dok) format"""

    def __init__(self, m, n, mat=None):
        super().__init__(m, n)
        self.dict = {}
        if mat:
            for i in range(m):
                for j in range(n):
                    if mat[i][j] != 0:
                        self.dict[(i, j)] = mat[i][j]
            self.nnz = len(self.dict)

    def __getitem__(self, key):
        """Implement M[i,j] and M[rmin:rmax, cmin:cmax]"""
        type_key = self.check_key(key)
        if type_key == int:
            return self.dict.get(key, 0)
        else:
            return self.__getslice(key[0].start, key[0].stop, key[1].start, key[1].stop)

    def __setitem__(self, key, value):
        """Implement assignments like M[i,j]=value"""
        type_key = self.check_key(key)
        if type_key == slice:
            raise ValueError("slices not supported in setitem")
        if value != 0:
            if self.dict.get(key, 0) == 0:
                self.nnz += 1
            self.dict[key] = value
        elif key in self.dict:
            self.dict.pop(key)
            self.nnz -= 1

    # TO BE COMPLETED
    def __add__(self, other):
        """
        Implement the addition of two sparse matrices.

        Parameters:
        - self: SparseDOK
            The current matrix
        - other: SparseMatrix
            The matrix to be added to the current matrix that could be either of type
            SparseDOK or SparseCOO.

        Returns:
        - SparseMatrix
            The result of the addition, i.e., a newly created SparseDOK matrix.

        Raises:
        - ValueError: If the dimensions of the matrices are incompatible or if the type
                      of the other matrix is not supported.
        """
        if issubclass(type(other), SparseMatrix):
            if self.m != other.m or self.n != other.n:
                raise ValueError("Addition: incompatible matrix dimensions")
            if type(other) is SparseDOK:
                return self.__add_with_dok(other)
            elif type(other) is SparseCOO:
                return self.__add_with_coo(other)
            else:
                raise ValueError(
                    "Addition of SparseDOK matrix with SparseMatrix"
                    f"of type {type(other)} is impossible"
                )
        else:
            raise ValueError(
                "Addition of SparseDOK matrix with object of"
                f"type {type(other)} is impossible"
            )

    def __addvalindok(self, key, value):
        if key in self.dict:
            newv = self.dict[key] + value
            if newv == 0:
                self.dict.pop(key)
                self.nnz -= 1
            else:
                self.dict[key] = newv
        else:
            self.dict[key] = value
            self.nnz += 1

    def __add_with_dok(self, other):
        newsm = SparseDOK(self.m, self.n)
        newsm.dict = self.dict.copy()
        for key, val in other.dict.items():
            newsm.__addvalindok(key, val)
        newsm.nnz = len(newsm.dict)
        return newsm
    
    def __add_with_coo(self, other):
        newsm = SparseDOK(self.m, self.n)
        newsm.dict = self.dict.copy()
        for i in range(other.nnz):
            newsm.__addvalindok(other.keys[i], other.values[i])
        newsm.nnz = len(newsm.dict)
        return newsm

    # TO BE COMPLETED
    def __getslice(self, rmin, rmax, cmin, cmax):
        """
        Extracts a submatrix from the current matrix. Equivalent to 
        M[rmin:rmax, cmin:cmax], with M the current matrix.

        Parameters:
        - rmin (int): The minimum row index.
        - rmax (int): The maximum row index (excluded).
        - cmin (int): The minimum column index.
        - cmax (int): The maximum column index (excluded).

        Returns:
        - submatrix: The extracted submatrix, a newly created SparseDOK matrix.

        """
        print("Not implemented yet")
        return None


class SparseCOO(SparseMatrix):
    """a class for sparse matrices using the (sorted) coordinate format"""

    def __init__(self, m, n, mat=None):
        super().__init__(m, n)
        self.keys = []  # ordered list of (i,j) positions of non zero values
        self.values = []  # the associated values
        if mat:
            for i in range(m):
                for j in range(n):
                    if mat[i][j] != 0:
                        self.keys.append((i, j))
                        self.values.append(mat[i][j])
            self.nnz = len(self.keys)

    def __getkeyindex(self, key):
        return bisect_left(
            self.keys, key[0] * self.n + key[1], key=lambda k: k[0] * self.n + k[1]
        )

    def __getitem__(self, key):
        """Implement M[i,j] and M[rmin:rmax, cmin:cmax]"""
        type_key = self.check_key(key)
        if type_key == int:
            i = self.__getkeyindex(key)
            if i != len(self.keys) and self.keys[i] == key:
                return self.values[i]
            else:
                return 0
        else:
            return self.__getslice(key[0].start, key[0].stop, key[1].start, key[1].stop)

    def __setitem__(self, key, value):
        """Implement assignations like M[i,j]=value"""
        type_key = self.check_key(key)
        if type_key == slice:
            raise ValueError("slices not supported in setitem")
        i = self.__getkeyindex(key)
        if i != len(self.keys) and self.keys[i] == key:
            if value == 0:
                self.keys.pop(i)
                self.values.pop(i)
                self.nnz -= 1
            else:
                self.values[i] = value
        else:
            self.keys.insert(i, key)
            self.values.insert(i, value)
            self.nnz += 1

    # TO BE COMPLETED
    def __add__(self, other):
        """
        Implement the addition of two sparse matrices.

        Parameters:
        - self: SparseCOO
            The current matrix
        - other: SparseMatrix
            The matrix to be added to the current matrix that could be either of type 
            SparseDOK or SparseCOO.

        Returns:
        - SparseMatrix
            The result of the addition, i.e., a newly created SparseCOO matrix.

        Raises:
        - ValueError: If the dimensions of the matrices are incompatible or if the type
                      of the other matrix is not supported.
        """
        print("Not implemented yet")
        return None

    # TO BE COMPLETED
    def __getslice(self, rmin, rmax, cmin, cmax):
        """
        Extracts a submatrix from the current matrix. Equivalent to 
        M[rmin:rmax, cmin:cmax], with M the current matrix.

        Parameters:
        - rmin (int): The minimum row index.
        - rmax (int): The maximum row index (excluded).
        - cmin (int): The minimum column index.
        - cmax (int): The maximum column index (excluded).

        Returns:
        - submatrix: The extracted submatrix, a newly created SparseCOO matrix.

        """
        print("Not implemented yet")
        return None


if __name__ == "__main__":
    print("Test DOK+DOK")
    print("------------")

    m1dok = SparseDOK(2, 3, [[4.0, 0.0, 1.0], [0.0, -3.5, -1.0]])
    m1dok[1, 1] = -3.5
    m2dok = SparseDOK(2, 3, [[0.0, -1.0, -1.0], [0.0, 4.0, 0.0]])

    print(f"m1dok (nnz={m1dok.getnbnz()}):")
    print(m1dok)

    print(f"m2dok (nnz={m2dok.getnbnz()}):")
    print(m2dok)

    sum12dok = m1dok + m2dok
    print(f"m1dok+m2dok (nnz={sum12dok.getnbnz()}):")
    print(sum12dok)

    sum21dok = m2dok + m1dok
    print(f"m2dok+m1dok (nnz={sum21dok.getnbnz()})")
    print(sum21dok)

    print("Test COO+COO")
    print("------------")

    m1coo = SparseCOO(2, 3, [[4.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
    m1coo[1, 1] = -3.5
    m2coo = SparseCOO(2, 3, [[0.0, -1.0, -1.0], [0.0, 4.0, 0.0]])

    print(f"m1coo (nnz={m1coo.getnbnz()}):")
    print(m1coo)

    print(f"m2coo (nnz={m2coo.getnbnz()}):")
    print(m2coo)

    sum12coo = m1coo + m2coo
    print(f"m1coo+m2coo (nnz={sum12coo.getnbnz()}):")
    print(sum12coo)

    sum21coo = m2coo + m1coo
    print(f"m2coo+m1coo (nnz={sum21coo.getnbnz()}):")
    print(sum21coo)

    print("Test DOK+COO")
    print("------------")

    print(f"m1dok (nnz={m1dok.getnbnz()}):")
    print(m1dok)
    print(f"m2coo (nnz={m2coo.getnbnz()}):")
    print(m2coo)

    sum12dokcoo = m1dok + m2coo
    print(f"m1dok+m2coo (nnz={sum12dokcoo.getnbnz()}):")
    print(sum12dokcoo)

    print(f"m2dok (nnz={m2dok.getnbnz()}):")
    print(m2dok)
    print(f"m1coo (nnz={m1coo.getnbnz()}):")
    print(m1coo)

    sum21dokcoo = m2dok + m1coo
    print(f"m2dok+m1coo (nnz={sum21dokcoo.getnbnz()}):")
    print(sum21dokcoo)

    print("Test COO+DOK")
    print("------------")

    print(f"m1coo (nnz={m1coo.getnbnz()}):")
    print(m1coo)
    print(f"m2dok (nnz={m2dok.getnbnz()}):")
    print(m2dok)

    sum12coodok = m1coo + m2dok
    print(f"m1coo+m2dok (nnz={sum12coodok.getnbnz()}):")
    print(sum12coodok)

    print(f"m2coo (nnz={m2coo.getnbnz()}):")
    print(m2coo)
    print(f"m1dok (nnz={m1dok.getnbnz()}):")
    print(m1dok)

    sum21coodok = m2coo + m1dok
    print(f"m2coo+m1dok (nnz={sum21coodok.getnbnz()}):")
    print(sum21coodok)

    print("Test getslice")
    print("-------------")

    mdok = SparseDOK(4, 4, [[0, 0, 1, 2], [1, 0, -2, 0], [0, 0, 0, 0], [-1, 0, 2, 0]])
    print("mdok:")
    print(mdok)

    print("mdok[1:4, 0:3]")
    print(mdok[1:4, 0:3])

    print("mdok[0:3, 1:4]")
    print(mdok[0:3, 2:4])

    mcoo = SparseCOO(4, 4, [[0, 0, 1, 2], [1, 0, -2, 0], [0, 0, 0, 0], [-1, 0, 2, 0]])
    print("mcoo:")
    print(mcoo)

    print("mcoo[1:4,0:3]")
    print(mcoo[1:4, 0:3])

    print("mcoo[0:3, 1:4]")
    print(mcoo[0:3, 2:4])
