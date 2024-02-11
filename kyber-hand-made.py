# playground from https://cryptographycaffe.sandboxaq.com/posts/kyber-01/

import numpy as np
from numpy.polynomial.polynomial import Polynomial

import random


def add_poly(a, b, q):
    """adds two polynomials modulo q

    Args:
        a (list): list of polynom a coefficients
        b (list): list of polynom b coefficients
        q (int): modulo

    Returns:
        list: a + b
    """
    result = [0] * max(len(a), len(b))
    for i in range(max(len(a), len(b))):
        if i < len(a):
            result[i] += a[i]
        if i < len(b):
            result[i] += b[i]
        result[i] %= q
    return result


def inv_poly(a, q):
    """modular negate of the polynom

    Args:
        a (list): list of polynom a coefficients
        q (int): modulo

    Returns:
        list: -a mod q
    """
    return list(map(lambda x: -x % q, a))


def sub_poly(a, b, q):
    """polynomal substraction

    Args:
        a (list): list of polynom a coefficients
        b (list): list of polynom b coefficients
        q (int): modulo

    Returns:
        list: a - b
    """
    return add_poly(a, inv_poly(b, q), q)


def mul_poly_simple(a, b, f, q):
    """schoolbook multiplication a * b

    Every component of a are multiplies b and we sum the components up in
    the end.
    For the final result, not only do we need to take modulo q for all the
    coefficients, we also need to take modulo f (X^n + 1) for the final
    polynomial such that we obtain a result with degree n-1

    Args:
        a (list): list of polynom a coefficients
        b (list): list of polynom b coefficients
        f (list): polynomial modulo
        q (int): integer modulo

    Returns:
        list: a * b
    """
    tmp = [0] * (
        len(a) * 2 - 1
    )  # the product of two degree n polynomial cannot exceed 2n

    # schoolbook multiplication
    for i in range(len(a)):
        # perform a_i * b
        for j in range(len(b)):
            tmp[i + j] += a[i] * b[j]

    # take polynomial modulo f
    # since Kyber's poly modulus is always x^n + 1,
    # we can efficiently compute the remainder
    degree_f = len(f) - 1
    for i in range(degree_f, len(tmp)):
        tmp[i - degree_f] -= tmp[i]
        tmp[i] = 0

    # take coefficients modulo q
    tmp = list(map(lambda x: x % q, tmp))
    return tmp[:degree_f]


def add_vec(v0, v1, q):
    """Vector of polynoms addition

    Args:
        v0 (list): list of polynoms
        v1 (list): list of polynoms
        q (int): modulo

    Returns:
        list: list of polynoms
    """
    assert len(v0) == len(v1)  # sizes need to be the same
    result = []
    for i in range(len(v0)):
        result.append(add_poly(v0[i], v1[i], q))
    return result


def mul_vec_simple(v0, v1, f, q):
    """Vector multiplication (inner product)

    Args:
        v0 (list): list of polynoms
        v1 (list): list of polynoms
        f (list): reduction polynom
        q (int): coefficients modulo

    Returns:
        list: list of polynoms
    """
    assert len(v0) == len(v1)  # sizes need to be the same

    degree_f = len(f) - 1
    result = [0 for i in range(degree_f - 1)]

    # textbook vector inner product
    for i in range(len(v0)):
        result = add_poly(result, mul_poly_simple(v0[i], v1[i], f, q), q)

    return result


def mul_mat_vec_simple(m, a, f, q):
    """Matrix multiplication

    This works exactly like the plain Matrix-Vector 
    Multiplication, where we perform Vector Multiplication
    for every row of the left-hand term with every 
    column of the right-hand term

    Args:
        m (list): matrix of polynoms : two dimensional array of lists
        a (list): list of polynoms
        f (list): reduction polynom
        q (int): coefficient modulo

    Returns:
        list: m*a
    """
    result = []
    # textbook matrix-vector multiplication
    for i in range(len(m)):
        result.append(mul_vec_simple(m[i], a, f, q))

    return result


def transpose(m):
    """Matrix transpose

    Args:
        m (list): matrix

    Returns:
        list: transposed matrix
    """
    result = [[None for i in range(len(m))] for j in range(len(m[0]))]

    for i in range(len(m)):
        for j in range(len(m[0])):
            result[j][i] = m[i][j]
    
    return result

def encrypt(A, t, m_b, f, q, r, e_1, e_2):
    """Encryption

    Args:
        A (list): Public key - Matrix A of polynoms of degree n
        t (list): Public key - Matrix t = A*s + e
        m_b (list): Input message, as a polynomial with coefs in [0,1]
        f (list): reduction polynom
        q (int): coefficient modulo
        r (list): blinding vector
        e_1 (list): error vector
        e_2 (list): error polynom

    Returns:
        _type_: _description_
    """
    half_q = int(q / 2 + 0.5)
    m = list(map(lambda x: x * half_q, m_b))

    u = add_vec(mul_mat_vec_simple(transpose(A), r, f, q), e_1, q)
    v = sub_poly(add_poly(mul_vec_simple(t, r, f, q), e_2, q), m, q)

    return u, v

def decrypt(s, u, v, f, q):
    """Decrypt (u, v) from private key s

    Args:
        s (list): private key, polynomial with small coefficients
        u (list): encrypted message vector
        v (list): encrypted message polynom
        f (list): reduction polynom
        q (int): coefficient modulo

    Returns:
        list: decrypted message
    """
    m_n = sub_poly(v, mul_vec_simple(s, u, f, q), q)

    half_q = int(q / 2 + 0.5)
    def round(val, center, bound):
        dist_center = np.abs(center - val)
        dist_bound = min(val, bound - val)
        return center if dist_center < dist_bound else 0

    m_n = list(map(lambda x: round(x, half_q, q), m_n))
    m_b = list(map(lambda x: x // half_q, m_n))
    
    return m_b

############################
############################
############################
            
np.random.seed(0xdeadbeef)

def sign_extend(poly, degree):
    if len(poly) >= degree:
        return poly

    return [0] * (degree - len(poly))

if 0:
    def test_mul_poly(N, f, q):
        degree_f = len(f) - 1

        for i in range(N):
            a = (np.random.random(degree_f) * q).astype(int)
            b = (np.random.random(degree_f) * q).astype(int)
            
            a_mul_b = mul_poly_simple(a, b, f, q)
            
            # NumPy reference poly multiplication
            # note that we need to convert the coefficients to int and extend the list to match the fixed size of our impl
            a_mul_b_ref = list(map(lambda x: int(x) % q, ((Polynomial(a) * Polynomial(b)) % Polynomial(f)).coef))
            a_mul_b_ref = sign_extend(a_mul_b_ref, degree_f)

            assert(a_mul_b == a_mul_b_ref)

    # test_mul_poly(100, [1, 0, 0, 0, 1], 17)

    def test_mul_vec(N, k, f, q):
        """_summary_

        Args:
            N (_type_): Iterations
            k (_type_): Dimension of the vectors
            f (_type_): Reduction vector
            q (_type_): coefficient modulo
        """
        degree_f = len(f) - 1

        for i in range(N):
            m = (np.random.random([k, k, degree_f]) * q).astype(int)
            a = (np.random.random([k, degree_f]) * q).astype(int)

            m_mul_a = mul_mat_vec_simple(m, a, f, q)

            m_poly = list(map(lambda x: list(map(Polynomial, x)), m))
            a_poly = list(map(Polynomial, a))
            prod = np.dot(m_poly, a_poly)
            m_mul_a_ref = list(map(lambda x: list(map(lambda y: int(y) % q, sign_extend((x % Polynomial(f)).coef, degree_f))), prod))

            assert(m_mul_a == m_mul_a_ref)

    test_mul_vec(1, 2, [1, 0, 0, 0, 1], 17)


def test_enc_dec(N, k, f, q):
    degree_f = len(f) - 1

    A = (np.random.random([k, k, degree_f]) * q).astype(int)
    s = (np.random.random([k, degree_f]) * 3).astype(int) - 1
    e = (np.random.random([k, degree_f]) * 3).astype(int) - 1
    t = add_vec(mul_mat_vec_simple(A, s, f, q), e, q)

    # print(f"A[{k}][{k}]: {A}")

    failed = 0

    for i in range(N):
        m_b = (np.random.random(degree_f) * 2).astype(int)

        r = (np.random.random([k, degree_f]) * 3).astype(int) - 1
        e_1 = (np.random.random([k, degree_f]) * 3).astype(int) - 1
        e_2 = (np.random.random([degree_f]) * 3).astype(int) - 1

        u, v = encrypt(A, t, m_b, f, q, r, e_1, e_2)
        m_b2 = decrypt(s, u, v, f, q)

        if m_b.tolist() != m_b2:
            failed += 1

    print(f"[k={k}, f={f}, q={q}] Test result: {failed}/{N} failed decryption!")

test_enc_dec(100, 2, [1, 0, 0, 0, 1], 17)
test_enc_dec(100, 2, [1, 0, 0, 0, 1], 37)
test_enc_dec(100, 2, [1, 0, 0, 0, 1], 67)
