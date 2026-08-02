## Description #############################################################################
#
# Miscellaneous functions.
#
############################################################################################

"""
    _get_degree_and_order(P::AbstractMatrix, n_max::Integer, m_max::Integer) -> Int, Int

Return the maximum degree and order to compute the associated Legendre functions given the
matrix `P` and the configuration values `n_max` and `m_max`. If `n_max` or `m_max` are
negative, the dimensions of `P` are used. The returned values are clamped so the
computation always fits the matrix `P`.

# Returns

- `Int`: Maximum degree that must be computed.
- `Int`: Maximum order that must be computed.
"""
function _get_degree_and_order(P::AbstractMatrix, n_max::Integer, m_max::Integer)
    # Convert the inputs to `Int` to ensure type stability.
    n_max = Int(n_max)
    m_max = Int(m_max)

    # Get the size of the matrix.
    rows, cols = size(P)

    # If the order or degree is less than 0, then the user wants to use all the available
    # memory.
    if n_max < 0
        n_max = rows - 1
    end

    if m_max < 0
        m_max = (cols <= rows) ? cols - 1 : n_max
    end

    # Make sure that the degree and order fits the matrix.
    if n_max > rows - 1
        n_max = rows - 1
    end

    if (m_max > cols - 1) || (m_max > n_max)
        m_max = min(cols - 1, n_max)
    end

    return n_max, m_max
end

"""
    _get_degree_and_order(dP::AbstractMatrix, P::AbstractMatrix, n_max::Integer, m_max::Integer) -> Int, Int

Return the maximum degree and order to compute the derivative of the associated Legendre
functions given the matrices `dP` and `P`, and the configuration values `n_max` and
`m_max`. If `n_max` or `m_max` are negative, the dimensions of `dP` are used. The returned
values are clamped so the computation always fits the matrix `dP`. This function throws an
`ArgumentError` if `P` does not have the dimensions required to compute the derivative
with the resulting degree and order.

# Returns

- `Int`: Maximum degree that must be computed.
- `Int`: Maximum order that must be computed.
"""
function _get_degree_and_order(
    dP::AbstractMatrix, P::AbstractMatrix, n_max::Integer, m_max::Integer
)
    # Obtain the maximum degree and order that fits the matrix `dP`.
    n_max, m_max = _get_degree_and_order(dP, n_max, m_max)

    # The algorithm that computes the derivative accesses the terms `P[n, m + 1]` when
    # `m < n`. Hence, `P` must have an additional column if `m_max < n_max`.
    P_rows, P_cols = size(P)
    req_cols = (m_max < n_max) ? m_max + 2 : m_max + 1

    if (P_rows < n_max + 1) || (P_cols < req_cols)
        throw(
            ArgumentError(
                "The matrix `P` must have at least $(n_max + 1) rows and $req_cols columns " *
                "to compute the derivative with degree $n_max and order $m_max, but it has " *
                "$P_rows rows and $P_cols columns.",
            ),
        )
    end

    return n_max, m_max
end
