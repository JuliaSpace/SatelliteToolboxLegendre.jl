# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Miscellaneous functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Return the maximum degree and order to compute the Legendre associated functions given the
# matrix `P` and the configuration values `n_max` and `m_max`.
function _get_degree_and_order(P::AbstractMatrix, n_max::Integer, m_max::Integer)
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

# Return the maximum degree and order to compute the Legendre associated functions given the
# matrices `dP`, `P`, and the configuration values `n_max` and `m_max`.
function _get_degree_and_order(dP, P, n_max, m_max)
    # Get the size of the matrices.
    Prows,  Pcols  = size(P)
    dProws, dPcols = size(dP)

    rows = min(Prows, dProws)
    cols = min(Pcols, dPcols)

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
