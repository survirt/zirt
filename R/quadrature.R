# Gauss-Hermite quadrature for a standard normal density, obtained from the
# eigenvalues of the symmetric tridiagonal Jacobi matrix (Golub and Welsch).
gh_nodes <- function(n) {
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- sqrt(i)
  J[cbind(i + 1L, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = e$vectors[1L, o]^2)
}

# Two-dimensional grid for a bivariate standard normal with correlation rho.
gh_grid2 <- function(n, rho) {
  g <- gh_nodes(n)
  z1 <- rep(g$x, times = n)
  z2 <- rep(g$x, each = n)
  list(t0 = z1, t1 = rho * z1 + sqrt(1 - rho^2) * z2,
       w = rep(g$w, times = n) * rep(g$w, each = n))
}
