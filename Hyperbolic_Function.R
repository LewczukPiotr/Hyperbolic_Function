
##### Confirmation of the H. Reiber's original approximation
erfc = function(x) {2*pnorm(x*sqrt(2), lower.tail=F)}

za = seq(0, 2, 0.01)
erfc_za = 0.5 * erfc(za)
zb = 1.5 * za
erfc_zb = 0.5 * erfc(zb)

## Reiber's original parameters (1994)
a=c=0.8736
b=0.4087

qb_fit = a/b * sqrt(za^2 + b^2) - c
plot(erfc_za, erfc_zb, type="l", xlab="Q_A", ylab="Q_B", lwd=2)
lines(za, qb_fit, col="red", lwd=2)


### The same with inverse erfc ###
erfcinv = function(x) {qnorm(x/2, lower.tail=F)/sqrt(2)}
x = seq(0, 0.5, 0.01)
y = 0.5 * erfc(1.5 * erfcinv(2*x))
plot(x,y, type="l", xlab="Q_A", ylab="Q_B", lwd=2)
qb_fit_inv = a/b * sqrt(x^2 + b^2) - c
lines(x, qb_fit_inv, col="red", lwd=2)



########## A Simulated Dataset with a constant CV = 10% of y ##########
set.seed(1234)

## True parameters:
a_true = 0.2
b_true = 0.1
c_true = 0.05


x = seq(1e-4, 0.7, length=5e4)
y = a_true/b_true*sqrt(x^2 + b_true^2) - c_true

## introducing Normal noise in the outcome y
y = rnorm(5e4, y, 0.1*y)
plot(x,y, xlab="Simulated x", ylab="Simulated y")


############### Parametrization with Levenberg-Marquardt ###############
a=0.001; b=0.001; c=0.001

converged = F
iter = 0
lbd = 1e-3
epsi = 100

resid = function(a,b,c) {return(y - a/b*sqrt(x^2 + b^2) + c)}

while(!converged) {
  iter = iter + 1
  
  ## elements of the Jacobian mx
  J1 = 1/b * sqrt(x^2 + b^2)
  J2 = -a*x^2*(b^2*sqrt(x^2 + b^2))^(-1)
  J3 = -1
  J = matrix(rbind(J1, J2, J3), byrow=T, ncol=3)
  
  ## updating term, h
  h = solve(t(J) %*% J + diag(lbd, 3)) %*% t(J) %*% resid(a,b,c)
  
  ## LM criterion
  if (norm(resid(a,b,c), "2") > norm(resid(a+h[1], b+h[2], c+h[3]), "2"))
  {lbd = lbd/epsi}
  else {
    while (norm(resid(a,b,c), "2") < norm(resid(a+h[1], b+h[2], c+h[3]), "2"))
    {lbd = lbd*epsi
    h = solve(t(J) %*% J + diag(lbd, 3)) %*% t(J) %*% resid(a,b,c)}}
  
  ## updating of the parameters
  a = a + h[1]
  b = b + h[2]
  c = c + h[3]
  
  ## check for convergence
  converged = norm(h, "2") < 1e-12}


RSS = sum(resid(a,b,c)^2)

## Covariance mx at the convergence
(Cov_mx = RSS/(length(x)-3+1)*solve(t(J) %*% J))

cat("Numb. of iterations:", iter)

## Point Estimates and the St. Errors:
cat("a_hat:", a, "\t se_a_hat:", sqrt(diag(Cov_mx))[1])
cat("b_hat:", b, "\t se_b_hat:", sqrt(diag(Cov_mx))[2])
cat("c_hat:", c, "\t se_c_hat:", sqrt(diag(Cov_mx))[3])

## plotting the estimated curve
xx = seq(1e-4, 0.7, length=1000)
lines(xx, (a/b*sqrt(xx^2 + b^2) - c), col="green", lwd=2)

