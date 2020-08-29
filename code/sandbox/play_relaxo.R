require('relaxo')
require('glmnet')
data(diabetes)
## Center and scale variables
x <- scale(diabetes$x)
y <- scale(diabetes$y)
## Compute "Relaxed Lasso" solution and plot results
object <- relaxo(x,y)
plot(object)
## Compute cross-validated solution with optimal
## predictive performance and print relaxation parameter phi and
## penalty parameter lambda of the found solution
cvobject <- cvrelaxo(x,y) #, foldid=c(rep(1:10,44), 1, 2))
print(cvobject$phi)
print(cvobject$lambda)
## Compute fitted values and plot them versus actual values
fitted.values <- predict(cvobject, newX = x)
plot(fitted.values, y)




# Problem: cannot do cvrelaxo b/c it doesn't allow us to use particular folds :/

# Alternate:
#  cv structure around relaxo? 
#   -OR-
require('glmnet')
#foldid=c(rep(1:10,44), 1, 2)
data(QuickStartExample)
fit <- glmnet(x, y, relax=TRUE)
cvfit$lambda.1se
relaxo(x, y, phi=seq(0,1, length=10))
