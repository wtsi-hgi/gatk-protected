import math
import numpy

## a python library for the fitting of linear models. No good python3 library exists for this. Lots of special-cases because
## it turns out that's what I need. Also because of speedocity.

def logbinom(k,N):
 """ Calculates the log-binomial coefficient of (N choose K) using stirling's approximation
 """
 if ( k >= N or k <= 0 ): 
  return 1.0
 return (0.5+N)*math.log(N) - (0.5+k)*math.log(k) - (N-k-0.5)*math.log(N-k) - 0.5*math.log(2*math.pi)

class OLS:

 class Fit:

  def __init__(self):
   self.residuals = None
   self.coefficients = None
   self.ssr = None

 def fit(response,predictors):
  """ An ordinary least squares fit using numpy least squares
     @response
     @predictors
  >>> import numpy as np
  >>> x = np.array([0, 1, 2, 3])
  >>> y = np.array([-1, 0.2, 0.9, 2.1])
  >>> A = np.vstack([x, np.ones(len(x))]).T
  >>> l = fit(y,A)
  >>> l.coefficients
  array([ 1.  , -0.95])
  >>> A = numpy.vstack([[0,1,2,3,4,5],[1,0,2,5,4,3],[3,1,2,4,2,4],[1,1,1,1,1,1]]).T
  >>> y = numpy.array([-1,0.2,0.9,2.1,1.8,-0.6])
  >>> l = fit(y,A)
  >>> l.coefficients
  array([-0.21666667,  0.95392157, -0.82647059,  0.92745098])
  """
  # ensure dimensions are good
  assert response.shape[0] == predictors.shape[0]
  regression = numpy.linalg.lstsq(predictors,response)
  beta = regression[0] 
  resid = response - numpy.inner(beta,predictors) 
  fitObj = OLS.Fit()
  fitObj.residuals = resid
  fitObj.coefficients = beta
  return fitObj

 def fitCoefficients(response,predictors):
  """ todo -- comment
      @response
      @predictors
      @reutrn
  """
  assert len(response) == len(predictors)
  return numpy.linalg.lstsq(predictors,response)[0]

class GLM:

 class Logistic:

  class Fit:

   def __init__(self):
    self.residuals = None
    self.coefficients = None
    self.predictedValues = None

   # package GLM.Logistic.Fit
   def gradientDescent(response,predictors,N):
    """ Fit logistic regression via simple gradient descent. This is incredibly slow, and is perhaps
        most useful for the testing of other fitting algorithms.
        @response -
        @predictors -
        @N -
    """
    STEP_SIZE_FIXED = 1/(math.sqrt(N)*math.log(N)*len(response)) # unused
    assert len(response) == len(predictors)
    print("Warning: fitting via gradient descent can be painfully slow for serious applications.")
    # initial guess
    beta = OLS.fitCoefficients(response,predictors)
    beta = beta/numpy.linalg.norm(beta)
    betaChangeSq = 1.0
    iter = 0
    prev_step_size = 1
    while ( betaChangeSq > 0.00001 and iter < 100 ):
     gradient = GLM.Logistic.gradient(response,predictors,beta,N)
     direction = gradient/numpy.linalg.norm(gradient)
     stepSize = GLM.Logistic.lineSearch(response,predictors,beta,N,direction,gradient,prev_step_size) # fixed small size for now rather than bisection (requires recalculating likelihood -- ugh)
     step = stepSize*direction
     beta_new = beta + step
     betaChangeSq = numpy.linalg.norm(step)**2/(len(beta)**2)
     beta = beta_new
     iter += 1
     prev_step_size = stepSize/2
    fitObj = GLM.Logistic.Fit()
    fitObj.coefficients = beta
    fitObj.predictedValues = GLM.Logistic.predict(predictors,beta,N)
    fitObj.residuals = response - fitObj.predictedValues
    return fitObj

   # package GLM.Logistic.Fit
   def simpleNewton(response,predictors,N):
    """ Uses unweighted Newton's method to fit the logistic. Not calculated in the IRLS form.
    """
    STEP_SIZE_FIXED = math.pow(0.95,math.log(response.size)/2)
    assert response.size == len(predictors)
    # initial guess
    beta = OLS.fitCoefficients(response,predictors)
    beta = beta/numpy.linalg.norm(beta)
    betaChangeSq = 1.0
    iter = 0
    while ( betaChangeSq > 0.00001 and iter < 100 ):
     gradient = GLM.Logistic.gradient(response,predictors,beta,N)
     hessian = GLM.Logistic.hessian(response,predictors,beta,N)
     step = -STEP_SIZE_FIXED * numpy.transpose(numpy.inner(numpy.linalg.inv(hessian),gradient))
     beta_new = beta - step
     betaChangeSq = numpy.linalg.norm(step)**2/len(beta)
     beta = beta_new
     iter += 1
    fitObj = GLM.Logistic.Fit()
    fitObj.coefficients = beta
    fitObj.predictedValues = GLM.Logistic.predict(predictors,beta,N)
    fitObj.residuals = response - fitObj.predictedValues 
    return fitObj

   # package GLM.Logistic.Fit
   def newton(response,predictors,N):
    """ Uses Iteratively Reweighted Least Squares to implement Newton's Method (equivalent for L2 norm)
    """
    # initial guess
    beta = numpy.linalg.lstsq(predictors,response)[0] 
    # reweight by the squared residuals
    betaChangeSq = 1.0
    iter = 0
    while ( betaChangeSq > 0.00001 and iter < 100 ):
     print(beta)
     pred = GLM.Logistic.predict(predictors,beta,N)
     pred_P = pred/N
     resid = (response-pred) 
     W = numpy.diag(N*pred_P*(1.0-pred_P)) 
     XW = numpy.inner(numpy.transpose(predictors),W)
     inverted = numpy.linalg.inv(XW*predictors)
     step = resid*predictors*inverted
     beta_new = beta + step
     betaChangeSq = numpy.linalg.norm(step)**2/len(beta)
     beta = beta_new
    fitObj = GLM.Logistic.Fit()
    fitObj.coefficients = beta
    fitObj.predictedValues = GLM.Logistic.predict(predictors,beta,N)
    fitObj.residuals = response - fitObj.predictedValues
    return fitObj

   # package GLM.Logistic.Fit
   def quasiNewton(response,predictors,N):
    """ Uses a BGFS approximation scheme to the Newton update (e.g. approximate hessian) to
        fit the logistic coefficients
    """
    assert false
    beta = OLS.fitCoefficients(response,predictors)
    beta = beta/numpy.linalg.norm(beta) 
    betaChangeSq = 1.0
    prevGrad = None
    prevStepSize = 1.0
    iter = 0
    approxHessianInv = numpy.diag(list(map(lambda i: 1.0,range(predictors[0].size))))
    while ( betaChangeSq > 0.00001 and iter < 1000 ):
     gradient = GLM.Logistic.gradient(response,predictors,beta,N)
     direction = numpy.transpose(numpy.inner(approxHessianInv,gradient))
     direction = direction/numpy.linalg.norm(direction)
     #stepSize = GLM.Logistic.lineSearch(response,predictors,beta,N,direction,gradient,prevStepSize)
     stepSize = 0.99*prevStepSize
     step = stepSize*direction
     beta_new = beta + step
     if ( iter > 0 ):
      gradDiff = gradient - prevGrad 
      stepGradDiffInner = numpy.inner(step,-gradDiff)
      approxHessianInv = approxHessianInv + ( numpy.outer(step,step) - numpy.inner(numpy.outer(step,gradDiff),approxHessianInv) - numpy.inner(approxHessianInv,numpy.outer(step,gradDiff)) )/stepGradDiffInner
     betaChangeSq = numpy.linalg.norm(step)**2/len(beta)
     beta = beta_new
     prevGrad = gradient
     prevStepSize = stepSize
     iter += 1
    fitObj = GLM.Logistic.Fit()
    fitObj.coefficients = beta
    fitObj.predictedValues = GLM.Logistic.predict(predictors,beta,N)
    fitObj.residuals = response - fitObj.predictedValues
    return fitObj
      
   # package GLM.Logistic.Fit
   def conjugateGradient(response,predictors,N):
    """ pretty specific (& hardcoded) method for fitting a logistic GLM
        @response - an nx1 numpy array of responses
        @predictors - an nxk numpy matrix of predictors
        @N - the number of trials (1 for strict logistic, >1 for a binomial logistic)
    """
    STEP_SIZE_FIXED = math.pow(0.95,math.log(response.size)/2)
    # first thing: assert the dimensions are are good
    assert len(response) == len(predictors)
    # second: generate an initial guess from OLS (this works well in practice for genotype data, but
    # the standard initial guess is to throw a p/(1-p) into the slot with the most support)
    beta = OLS.fitCoefficients(response,predictors)
    beta = beta/numpy.linalg.norm(beta)
    # calculate the (variable part) of the likelihood under this guess
    likelihood = GLM.Logistic.likelihood(response,predictors,beta,N)
    # instantiate the "previous" values to defaults. The "old" gradient is just 0.
    previousGradient = numpy.array(list(map(lambda x: 0.0,range(beta.size))))
    previousDirection = numpy.matrix(list(map(lambda x: 1.0,range(beta.size))))
    iter = 0
    betaChangeSq = 1.0
    while ( betaChangeSq > 0.00001 and iter < 100 ): # todo -- maybe this is a parameter?
     gradient = GLM.Logistic.gradient(response,predictors,beta,N)
     if ( iter == 0 ):
      # take a newton step
      gradient = GLM.Logistic.gradient(response,predictors,beta,N)
      hessian = GLM.Logistic.hessian(response,predictors,beta,N)
      step = -STEP_SIZE_FIXED*numpy.transpose(numpy.inner(numpy.linalg.inv(hessian),gradient))
      beta_new = beta - step
      nrm = numpy.linalg.norm(step)
      betaChangeSq = nrm ** 2 
      direction = step/nrm
     else:
      # take a CG step
      conjugateFactor = GLM.Logistic.conjugateFactor(gradient,previousGradient,previousDirection)
      hessian = GLM.Logistic.hessian(response,predictors,beta,N)
      direction = gradient - conjugateFactor*previousDirection
      # normalize the direction
      direction = direction/numpy.linalg.norm(direction)
      # quickly calculate the inner product numpy.transpose(direction)*hessian*direction
      dHd = 0.0
      for predictorVector in predictors:
       dHd += (numpy.inner(direction,predictorVector)**2)*N*GLM.Logistic.linkDerivative(predictorVector,beta)
      delta = numpy.inner(gradient,direction)/(0.0-dHd)
      step = 1.0 * delta * direction
      beta_new = beta - step
      betaChangeSq = numpy.linalg.norm(step)**2
     likelihood_new = GLM.Logistic.likelihood(response,predictors,beta_new,N)[0]
     likelihood = likelihood_new
     beta = beta_new
     previousGradient = gradient
     previousDirection = direction
     iter += 1
    fitObj = GLM.Logistic.Fit()
    fitObj.coefficients = beta
    fitObj.predictedValues = GLM.Logistic.predict(predictors,beta,N)
    fitObj.residuals = response - fitObj.predictedValues
    return fitObj

  # package GLM.Logistic
  def conjugateFactor(grad,prevGrad,prevDir):
   """ Calculates the conjugate factor for conjugate gradient descent via the Hestenes-Steifel formula
       The idea behind conjugate gradient descent is that your next guess is a (normalized) projection of the gradient
       Onto a direction that is "more orthogonal" to the last step, that is
       dir_new = grad - dir_old*FACTOR
       this method calculates FACTOR
       @grad - the gradient of the likelihood function just calculated
       @prevGrad - the previous gradient of the likelihood function
       @prevDir - the previous direction of conjugate gradient descent
       @return - the Hestenes-Steifel formula for the step factor (inner(g,g-prev_g)/inner(prev_d,g-prev_g))
   """
   dif = grad - prevGrad
   return numpy.inner(grad,dif)/numpy.inner(prevDir,dif)

  # package GLM.Logistic
  def lineSearch(y,x,beta,N,direction,gradient,init_step):
   """ Performs a line search of the logistic likelihood FOR THE PARAMETER BETA in the direction direction
       using the Weak Wolfe conditions. The gradient at beta is passed in (pre-computed) and assumed to be true.
   """
   c1 = 0.0001
   c2 = 0.8
   lik0 = GLM.Logistic.likelihood(y,x,beta,N)
   grad0 = numpy.inner(gradient,direction)
   p = init_step
   q = 0
   r = None
   iter = 0
   while ( iter < 5 ):
    testBeta = beta + p*direction
    testLik = GLM.Logistic.likelihood(y,x,testBeta,N)
    if ( testLik < lik0 - p*c1*grad0 ):
     r = p
    else:
     testGrad = GLM.Logistic.gradient(y,x,testBeta,N)
     if ( abs(numpy.inner(testGrad,direction)) > abs(c2*grad0) ):
      q = p
     else:
      break
    if ( r != None ):
     p = (q+r)/2
    else:
     p = 2*p
    iter += 1
   return p

  # package GLM.Logistic
  def likelihood(y,x,beta,N):
   """ Computes the (un-optimizable) part of the logistic likelihood.
       The likelihood is given by
       P[y|x] = prod_i (N choose y_i) (p(x_i))^y_i(1-p(x_i))^(1-y_i)
       taking a log and ignoring the N choose y_i term
       L[y|x] = sum_i y_iLog[p(x_i)] + (N-y_i)Log[1-p(x_i)]
       which plugging in for the logit expression simplifies to
       L[y|x] = sum_i y_i dot(x,beta) - N log(1+exp(dot(x,beta)))
       @y - the response
       @x - the predictors
       @N - the number of trials (for binomial LR)
       @beta - the coefficients
   """
   # todo - there's probably a slightly faster way to do this using the fact that for large z, log(1+e^z) ~ z - c for some c
   # note that math.log and math.exp are significantly faster than numpy.log and numpy.exp
   loglik = 0.0
   for i in range(len(y)):
    xi_times_beta = numpy.inner(x[i],beta)
    loglik += y[i]*xi_times_beta
    loglik -= N*math.log(1+math.exp(xi_times_beta))
   return loglik[0]

  # package GLM.Logistic
  def gradient(y,x,beta,N):
   """ Calculates the gradient of the likelihood function given the data (y,x) and coefficients beta.
       The gradient of the logistic simplifies down quite a bit. It's just the sum of the inner product
       of the ith residuals with the ith predictors. The ith residuals is just a single scalar.
       @y - the response
       @x - the predictors
       @beta - the coefficients
       @N - the number of trials (for binomial LR)
   """
   grad = numpy.matrix(list(map(lambda x: 0.0,range(beta.size))))
   for i in range(len(y)):
    grad += (y[i]-N*GLM.Logistic.logistic(x[i],beta))*x[i]
   return grad

  # package GLM.Logistic
  def hessian(y,x,beta,N):
   """ Calculates the hessian of the likelihood function given the data (y,x) and coefficients beta.
       The hessian simplifies like the gradient quite nicely to the sum of the data over:
       -N*pi(x)(1-pi(x))*outer(x,x)
       @y - the response
       @x - the predictors
       @beta - the coefficients
       @N - the number of trials (for binomial LR)
   """
   hessian = numpy.matrix(numpy.zeros(shape=(beta.size,beta.size)))
   for i in range(len(y)):
    prob_yi = GLM.Logistic.logistic(x[i],beta)
    hessian += N*prob_yi*(1-prob_yi)*numpy.matrix(numpy.outer(x[i],x[i]))
   return hessian
  
  # package GLM.Logistic
  def logistic(x,b):
   """ Calculates the (simple) logistic equation exp(dot(x,b))/(1+exp(dot(x,b)))
       This can likely be significantly sped up with the use of a lookup table and pade approximation. Would really only need points at 0,5,-5, and thresholds after which return {1,0}.
       For now, do the potentially more expensive thing.
   """
   eToProd = math.exp(numpy.inner(x,b))
   return eToProd/(1+eToProd)

  # package GLM.Logistic
  def linkDerivative(x,b):
   """ Calculates the derivative of the link function at x,b. For logistic, this happens to be p'(x,b) = p(x,b)(1-p(x,b)).
   """
   p = GLM.Logistic.logistic(x,b)
   return p*(1.0-p)

  # package GLM.Logistic
  def predict(x,beta,N):
   """ Predicts what the value of y should be given x, beta, and N. Basically just N*pi(x).
       @x - predictor variables. A matrix.
       @beta - coefficient vector.
       @N - number of trials (for binomial logistic).
   """
   return numpy.array(list(map(lambda u: N*GLM.Logistic.logistic(u,beta),x)))
