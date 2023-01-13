import numpy as np

# === Deconvolution with gradient descent ===

def dotIm(im1, im2):
  d = np.sum(np.ndarray.flatten(im1*im2))
  if d==0: return 1e-6
  return d

def duplicate(M):
  efunc = lambda row: np.array([[x,x,x] for x in row])
  return np.array([efunc(row) for row in M])


def applyKernel(im, kernel):
  ''' return Mx, where x is im '''
  sx,sy = kernel.shape
  s1,s2=int((sx-1)/2), int((sy-1)/2)
  padded_im = np.pad(im, ((s1,s1), (s2,s2), (0,0)),mode='edge')
  M1 = duplicate(kernel)
  result = np.zeros_like(im)
  for i in range(im.shape[0]):
      for j in range(im.shape[1]):
          result[i][j] = np.sum(padded_im[i:i+sx,j:j+sy]*M1,axis=(0,1))
  return result

def applyConjugatedKernel(im, kernel):
  ''' return M^T x, where x is im '''
  return applyKernel(im,np.transpose(kernel))

def computeResidual(kernel, x, y):
  ''' return Mx - y '''
  return y-applyKernel(x,kernel)


def computeStepSize(r, kernel):
  mtmr = applyConjugatedKernel(applyKernel(r,kernel),kernel)
  rmtmr = dotIm(r, mtmr)
  rmtmr[rmtmr==0]=1e-6
  return dotIm(r,r)/rmtmr

def deconvGradDescent(im_blur, kernel, niter=10):
  ''' return deblurred image '''
  x = np.zeros_like(im_blur)
  for _ in range(niter):
    r = applyConjugatedKernel(computeResidual(kernel,x,im_blur),kernel)
    alpha = computeStepSize(r, kernel)
    x = x + alpha*r

  return x

# === Deconvolution with conjugate gradient ===

def computeGradientStepSize(r, d, kernel):
  return dotIm(r,d)/(10**(-6)+dotIm(d,applyConjugatedKernel(applyKernel(d,kernel),kernel)))

def computeConjugateDirectionStepSize(old_r, new_r):
  return dotIm(new_r,new_r)/(10**(-6)+dotIm(old_r,old_r))

def deconvCG(im_blur, kernel, niter=10):
  ''' return deblurred image '''
  x = np.zeros_like(im_blur)
  r = applyConjugatedKernel(computeResidual(kernel,x,im_blur),kernel)
  d = np.array(r)
  for _ in range(niter):
    Ad = applyConjugatedKernel(applyKernel(d,kernel),kernel)
    r_dot = dotIm(r,r)
    alpha = r_dot/dotIm(d,Ad)  #computeGradientStepSize(r, d, kernel)
    x = x + alpha*d
    r = r - alpha*Ad
    beta = dotIm(r,r)/r_dot
    d = r + beta*d
  return x

def laplacianKernel():
  ''' a 3-by-3 array '''
  return np.array([
    [0,-1,0],
    [-1,4,-1],
    [0,-1,0]
  ])

def applyLaplacian(im):
  ''' return Lx (x is im)'''
  return applyKernel(im,laplacianKernel())

def applyAMatrix(im, kernel):
  ''' return Ax, where A = M^TM'''
  return applyConjugatedKernel(applyKernel(im,kernel),kernel)

def applyRegularizedOperator(im, kernel, lamb):
  ''' (A + lambda L )x'''
  return applyAMatrix(im,kernel) + lamb*applyLaplacian(im)


def computeGradientStepSize_reg(grad, p, kernel, lamb):
  return dotIm(grad,p)/(10**(-6)+dotIm(p,applyRegularizedOperator(p,kernel,lamb)))

def deconvCG_reg(im_blur, kernel, lamb=0.05, niter=10):
  ''' return deblurred and regularized im '''
  x = np.zeros_like(im_blur)
  r = applyConjugatedKernel(im_blur,kernel) - applyRegularizedOperator(x,kernel,lamb)
  d = np.array(r)
  for _ in range(niter):
    Ad = applyRegularizedOperator(d,kernel,lamb)
    r_dot = dotIm(r,r)
    alpha = r_dot/dotIm(d,Ad) 
    x = x + alpha*d
    r = r - alpha*Ad
    beta = dotIm(r,r)/r_dot
    d = r + beta*d
  return x

    
def naiveComposite(bg, fg, mask, y, x):
  ''' naive composition'''
  res = np.array(bg)
  sy,sx,_ = fg.shape
  res[y:y+sy,x:x+sx] = mask*fg + (1-mask)*bg[y:y+sy,x:x+sx]
  return res

def Poisson(bg, fg, mask, niter=200):
  ''' Poisson editing using gradient descent'''
  x = (1-mask)*np.array(bg)
  b = applyLaplacian(fg)
  for _ in range(niter):
    r = b-applyLaplacian(x)
    alpha = dotIm(r,r)/(10**(-6)+ dotIm(r, applyLaplacian(r)))
    x = x + alpha*mask*r
            
  return x



def PoissonCG(bg, fg, mask, niter=200):
  ''' Poison editing using conjugate gradient '''
  x = (1-mask)*np.array(bg)
  r = applyLaplacian(fg)-applyLaplacian(x)
  d = np.array(r)
  for _ in range(niter):
    Ad = applyLaplacian(d)
    r_dot = dotIm(r,r)
    alpha = r_dot/dotIm(d,Ad)  
    x = x + alpha*d
    r = mask*(r - alpha*Ad)
    beta = dotIm(r,r)/r_dot
    d = r + beta*d
  return x 
  
  

#==== Helpers. Use them as possible. ==== 

def convolve3(im, kernel):
  from scipy import ndimage
  center=(0,0)
  r=ndimage.filters.convolve(im[:,:,0], kernel, mode='reflect', origin=center) 
  g=ndimage.filters.convolve(im[:,:,1], kernel, mode='reflect', origin=center) 
  b=ndimage.filters.convolve(im[:,:,2], kernel, mode='reflect', origin=center) 
  return (np.dstack([r,g,b]))

def gauss2D(sigma=2, truncate=3):
  kernel=horiGaussKernel(sigma, truncate);
  kerker=np.dot(kernel.transpose(), kernel)
  return kerker/sum(kerker.flatten())

def horiGaussKernel(sigma, truncate=3):
  from scipy import signal
  sig=signal.gaussian(2*int(sigma*truncate)+1,sigma)
  return np.array([sig/sum(sig)])


if __name__ == '__main__':
  pass

