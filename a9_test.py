import cv2
import a9
import numpy as np 

class ImageIO:
  def __init__(self, path:str) -> None:
    self.path = path
  def imread(self,fname:str)->np.ndarray:
    im = cv2.imread(self.path+'Input/'+fname)
    return np.array(im, dtype='float')/255
  def imwrite(self,im:np.ndarray, fname:str)->None:
    cv2.imwrite(self.path+'Output/'+fname, 255*im)

io = ImageIO('/Users/kartikeshmishra/Kartikesh/68370/a10/asst/')


def test_grad_descent():
  im=io.imread('pru.png')
  kernel=a9.gauss2D(1)
  im_blur=a9.convolve3(im, kernel)
  io.imwrite(im_blur, 'pru_blur.png')
  im_sharp=a9.deconvGradDescent(im_blur, kernel)
  io.imwrite(im_sharp, 'pru_sharp.png')

def test_conjugate_grad_descent():
  im=io.imread('pru.png')
  kernel=a9.gauss2D(1)
  im_blur=a9.convolve3(im, kernel)

  io.imwrite(im_blur, 'pru_blur.png')
  im_sharp=a9.deconvCG(im_blur, kernel,20)
  io.imwrite(im_sharp, 'pru_sharp_CG.png')

def test_real_psf():
  im=io.imread('pru.png')
  f=open('psf', 'r')
  psf=[map(float, line.split(',')) for line in f ]
  kernel=np.array(psf)
  im_blur=a9.convolve3(im, kernel)
  #kernel=kernel[::-1, ::-1]
  io.imwrite(im_blur, 'pru_blur_real.png')
  io.imwriteGrey(kernel/np.max(kernel), 'psf.png')
  im_sharp=a9.deconvCG(im_blur, kernel, 4)
  io.imwrite(im_sharp, 'pru_sharp_CG_real.png')



def test_conjugate_grad_descent_reg():
  im=io.imread('pru.png')
  kernel=a9.gauss2D(1)
  im_blur=a9.convolve3(im, kernel)
  noise=np.random.random(im_blur.shape)-0.5
  im_blur_noisy=im_blur+0.05*noise

  io.imwrite(im_blur_noisy, 'pru_blur_noise.png')
  im_sharp=a9.deconvCG_reg(im_blur_noisy, kernel,10)
  im_sharp_wo_reg=a9.deconvCG(im_blur_noisy, kernel,20)
  io.imwrite(im_sharp, 'pru_sharp_CG_reg.png')
  io.imwrite(im_sharp_wo_reg, 'pru_sharp_CG_wo_reg.png')


def test_naive_composite():
  fg=io.imread('bear.png')
  bg=io.imread('waterpool.png')
  mask=io.imread('mask.png')
  out=a9.naiveComposite(bg, fg, mask, 50, 1)
  io.imwrite(out, 'naive_composite.png')

def test_Poisson():

  y=50
  x=10

  fg=io.imread('bear.png')
  bgn=io.imread('waterpool.png')
  mask=io.imread('mask.png')

  mask[mask>0.5]=1.0
  mask[mask<0.6]=0.0
  mask = a9.duplicate(mask)
  bg = bgn[50:210,10:303]

  bg[bg==0]=1e-4
  fg[fg==0]=1e-4
  bg=np.log(bg)+3
  fg=np.log(fg)+3

  tmp=a9.Poisson(bg, fg, mask)
  w ,h , _ = fg.shape
  out = np.array(bgn)
  out[y:y+h, x:x+w]=np.exp(tmp-3)
  out[y:y+h, x:x+w]=tmp

  io.imwrite(out, 'poisson.png')

def test_PoissonCG():

  y=50
  x=10

  fg=io.imread('bear.png')
  bgn=io.imread('waterpool.png')
  mask=io.imread('mask.png')

  mask[mask>0.5]=1.0
  mask[mask<0.6]=0.0
  mask = a9.duplicate(mask)
  bg = bgn[50:210,10:303]

  bg[bg==0]=1e-4
  fg[fg==0]=1e-4
  bg=np.log(bg)+3
  fg=np.log(fg)+3

  tmp=a9.PoissonCG(bg, fg, mask)
  w ,h , _ = fg.shape
  out = np.array(bgn)
  out[y:y+h, x:x+w]=np.exp(tmp-3)
  out[y:y+h, x:x+w]=tmp

  io.imwrite(out, 'poisson_CG.png')



# test_grad_descent()
# test_conjugate_grad_descent()
# test_real_psf()
test_conjugate_grad_descent_reg()
test_naive_composite()
# test_Poisson()
# test_PoissonCG()

