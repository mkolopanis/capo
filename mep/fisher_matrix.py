import aipy as a, numpy as n, pylab as p
import useful_functions as uf
import matplotlib as mpl


def gaussian_model((A,nu0,sigma),nuvec):
    return -A*n.exp(-(nuvec-nu0)**2/sigma**2/2.)

def gaussian_model_derivs((A,nu0,sigma),nuvec):
    Aderiv = -n.exp(-0.5*(nuvec-nu0)**2/sigma**2)
    nu0deriv = (A/sigma**2)*(nuvec-nu0)*n.exp(-0.5*(nuvec-nu0)**2/sigma**2)
    sigderiv = (A/sigma**3)*(nuvec-nu0)**2*n.exp(-0.5*(nuvec-nu0)**2/sigma**2)
    return Aderiv, nu0deriv, sigderiv

def compute_bias((A,nu0,sigma),nuvec,Cinv,expBias):
    bias = n.zeros(3)
    derivs = gaussian_model_derivs((A,nu0,sigma),nuvec)
    Finv = n.linalg.inv(fisher_matrix((A,nu0,sigma),nuvec,Cinv))
    for ii in range(3):
        bias[ii] = n.dot(derivs[ii],n.dot(Cinv,expBias))
    bias = n.dot(Finv,bias)
    return bias

def fisher_matrix((A,nu0,sigma),nuvec,Cinv):
    fisher = n.zeros((3,3))
    derivs = gaussian_model_derivs((A,nu0,sigma),nuvec)
    for ii in range(3):
        for jj in range(3):
            fisher[ii,jj] = n.dot(n.transpose(derivs[ii]),uf.vdot(Cinv,derivs[jj]))
    print 'fisher ',fisher 
    return fisher 

def fisher_select_pair(F,si,sj):
    Finv = n.linalg.inv(F) #uf.pseudo_inverse(F) 
    for kk in range(Finv.shape[0])[::-1]:
        if kk==si or kk==sj:
            continue
        Finv = n.delete(Finv,kk,axis=0)
        Finv = n.delete(Finv,kk,axis=1)
    #print Finv 
    Fnew = n.linalg.inv(Finv) #uf.pseudo_inverse(Finv)
    return Fnew 

def plot_pairwise_contours(theta,nuvec,Cinv,lvls=(-2.291,-6.158,-11.618)):
    """
    > theta is a (3,) vector that contains the model parameters
    > thetavecs is a (n,3) matrix that contains the values of the parameters 
        that will be plotted over
    """
    labels = ['A','nu0','sigma']
    fisher = fisher_matrix(theta,nuvec,Cinv)
    Finv = n.linalg.inv(fisher)
    print "F-inverse"
    print Finv
    thetavecs = n.zeros((50,theta.shape[0]))
    for ii in range(theta.shape[0]):
        thetavecs[:,ii] = n.linspace(theta[ii]-5*n.sqrt(Finv[ii,ii]),theta[ii]+5*n.sqrt(Finv[ii,ii]),num=50)
    #print thetavecs
    for ii,jj in ((0,1),(0,2),(1,2)):
        print ii,jj
        ts = thetavecs[:,[ii,jj]]
        #print thetavecs.shape
        #print ts.shape
        fs = fisher_select_pair(fisher,ii,jj)
        #print fs.shape
        t0,t1 = n.meshgrid(ts[:,0],ts[:,1])
        #print t0.shape

        Z = fs[0,0]*(t0-theta[ii])*(t0-theta[ii]) + (fs[0,1]+fs[1,0])*(t0-theta[ii])*(t1-theta[jj]) + fs[1,1]*(t1-theta[jj])*(t1-theta[jj])
        Z *= -1

        p.pcolor(t0,t1,Z)
        p.colorbar()
        CS = p.contour(t0,t1,Z,levels=lvls) #levels=lvls
        p.clabel(CS, inline=1, fontsize=10)
        #p.contour(t0,t1,Z,lvls)
        p.xlabel(labels[ii])
        p.ylabel(labels[jj])
        p.savefig('./figures/fisher/contours_{0}_{1}.pdf'.format(labels[ii],labels[jj]))
        p.clf()

if __name__=='__main__':
    theta = n.array([2,50.,1])
    nuvec = n.arange(30,100)
    #thetavecs = n.vstack((n.linspace(1,20,num=50),n.linspace(-5,5,num=50),n.linspace(0.1,1,num=50)))
    #Cinv = n.dot(nuvec.reshape((nuvec.shape[0],1)),nuvec.reshape(1,nuvec.shape[0]))
    Cinv = n.eye(70)#uf.invertible_matrix(70)
    plot_pairwise_contours(theta,nuvec,Cinv)




