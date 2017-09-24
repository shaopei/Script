#http://www.channelgrubb.com/blog/2015/2/27/beta-binomial-in-python

########### PARAMETERS ################
# alpha: scalar, the alpha parameter of the prior beta distribution
# beta: scalar, the beta parameter of the prior beta distribution
# n: scalar, the number of binomial trials
# k: array/vector-like, the number of successes

###########BETA-BINOMIAL: LOG TRANSFORM OF GAMMA FUNCTION ############################
from scipy import special
import numpy as np

#in R dbetabinom.ab(x, size, shape1, shape2, log = FALSE, Inf.shape = 1e6)
def dbetabinom_ab (k, n, alpha,beta):
    #return the pdf
    part_1 = special.comb(n,k)
    part_2 = special.betaln(k+alpha,n-k+beta)
    part_3 = special.betaln(alpha,beta)
    result = (np.log(part_1) + part_2)- part_3
    return np.exp(result)

def pbetabinom_ab (k, n, alpha,beta):
    assert k <= n, str(k)+' is too large'
    p = 0
    for i in xrange(k+1):
        p += dbetabinom_ab(i, n, alpha,beta)
    return min(round(p,10),1)
    
    
#def betaBinom_test_ab(k, n, alpha,beta):
#    assert k <= n, str(k)+' is too large'
#    p = 0
#    for i in xrange(k+1):
#        p += dbetabinom_ab(i, n, alpha,beta)
#    return min(round((min(p, 1-p))*2,10),1)


def betaBinom_test_ab(k, n, alpha,beta):
    assert k <= n, str(k)+' is too large'
    temp = pbetabinom_ab (k, n, alpha,beta)
    return 2*min(temp,1-temp)




def binomtest_modified_to_betaBinomTest(x, n, p):
    #return (scipy.stats.binom_test(x, n, p), normal_approx(x, n, p))
    return betaBinom_test_ab (x, n, 2.019315,2.031387)

def binomtest(x, n, p):
    return binomtest_modified_to_betaBinomTest(x, n, p)



for x in xrange(201):
    print x, betaBinom_test_ab (x, 200, 2.019315,2.031387)

for x in xrange(10):
    print x, betaBinom_test_ab (x, 10, 600,400)



    

#def binomtest(k, n, p):
#    assert k <= n, str(k)+' is too large'
#    k = min(k, n-k)
#    p_out = 0
#    for i in xrange(k+1):
#        p_out += scipy.stats.binom.pmf(i, n, p)
#        print i, scipy.stats.binom.pmf(i, n, p)
#    print p_out, 1-p_out
#    return min(round(p_out*2,10),1)


#def binomtest(k, n, p):
#    k = min(k, n-k)
#    p_out = scipy.stats.binom.cdf(k,n,p)
#    return min(round(p_out*2,10),1)




import matplotlib.pyplot as plt
xx = range(501)
xx_pvalue=[]
for i in xx:
    xx_pvalue.append(betaBinom_test_ab (i, 500, 2.019315,2.031387))
plt.scatter(xx, xx_pvalue)
plt.title('Scatter plot pythonspot.com')
plt.xlabel('x')
plt.ylabel('y')
plt.show()