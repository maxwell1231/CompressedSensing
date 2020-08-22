"""An iterative algorithm for Set Query.

Author: Maxwell Pang (maxwell.amao@gmail.com)
"""

import math
import random
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def getH(B):
  """Returns a hash function used in creating phi.

  Args:
    B: number of buckets.

  Returns:
    A pair-wise independent hash function mapping integers to range(B).
  """
  print("B is %s" % B)
  a = random.randint(0, B - 1) * 2 + 1
  b = random.randint(0, B - 1)
  return lambda x: (a * x + b) % B

def testH():
  """Tests a hash function used in creating phi.

  Args:
    None.

  Returns:
    None.
  """
  h = getH(2 ** 10)
  results = [0] * 1024
  for i in range(2 ** 10):
    results[h(i)] += 1
  total = 0
  for i in results:
    if i >= 2:
      total += 1
  assert total <= 2
  hprime = getH(2 ** 10)
  collisions = 0
  for i in range(2 ** 10):
    if h(i) == hprime(i):
      collisions += 1
  assert collisions <= 3

  
def getSigma(n):
  """Returns a hash function used in creating phi.

  Args:
    n: length of x.

  Returns:
    A hash function mapping range(n) to 1 or -1.
  """
  k = random.randint(0, n - 1)
  return lambda x : (hash(str((x + k) % n)) % 2) * 2 - 1


def testSigma():
  """Tests a hash function used in creating phi.

  Args:
    None.

  Returns:
    None.
  """
  sigma = getSigma(2 ** 10)
  numOnes = 0
  for i in range(2 ** 10):
    if (sigma(i) == 1):
      numOnes += 1
  assert numOnes > 2**9 - 64
  assert numOnes < 2**9 + 64

# Not used, invalid syntax  
def createPhi(h, sigma, n, B):
  H = [h(i) for i in range(n)]
  G = [sigma(i) for i in range(n)]
  phi = np.zeros((n, B))
  phi[np.arrange(n), H] = 1
  phi *= G
  return phi


class IterativeSetQuery:
  """Creates a class to set up and solve compressed sensing."""
  def __init__(self, eps, k, n):
    """Initializes the parameters.

    Args:
      eps: A variable relating to error bound.
      k: The number of non-zero indices.
      n: The length of x.

    Returns:
      Nothing.
    """
    self.c = 20
    self.gamma = 1./600
    self.blist = []
    self.num_blocks = int(math.log(k, 10))
    for i in range(1, self.num_blocks + 1):
      eps *= 10
      self.blist.append(int(self.c * k / eps) + 1)
    self.m = 0
    for i in self.blist:
      self.m += i
    self.hlist = [] # Keeps track of all the hash functions
    self.sigmalist = [] # Keeps track of all the other hash functions
    self.philist = [] # Keeps track of all the arrays made by hash function pairs
    for i in range(self.num_blocks):
      h = getH(self.blist[i]) # Creates a hash function
      self.hlist.append(h)
      sigma = getSigma(n)
      self.sigmalist.append(sigma)
    
  def measure(self, x, start_block=0):
    """Creates y from x.

    Args:
      x: A k-sparse vector.

    Returns:
      The product of phi and x.
    """
    y = []
    for j in range(start_block, self.num_blocks):
      t = [0] * self.blist[j]
      for i in x:
        t[self.hlist[j](i)] += x[i] * self.sigmalist[j](i)
      y.extend(t)
    return np.array(y)

  def query(self, y, S):
    """Finds x from y and the list of non-zero indices.

    Args:
      y: The product of phi and x.
      S: The list of non-zero indices.

    Returns:
      An approximation for x based on y and S.
    """
    S = set(S)
    xprime = {}
    print("blist=%s" % self.blist)
    for j in range(self.num_blocks):
      counter = {} # Counts number of collisions for each bucket
      hprime = self.hlist[j] # Hash function for this iteration, should be fixed with above
      B = self.blist[j]
      xhat = {}
      removable = []
      for i in S:
        hashv = hprime(i)
        if hashv in counter:
          counter[hashv] = False
        else:
          counter[hashv] = True
      for i in S:
        hashv = hprime(i)
        if counter.get(hashv, False):
          xprime[i] = xhat[i] = y[hprime(i)] * self.sigmalist[j](i)
          removable.append(i)
      for i in removable:
        S.remove(i)
      y = y[B:] - self.measure(xhat, start_block=j+1)
    return xprime

  def createRep(self, S):
    Sp = list(S)
    Sp.sort()
    print("Sp is %s" % Sp)
    phiArr = []
    for i in range(self.m):
      phiArr.append([0] * len(Sp))
    print("The number of blocks is %s" % self.num_blocks)
    print("Phi array looks like %s" % phiArr)
    l = 0
    for i in range(self.num_blocks):
      for j in range(len(Sp)):
        s_j = Sp[j]
        hashv = self.hlist[i](s_j)
        print("%s, %s" % (j, l + hashv))
        phiArr[l + hashv][j] = self.sigmalist[i](s_j)
      l += self.blist[i]
    return phiArr
    


  # How can I make it so that S is both used in constructing x but also in retrieving x; does it have to appear both in the function and class??


def createX(n, S):
  """Creates an x.

  Args:
    n: The length of x.
    S: The list of non-zero indices.

  Returns:
    A k-sparse vector of length n.
  """
  return {i: 1. for i in S}


def makePerm(n, k):
  """Returns a random list of k integers from 0 to n - 1.

  Args:
    n: The range to choose from.
    k: The length of the list.

  Returns:
    A random list of k integers from 0 to n - 1.
  """
  result = np.random.permutation(n)
  return set(result[:k])

n = 10000000
k = 100
S = makePerm(n, k)
x = createX(n, S)
print("X = %s\nS = %s" % (x, S))

isq = IterativeSetQuery(eps=.1, k=k, n=n)
y = isq.measure(x)
print("Y = %s" % y)

tic = time.perf_counter()
xprime = isq.query(y, S)
toc = time.perf_counter()
print("Approximate of X = %s" % xprime)
print("Amount of time = %s" % ((toc - tic) * 1e6))

phirep = isq.createRep(S)
print("Phi array looks like %s" % phirep)

data = np.array(phirep)
cmap = plt.cm.bwr
norm = plt.Normalize(vmin=data.min(), vmax=data.max())
image = cmap(norm(data))
plt.imsave('isqphi.png', image)

# Notes for next time
# Plot phi, results after each execution
# Try half 1, half -1 instead of all 1
# S could just have intersection with some of the non-zero indices
# Check edge cases
# Measure time and space used, use timer.start??
# Make separate documents for inputs and outputs
