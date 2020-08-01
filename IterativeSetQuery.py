"""An iterative algorithm for Set Query.

Author: Maxwell Pang (maxwell.amao@gmail.com)
"""

import math
import random
import numpy as np
from itertools import permutations


def getH(B):
  """Returns a hash function used in creating phi.

  Args:
    B: number of buckets.

  Returns:
    A pair-wise independent hash function mapping integers to range(B).
  """
  a = random.randint(0, B//2 - 1) * 2 + 1
  b = random.randint(0, B - 1)
  return lambda x: (a * x + b) % B

def testH():
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
  k = random.randint(0, 2 ** n - 1)
  return lambda x : (((k // (2 ** x)) % 2) - .5) * 2


def testSigma():
  sigma = getSigma(2 ** 10)
  numOnes = 0
  for i in range(2 ** 10):
    if (sigma(i) == 1):
      numOnes += 1
  assert numOnes > 2**9 - 64
  assert numOnes < 2**9 + 64
  
def createPhi(h, sigma, n, B):
  H = [h(i) for i in range(n)]
  G = [sigma(i) for i in range(n)]
  phi = np.zeros((n, B))
  phi[np.arrange(n), H] = 1
  phi *= G
  return phi


class IterativeSetQuery:
  def __init__(self, eps, k, n, B):
    self.c = 20
    self.gamma = 1/600
    self.klist = []
    self.epslist = []
    self.blist = []
    self.num_blocks = math.log(k, 2)
    for i in range(1, self.num_blocks + 1):
      self.klist.append(k * (self.gamma ** i))
      self.epslist.append(self.eps * ((10 * self.gamma) ** i))
      self.blist.append(self.c * self.klist[-1] / self.epslist[-1])
    self.m = 0
    for i in self.blist:
      self.m += i
    self.hlist = [] # Keeps track of all the hash functions
    self.sigmalist = [] # Keeps track of all the other hash functions
    self.philist = [] # Keeps track of all the arrays made by hash function pairs
    for i in range(self.num_blocks):
      h = getH(B) # Creates a hash function
      self.hlist.append(h)
      sigma = getSigma(n)
      self.sigmalist.append(sigma)
      phi = createPhi(h, sigma, n, B)
      self.philist.append(phi)
    
  def measure(self, x):
    y = [0] * self.m
    for i in range(self.num_blocks):
      helper = np.dot(self.philist[i], x)
      for j in range(self.m):
        y[m] += helper[m]
    return np.array(y)

  def query(self, y, S):
    xprime = [0] * n # 
    for j in range(1, self.num_blocks + 1):
      T = [] # T is the list of indices that don't have any collisions
      goodbucks = [] # The list of good buckets
      counter = {} # Counts number of collisions for each bucket
      hprime = self.hlist[j-1] # Hash function for this iteration, should be fixed with above
      for i in range(B):
        counter[i] = 0
      for i in S:
        counter[hprime(i)] += 1
      for i in range(B):
        if counter[i] == 1:
          goodbucks.append(i)
      for i in S:
        if hprime(i) in goodbucks:
          T.append(i)
        l = 0 # Change this to initialize before for loop
        for k in range(j-2):
          l += blist[k]
          xhat = [0] * n
      for i in T:
        xhat[i] = y[l + hprime(i)] # Need to use sigma here
        S.remove(i)
        y = y - np.dot(phi, np.array(xhat)) # Need to make this faster, matrix multiplication bad
    xprime[i for i in T] = xhat[i]
    return xprime


  # How can I make it so that S is both used in constructing x but also in retrieving x; does it have to appear both in the function and class??

# Make separate function to create x and S
def create(x, philist, S, num_blocks, m):
  xh = [0] * n # Helper for x
  for i in S: # Make the values the same
    xh[i] = 1
  x = np.array(xh)
  y = [0] * m
  for i in range(num_blocks):
    prod = np.dot(philist[i], x)
    for j in range(m):
      y[j] += prod[j]
  realy = np.array(y)
  return realy

def makePerm(n, k):
  perm = permutations(range(n))
  return perm[:k]
