import math
import random
import numpy as np
from itertools import permutations

"""
Input: Number of buckets
Output: Hash function
Maxwell
"""
def getH(B):
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
      

def iterativeSetQuery(eps, k, n, B):
  # Start of initialization
  c = 20
  gamma = 1/600
  klist = []
  epslist = []
  blist = []
  for i in range(1, math.log(k, 2) + 1):
    klist.append(k * (gamma ** i))
    epslist.append(eps * ((10 * gamma) ** i))
    blist.append(c * klist[-1] / epslist[-1])
  m = 0
  for i in blist:
    m += i
  h = getH(B) # Create multiple hash functions and store in a list
  sigma = getSigma(n)
  phi = createPhi(h, sigma, n, B)
  # Make separate function to create x and S
  xh = [0] * n # Helper for x
  perm = permutations(range(1, n+1)) # Decide which k indices are non-zero
  for i in range(k): # Make the values the same
    xh[perm[i]] = 1
  x = np.array(xh)
  S = perm[:k] # The list of non-zero indices
  y = np.dot(phi, x)
  # End of initialization
  
  xprime = [0] * n # 
  for j in range(1, math.log(k, 2) + 1):
    T = [] # T is the list of indices that don't have any collisions
    goodbucks = [] # The list of good buckets
    counter = {} # Counts number of collisions for each bucket
    hprime = getH(B) # Hash function for this iteration, should be fixed with above
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
    
    
  

testH()
testSigma()
