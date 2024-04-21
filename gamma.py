from math import gamma

D = 10 / 2 

def mine(val, p): 
  a = 1
  return a * (((1 - p) ** (-1 / a)) - 1)

def real(val):
  return gamma(val)

print(real(3))
print(mine(3, D))