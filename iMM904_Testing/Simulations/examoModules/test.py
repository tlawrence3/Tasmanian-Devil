@profile
def tester(a, b, c, i):
   c = c + a + b + i 
   print c
   return c

if __name__ == '__main__':
    i = 0
    a = 2
    b = 3
    c = 0
    while c < 100:
        i += 1
        c = tester(a, b, c, i)
	
