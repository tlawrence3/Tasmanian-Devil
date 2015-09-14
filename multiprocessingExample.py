import multiprocessing as mp

def chunks(l, n):
   n = max(1, n)
   return [l[i:i + n] for i in xrange(0, len(l), n)]

def prune_rxn(foo, something, boo):
  do something

processes = []
for worker in xrange(num_proc):
  p = Process(target = prune_rxn, args = (foo, something, boo))
  p.start()
  processes.append(p)

for p in processes:
    p.join()
