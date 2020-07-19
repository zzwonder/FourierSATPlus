import time

start = time.time()
number = 0
for i in range(200000000):
    number += 5
print(time.time()-start)
