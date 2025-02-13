def fib(n):
    if n == 1:
        return 1
    elif n == 2:
        return 1      
    else:
        return fib(n-1) + fib(n-2)
    

# print([fib(n) for n in range(10)])
# print([fib(t) for t in range(10)])
for t in range(2,1000):
    print(fib(t)/fib(t-1))