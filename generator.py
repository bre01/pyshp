def h():
    for i in range(19):
        yield i

m=h()
print(m.__next__())
print(m)