class Squares(object):
    def __init__(self, start, stop):
       self.start = start
       self.stop = stop

    def __iter__(self): 
        return self

    def __next__(self): # next in Python 2
       if self.start >= self.stop:
           raise StopIteration
       current = self.start * self.start
       self.start += 1
       return current
    def current(self):
        return  self.start


iterator = Squares(1, 10)

class MyI(object):
    def __init__(self,start,stop):
        self.my_list=[i for i in range(20)]
        self.start=start
        self.stop=stop
        
    def __iter__(self):
        print("iter")
        return self
    
    def __next__(self):
        print("next")
        if self.start >=self.stop:
            raise StopIteration
        current=self.start
        self.start+=2
        return current


myI=MyI(1,9)
a=list(myI)
print(a)