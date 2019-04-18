import fileinput
import numpy

from matplotlib import pylab

"""
Plotting a zero-mean cumulative sum of numbers.
"""

def plot(filename):
    """
    Read single-column numbers in filename and plot zero-mean
    cumulative sum.
    """
    numbers = []
    for line in fileinput.input(filename): # (2)
        numbers.append(float(line.split(’\n’)[0]))
    mean = numpy.mean(numbers) # (3)
    cumulative_sum = numpy.cumsum([number -
                                   mean for number in numbers])
    pylab.plot(cumulative_sum[0::10],’k-’) # (4)
    pylab.xlabel(’i’)
    pylab.title(’Zero Mean Cumulative Sum’)
    pylab.savefig(filename+’.png’) # (5)
    pylab.show()
