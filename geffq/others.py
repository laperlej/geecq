"""
Created on May 13, 2014

@author: Jonathan Laperle(jonathan.laperle@usherbrooke.ca)
"""
def mean(container):
    """Averages all elements in container
    
    Args:
        container: a container of sum()able elements
        
    Returns:
        The mean value of all elements
    """
    return sum(container)/len(container)

def argmax(container):
    """Finds the index of maximal element
    
    Args:
        a container of comparable elements
        
    Returns:
        the index of the maximal element
    """
    maxIndex = 0
    for i in range(len(container)):
        if container[maxIndex] < container[i]:
            maxIndex = i
    return maxIndex

def writeFile(fileName, content):
    """Writes content into a file
    """
    f = open(fileName, 'w')
    f.write(content)
    f.close()
        
def appendFile(fileName, content):
    """Appends content into a file
    """
    f = open(fileName, 'a')
    f.write(content)
    f.close()