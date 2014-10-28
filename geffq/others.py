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
    max_index = 0
    for i in range(len(container)):
        if container[max_index] < container[i]:
            max_index = i
    return max_index

def write_file(file_name, content):
    """Writes content into a file
    """
    f = open(file_name, 'w')
    f.write(content)
    f.close()

def append_file(file_name, content):
    """Appends content into a file
    """
    f = open(file_name, 'a')
    f.write(content)
    f.close()
