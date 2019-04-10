import numpy as np
from sklearn.metrics import pairwise_distances
import os
import psutil
import sys


f = open("data/drugname568.txt", "r")
NAMES = []
for item in f.readlines():
    NAMES.append(item.strip())
f.close()

CNT = 568
T = 100

def memory():
    process = psutil.Process(os.getpid())
    num = process.memory_info().rss
    uk = 1024
    um = uk*1024
    ug = um*1024
    g = num // ug
    m = (num%ug) // um
    k = num%uk
    print("{}G {}M {}K".format(g, m, k))

def save(order, pairs, threshold, filename):
    f = open(filename, "w")
    for i in range(threshold):
        ids = order[i]
        f.write(",".join(pairs[ids]) + "\n")
        print(",".join(pairs[ids]))
    f.close()


def similar2(array, sufix):
    dist = pairwise_distances(array)
    values = []
    pairs = []
    for i in range(CNT):
        for j in range(i+1, CNT):
            values.append(dist[i, j])
            pairs.append([NAMES[i], NAMES[j]])

    order = np.argsort(values)
    save(order, pairs, T, "results/doulbe_pairs_" + str(sufix) +".txt")

def similar3(array, sufix):
    values = []
    pairs = []
    for i in range(CNT):
        print(i)
        base = CNT - i
        h = (base - 1) * (base - 2) // 2
        m = np.zeros(shape=(h, CNT))
        row = 0
        for j in range(i+1, CNT):
            for k in range(j+1, CNT):
                m[row][i] = 1
                m[row][j] = 1
                m[row][k] = -1
                row += 1
                pairs.append([NAMES[i], NAMES[j], NAMES[k]])
        temp = np.matmul(m, array)**2
        temp = temp.sum(axis=1)
        values += list(temp)

    order = np.argsort(values)
    save(order, pairs, T, "results/triple_pairs_" + str(sufix) +".txt")

def savedouble():
    filename = "results/99.npy"
    array = np.load(filename)
    
    row_dict = {}
    col_dict = {}
    index = 0
    length = CNT*(CNT-1)//2
    print(length)
    partial = np.zeros(shape=(length, 100), dtype=np.float64)
    for row in range(0, CNT):
        print(row)
        for col in range(row+1, CNT):
            partial[index] = array[row] + array[col]
            row_dict[index] = row
            col_dict[index] = col
            index += 1
    print(index)
    np.save("results/doublearray.npy", partial)
    np.save("results/rowdict.npy", row_dict)
    np.save("results/coldict.npy", col_dict)

def similiar4(start, end):
    double = np.load("results/doublearray.npy")
    row_dict = np.load("results/rowdict.npy").item()
    col_dict = np.load("results/coldict.npy").item()
    step = 1000
    for i in range(start, end, step):
        print(i)
        similiar4partial(i, step, double, row_dict, col_dict)

def similiar4partial(index, step, double, row_dict, col_dict):
    length = len(row_dict)
    step = min(step, length - index)
    width = length - index
    
    X = double[list(range(index, index+step))]
    Y = double[list(range(index, length))]
    dist = pairwise_distances(X, Y)
    cnt = 0

    size = (width + width - step + 1) * step // 2
    values = np.zeros(shape=(size))
    pairs = np.zeros(shape=(size, 4))
    for i in range(step):
        for j in range(i+1, width):
            idi = i + index
            idj = j + index
            a1 = row_dict[idi]
            a2 = col_dict[idi]
            a3 = row_dict[idj]
            a4 = col_dict[idj]
            values[cnt] = dist[i, j]
            pairs[cnt] = np.array([a1, a2, a3, a4])
            cnt += 1
    order = np.argsort(values)
    
    tempV = []
    tempP = []

    for i in range(size):
        ids = order[i]
        if len(set(pairs[ids])) == 4:
            tempV.append(values[ids])
            tempP.append(pairs[ids])
        if len(tempV) == 100:
            break
    tempV = np.array(tempV)
    tempP = np.array(tempP)
    d = {}
    d["v"] = tempV
    d["p"] = tempP
    print(tempP)
    np.save("results/temp_" + str(index) + ".npy", d)

    del values
    del pairs
    del order
    del dist
    del tempV
    del tempP
    
def merge4():
    values = []
    pairs = []
    for i in range(0, 161028, 1000):
        filename = "results/temp_" + str(i) + ".npy"
        d = np.load(filename).item()
        values.extend(list(d["v"]))
        pairs.extend(list(d["p"]))
        print(i)
    order = np.argsort(values)
    str_pairs = []
    for i in range(len(pairs)):
        a, b, c, d = pairs[i]
        a = int(a)
        b = int(b)
        c = int(c)
        d = int(d)
        str_pairs.append([NAMES[a], NAMES[b], NAMES[c], NAMES[d]])

    save(order, str_pairs, T, "results/quadra_pairs.txt")

if __name__ == "__main__":
    #savedouble()   
    #similiar4partial(0, 1000)
    merge4()
    exit()

    if len(sys.argv) != 3:
        print("wrong number of para")
        exit()
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    similiar4(start, end)
