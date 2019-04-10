import numpy as np
import pubchempy as pcp


def getId():
    f = open("data/drugname568.txt", "r")
    names = []
    for item in f.readlines():
        names.append(item.strip())
    f.close()
    return names


def writeSmileFiles():
    drugnames = getId()
    f = open("data/drugbanksmile.csv", "r")
    cid_name_dict = {}
    smile_name_dict = {}
    cid_cid_dict = {}
    smile_cid_dict = {}
    for line in f.readlines():
        items = line.strip().split(",")
        cid_name_dict[items[0]] = items[1]
        smile_name_dict[items[0]] = items[2]
        cid_cid_dict[items[0]] = items[3]
        smile_cid_dict[items[0]] = items[4]
    f.close()
    fw = open("data/SmileByName.txt", "w")
    for name in drugnames:
        fw.write(name +',\t'+ cid_name_dict[name] +',\t'+ smile_name_dict[name] + "\n")
    fw.close()

    fw = open("data/SmileByCid.txt", "w")
    for name in drugnames:
        fw.write(name +',\t'+ cid_cid_dict[name] +',\t'+ smile_cid_dict[name] + "\n")
    fw.close()

def downloadPNG():
    f = open("data/SmileByName.txt", 'r')
    for line in f.readlines():
        items = line.strip().split(",")
        print(items[0])
        if items[1].strip() == '':
            continue
        pcp.download("PNG", "data/ByName/" + items[0].strip(), items[1].strip(), 'cid')
    f.close()
    print("ByName finished")
    f = open("data/SmileByCid.txt", 'r')
    for line in f.readlines():
        items = line.strip().split(",")
        print(items[0])
        if items[1].strip() == '':
            continue
        pcp.download("PNG", "data/ByCid/" + items[0].strip(), items[1].strip(), 'cid')
    f.close()
    print("ByCid finished")

def createGraph():
    

if __name__ == "__main__":
    downloadPNG()
