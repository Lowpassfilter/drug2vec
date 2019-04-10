import numpy as np

def fullDrug():
    full_dict = np.load("full_dict.npy").item()
    drug_list = full_dict["drug"]
    return drug_list

def sampleDrug():
    f = open("sample_drugs.txt")
    drug_list = []
    import json
    for line in f.readlines():
        drug = json.loads(line.strip())
        drug_list.append(drug)
    return drug_list

def getID(drug):
    ids = drug["drugbank-id"]
    id_list = []
    for i in ids:
        if(len(i) == 2):
            id_list.insert(0, i["$"])
        else:
            id_list.append(i)
    return id_list

def bioDrug(drug):
    # id, name, pubchem sid
    name = drug["name"]
    id_list = getID(drug)
    name = drug["name"]
    try:
        ext_ids = drug["external-identifiers"]["external-identifier"]  # which is a list now
    except:
        #print("bio error on " + id_list[0])
        return ["bio err", id_list[0], name]
    pubchem_sid = ""
    for ids in ext_ids:
        resource = ids["resource"]
        identifier = ids["identifier"]
        if "PubChem Substance" in resource:
            pubchem_sid = identifier
        if pubchem_sid == "":
            print("error on " + id_list[0])
    return ["bio", name,  id_list[0], pubchem_sid, ""]

def smDrug(drug):
    # type, id, name, cid, sid
    name = drug["name"]
    id_list = getID(drug)
    name = drug["name"]
    try:
        ext_ids = drug["external-identifiers"]["external-identifier"]  # which is a list now
    except:
        #print("sm error on " + id_list[0])
        return ["sm err ", id_list[0], name,"" , ""]
    pubchem_sid = ""
    pubchem_cid = ""
    for ids in ext_ids:
        resource = ids["resource"]
        identifier = ids["identifier"]
        if "PubChem Substance" in resource:
            pubchem_sid = identifier
        if "PubChem Compound" in resource:
            pubchem_cid = identifier
        if pubchem_sid == "" and pubchem_cid == "":
            print("error on " + id_list[0])
    return ["sm ", id_list[0],name, pubchem_cid, pubchem_sid]

def getSMILES():
    f_r = open("drug_list.txt", "r")
    f_e = open("drug_miss.txt", "w")
    import pubchempy as pcp
    for line in f_r.readlines():
        items = line.strip().split("#")
        error = items[1]
        print("downloading " + items[1] + "  " + items[2])
        try:
            pcp.download('CSV', "temp/" + items[1] + "-name.csv", items[2], 'name', operation='property/CanonicalSMILES') 
        except:
            print("name " +  items[2])
            error += " name " + items[2]
        try:
            pcp.download('CSV', "temp/" + items[1] + "-cid.csv", [int(items[3])],  operation='property/CanonicalSMILES') 
        except:
            print("cid " + items[3])
            error += ", \t\t cid " + items[3]
        if error != items[1]:
            f_e.write(error + "\n")
    f_e.close()
    f_r.close()

def saveDrugList(drug_list):
    f = open("drug_list.txt", "w")
    for drug in drug_list:
        if drug["@type"] == "biotech":
            #l = bioDrug(drug)
            pass
        else:
            l = smDrug(drug)
            s = "#".join(l)
            f.write(s + "\n")
    f.close()

def saveInteractions(drug_list):
    f = open("drug_interaction.csv", "w")
    for drug in drug_list:
        id_list = [getID(drug)[0]]
        try:
            interactions = drug["drug-interactions"]["drug-interaction"]
            for inters in interactions:
                neibour = inters["drugbank-id"]
                desc = inters["description"]
                if "increase" in desc and "decrease" not in desc and "reduce" not in desc:
                    value = 1
                elif "increase" not in desc and ("decrease" in desc or "reduce" in desc):
                    value =-1
                else:
                    print(desc)
                id_list.append(inters["drugbank-id"])
                id_list.append(str(value))
        except:
            pass
        s = ",".join(id_list)
        f.write(s + "\n")

def mergeSMILES():
    import os
    files = os.listdir("temp/")
    SmileByName = {}
    SmileByCid = {}
    CidByName = {}
    CidByCid = {}
    for path in files:
        dbid = path.split("-")[0]
        f = open("temp/" + path, "r")
        f.readline()
        items = f.readline().strip().split(",")
        cid = items[0]
        smile = items[1].strip("\"")
        if "name" in path:
            CidByName[dbid] = cid
            SmileByName[dbid] = smile
        if "cid" in path:
            CidByCid[dbid] = cid
            SmileByCid[dbid] = smile
        f.close()
    nameKeys = set(SmileByName.keys())
    cidKeys = set(SmileByCid.keys())
    wholeKeys = list(nameKeys.union(cidKeys))
    wholeKeys.sort()
    fw = open("smile_download.csv", "w")
    for k in wholeKeys:
        c1 = ""
        s1 = ""
        c2 = ""
        s2 = ""
        if k in nameKeys:
            c1 = CidByName[k]
            s1 = SmileByName[k]
        if k in cidKeys:
            c2 = CidByCid[k]
            s2 = SmileByCid[k]
        l = [k, c1, s1, c2, s2]
        s = ",".join(l) + "\n"
        fw.write(s)
    fw.close()

    fw = open("smile_inconsistant.csv", "w")
    commonKeys = nameKeys.intersection(cidKeys)
    cnt = 0 
    for k in commonKeys:
        if SmileByName[k] != SmileByCid[k]:
            fw.write(k + "," + CidByName[k] + "," + CidByCid[k] + "\n")
            cnt += 1
    print(cnt)
    fw.close()

def getPubChemKey():
    f = open("smile_inconsistant.csv", "r")
    dbid = []
    cid1 = []
    cid2 = []
    for line in f.readlines():
        items = line.strip().split(",")
        dbid.append(items[0])
        cid1.append(items[1].strip("\""))
        cid2.append(items[2].strip("\""))
    import pubchempy as pcp
    pcp.download('CSV', 'name_keys.csv', cid1, operation='property/InChIKey') 
    pcp.download('CSV', 'cid_keys.csv', cid2, operation='property/InChIKey') 
    

def getDBKey(drug_list):
    f = open("smile_inconsistant.csv", "r")
    dbid = []
    for line in f.readlines():
        items = line.strip().split(",")
        dbid.append(items[0])
    drugs = []
    for drug in drug_list:
        drugID = getID(drug)[0]
        if drugID in dbid:
            drugs.append(drug)
    key = []
    for db in dbid:
        for drug in drugs:
            drugID = getID(drug)[0]
            if db == drugID:
                ps = drug['calculated-properties']['property']
                for p in ps:
                    if p["kind"] == "InChIKey":
                        key.append(p["value"])
    return(key)
                
def reduceIncon(drug_list):
    f = open("smile_inconsistant.csv", "r")
    dbid = []
    cid1 = []
    cid2 = []
    for line in f.readlines():
        items = line.strip().split(",")
        dbid.append(items[0])
        cid1.append(items[1].strip("\""))
        cid2.append(items[2].strip("\""))
    f.close()
    dbkey = getDBKey(drug_list)
    namekey = []
    cidkey = []

    f = open("name_keys.csv", "r")
    f.readline()
    for line in f.readlines():
        items = line.strip().split(",")
        namekey.append(items[1].strip("\""))
    f.close()
    
    f = open("cid_keys.csv", "r")
    f.readline()
    for line in f.readlines():
        items = line.strip().split(",")
        cidkey.append(items[1].strip("\""))
    f.close()
    
    fw = open("smile_inconsistent2.csv", "w")
    for i in range(len(dbid)):
        a = dbid[i]
        b = cid1[i]
        c = cid2[i]
        if dbkey[i] != namekey[i] and dbkey[i] == cidkey[i]:
            fw.write(a + "," + c + "\n")
        else:
            fw.write(a + "," + b + "," + c + "\n")
    fw.close()
 
def interactionStatics():
    f = open("drug_interaction.csv")
    drug = set([])
    cntP = 0
    cntN = 0
    cnt = 0
    for line in f.readlines():
        items = line.strip().split(",")
        d = items.pop(0)
        drug.add(d)
        flagP = False
        flagN = False
        length = int(len(items) / 2)
        for i in range(length):
            drug.add(items[i*2])
            if int(items[i*2 + 1]) == -1:
                cntN += 1
                flagP = True
            if int(items[i*2 + 1]) == 1:
                cntP += 1
                flagN = True
        if flagP and flagN:
            cnt += 1
    
    print(cnt)
    print(len(drug))
    print(cntP, cntN)
    print(cntP/(len(drug)**2))
    print(cntN/(len(drug)**2))

def BMCtest():
    f = open("drug_name_568.txt", "r")
    drug_list = []
    for line in f.readlines():
        drug_list.append(line.strip())
    f.close()
    
    f = open("smile_download.csv", "r")
    drug_dict = {}
    for line in f.readlines():
        items = line.strip().split()
        drug_dict[items[0].strip()] = "abc"
    f.close()

    for drug in drug_list:
        if not drug not in drug_dict.keys():
            print(drug)


if __name__ == "__main__":
    #drug_list = sampleDrug()
    #ps = drug['calculated-properties']['property']
    #for p in ps:
    #    if p["kind"] == "InChIKey":
    #        print(p["value"])

    #interactionStatics()
    #exit()
    #drug_list = fullDrug()
    
    #saveDrugList(drug_list)
    #getSMILES()
    #mergeSMILES()
    #reduceIncon(drug_list)
    BMCtest()
