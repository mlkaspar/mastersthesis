#! python2
#use this on windows cmd line: FOR %i IN (*.txt) DO python biomass.py %i nh
#for bash, use: for i in *.txt; do python getBiomass.py %i nh; done
import sys
import csv
import math  

if len(sys.argv) < 3:
    print("1. argument = .txt files, 2. argument = h for header, nh without header")
    sys.exit()

#functions
def mean(values):
    result = 0
    for i in values:
        result += i
    result = result / (len(values))
    return result

def stats(value, name):
    result = ("Min and Max and Mean value of " + name + ":\n"
              + str(min(value))
              + ", "
              + str(max(value))
              + ", "
              + str(mean(value)))
    return result

#Calculate width, length, volume, carbon and save as new array

with open(sys.argv[1]) as f:
    print(sys.argv[1])
    #skip first 11 lines!
    f = csv.reader(f.readlines()[11:], delimiter="#")
    calculated_width = []
    calculated_length = []
    calculated_volume = []
    calculated_carbon = []

    for i, line in enumerate(f):
        area = float(line[1])
        voleqcylinder = float(line[2])
        perimeter = float(line[3])
        length = float(line[4])
        width = float(line[5])
        maxferet = float(line[6])
        minferet = float(line[7])

        #length
        calc_length = max(length, maxferet)
        if calc_length <= 0:
            print("invalid values at line", i)
            continue
        calculated_length.append(calc_length)

        #width
        if length > maxferet:
            calc_width = width
        else:
            calc_width = minferet

        if calc_width <= 0:
            print("invalid values at line", i)
            continue
        calculated_width.append(calc_width)

        #volume
        #formula of cylindrical decomposition
        calc_vol = (((calc_width**2) * math.pi) / 4) * (calc_length - (calc_width / 3))
        if calc_vol <= 0:
            print("invalid values at line", i)
            continue
        calculated_volume.append(calc_vol)

        #carbon
        calc_carbon = (calc_vol**0.86) * 218
        calculated_carbon.append(calc_carbon)

#combine calculated length, width and volume to one list
combined_calcs = zip(calculated_length, calculated_width,
                     calculated_volume, calculated_carbon)

#write output file
filename = sys.argv[1]
filename2 = filename.replace(".txt", "")
with open(filename2 + ".csv", "wb") as out:
    w = csv.writer(out, delimiter=";", dialect="excel")
    if sys.argv[2] == 'h':
        w.writerow(["Name"] + ["Calculated_Length"] + ["Calculated_Width"] +
                   ["Calculated_Volume"] + ["Calculated_Carbon"] + ["Light"] + ["Medium"] + ["Triplicate"])
        for row in combined_calcs:
          #w.writerow([filename2, row[0], row[1], row[2], row[3]])
          w.writerow([filename2, row[0], row[1], row[2], row[3],
                      filename2[0], filename2[1], filename2[3]])
    elif sys.argv[2] == 'nh':
        for row in combined_calcs:
          w.writerow([filename2, row[0], row[1], row[2], row[3],
                      filename2[0], filename2[1], filename2[3]])
    else:
        print("header option not defined, writing csv file with header...")
        w.writerow(["Name"] + ["Calculated_Length"] + ["Calculated_Width"] +
                   ["Calculated_Volume"] + ["Calculated_Carbon"] + ["Light"] + ["Medium"] + ["Triplicate"])
        for row in combined_calcs:
          w.writerow([filename2, row[0], row[1], row[2], row[3],
                      filename2[0], filename2[1], filename2[3]])


print("Total number of values: " + str(len(calculated_volume)))

print(stats(calculated_length, "calc. length"))
print(stats(calculated_width, "calc. width"))
print(stats(calculated_volume, "calc. volume"))
print(stats(calculated_carbon, "calc. carbon"))

#combine csv on commandline windows: copy *1*.csv strain1.csv
#on bash, create file with header only first: head -1 DH1B.csv > final.csv
#then: for filename in $(ls *1*.csv); do sed 1d $filename >> final.csv; done
