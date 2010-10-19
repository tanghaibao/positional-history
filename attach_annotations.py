#!/usr/bin/env python
# -*- coding: UTF-8 -*-


import csv

def main():
    # read in the annotations first
    reader = csv.DictReader(file("TAIR9_functional_descriptions"), delimiter="\t")
    annotations = {}
    for row in reader:
        gene = row["Model_name"]
        gene = gene.rsplit(".", 1)[0]
        annotations[gene] = [row["Type"], row["Short_description"]]

    reader = csv.reader(file("master_list.tab"), delimiter="\t")
    writer = csv.writer(file("master_list-2010-06-22.tab", "w"), delimiter="\t")
    row = reader.next()
    new_row = [row[0]] + ["type", "description"] + row[1:]
    writer.writerow(new_row)
    for row in reader:
        gene = row[0].split("|")[0]
        new_row = [row[0]] + annotations[gene] + row[1:]
        writer.writerow(new_row)


if __name__ == '__main__':
    main()
