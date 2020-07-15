#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:22:21 2019

@author: xander
"""
import matplotlib.pyplot as plt, os
os.path.join(os.path.dirname(__file__))

treated_Files = ['Treated_ChIP.xls','Treated_ATAC_Broad.xls']
wild_Files = ['Wild_ChIP.xls','Wild_ATAC_Broad.xls']

def Get_Profile (chip_File, atac_File):
    bound_Regions = Parse_PeakCaller(chip_File)
    open_Regions = Parse_PeakCaller(atac_File)
    distances = []
    indices = []
    current_Dist = 100
    for i in bound_Regions.keys():
        if (i not in open_Regions.keys()):
            print (i + ' contained in CHiP data but not in ATAC data, the data will be excluded from the analysis.')
        else:    
            for j in range(len(bound_Regions[i])):
                shortest_Dist = 250000000
                bound_Start = int(bound_Regions[i][j][0])
                bound_End = int(bound_Regions[i][j][1])
                for k in range(len(open_Regions[i])):
                    open_Start = int(open_Regions[i][k][0])
                    open_End = int(open_Regions[i][k][1])
                    if (bound_Start >= open_Start) and (bound_Start <= open_End):
                        shortest_Dist = 0
                        break
                    elif (bound_Start < open_Start) and (bound_End > open_Start) and (bound_End < open_End):
                        shortest_Dist = (-1 + (bound_End - open_Start) / (bound_End - bound_Start))
                        break
                    elif (bound_Start > open_Start) and (bound_Start < open_End) and (bound_End > open_End):
                        shortest_Dist = (1 - (open_End - open_Start) / (bound_End - bound_Start))
                        break
                    elif (bound_Start < open_Start) and (bound_End < open_Start):
                        current_Dist = bound_End - open_Start
                    elif (bound_Start > open_End):
                        current_Dist = bound_Start - open_End
                    if (abs(current_Dist) < shortest_Dist):
                        shortest_Dist = current_Dist
                indices.append([i, bound_Start, bound_End])
                distances.append(shortest_Dist)
    return [indices, distances]

def Parse_PeakCaller (peak_File,p_Value = 1.3, q_Value = 1.3):
    dict_Coord = {}
    if 'Broad' in peak_File:
        p_Loc = 5
    else:
        p_Loc = 6
    with open (peak_File) as file:
        text = file.readlines()
        for i in range(0,len(text)):
            line = text[i].strip().split('\t')
            if (len(line) > 1):
                start_Loc = i + 1
                break
        for i in range(start_Loc,len(text)):
            line = text[i].strip().split('\t')
            if (line[0] not in dict_Coord.keys()):
                dict_Coord[line[0]] = []
            if (float(line[p_Loc]) > p_Value) and (float(line[8]) > q_Value):
                dict_Coord[line[0]].append([line[1],line[2]])
        file.close()
    return dict_Coord

def Export_Data (file_List, cell_Status):
    if (cell_Status == 0):
        title = 'Wild_Type'
    else:
        title = 'Treated'
    with open (title + '_Broad_Data.csv', 'w') as file:
        for i in range(len(file_List[0])):
            file.write(file_List[0][i][0] + ',' + str(file_List[0][i][1]) + ',' + str(file_List[0][i][2]) + ',' + str(file_List[1][i]) + ',')
            if (file_List[1][i] <= -1):
                label = 'A'
            elif (file_List[1][i] > -1) and (file_List[1][i] < 0):
                label = 'B'
            elif (file_List[1][i] == 0):
                label = 'C'
            elif (file_List[1][i] > 0) and (file_List[1][i] < 1):
                label = 'D'
            else:
                label= 'E'
            file.write(label + '\n')
        file.close()

def Graph_Data (graph_Data, cell_Status):
    if (cell_Status == 0):
        title = 'Wild_Type'
    else:
        title = 'Treated'
    grouped_Dist = {}
    with open (graph_Data) as file:
        text = file.readlines()
        for line in text:
            line = line.strip().split(',')
            if (line[4] not in grouped_Dist.keys()):
                grouped_Dist[line[4]] = []
            grouped_Dist[line[4]].append(float(line[3]))
        file.close()
    labels = ['A','B','C','D','E']
    y_Data = [0,0,0,0,0]
    for i in grouped_Dist.keys():
        y_Data[labels.index(i)] = len(grouped_Dist[i])
    plt.xlabel('Category')
    plt.ylabel('Counts')
    plt.bar(labels,y_Data)
    plt.savefig(title + '_Broad_Profile.png', )
    plt.close()

wild_Profile = Get_Profile(wild_Files[0], wild_Files[1])
Export_Data(wild_Profile, 0)
Graph_Data ('Wild_Type_Broad_Data.csv', 0)

treat_Profile = Get_Profile(treated_Files[0],treated_Files[1])
Export_Data(treat_Profile, 1)
Graph_Data ('Treated_Broad_Data.csv', 1)