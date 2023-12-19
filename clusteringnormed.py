#!/bin/python3
import time
import re
import numpy as np
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
start_time = time.time()
def transform_string(input_string):
    # Replace 1/1 with 2
    transformed_string = re.sub(r'0/0', 'a', input_string)
    transformed_string = re.sub(r'0/1', 'b', transformed_string)
    transformed_string = re.sub(r'1/0', 'b', transformed_string)
    transformed_string = re.sub(r'1/1', 'B', transformed_string)
    transformed_string = re.sub(r'\./\.:.', 'xoor', transformed_string)
    transformed_string = re.sub(r'[0-9,:\s./]', '', transformed_string)
    pattern = re.compile(r'\s')
    result = re.sub(pattern, '', transformed_string)
    return result
mapping = {'a': 0, 'b': 1, 'B': 2 ,'x': 0 ,}
result_array = np.array([])
result_list=[]
all_list=[]
h=0
counter=0
with open ("firstmil.vcf","r") as vcf:
    for line in vcf:
        if line.startswith('#'):
            h+=1
            continue  
        #pattern = re.compile(r'^(\S+)\s+(\S+).+GT:.+\s+(.+)$')
        pattern = re.compile(r'^(\S+)\s+(\d+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+GT:PL\s(.+)$')
        match = pattern.match(line)
        if match:
            chromosome = match.group(1)
            position = match.group(2)
            rest_of_line = match.group(3)
            string = transform_string(rest_of_line)
            #if "x" in string:
            #    continue
            #counter+=1
            #print(counter)
            result_list = [mapping[char] for char in string]
            if not len(result_list)==476:
                print("error")
                print(len(result_list))
                print(string)
            else:
                all_list.append(result_list)
    #print(all_list)
    numpy_array = np.array(all_list)
    print("Transforming array")
    transposed_array = np.transpose(numpy_array)
    print("plotting")
    # Perform hierarchical clustering
    linkage_matrix = linkage(transposed_array, method='ward')

    # Truncate dendrogram
    plt.figure(figsize=(60, 40))  
    dendrogram(linkage_matrix, orientation='right',show_leaf_counts=True)
    plt.title('Clustering Dendrogram')
    plt.xlabel('Distance')
    plt.ylabel('Data Points')
    plt.savefig('dendrogram_norm_x.pdf')
    plt.clf()
    print("dendrogram done")
    distance_threshold = 10
    clusters = fcluster(linkage_matrix, distance_threshold, criterion='distance')
    print("second plot")
    # Scatter plot of important data points and save the plot
    plt.figure(figsize=(10, 8))
    plt.scatter(range(len(clusters)), clusters, c=clusters, cmap='viridis', marker='o', s=50)
    plt.title('Scatter Plot of Important Data Points')
    plt.xlabel('Data Points')
    plt.ylabel('Cluster Label')
    plt.savefig('scatter_plot_norm_x_x.pdf')  # Save the plot
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Time taken: {elapsed_time} seconds")
print("end")
