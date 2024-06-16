import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import random

def calculate_cophenetic_correlation(distance_matrix, cophenetic_matrix):
    n = len(distance_matrix)
    
    dist_vector = []
    cophenetic_vector = []
    
    for i in range(n):
        for j in range(i + 1, n):
            dist_vector.append(distance_matrix[i, j])
            cophenetic_vector.append(cophenetic_matrix[i, j])
    
    mean_dist = sum(dist_vector) / len(dist_vector)
    mean_coph = sum(cophenetic_vector) / len(cophenetic_vector)
    
    covariance = sum((dist_vector[i] - mean_dist) * (cophenetic_vector[i] - mean_coph) for i in range(len(dist_vector))) / len(dist_vector)
    
    std_dist = (sum((dist_vector[i] - mean_dist) ** 2 for i in range(len(dist_vector))) / len(dist_vector)) ** 0.5
    std_coph = (sum((cophenetic_vector[i] - mean_coph) ** 2 for i in range(len(cophenetic_vector))) / len(cophenetic_vector)) ** 0.5
    
    ccc = covariance / (std_dist * std_coph)
    
    return ccc

def get_clusters_at_rank(labels_, distance_matrix_, mode, result_file):
    labels = labels_[:]
    distance_matrix = np.copy(distance_matrix_)
    clusters = [[i] for i in range(len(labels))]


    cophenetic_matrix = np.zeros((len(labels), len(labels)))

    f = open(result_file, "w")
    print(distance_matrix)
    my_method = 'average'
    if mode == "MAX":
        my_method = 'complete'
    elif mode == "MIN":
        my_method = 'single'
    upper_triangular_indices = np.triu_indices_from(distance_matrix, k=1)
    upper_triangular_values = distance_matrix[upper_triangular_indices]

    # Step 4: Condense the matrix into a 1D array
    condensed_array = upper_triangular_values

    Z = linkage(condensed_array, method=my_method)
    original_labels = labels[:]

    f.write("Matriz de Distancias Original:\n")
    f.write("-\t")
    for l in labels:
        f.write("{}\t".format(l))
    f.write('\n')
    for i in range(len(labels)):
        f.write("{}\t".format(labels[i]))
        for j in range(len(labels)):
            f.write("{:.2f}\t".format(distance_matrix[i, j]))
        f.write('\n')
        
    f.write('\n')

    for k in range(len(original_labels) - 1):
        min_dist = sys.float_info.max
        min_row = 0
        min_col = 0
        for i in range(len(labels)):
            for j in range(i):
                if distance_matrix[i, j] < min_dist:
                    min_dist = distance_matrix[i, j]
                    min_row = i;
                    min_col = j;
        print(min_dist)
        print(min_row, " ", min_col)

        new_labels = labels[:]
        new_labels.pop(min_row)
        new_labels[min_col] = new_labels[min_col] + labels[min_row]

        for i in range(len(clusters[min_row])):
            for j in range(len(clusters[min_col])):
                cophenetic_matrix[clusters[min_row][i], clusters[min_col][j]] = min_dist
                cophenetic_matrix[clusters[min_col][j], clusters[min_row][i]] = min_dist

        new_clusters = []
        for cluster in clusters:
            new_clusters.append(cluster[:])

        new_clusters.pop(min_row)
        new_clusters[min_col] = new_clusters[min_col] + clusters[min_row] 

        new_distance_matrix = np.zeros((len(new_labels), len(new_labels)))

        for i in range(len(new_labels)):
            for j in range(i):
                ii = 0
                jj = 0
                if i >= min_row:
                    ii = i + 1
                else:
                    ii = i
                if j >= min_row:
                    jj = j + 1
                else:
                    jj = j

                if i == min_col:
                    if mode == "MIN":
                        new_distance_matrix[i, j] = min(distance_matrix[ii, jj], distance_matrix[min_row, jj])
                    elif mode == "MAX":
                        new_distance_matrix[i, j] = max(distance_matrix[ii, jj], distance_matrix[min_row, jj])
                    else:
                        new_distance_matrix[i, j] = (distance_matrix[ii, jj] + distance_matrix[min_row, jj]) / 2
                elif j == min_col:
                    if mode == "MIN":
                        new_distance_matrix[i, j] = min(distance_matrix[ii, jj], distance_matrix[ii, min_row])
                    elif mode == "MAX":
                        new_distance_matrix[i, j] = max(distance_matrix[ii, jj], distance_matrix[ii, min_row])
                    else:
                        new_distance_matrix[i, j] = (distance_matrix[ii, jj] + distance_matrix[ii, min_row]) / 2
                else:
                    new_distance_matrix[i, j] = distance_matrix[ii, jj]
                new_distance_matrix[j, i] = new_distance_matrix[i, j]

        f.write("Paso {}:\n".format(k + 1))
        f.write("Se unen las secuencias {} + {} -> {} con el valor {:.2f}\n\n".format(labels[min_col], labels[min_row], new_labels[min_col], min_dist))
        labels = new_labels[:]
        clusters = new_clusters[:]
        distance_matrix = np.copy(new_distance_matrix)

        f.write("Matriz de Distancias despues del {} paso:\n".format(k + 1))
        f.write("-\t")
        for l in labels:
            f.write("{}\t".format(l))
        f.write('\n')
        for i in range(len(labels)):
            f.write("{}\t".format(labels[i]))
            for j in range(len(labels)):
                f.write("{:.2f}\t".format(distance_matrix[i, j]))
            f.write('\n')
        f.write('\n')


    for label in labels:
        print(label, end=' ')
    print("")   

    for i in range(len(labels)):
        for j in range(len(labels)):
            print(distance_matrix[i, j], end=' ')
        print("")


    f.write("Matriz Cofenetica:\n")
    f.write("-\t")
    for l in original_labels:
        f.write("{}\t".format(l))
    f.write('\n')
    for i in range(len(cophenetic_matrix)):
        f.write("{}\t".format(original_labels[i]))
        for j in range(len(cophenetic_matrix)):
            f.write("{:.2f}\t".format(cophenetic_matrix[i, j]))
        f.write('\n')
    f.write('\n')


    f.write("CCC: {}\n".format(calculate_cophenetic_correlation(distance_matrix_, cophenetic_matrix)))
    f.close()
    plt.figure(figsize=(8,7))   
    dendrogram(Z, labels=original_labels)
    plt.title('Dendrograma')
    plt.xlabel('Indice del cluster')
    plt.ylabel('Distancias')
    plt.savefig("dendrogram3_{}".format(mode))
    plt.show()
    print(calculate_cophenetic_correlation(distance_matrix_, cophenetic_matrix))

    return calculate_cophenetic_correlation(distance_matrix_, cophenetic_matrix)


def create_random_dist_matrix(size):
    file = open("input.txt", "w")
    matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(i):
            matrix[i ,j] = random.uniform(0.0, 5.0)
            matrix[j, i] = matrix[i, j]


    for i in range(size):
        for j in range(size):
            file.write("{} ".format(matrix[i,j]))
        file.write('\n')
    return matrix


if __name__ == "__main__":
    dist_matrix = np.loadtxt("new_distancia.txt")
    labels = ["A","B","C","D","E",'F']
    ccc1 = get_clusters_at_rank(labels, dist_matrix, 'MIN', "result_min.txt")
    ccc2 = get_clusters_at_rank(labels, dist_matrix, 'MAX', "result_max.txt")
    ccc3 = get_clusters_at_rank(labels, dist_matrix, 'AVG', "result_avg.txt")

    best = None

    if ccc1 >= ccc2 and ccc1 >= ccc3:
        best = "distancia minima"
    elif ccc2 >= ccc1 and ccc2 >= ccc3:
        best = "distancia maxima"
    else:
        best = "distancia promedio" 

    file_min = open("result_min.txt","a")
    file_max = open("result_max.txt","a")
    file_avg = open("result_avg.txt","a")

    file_min.write("La mejor estrategia es de {}.\n".format(best))
    file_max.write("La mejor estrategia es de {}.\n".format(best))
    file_avg.write("La mejor estrategia es de {}.\n".format(best))

