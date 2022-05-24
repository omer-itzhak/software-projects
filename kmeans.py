import sys


def k_mean(k, max_iter, input_filename, output_filename):
    epsilon = 0.001
    data_points = read_file(input_filename)
    if k > len(data_points):
        invalid_input()
    more_than_epsilon = True
    iteration = 0
    centroids = [data_points[i] for i in range(k)]
    while more_than_epsilon and iteration < max_iter:
        iteration += 1
        clusters = [[] for i in range(k)]
        for data_point in data_points:
            cluster_i = 0
            min_euclidean_dist = float('inf')
            for i, centroid in enumerate(centroids):
                euclidean_dist = euclidean_norm(sub_vector(data_point, centroid))
                if euclidean_dist < min_euclidean_dist:
                    min_euclidean_dist = euclidean_dist
                    cluster_i = i
            clusters[cluster_i].append(data_point)
        new_centroids = []
        for cluster in clusters:
            new_centroids.append(avg(cluster))
        change_vector = []
        for i in range(k):
            sub = sub_vector(centroids[i], new_centroids[i])
            norm = euclidean_norm(sub)
            change_vector.append(norm)
        distance = euclidean_norm(change_vector)
        more_than_epsilon = distance > epsilon
        centroids = new_centroids
    f = open(output_filename,"w")
    for centroid in centroids:
        for i, x in enumerate(centroid):
            f.write('%.4f' % x)
            if i != len(centroid) - 1:
                f.write(',')
        f.write('\n')
    f.close()
    return 0


def read_file(input_filename):
    f = open(input_filename)
    assert f != None, "An Error Has Accured!"
    x = f.read()
    temp = x.split('\n')
    temp = temp[:len(temp) - 1:]
    data_points = []
    for vector in temp:
        data_points.append(vector.split(','))
        for i in range(len(data_points[-1])):
            data_points[-1][i] = float(data_points[-1][i])
    f.close()
    return data_points


def avg(cluster):
    n = len(cluster[0])
    centroid = []
    for i in range(n):
        temp = 0
        for vector in cluster:
            temp += vector[i]
        centroid.append(temp / len(cluster))
    return centroid


def sub_vector(old, new):
    return [old[i] - new[i] for i in range(len(old))]


def euclidean_norm(vector):
    res = 0
    for x in vector:
        res += x ** 2
    res = res ** 0.5
    return res


def main():
   return args_parsing()


def args_parsing():
    args = sys.argv
    if len(args) < 4 or len(args) > 5:
        invalid_input()
    k = args[1]
    if not k.isnumeric():
        invalid_input()
    k = int(k)
    if k <= 0:
        invalid_input()
    if len(args) == 4:
        max_iter = 200
        idx_of_file = 2
    else:
        idx_of_file = 3
        max_iter = args[2]
        if not max_iter.isnumeric():
            invalid_input()
        max_iter = int(max_iter)
        if max_iter <= 0:
            invalid_input()
    input_filename = args[idx_of_file]
    output_filename = args[idx_of_file + 1]
    try:
        f = open(output_filename, 'r')
    except OSError:
        # invalid_input()
        pass
    return k_mean(k, max_iter, input_filename, output_filename)


def invalid_input():
    print("Invalid Input!")
    exit(1)


if __name__ == '__main__':
    main()
