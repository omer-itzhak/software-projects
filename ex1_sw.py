import sys


def k_mean(k, max_iter, input_filename, output_filename):
    if k <= 0 or max_iter <= 0:
        invalid_input()
        return
    try:
        f = open(output_filename, 'w')
    except OSError:
        invalid_input()
        return
    epsilon = 0.001
    data_points = read_file(input_filename)
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
    for centroid in centroids:
        for i, x in enumerate(centroid):
            f.write('%.4f' % x)
            if i != len(centroid) - 1:
                f.write(',')
        f.write('\n')
    f.close()


def read_file(input_filename):
    f = open(input_filename)
    x = f.read()
    temp = x.split('\n')
    temp = temp[:len(temp) - 1:]
    res = []
    for vector in temp:
        res.append(vector.split(','))
        for i in range(len(res[-1])):
            res[-1][i] = float(res[-1][i])
    f.close()
    return res


def avg(cluster):
    res = 0
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
    args_parsing()


def args_parsing():
    max_iter = 200
    args = sys.argv
    if len(args) < 4 or len(args) > 5:
        invalid_input()
        return
    k = args[1]
    if not k.isnumeric():
        invalid_input()
        return
    k = int(k)
    possible_max_iter = args[2]
    if possible_max_iter.isnumeric():
        max_iter = int(max_iter)
        input_filename = args[3]
        output_filename = args[4]
    else:
        input_filename = args[2]
        output_filename = args[3]
    return k_mean(k, max_iter, input_filename, output_filename)


def invalid_input():
    print("Invalid Input!")
    return


if __name__ == '__main__':
    main()
