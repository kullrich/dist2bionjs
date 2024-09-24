#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <getopt.h>
#include "njs.hpp"
#include "bionjs.hpp"

// Structure to hold the distance matrix
struct DistMatrix {
    std::vector<std::vector<double>> matrix;
    std::vector<std::string> labels;
};

// Structure to hold the phylogenetic tree
struct Phylo {
    std::vector<std::pair<int, int>> edge;
    std::vector<double> edge_length;
    std::vector<std::string> tip_label;
    int Nnode;
};

Phylo bionjs(const std::vector<std::vector<double>>& X, int fs = 15) {
    if (fs < 1) {
        throw std::invalid_argument("argument 'fs' must be a non-zero positive integer");
    }
        

    int N = X.size();
    if (N < 3) {
        throw std::invalid_argument("cannot build a BIONJ* tree with less than 3 observations");
    }

    std::vector<double> distMatrix((N*(N-1))/2);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double value = X[i][j];
            if (std::isnan(value) || value < 0) value = -1;
            //std::cout << "i: " << i+1 << " j: " << j+1 << " value: " << value << " index: " << give_index(i+1, j+1, N) << ";" << std::endl;
            distMatrix[give_index(i+1, j+1, N)] = value;
        }
    }

    std::vector<int> edge1(2 * N - 3);
    std::vector<int> edge2(2 * N - 3);
    std::vector<double> edge_length(2 * N - 3);
    
    cpp_bionjs(distMatrix, N, edge1, edge2, edge_length, fs);

    std::vector<std::string> labels(N);
    for (int i = 0; i < N; ++i) {
        labels[i] = std::to_string(i + 1);
    }

    Phylo obj;
    for (int i = 0; i < 2 * N - 3; ++i) {
        obj.edge.emplace_back(edge1[i], edge2[i]);
    }
    obj.edge_length = edge_length;
    obj.tip_label = labels;
    obj.Nnode = N - 2;

    // Reorder function is skipped as it is specific to how the ordering is needed
    return obj;
}

Phylo njs(const std::vector<std::vector<double>>& X, int fs = 15) {
    if (fs < 1) {
        throw std::invalid_argument("argument 'fs' must be a non-zero positive integer");
    }
        

    int N = X.size();
    if (N < 3) {
        throw std::invalid_argument("cannot build a BIONJ* tree with less than 3 observations");
    }

    std::vector<double> distMatrix((N*(N-1))/2);
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double value = X[i][j];
            if (std::isnan(value) || value < 0) value = -1;
            //std::cout << "i: " << i+1 << " j: " << j+1 << " value: " << value << " index: " << give_index(i+1, j+1, N) << ";" << std::endl;
            distMatrix[give_index(i+1, j+1, N)] = value;
        }
    }

    std::vector<int> edge1(2 * N - 3);
    std::vector<int> edge2(2 * N - 3);
    std::vector<double> edge_length(2 * N - 3);
    
    cpp_njs(distMatrix, N, edge1, edge2, edge_length, fs);

    std::vector<std::string> labels(N);
    for (int i = 0; i < N; ++i) {
        labels[i] = std::to_string(i + 1);
    }

    Phylo obj;
    for (int i = 0; i < 2 * N - 3; ++i) {
        obj.edge.emplace_back(edge1[i], edge2[i]);
    }
    obj.edge_length = edge_length;
    obj.tip_label = labels;
    obj.Nnode = N - 2;

    // Reorder function is skipped as it is specific to how the ordering is needed
    return obj;
}

DistMatrix read_dist_matrix(std::istream& input, int skip_lines = 0) {
    std::string line;
    for (int i = 0; i < skip_lines; ++i) {
        std::getline(input, line); // Skip the specified number of lines
    }

    std::vector<std::string> labels;
    std::vector<std::vector<double>> matrix;

    while (std::getline(input, line)) {
        std::stringstream ss(line);
        std::string item;

        std::getline(ss, item, ','); // First item is the label
        labels.push_back(item);

        std::vector<double> row;
        while (std::getline(ss, item, ',')) {
            row.push_back(std::stod(item));
        }
        matrix.push_back(row);
    }
    // Create and return the DistMatrix struct
    DistMatrix distMatrix;
    distMatrix.matrix = matrix;
    distMatrix.labels = labels;

    return distMatrix;
}

DistMatrix read_dist_matrix_from_stdin(int skip_lines = 0) {
    return read_dist_matrix(std::cin, skip_lines);
}

DistMatrix read_dist_matrix_from_file(const std::string& filename, int skip_lines = 0) {
    std::ifstream input(filename);
    if (input.is_open()) {
        return read_dist_matrix(input, skip_lines);
    } else {
        std::cerr << "Error: unable to open file " << filename << std::endl;
        return DistMatrix(); // Return an empty DistMatrix
    }
}

int main(int argc, char** argv) {
    int opt;
    int skip_lines = 0;
    std::string dist_file;
    std::string mode;

    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"mode", required_argument, 0, 'm'},
        {"skip", required_argument, 0, 's'},
        {0, 0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "i:m:s:", long_options, NULL)) != -1) {
        switch (opt) {
        case 'i':
            dist_file = optarg;
            break;
        case 'm':
            mode = optarg;
            break;
        case 's':
            skip_lines = std::stoi(optarg);
            break;
        default:
            std::cerr << "Invalid option: " << opt << std::endl;
            return 1;
        }
    }

    // Check if dist file
    if (dist_file.empty()) {
        std::cerr << "Error: dist file not specified" << std::endl;
        return 1;
    }

    // Check if mode set
    if (mode.empty()) {
        std::cerr << "Error: mode not specified" << std::endl;
        return 1;
    }

    // Reading the distance matrix from standard input
    DistMatrix X;
    if (dist_file == "-") {
        X = read_dist_matrix_from_stdin(skip_lines);
    } else {
        try {
            X = read_dist_matrix_from_file(dist_file, skip_lines);
        } catch (const std::exception& e) {
            std::cerr << "Error reading dist file: " << e.what() << std::endl;
            return 1;
        }
    }

    try {
        Phylo tree;
        if (mode == "njs") {
            tree = njs(X.matrix, 15);
        } else if (mode == "bionjs") {
            tree = njs(X.matrix, 15);
        } else {
            std::cerr << "Error: invalid mode" << std::endl;
            return 1;
        }
        std::cout << "Tree constructed successfully!" << std::endl;
        // Output tree information
        for (size_t i = 0; i < tree.edge.size(); ++i) {
            std::cout << "Edge: " << tree.edge[i].first << " - " << tree.edge[i].second
                      << ", Length: " << tree.edge_length[i] << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
