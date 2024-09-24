/* Forked from https://github.com/emmanuelparadis/ape/src/njs.c */

/* bionjs.c    2014-03-21 */

/* Copyright 2011-2014 Andrei-Alin Popescu */

/* Changes made by Kristian Ullrich 2024-06-10 */

#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include "njs.hpp"

int give_index(
    int i,
    int j,
    int n) {
    if (i>j) {
        return n*(j - 1) - j*(j - 1)/2 + i - j - 1;
    } else {
        return n*(i - 1) - i*(i - 1)/2 + j - i - 1;
    }
    
}

void error(
    const std::string& message,
    int OTU1,
    int OTU2) {
    
    std::cerr << "Error: " << message << " " << OTU1 << " " << OTU2 << std::endl;
    throw std::runtime_error(message);
}

int H(
    double t) {
    
    if (t >= 0 - 1e-10) {
        return 1;
    }
    return 0;
}

double cnxy(
    int x,
    int y,
    int n,
    std::vector<double>& D) {
    
    double nMeanXY = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            if ((i == x && j == y) || (j == x && i == y)) {
                continue;
            }
            double n1 = 0;
            double n2 = 0;
            if (i != x) {
                n1 = D[give_index(i, x, n)]; // Accessing D using give_index function
            }
            if (j != y) {
                n2 = D[give_index(j, y, n)]; // Accessing D using give_index function
            }
            if (n1 == -1 || n2 == -1 || D[give_index(i, j, n)] == -1) {
                continue; // Accessing D using give_index function
            }
            nMeanXY += (n1 + n2 - D[give_index(x, y, n)] - D[give_index(i, j, n)]); // Accessing D using give_index function
        }
    }
    return nMeanXY;
}

int mxy(
    int x,
    int y,
    int n,
    std::vector<double>& D) {

    std::vector<int> mx(n + 1, 0);
    std::vector<int> my(n + 1, 0);
    int xmy=0;
    int ymx=0;

    for (int i = 1; i <= n; i++) {
        if (i != x && D[give_index(x, i, n)] == -1) {
            mx[i] = 1;
        }
        if (i != y && D[give_index(y, i, n)] == -1) {
            my[i] = 1;
        }
    }
    for (int i = 1; i <= n;i++) {
        if (i != x && mx[i] == 1 && my[i] == 0) {
            xmy++;
        }
        if (i != y && my[i] == 1 && mx[i] == 0) {
            ymx++;
        }
    }
    return xmy+ymx;
}

double nxy(
    int x,
    int y,
    int n,
    std::vector<double>& D) {

    int sCXY = 0;
    double nMeanXY = 0;
    
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            if ((i == x && j == y) || (j == x && i == y)) {
                continue;
            }
            double n1 = 0;
            double n2 = 0;
            if (i != x) {
                n1 = D[give_index(i, x, n)];
            }
            if (j != y) {
                n2 = D[give_index(j, y, n)];
            }
            if (n1 == -1 || n2 == -1 || D[give_index(i, j, n)] == -1) {
                continue;
            }
            sCXY++;
            nMeanXY += H(n1 + n2 - D[give_index(x, y, n)] - D[give_index(i, j, n)]);
        }
    }
    if (sCXY == 0) {
        return 0;
    }
    return nMeanXY / sCXY;
}

int cxy(
    int x,
    int y,
    int n,
    std::vector<double>& D) {

    int sCXY = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                continue;
            }
            if ((i == x && j == y) || (j == x && i == y)) {
                continue;
            }
            double n1 = 0;
            double n2 = 0;
            if (i != x) {
                n1 = D[give_index(i, x, n)];
            }
            if (j != y) {
                n2 = D[give_index(j, y, n)];
            }
            if (n1 == -1 || n2 == -1 || D[give_index(i, j, n)] == -1) {
                continue;
            }
            sCXY++;
        }
    }
    return sCXY;
}

void cpp_choosePair(std::vector<double>& D, int n, std::vector<double>& R, std::vector<int>& s, int& sw, int& x, int& y, int fS) {
    int i = 0;
    int j = 0;
    int k = 0;
    int sww = 0;
    std::vector<double> cFS(fS, -1e50);
    std::vector<int> iFS(fS, 0);
    std::vector<int> jFS(fS, 0);

    double max = -1e50;

    for (i = 1; i < n; i++) {
        for (j = i + 1; j <= n; j++) {
            if (D[give_index(i, j, n)] == -1) {
                sww = 1;
                continue;
            }
            if (s[give_index(i, j, n)] <= 2) {
                continue;
            }

            int tr = 0;
            double numb = (R[give_index(i, j, n)] / (s[give_index(i, j, n)] - 2)) - D[give_index(i, j, n)];

            for (k = 0; k < fS && cFS[k] > numb; k++);

            for (tr = fS - 1; tr > k; tr--) {
                cFS[tr] = cFS[tr - 1];
                iFS[tr] = iFS[tr - 1];
                jFS[tr] = jFS[tr - 1];
            }

            if (k < fS) {
                cFS[k] = numb;
                iFS[k] = i;
                jFS[k] = j;
            }
        }
    }
    // No missing distances, just return the one with maximal Q-criterion
    if (sww == 0) {
        x = iFS[0];
        y = jFS[0];
        sw = 0;
        return;
    }
    //calculate N*(x,y)
    for (i = 0; i < fS; i++) {
        if (iFS[i] == 0 || jFS[i] == 0) {
            continue;
        }
        double nb = nxy(iFS[i], jFS[i], n, D);
        if (nb > max) {
            max = nb;
        }
        cFS[i] = nb;
    }
    int dk=0;
    //shift the max N*xy to the front of the array
    for (i = 0; i < fS; i++) {
        if (cFS[i] == max) {
            cFS[dk] = cFS[i];
            iFS[dk] = iFS[i];
            jFS[dk] = jFS[i];
            dk++;
        }
    }
    //if just one pair realises max N*xy, return it
    if (dk == 1) {
        x = iFS[0];
        y = jFS[0];
        return;
    }
    fS = dk;
    max = -1e50;
    //on the front of the array containing max N*xy values compute cxy
    for (i = 0; i < fS; i++) {
        if (iFS[i] == 0 || jFS[i] == 0) {
            continue;
        }
        double nb = cxy(iFS[i], jFS[i], n, D);
        if (nb > max) {
            max = nb;
        }
        cFS[i] = nb;
    }
    //and again shift maximal C*xy values at the fron of the array
    dk = 0;
    for (i = 0; i < fS; i++) {
        if (cFS[i] == max) {
            cFS[dk] = cFS[i];
            iFS[dk] = iFS[i];
            jFS[dk] = jFS[i];
            dk++;
        }
    }
    //if just one C*xy with maximal value, return pair realising it
    if (dk == 1) {
        x = iFS[0];
        y = jFS[0];
        return;
    }
    fS = dk;
    max = -1e50;
    //on the front of the array containing max C*xy compute m*xy
    for (i = 0; i < fS; i++) {
        if (iFS[i] == 0 || jFS[i] == 0) {
            continue;
        }
        double nb = mxy(iFS[i], jFS[i], n, D);
        if (nb > max) {
            max = nb;
        }
        cFS[i] = nb;
    }
    //again shift maximal m*xy values to the fron of the array
    dk = 0;
    for (i = 0; i < fS; i++) {
        if (cFS[i] == max) {
            cFS[dk] = cFS[i];
            iFS[dk] = iFS[i];
            jFS[dk] = jFS[i];
            dk++;
        }
    }
    //if just one maximal value for m*xy return the pair realising it, found at 0
    if (dk == 1) {
        x = iFS[0];
        y = jFS[0];
        return;
    }
    fS = dk;
    //and calculate cnxy on these values, but this time we do not shift, but simply
    //return the pair realising the maximum, stored at iPos
    max = -1e50;
    int iPos = 0;
    for (i = 0; i < fS; i++) {
        if (iFS[i] == 0 || jFS[i] == 0) {
            continue;
        }
        double nb = cnxy(iFS[i], jFS[i], n, D);
        if (nb > max) {
            max = nb;
            iPos = i;
        }
        cFS[i] = nb;
    }
    if (iFS[iPos] == 0 || jFS[iPos] == 0) {
        error("distance information insufficient to construct a tree, cannot calculate agglomeration criterion", iPos, iPos);
    }
    x = iFS[iPos];
    y = jFS[iPos];
}

void cpp_njs(
    std::vector<double>& D,
    int N,
    std::vector<int>& edge1,
    std::vector<int>& edge2,
    std::vector<double>& edge_length,
    int fsS) {
    //double Sdist;
    //double Ndist;
    double A;
    double B;
    double smallest_S;

    //int i;
    //int j;
    //int ij;
    int OTU1;
    int OTU2;
    //int o_l;

    int k = 0;
    int sw = 0;
    int n = N;
    int cur_nod = 2 * n - 2;

    std::vector<double> R(n * (n - 1) / 2, 0);
    std::vector<double> v(n * (n - 1) / 2, 0);
    std::vector<double> new_v(n * (n - 1) / 2, 0);
    std::vector<double> S(n + 1, 0);
    std::vector<double> newR(n * (n - 1) / 2, 0);
    std::vector<double> new_dist(n * (n - 1) / 2, 0);
    std::vector<int> otu_label(n + 1);
    std::vector<int> s(n * (n - 1) / 2, 0);//s contains |Sxy|, which is all we need for agglomeration
    std::vector<int> newS(n * (n - 1) / 2, 0);

    for (int i = 1; i <= n; ++i) {
        otu_label[i] = i; // otu_label[0] is not used
    }

    // Compute Sxy and Rxy
    //std::cout << "Sxy and Rxy:" << std::endl;
    for (int i = 1; i < n; ++i) {
        for (int j = i + 1; j <= n; ++j) {//algorithm assumes i,j /in Sij, so skip pair if it is not known
            if (D[give_index(i, j, n)] == -1) {
                continue;
            }
            //std::cout << "D: " << D[give_index(i, j, n)] << " i: " << i << " j: " << j << " n: " << n << std::endl;
            for (int k = 1; k <= n; ++k) {
                //ij is the pair for which we compute
                //skip k if we do not know the distances between it and i AND j
                if (k == i || k == j) {
                    if (i != k) {
                        R[give_index(i, j, n)] += D[give_index(i, k, n)];
                        //std::cout << "i != k:" << std::endl;
                        //std::cout << "R: " << R[give_index(i, j, n)] << " i: " << i << " j: " << j << " n: " << n << " k: " << k << std::endl;
                    }
                    if (j != k) {
                        R[give_index(i, j, n)] += D[give_index(j, k, n)];
                        //std::cout << "j != k:" << std::endl;
                        //std::cout << "R: " << R[give_index(i, j, n)] << " i: " << i << " j: " << j << " n: " << n << " k: " << k << std::endl;
                    }
                    s[give_index(i, j, n)]++;
                    continue;
                }
                if (D[give_index(i, k, n)] == -1 || D[give_index(j, k, n)] == -1) {
                    continue;
                }
                s[give_index(i, j, n)]++;
                R[give_index(i, j, n)] += D[give_index(i, k, n)];
                R[give_index(i, j, n)] += D[give_index(j, k, n)];
            }
        }
    }

    k = 0;
    sw = 1;

    while (n > 3) {
        int ij = 0;

        for (int i = 1; i < n; ++i) {
            for (int j = i + 1; j <= n; ++j) {
                newR[give_index(i, j, n)] = 0;
                newS[give_index(i, j, n)] = 0;
            }
        }

        smallest_S = -1e50;

        if (sw == 0) {
            for (int i = 1; i <= n; ++i) {
                S[i] = 0;
            }
        }

        B = n - 2;

        if (sw == 1) {
            cpp_choosePair(D, n, R, s, sw, OTU1, OTU2, fsS);
        } else {
            for (int i = 1; i <= n; ++i) {
                for (int j = 1; j <= n; ++j) {
                    if (i == j) {
                        continue;
                    }
                    S[i] += D[give_index(i, j, n)];
                }
            }

            B = n - 2;

            for (int i = 1; i < n; ++i) {
                for (int j = i + 1; j <= n; ++j) {
                    double A = S[i] + S[j] - B * D[give_index(i, j, n)];
                    if (A > smallest_S) {
                        OTU1 = i;
                        OTU2 = j;
                        //std::cout << "OTU1: " << OTU1 << " OTU2: " << OTU2 << std::endl;
                        smallest_S = A;
                    }
                    ij++;
                }
            }
        }
        //update Rxy and Sxy, only if matrix still incomplete
        if (sw == 1) {
            for (int i = 1; i < n; ++i) {
                if (i == OTU1 || i == OTU2) {
                    continue;
                }
                for (int j = i + 1; j <= n; ++j) {
                    if (j == OTU1 || j == OTU2) {
                        continue;
                    }
                    if (D[give_index(i, j, n)] == -1) {
                        continue;
                    }
                    if (D[give_index(i, OTU1, n)] != -1 && D[give_index(j, OTU1, n)] != -1) {
                        R[give_index(i, j, n)] -= (D[give_index(i, OTU1, n)] + D[give_index(j, OTU1, n)]);
                        s[give_index(i, j, n)]--;
                    }
                    if (D[give_index(i, OTU2, n)] != -1 && D[give_index(j, OTU2, n)] != -1) {
                        R[give_index(i, j, n)] -= (D[give_index(i, OTU2, n)] + D[give_index(j, OTU2, n)]);
                        s[give_index(i, j, n)]--;
                    }
                }
            }
        }
        //std::cout << "k: " << k << " OTU1: " << OTU1 << " OTU2: " << OTU2 << " cur_nod: " << cur_nod << std::endl;
        edge2[k] = otu_label[OTU1];
        edge2[k + 1] = otu_label[OTU2];
        edge1[k] = edge1[k + 1] = cur_nod;

        // get the distances between all OTUs but the 2 selected ones
        // and the latter:
        // a) get the sum for both
        // b) compute the distances for the new OTU
        double sum = 0;

        for (int i = 1; i <= n; ++i) {
            if (i == OTU1 || i == OTU2) {
                continue;
            }
            if (D[give_index(OTU1, i, n)] == -1 || D[give_index(OTU2, i, n)] == -1) {
                continue;
            }
            sum += (D[give_index(OTU1, i, n)] - D[give_index(OTU2, i, n)]);
        }
        // although s was updated above, s[otu1,otu2] has remained unchanged
        // so it is safe to use it here
        // if complete distanes, use N-2, else use S
        int down = B;
        if (sw == 1) {
            down = s[give_index(OTU1, OTU2, n)] - 2;
        }
        if (down <= 0) {
            throw std::runtime_error("Distance information insufficient to construct a tree, leaves " + std::to_string(OTU1) + " and " + std::to_string(OTU2) + " isolated from tree");
        }
        sum *= (1.0 / (2 * down));
        double dxy = D[give_index(OTU1, OTU2, n)] / 2;

        //std::cout << "k: " << k << " OTU1: " << OTU1 << " OTU2: " << OTU2 << " cur_nod: " << cur_nod << std::endl;
        edge_length[k] = dxy + sum; // OTU1
        edge_length[k + 1] = dxy - sum; // OTU2
        // no need to change distance matrix update for complete distance
        // case, as pairs will automatically fall in the right category
        A = D[give_index(OTU1, OTU2, n)];
        ij = 0;
        for (int i = 1; i <= n; i++) {
            if (i == OTU1 || i == OTU2) {
                continue;
            }
            if (D[give_index(OTU1, i, n)] != -1 && D[give_index(OTU2, i, n)] != -1) {
                new_dist[ij] = 0.5 * (D[give_index(OTU1, i, n)] - edge_length[k] + D[give_index(OTU2, i, n)] - edge_length[k + 1]);
            } else {
                if (D[give_index(OTU1, i, n)] != -1) {
                    new_dist[ij] = D[give_index(OTU1, i, n)] - edge_length[k];
                } else {
                    if (D[give_index(OTU2, i, n)] != -1) {
                        new_dist[ij] = D[give_index(OTU2, i, n)] - edge_length[k + 1];
                    } else {
                        new_dist[ij] = -1;
                    }
                }
            }
            ij++;
        }
        for (int i = 1; i < n; i++) {
            if (i == OTU1 || i == OTU2) {
                continue;
            }
            for (int j = i + 1; j <= n; j++) {
                if (j == OTU1 || j == OTU2) {
                    continue;
                }
                new_dist[ij] = D[give_index(i, j, n)];
                ij++;
            }
        }
        // compute Rui, only if the distance matrix is still incomplete
        ij = 0;
        if (sw == 1) {
            for (int i = 2; i < n; i++) {
                ij++;
                if (new_dist[give_index(i, 1, n - 1)] == -1) {
                    continue;
                }
                for (int j = 1; j < n; j++) {
                    if (j == 1 || j == i) {
                        if (i != j) {
                            newR[give_index(1, i, n - 1)] += new_dist[give_index(i, j, n - 1)];
                        }
                        if (j != 1) {
                            newR[give_index(1, i, n - 1)] += new_dist[give_index(1, j, n - 1)];
                        }
                        newS[give_index(1, i, n - 1)]++;
                        continue;
                    }
                    if (new_dist[give_index(i, j, n - 1)] != -1 && new_dist[give_index(1, j, n - 1)] != -1) {
                        newS[give_index(1, i, n - 1)]++;
                        newR[give_index(1, i, n - 1)] += new_dist[give_index(i, j, n - 1)];
                        newR[give_index(1, i, n - 1)] += new_dist[give_index(1, j, n - 1)];
                    }
                }
            }          
        }
        // fill in the rest of R and S, again only if the distance matrix still
        // incomplete
        if (sw == 1) {
            for (int i = 1; i < n; i++) {
                if (i == OTU1 || i == OTU2) {
                    continue;
                }
                for (int j = i + 1; j <= n; j++) {
                    if (j == OTU1 || j == OTU2) {
                        continue;
                    }
                    newR[ij] = R[give_index(i, j, n)];
                    newS[ij] = s[give_index(i, j, n)];
                    ij++;
                }
            }          
        }
        // update newR and newS with the new taxa, again only if the distance
        // matrix is still incomplete
        if (sw == 1) {
            for (int i = 2; i < n - 1; i++) {
                if (new_dist[give_index(1, i, n - 1)] == -1) {
                    continue;
                }
                for (int j = i + 1; j <= n - 1; j++) {
                    if (new_dist[give_index(1, j, n - 1)] == -1) {
                        continue;
                    }
                    if (new_dist[give_index(i, j, n - 1)] == -1) {
                        continue;
                    }
                    newR[give_index(i, j, n - 1)] += (new_dist[give_index(1, i, n - 1)] + new_dist[give_index(1, j, n - 1)]);
                    newS[give_index(i, j, n - 1)]++;
                }
            }          
        }
        /* update before the next loop
        (we are sure that OTU1 < OTU2) */
        if (OTU1 != 1) {
            for (int i = OTU1; i > 1; i--) {
                //std::cout << "k: " << k << " OTU1: " << OTU1 << " OTU2: " << OTU2 << " cur_nod: " << cur_nod << std::endl;
                otu_label[i] = otu_label[i - 1];
            }
        }
        if (OTU2 != n) {
            for (int i = OTU2; i < n; i++) {
                //std::cout << "k: " << k << " OTU1: " << OTU1 << " OTU2: " << OTU2 << " cur_nod: " << cur_nod << std::endl;
                otu_label[i] = otu_label[i + 1];
            }
        }
        otu_label[1] = cur_nod;

        n--;
        for (int i = 0; i < n * (n - 1) / 2; i++) {
            D[i] = new_dist[i];
            if (sw == 1) {
                R[i] = newR[i];
                s[i] = newS[i];
            }
        }
        //std::cout << "n: " << n << " k: " << k << " OTU1: " << OTU1 << " OTU2: " << OTU2 << " cur_nod: " << cur_nod << std::endl;
        cur_nod--;
        k = k + 2;
        //std::cout << "n: " << n << " k: " << k << " OTU1: " << OTU1 << " OTU2: " << OTU2 << " cur_nod: " << cur_nod << std::endl;
    }
    int dK = 0; // number of known distances in the final distance matrix
    int iUK = -1; // index of unknown distance, if we have one missing distance
    int iK = -1; // index of the only known distance, only needed if dK==1
    for (int i = 0; i < 3; i++) {
        edge1[N * 2 - 4 - i] = cur_nod;
        edge2[N * 2 - 4 - i] = otu_label[i + 1];
        if (D[i] != -1) {
            dK++;
            iK = i;
        } else {
            iUK = i;
        }
    }
    if (dK == 2) {
        // if two distances are known: assume our leaves are x,y,z, d(x,z) unknown
        // and edge weights of three edges are a,b,c, then any b,c>0 that
        // satisfy c-b=d(y,z)-d(x,y) a+c=d(y,z) are good edge weights, but for
        // simplicity we assume a=c if d(yz)<d(xy) a=b otherwise, and after some
        // algebra we get that we can set the missing distance equal to the
        // maximum of the already present distances
        double max = -1e50;
        for (int i = 0; i < 3; i++) {
            if (i == iUK) {
                continue;
            }
            if (D[i] > max) {
                max = D[i];
            }
        }
        D[iUK] = max;
    }
    if (dK == 1) {
        // through similar motivation as above, if we have just one known distance
        // we set the other two distances equal to it
        for (int i = 0; i < 3; i++) {
            if (i == iK) {
                continue;
            }
            D[i] = D[iK];
        }
    }
    if (dK == 0) {
        // no distances are known, we just set them to 1
        for (int i = 0; i < 3; i++) {
            D[i] = 1;
        }
    }
    edge_length[N * 2 - 4] = (D[0] + D[1] - D[2]) / 2;
    edge_length[N * 2 - 5] = (D[0] + D[2] - D[1]) / 2;
    edge_length[N * 2 - 6] = (D[2] + D[1] - D[0]) / 2;
}
