
/*
Input:
n - matrix size
A - symmetric matrix (n x n)
n_iterations - number of Jacobi rotations

Output:
epsilon - off-diagonal norm after the rotations
lambda - approximate eigenvalues
H - matrix of approximate eigenvectors

Constraints:
n >= 2
A must be symmetric
n_iterations >= 0

Example input:
3
4 1 1
1 3 1
1 1 2
3
*/

#include <iostream>
#include <vector>
#include <cmath>

int main() {
    int n, i_max = 0, j_max = 0, n_iterations;
    double max_el, alpha = 0, new_el_b, new_el_a_new, new_el_h, epsilon = 0;

    std::cin >> n;

    std::vector<std::vector<double>> a(n, std::vector<double> (n));
    std::vector<std::vector<double>> a_new(n, std::vector<double> (n));
    std::vector<std::vector<double>> I(n, std::vector<double> (n));
    std::vector<std::vector<double>> u(n, std::vector<double> (n));
    std::vector<std::vector<double>> u_t(n, std::vector<double> (n));
    std::vector<std::vector<double>> b(n, std::vector<double> (n));
    std::vector<std::vector<double>> h(n, std::vector<double> (n));
    std::vector<std::vector<double>> h_new(n, std::vector<double> (n));

    for(int i = 0; i < n; i++) {
        I[i][i] = 1;
    }

    h = I;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
            std::cin >> a[i][j];
        }
    }

    std::cout << "Enter the required number of iterations: ";
    std::cin >> n_iterations;

    for(int i = 0; i < n_iterations; i++) {

        max_el = 0;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i < j && std::abs(a[i][j]) > max_el) {
                    max_el = std::abs(a[i][j]);
                    i_max = i;
                    j_max = j;
                }
            }
        }

        if (max_el == 0) {
                    break;
                }

        alpha = 0.5 * std::atan2(-2 * a[i_max][j_max], a[i_max][i_max] - a[j_max][j_max]);

        u = I;

        u[i_max][i_max] = std::cos(alpha);
        u[j_max][j_max] = std::cos(alpha);
        u[i_max][j_max] = -1 * std::sin(alpha);
        u[j_max][i_max] = std::sin(alpha);

        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                new_el_h = 0;
                for (int r = 0; r < n; r++) {
                    new_el_h += h[i][r] * u[r][j];
                }
                h_new[i][j] = new_el_h;
            }
        }

        h = h_new;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                u_t[j][i] = u[i][j];
            }
        }

        /* 
        B = A * U
        A_new = u_t * B
        */

        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                new_el_b = 0;
                for (int r = 0; r < n; r++) {
                    new_el_b += a[i][r] * u[r][j];
                }
                b[i][j] = new_el_b;
            }
        }

        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                new_el_a_new = 0;
                for (int r = 0; r < n; r++) {
                    new_el_a_new += u_t[i][r] * b[r][j];
                }
                a_new[i][j] = new_el_a_new;
            }
        }
        
        a = a_new;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            epsilon += a[i][j] * a[i][j];
        }
    }

    epsilon = std::sqrt(epsilon);

    std::cout << "epsilon = " << epsilon << "\n";

    for (int i = 0; i < n; i++) {
        std::cout << "lambda " << i + 1 << " = " << a[i][i] << "\n";
    }

    std::cout << "H:\n";

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << h[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    return 0;
}