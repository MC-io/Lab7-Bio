#include <iostream>
#include <string>
#include <fstream>
#include <vector>


void create_new_distance_matrix(std::string filename, int start, int end)
{
    std::ifstream file(filename.c_str());
    std::vector<std::vector<float>> matrix(end, std::vector<float>(end, 0)); 
    for (int i = 0; i < end; i++)
    {
        for (int j = 0; j < i; j++)
        {   
            file >> matrix[i][j];
            matrix[j][i] = matrix[i][j];
        }
    }
    file.close();

    std::ofstream out("new_distancia.txt");

    for (int i = start - 1; i < end; i++)
    {
        for (int j = start - 1; j < end; j++)
        {   
            out << matrix[i][j] << '\t';
        }
        out << '\n';
    }

    out.close();
}


int main()
{
    create_new_distance_matrix("distancia.txt", 1, 6);

    return 0;
}