#include "MyMatrix.hpp"

using namespace algebra;
template <typename T>
using map = MatrixTraits<T>::map;

int main(int argc, char *argv[]){
    map<double> m{ {{1,3}, 2},
            {{0,3}, 4},
            {{2,5}, 1}                  
    };

    MyMatrix<double> mat(m);

    std::cout << "nrows: " << mat.nrows() << "\nncols: " << mat.ncols();
    std::cout << std::endl;

    return 0;
}