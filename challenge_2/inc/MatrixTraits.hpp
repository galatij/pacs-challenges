#ifndef _MATRIX_TRAITS_
#define _MATRIX_TRAITS_
#include <map>
#include <vector>
#include <array>

namespace algebra{
    template <class T>
    struct MatrixTraits{
        using map = std::map<std::array<std::size_t, 2>, T>;
        using idx_t = std::map<std::array<std::size_t, 2>, T>::key_type;
        using value_t = std::map<std::array<std::size_t, 2>, T>::value_type;
        using vec = std::vector<T>;
    };

    enum class StorageOrder
    {
        ROWMAJOR,
        COLMAJOR
    };
}



#endif // _MATRIX_TRAITS_