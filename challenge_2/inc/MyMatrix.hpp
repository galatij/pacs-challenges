#ifndef _MY_MATRIX_HPP_
#define _MY_MATRIX_HPP_

#include <iostream>
#include <type_traits>
#include "MatrixTraits.hpp"

namespace algebra{
    template<typename T>
    using T_map_t = typename MatrixTraits<T>::map;

    template <typename T, StorageOrder O = StorageOrder::ROWMAJOR>
    class MyMatrix: public MatrixTraits<T>
    {
    public:
    
        MyMatrix(const T_map_t<T> & data) : m_data(data){
            static_assert(std::is_empty_v<T_map_t<T>> == false, "Empty data");
            if constexpr (O == StorageOrder::ROWMAJOR)
            {
                m_nrows = data.crbegin()->first[0] + 1;
                auto max_col_iter = std::max_element(m_data.cbegin(), m_data.cend(),
                            [](const auto& a, const auto& b) {return a.first[1] < b.first[1];}
                );
                std::cout << max_col_iter->second;
                m_ncols = max_col_iter->first[1] + 1;
            }
            else
            {
                m_ncols = data.crbegin()->first[0] + 1;
                auto max_row_iter = std::max_element(m_data.cbegin(), m_data.cend(),
                            [](const auto& a, const auto& b) {return a.first[1] < b.first[1];}
                );
                m_nrows = max_row_iter->first[1] + 1;
            }
        };
        //value_t operator=()(size_t i, size_t j) const;
        

        size_t nrows() const
        {
            return m_nrows;
        }

        size_t ncols() const
        {
            return m_ncols;
        }


    private:
        size_t m_nrows;
        size_t m_ncols;
        T_map_t<T> m_data;

    };
}






#endif // _MY_MATRIX_HPP_