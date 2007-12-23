/*
 * File:        vector_template_instantiation.cc
 * Copyright:   (c) 2005-2008 Kevin T. Chu and Masa Prodanovic
 * Revision:    $Revision: 1.1 $
 * Modified:    $Date: 2006/02/16 14:52:08 $
 * Description: Explicit template instantiation for LSMLIB 
 */

#include <vector>

template class std::vector<bool>;
template class std::vector<char>;
template class std::vector<int>;
template class std::vector<float>;
template class std::vector<double>;
