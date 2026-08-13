/*
 * A2AVectorsHigh.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#include <Hadrons/Modules/MSolver/A2AVectorsHigh.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MSolver;

template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 48>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 64>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 96>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 128>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 144>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 192>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<FIMPL, 256>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 48>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 64>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 96>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 128>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 144>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 192>;
template class HADRONS_NAMESPACE::MSolver::TA2AVectorsHigh<ZFIMPL, 256>;
