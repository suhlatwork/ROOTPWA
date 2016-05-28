///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2009 Sebastian Neubert
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
#ifdef __CINT__


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;


#pragma link C++ class rpwa::fitResult+;
#pragma link C++ class rpwa::complexMatrix-;


#pragma read                                                                                        \
        sourceClass="TCMatrix"                                                                      \
        source="TMatrixD _re; TMatrixD _im"                                                         \
        version="[1-]"                                                                              \
        targetClass="rpwa::complexMatrix"                                                           \
        target="_size1, _size2, _nmbDataElements, _data"                                            \
        include="TMatrixD.h"                                                                        \
        code="{                                                                                     \
{                                                                                                   \
    _size1 = onfile._re.GetNrows();                                                                 \
    _size2 = onfile._re.GetNcols();                                                                 \
    _nmbDataElements = _size1 * _size2;                                                             \
    _data = new std::complex<double>[_nmbDataElements];                                             \
    unsigned int dataIndex = 0;                                                                     \
    for (unsigned int row = 0; row < _size1; ++row) {                                               \
        for (unsigned int col = 0; col < _size2; ++col) {                                           \
            _data[dataIndex++] = std::complex<double>(onfile._re(row, col), onfile._im(row, col));  \
        }                                                                                           \
    }                                                                                               \
    newObj->readMatrix();                                                                           \
}                                                                                                   \
              }";
#pragma read                                                                                        \
        sourceClass="rpwa::complexMatrix"                                                           \
        source="TMatrixD _re; TMatrixD _im"                                                         \
        version="[1]"                                                                               \
        targetClass="rpwa::complexMatrix"                                                           \
        target="_size1, _size2, _nmbDataElements, _data"                                            \
        include="TMatrixD.h"                                                                        \
        code="{                                                                                     \
{                                                                                                   \
    _size1 = onfile._re.GetNrows();                                                                 \
    _size2 = onfile._re.GetNcols();                                                                 \
    _nmbDataElements = _size1 * _size2;                                                             \
    _data = new std::complex<double>[_nmbDataElements];                                             \
    unsigned int dataIndex = 0;                                                                     \
    for (unsigned int row = 0; row < _size1; ++row)                                                 \
        for (unsigned int col = 0; col < _size2; ++col)                                             \
            _data[dataIndex++] = std::complex<double>(onfile._re(row, col), onfile._im(row, col));  \
    newObj->readMatrix();                                                                           \
}                                                                                                   \
              }";
#pragma read                                                                                        \
        sourceClass="rpwa::fitResult"                                                               \
        source="std::vector<std::string> _prodAmpNames; std::vector<std::string> _waveNames"        \
        version="[-6]"                                                                              \
        targetClass="rpwa::fitResult"                                                               \
        target="_prodAmpRanks, _prodAmpWaveIndices, _waveProdAmpIndices, _waveRefls"                \
        include="TMatrixD.h, partialWaveFitHelper.h"                                                \
        code="{                                                                                     \
{                                                                                                   \
    _prodAmpRanks.resize(onfile._prodAmpNames.size());                                              \
    _prodAmpWaveIndices.resize(onfile._prodAmpNames.size());                                        \
    _waveProdAmpIndices.assign(onfile._waveNames.size(), std::vector<unsigned int>());              \
    _waveRefls.resize(onfile._waveNames.size());                                                    \
    std::vector<std::string> waveNames;                                                             \
    for (size_t i = 0; i < onfile._prodAmpNames.size(); ++i) {                                      \
        const std::string& prodAmpName = onfile._prodAmpNames[i];                                   \
        if (prodAmpName.length() == 0 or prodAmpName[0] != 'V'                                      \
            or prodAmpName.find('_') == std::string::npos) {                                        \
            printErr << \"production amplitude name '\" << prodAmpName                              \
                     << \"' does not follow the naming convention. \"                               \
                     << \"cannot deduce corresponding wave name. Aborting...\" << std::endl;        \
            throw;                                                                                  \
        }                                                                                           \
        const std::string waveName = prodAmpName.substr(prodAmpName.find('_') + 1);                 \
        const int         rank     = atoi(&prodAmpName.c_str()[1]);                                 \
        const int         refl     = rpwa::partialWaveFitHelper::getReflectivity(waveName);         \
        const size_t      waveIdx  = std::find(waveNames.begin(), waveNames.end(), waveName)        \
                                     - waveNames.begin();                                           \
        if (waveIdx == waveNames.size())                                                            \
            waveNames.push_back(waveName);                                                          \
        if (waveName != onfile._waveNames[waveIdx]) {                                               \
            printErr << \"order of wave names differs from order extracted by production \"         \
                     << \"amplitude names. Aborting...\" << std::endl;                              \
            throw;                                                                                  \
        }                                                                                           \
        _prodAmpRanks[i]       = rank;                                                              \
        _prodAmpWaveIndices[i] = waveIdx;                                                           \
        _waveProdAmpIndices[waveIdx].push_back(i);                                                  \
        _waveRefls[waveIdx]    = refl;                                                              \
    }                                                                                               \
    if (waveNames.size() != onfile._waveNames.size()) {                                             \
        printErr << \"number of wave names differs from read fit result. Aborting...\" << std::endl;\
        throw;                                                                                      \
    }                                                                                               \
}                                                                                                   \
              }";


#endif
