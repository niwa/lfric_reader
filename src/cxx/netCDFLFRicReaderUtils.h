#ifndef netCDFLFRicReaderUtils_h
#define netCDFLFRicReaderUtils_h

#include <vector>

void resolvePeriodicGrid(std::vector<double> & nodeCoordsX,
                         std::vector<double> & nodeCoordsY,
                         std::vector<long long> & faceNodeConnectivity,
                         const size_t numFaces,
                         const size_t numVertsPerFace);

#endif
