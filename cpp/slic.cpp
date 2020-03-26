/****************************************************************************
**
** Copyright (C) 2020 Harlen Batagelo <hbatagelo@gmail.com> and
**                    João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program. If not, see <http://www.gnu.org/licenses/>.
**
****************************************************************************/

#include <QtConcurrent>

#include "slic.h"
#include "util.h"

SLIC::SLIC(QImage *image,
           int superpixelSize,
           double compactness,
           int numIterations) :
    m_image(image),
    m_superpixelSize(superpixelSize),
    m_compactness(compactness),
    m_numIterations(numIterations)
{
    m_width = m_image->width();
    m_height = m_image->height();
    m_pixelCount = m_width * m_height;

    ConvertImageToSLICVec();
}

void SLIC::ConvertImageToSLICVec()
{
    m_imageVec.clear();
    m_imageVec.resize(m_pixelCount);

    auto st = reinterpret_cast<QRgb *>(m_image->bits());
    for (int p = 0; p < m_pixelCount; ++p)
        m_imageVec[p] = SLICCenter(QColor(st[p]));

    auto future = QtConcurrent::map(m_imageVec.begin(), m_imageVec.end(),
        [&] (SLICCenter &sc) { ConvertImageToSLICVecKernel(sc); });
    future.waitForFinished();
}

void SLIC::ConvertImageToSLICVecKernel(SLICCenter &sc)
{
    QColor qColor;
    qColor.setRgbF(sc.m_d[0], sc.m_d[1], sc.m_d[2]);

    auto color = Color(qColor);
    sc.m_d[0] = color.getLab_L();
    sc.m_d[1] = color.getLab_a();
    sc.m_d[2] = color.getLab_b();

    sc.m_nd[0] = sc.m_d[0] / 100.0;
    sc.m_nd[1] = (sc.m_d[1] + 100.0) / 200.0;
    sc.m_nd[2] = (sc.m_d[2] + 100.0) / 200.0;
}

void SLIC::ComputeLabels()
{
    DetectLabEdges();
    ComputeHexSeeds();
    ComputeSuperpixels();    
    UpdateConnectivity();
}

void SLIC::SetPixelsLabel(int labelId, const std::vector<QPoint> &pixels)
{
    for (const auto& pt : pixels)
    {
        m_pixelLabels[pt.y() * m_width + pt.x()] = labelId;
    }
}

void SLIC::ComputeHexSeeds()
{
    const int neighbors = 8;
    const int dx8[8] = { -1, -1,  0,  1, 1, 1, 0, -1 };
    const int dy8[8] = {  0, -1, -1, -1, 0, 1, 1,  1 };

    m_step = static_cast<int>(std::sqrt(m_superpixelSize) + 0.5);

    auto xStrips = static_cast<int>(0.5 + double(m_width) / m_step);
    auto yStrips = static_cast<int>(0.5 + double(m_height) / m_step);

    auto xErr = m_width - m_step * xStrips;
    if (xErr < 0)
        xErr = m_width - m_step * --xStrips;

    auto yErr = m_height - m_step * yStrips;
    if (yErr < 0)
        yErr = m_height - m_step * --yStrips;

    auto xErrPerStrip = double(xErr) / xStrips;
    auto yErrPerStrip = double(yErr) / yStrips;

    int xOffset = m_step / 2;
    int yOffset = m_step / 2;

    m_seeds.clear();
    m_seeds.reserve(xStrips * yStrips);

    for (auto y = 0; y < yStrips; ++y)
    {
        auto ye = static_cast<int>(y * yErrPerStrip);
        for (auto x = 0; x < xStrips; ++x)
        {
            auto xe = static_cast<int>(x * xErrPerStrip);
            auto seedx = std::min(x * m_step + (xOffset << (y & 0x1)) + xe, m_width - 1);
            auto seedy = (y * m_step + yOffset + ye);

            SLICCenter sc = m_imageVec[seedy * m_width + seedx];
            sc.m_d[3] = seedx;
            sc.m_d[4] = seedy;

            m_seeds.push_back(sc);
        }
    }

    for (auto &seed : m_seeds)
    {
        auto originalPos = QPoint(seed.m_d[3], seed.m_d[4]);
        auto originalIndex = originalPos.y() * m_width + originalPos.x();

        auto savedIndex = originalIndex;
        for (auto i = 0; i < neighbors; ++i)
        {
            auto newPos = originalPos + QPoint(dx8[i], dy8[i]);

            if (inBounds(newPos, QPoint(0, 0), QPoint(m_width, m_height)))
            {
                auto newIndex = newPos.y() * m_width + newPos.x();
                if (m_labEdges[newIndex] < m_labEdges[savedIndex])
                    savedIndex = newIndex;
            }
        }
        if (savedIndex != originalIndex)
        {
            seed.m_d = {
                m_imageVec[savedIndex].m_d[0],
                m_imageVec[savedIndex].m_d[1],
                m_imageVec[savedIndex].m_d[2],
                double(savedIndex % m_width),
                double(savedIndex / m_width)
            };
        }
    }
}

void SLIC::DetectLabEdges()
{
    m_labEdges.clear();
    m_labEdges.resize(m_width * m_height, -1);

    for (auto i = 1; i < m_height - 1; ++i)
    {
        auto offset = i * m_width;
        for (auto j = 1; j < m_width - 1; ++j)
        {
            auto tmp = offset + j;
            m_labEdges[tmp] = tmp;
        }
    }

    auto future = QtConcurrent::map(m_labEdges.begin(), m_labEdges.end(),
        [&] (double &ind) { DetectLabEdgesKernel(ind); });
    future.waitForFinished();
}

void SLIC::DetectLabEdgesKernel(double &index)
{
    if (index < 0)
        return;

    int i = index;

    auto im1 = i - 1;
    auto ip1 = i + 1;
    auto imw = i - m_width;
    auto ipw = i + m_width;

    auto dx = 0.0, dy = 0.0;
    for (auto j = 0; j < 3; ++j)
    {
        auto vm1 = m_imageVec[im1].m_d[j];
        auto vp1 = m_imageVec[ip1].m_d[j];
        auto vmw = m_imageVec[imw].m_d[j];
        auto vpw = m_imageVec[ipw].m_d[j];
        dx += (vm1 - vp1) * (vm1 - vp1);
        dy += (vmw - vpw) * (vmw - vpw);
    }

    index = dx * dx + dy * dy;
}

void SLIC::SLICThread(int startK, int numK)
{
    auto stepByComp = m_step / m_compactness;
    auto invWeight = 1.0 / (stepByComp * stepByComp);

    for (auto n = startK; n < numK; ++n)
    {
        int y1 = std::max(0.0, m_seeds[n].m_d[4] - m_step);
        int y2 = std::min(double(m_height), m_seeds[n].m_d[4] + m_step);
        int x1 = std::max(0.0, m_seeds[n].m_d[3] - m_step);
        int x2 = std::min(double(m_width),  m_seeds[n].m_d[3] + m_step);

        for (auto y = y1; y < y2; ++y)
        {
            int offset = y * m_width;
            for (auto x = x1; x < x2; ++x)
            {
                auto i = offset + x;

                auto l = m_imageVec[i].m_d[0];
                auto a = m_imageVec[i].m_d[1];
                auto b = m_imageVec[i].m_d[2];

                auto ls = l - m_seeds[n].m_d[0];
                auto as = a - m_seeds[n].m_d[1];
                auto bs = b - m_seeds[n].m_d[2];
                auto xs = x - m_seeds[n].m_d[3];
                auto ys = y - m_seeds[n].m_d[4];

                auto dist = std::sqrt((ls * ls) + (as * as) + (bs * bs)
                                      + (xs * xs) + (ys * ys)
                                      * invWeight);

                if (dist < m_distVec[i])
                {
                    m_distVec[i] = dist;
                    m_pixelLabels[i] = n;
                }
            }
        }
    }
}

void SLIC::ComputeSuperpixels()
{   
    int kMeansSteps = m_seeds.size();

    std::vector<double> edgeSum(kMeansSteps, 0);
    std::vector<double> clusterSize(kMeansSteps, 0);
    std::vector<double> inv(kMeansSteps, 0);

    std::vector<SLICCenter> sigma(kMeansSteps);
    m_distVec = std::vector<double>(m_pixelCount, std::numeric_limits<double>::max());

    auto maxThreads = QThread::idealThreadCount();
    auto centersPerThread = static_cast<int>(std::ceil(kMeansSteps / double(maxThreads)));

    m_pixelLabels.clear();
    m_pixelLabels.resize(m_pixelCount, 0);

    for (auto iteration = 0; iteration < m_numIterations; ++iteration)
    {
        m_distVec.assign(m_pixelCount, std::numeric_limits<qreal>::max());

        auto future = std::make_unique<QFuture<void>[]>(static_cast<size_t>(maxThreads));

        auto centersTotal = kMeansSteps;
        auto startingK = 0;
        auto threadCount = 0;
        while (centersTotal > 0)
        {
            // Distribute batch over the thread pool
            auto centersThisThread = centersPerThread;
            if (centersTotal - centersPerThread < 0)
                centersThisThread = centersTotal;

            centersTotal -= centersPerThread;

            future[threadCount] = QtConcurrent::run(this, &SLIC::SLICThread,
                                                    startingK,
                                                    startingK + centersThisThread);
            startingK += centersThisThread;

            threadCount++;

            if (threadCount == maxThreads || centersTotal <= 0)
            {
                // Wait until this batch finishes
                for (auto i = 0; i < threadCount; ++i)
                    future[i].waitForFinished();

                threadCount = 0;
            }
        }

        // Recalculate the centroid and store in the seed values.
        // Instead of reassigning memory on each iteration, just reset.
        sigma.assign(kMeansSteps, SLICCenter());
        clusterSize.assign(kMeansSteps, 0);
        edgeSum.assign(kMeansSteps, 0);

        auto index = 0;
        for (auto i = 0; i < m_height; ++i)
        {
            for (auto j = 0; j < m_width; ++j)
            {
                auto label = m_pixelLabels[index];

                sigma[label].m_d[0] += m_imageVec[index].m_d[0];
                sigma[label].m_d[1] += m_imageVec[index].m_d[1];
                sigma[label].m_d[2] += m_imageVec[index].m_d[2];
                sigma[label].m_d[3] += j;
                sigma[label].m_d[4] += i;

                edgeSum[label] += m_labEdges[index];

                clusterSize[label] += 1.0;
                index++;
            }
        }

        for (auto k = 0; k < kMeansSteps; ++k)
        {
            if(clusterSize[k] <= 0)
                clusterSize[k] = 1;

            inv[k] = 1.0 / clusterSize[k];
        }

        for (auto k = 0; k < kMeansSteps; ++k)
        {
            for (auto l = 0; l < 5; ++l)
                m_seeds[k].m_d[l] = sigma[k].m_d[l] * inv[k];

            edgeSum[k] *= inv[k];
        }
    }
}

void SLIC::ComputeGraph(const QImage *seedsImage,
                        const QColor fgColor,
                        const QColor bgColor)
{
    const int neighbors = 4;
    const int dx4[neighbors] = { -1,  0,  1,  0 };
    const int dy4[neighbors] = {  0, -1,  0,  1 };

    m_labelNodes.clear();
    m_labelNodes.reserve(m_numLabels);
    for (auto i = 0; i < m_numLabels; ++i)
    {
        m_labelNodes.emplace_back(SLICSuperpixelNode(i));
    }

    for (auto i = 0; i < m_height; ++i)
    {
        for (auto j = 0; j < m_width; ++j)
        {
            auto index = i * m_width + j;
            auto labelId = m_pixelLabels[index];
            SLICSuperpixelNode &node = m_labelNodes[labelId];
            node.AddPixel(QPoint(j, i));

            // Update number of seeded pixels
            QColor pixelColor = seedsImage->pixelColor(j, i);
            if (pixelColor == fgColor)
                node.AddSeededPixels(true);
            else if (pixelColor == bgColor)
                node.AddSeededPixels(false);

            // Find adjacent label
            for (auto n = 0; n < neighbors; ++n)
            {
                QPoint pos = QPoint(j, i) + QPoint(dx4[n], dy4[n]);
                if (inBounds(pos, QPoint(0, 0), QPoint(m_width, m_height)))
                {
                    auto nLabelId = m_pixelLabels[pos.y() * m_width + pos.x()];
                    if (nLabelId != labelId)
                    {
                        // A neighbor has been found. Update both.
                        node.AddNeighbor(nLabelId);
                        m_labelNodes[nLabelId].AddNeighbor(labelId);
                    }
                }
            }          
        }
    }

    m_minVertexValency = std::numeric_limits<int>::max();
    m_maxVertexValency = 0;

    // Update center of mass and valency
    for (auto &labelNode : m_labelNodes)
    {
        labelNode.UpdateCenterOfMass();

        auto numNeighbors = labelNode.GetNeighbors().size();
        if (numNeighbors > m_maxVertexValency)
            m_maxVertexValency = numNeighbors;
        if (numNeighbors < m_minVertexValency)
            m_minVertexValency = numNeighbors;
    }

    // Simple voting scheme for deciding whether a superpixel is seeded
    // as foreground or background
    for (auto &labelNode : m_labelNodes)
    {
        auto fgPixels = labelNode.GetNumberOfForegroundPixels();
        auto bgPixels = labelNode.GetNumberOfBackgroundPixels();

        if (fgPixels + bgPixels > 0)
        {
            // This is a foreground superpixel if the number of internal
            // foreground pixels is greater than the number of background
            // pixels (and vice-versa)
            if (fgPixels > bgPixels)
                labelNode.SetSeedType(1); // Foreground
            else
                labelNode.SetSeedType(-1); // Background
        }
    }
}

void SLIC::ComputeFeatures()
{
    for (auto &labelNode : m_labelNodes)
    {
        labelNode.m_featureVector.setZero(3);

        auto sumR = 0, sumG = 0, sumB = 0;
        for (const auto &pt : labelNode.GetPixels())
        {
            QRgb qRGB = m_image->pixel(pt.x(), pt.y());
            sumR += qRed(qRGB);
            sumG += qGreen(qRGB);
            sumB += qBlue(qRGB);
        }
        auto numPixels = labelNode.GetPixels().size();
        sumR /= numPixels;
        sumG /= numPixels;
        sumB /= numPixels;

        labelNode.m_featureVector(0) = sumR / 255.0;
        labelNode.m_featureVector(1) = sumG / 255.0;
        labelNode.m_featureVector(2) = sumB / 255.0;
    }
}

void SLIC::UpdateConnectivity()
{   
    const int neighbors = 4;
    const int dx4[neighbors] = { -1,  0,  1,  0 };
    const int dy4[neighbors] = {  0, -1,  0,  1 };

    auto numberOfSuperpixels = static_cast<int>(double(m_pixelCount) / double(m_step * m_step));
    auto segmentSizeLimit = m_pixelCount / numberOfSuperpixels / neighbors;

    std::vector<int> newLabels;
    newLabels.resize(m_pixelCount, -1);

    auto newNumLabels = 0;
    auto posVec = std::make_unique<QPoint[]>(static_cast<size_t>(m_pixelCount));
    auto adjacentLabel = 0;
    for (auto pixelIndex = 0; pixelIndex < m_pixelCount; ++pixelIndex)
    {
        if (newLabels[pixelIndex] >= 0)
            continue;

        int i = pixelIndex / m_width;
        int j = pixelIndex % m_width;
        posVec[0] = QPoint(j, i);
        newLabels[pixelIndex] = newNumLabels;

        for (auto n = 0; n < neighbors; ++n)
        {
            auto pos = posVec[0] + QPoint(dx4[n], dy4[n]);

            if (inBounds(pos, QPoint(0, 0), QPoint(m_width, m_height)))
            {
                auto neighborIndex = pos.y() * m_width + pos.x();
                if (newLabels[neighborIndex] >= 0)
                    adjacentLabel = newLabels[neighborIndex];
            }
        }

        auto labelSize = 1;
        for (auto k = 0; k < labelSize; ++k)
        {
            for (auto n = 0; n < neighbors; ++n)
            {
                auto pos = posVec[k] + QPoint(dx4[n], dy4[n]);

                if (inBounds(pos, QPoint(0, 0), QPoint(m_width, m_height)))
                {
                    auto neighborIndex = pos.y() * m_width + pos.x();
                    if (newLabels[neighborIndex] >= 0)
                        continue;

                    if (m_pixelLabels[pixelIndex] == m_pixelLabels[neighborIndex])
                    {
                        posVec[labelSize++] = pos;
                        newLabels[neighborIndex] = newNumLabels;
                    }
                }
            }
        }

        if (labelSize <= segmentSizeLimit)
        {
            for (auto k = 0; k < labelSize; ++k)
            {
                newLabels[posVec[k].y() * m_width + posVec[k].x()] = adjacentLabel;
            }
            newNumLabels--;
        }
        newNumLabels++;
    }

    m_numLabels = newNumLabels;
    m_pixelLabels = newLabels;
}
