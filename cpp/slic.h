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

#ifndef SLIC_H
#define SLIC_H

#include <QImage>
#include <QPoint>
#include <Eigen/Core>

#include "color.h"

class SLIC;

// A superpixel as a graph node. The node is identified by its
// label id and a list of neighboring superpixel label ids
class SLICSuperpixelNode
{
    friend class SLIC; // Grant private access to SLIC

public:
    SLICSuperpixelNode(const int& id) : m_id(id)
    {
        m_innerPixels.clear();
        m_centerOfMass = QPointF(0, 0);
        ResetSeededPixels();
    }

    // Add a (x, y) pixel coordinates to the list of inner pixels
    void AddPixel(const QPoint& pt)
    {
        m_innerPixels.push_back(pt);
    }

    // Add a neighbor label if not a neighbor already
    bool AddNeighbor(const int& id)
    {
        if (!IsNeighbor(id))
        {
            m_neighbors.push_back(id);
            return true;
        }
        return false;
    }

    // Is id a neighbor label?
    bool IsNeighbor(const int& id) const
    {
        auto res = std::find(std::begin(m_neighbors), std::end(m_neighbors), id);
        return res != std::end(m_neighbors);
    }

    // Compute the center of mass of the inner pixels and reserves
    // memory for the neighbors' segment pixels
    void UpdateCenterOfMass()
    {
        if (m_innerPixels.size() == 0)
            return;

        double sumX = 0, sumY = 0;
        for (const QPoint& pt : m_innerPixels)
        {
            sumX += pt.x();
            sumY += pt.y();
        }
        sumX /= static_cast<double>(m_innerPixels.size());
        sumY /= static_cast<double>(m_innerPixels.size());
        m_centerOfMass = QPointF(sumX, sumY);
    }

    void ResetSeededPixels() { m_seededPixels[0] = m_seededPixels[1] = 0; m_seedType = 0; }
    void AddSeededPixels(bool addToForeground) { m_seededPixels[addToForeground]++; }
    void SetSeedType(int seedType) { m_seedType = seedType; }

    // Getters
    const int& GetId() const { return m_id; }
    const int& GetSeedType() const { return m_seedType; }
    const int& GetNumberOfBackgroundPixels() const { return m_seededPixels[0]; }
    const int& GetNumberOfForegroundPixels() const { return m_seededPixels[1]; }
    const std::vector<int>& GetNeighbors() const { return m_neighbors; }
    const std::vector<QPoint>& GetPixels() const { return m_innerPixels; }
    const QPointF& GetCenterOfMass() const { return m_centerOfMass; }
    const Eigen::VectorXd& GetFeatureVector() const { return m_featureVector; }

private:
    int m_id; // Superpixel label id
    std::vector<int> m_neighbors; // List of neighboring superpixel label ids

    // Number of seeded pixels in this superpixel
    // (used for seeded segmentation)
    int m_seededPixels[2]; // [0] background pixels, [1] foreground pixels

    // Whether this superpixel is a foreground (1), background (-1) or not
    // seeded (0) superpixel
    int m_seedType;   

    // List of (x, y) pixel coordinates of the pixels contained in this
    // superpixel
    std::vector<QPoint> m_innerPixels;

    QPointF m_centerOfMass; // Center of mass with respect to m_innerPixels

    // This is the feature vector. Its size depend on the flags used with
    // SLIC::ComputeSuperpixelFeatures()
    Eigen::VectorXd m_featureVector; // Mean RGB
};

class SLICCenter
{
public:
    std::array<double, 5> m_d; // Color in CIE Lab format + 2D pixel position
    std::array<double, 3> m_nd; // Normalized Lab

    SLICCenter() : m_d{0} { }
    // Workaround constructor to store a QColor instead
    SLICCenter(const QColor &qColor)
    {
        m_d = { qColor.redF(), qColor.greenF(), qColor.blueF() };
    }
    // Constructor for CIE Lab color and xy pixel position
    SLICCenter(const Color &color, const QPoint &point = QPoint(0, 0))
    {
        m_d = { static_cast<double>(color.getLab_L()),
                static_cast<double>(color.getLab_a()),
                static_cast<double>(color.getLab_b()),
                static_cast<double>(point.x()),
                static_cast<double>(point.y()) };

        m_nd[0] = m_d[0] / 100.0;
        m_nd[1] = (m_d[1] + 100.0) / 200.0;
        m_nd[2] = (m_d[2] + 100.0) / 200.0;
    }
};

class SLIC
{
public:
    SLIC(QImage *image,
         int superpixelSize = 100,
         double compactness = 10.0,
         int numIterations = 10);

    void ComputeLabels();
    void ComputeGraph(const QImage *seedsImage,
                      const QColor fgColor = Qt::red,
                      const QColor bgColor = Qt::blue);
    void ComputeFeatures();
    int GetMinVertexValency() const { return m_minVertexValency; }
    int GetMaxVertexValency() const { return m_maxVertexValency; }

    int GetNumLabels() const { return m_numLabels; }
    std::vector<SLICSuperpixelNode>& GetLabelNodes() { return m_labelNodes; }

private:
    QImage *m_image;

    int m_width;
    int m_height;
    int m_pixelCount; // width x height

    int m_superpixelSize;
    double m_compactness;
    int m_numIterations;

    // Number of superpixels
    int m_numLabels;

    // List of superpixel ids
    std::vector<int> m_pixelLabels;

    // Magnitude of gradient of Lab image
    std::vector<double> m_labEdges;

    // List of superpixel nodes (superpixel label graph)
    std::vector<SLICSuperpixelNode> m_labelNodes;

    int m_minVertexValency = 0; // Minimum vertex valency of the superpixel graph
    int m_maxVertexValency = 0; // Maximum vertex valency of the superpixel graph

    std::vector<SLICCenter> m_imageVec;
    std::vector<SLICCenter> m_seeds;
    std::vector<double> m_distVec;
    int m_step;

    void SetPixelsLabel(int labelId, const std::vector<QPoint> &pixels);
    void ComputeHexSeeds();
    void ConvertImageToSLICVec();
    void ConvertImageToSLICVecKernel(SLICCenter &sc);
    void DetectLabEdges();
    void DetectLabEdgesKernel(double &ind);
    void SLICThread(int startK, int numK);
    void ComputeSuperpixels();
    void UpdateConnectivity();
};

#endif // SLIC_H
