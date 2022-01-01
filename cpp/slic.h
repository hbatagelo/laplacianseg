/****************************************************************************
**
** Copyright (C) 2020-2022 Harlen Batagelo <hbatagelo@gmail.com>,
**                         João Paulo Gois <jpgois@gmail.com>.
**
** This file is part of the implementation of the paper
** 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
** by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
** Luis Gustavo Nonato. DOI 10.1109/TPAMI.2020.2974475.
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

#include <Eigen/Core>
#include <QColor>
#include <QImage>

#include "util.h"

namespace lc {

enum class SeedType { None = 0, Background = -1, Foreground = 1 };

// Superpixel graph node
class SLICNode {
  friend class SLIC;

  SeedType m_seedType{SeedType::None};
  int m_numBackgroundPixels{};
  int m_numForegroundPixels{};
  std::vector<Eigen::Vector2i> m_pixels; // Pixel coordinates
  std::vector<int> m_neighbors;          // Index of adjacent superpixels
  Eigen::Vector2d m_centerOfMass;        // Center of mass w.r.t. m_pixels
  Eigen::Vector3d m_featureVector{Eigen::Vector3d::Zero()}; // Mean RGB

public:
  // Getters
  [[nodiscard]] SeedType GetSeedType() const { return m_seedType; }
  [[nodiscard]] int GetNumBackgroundPixels() const {
    return m_numBackgroundPixels;
  }
  [[nodiscard]] int GetNumForegroundPixels() const {
    return m_numForegroundPixels;
  }
  [[nodiscard]] const std::vector<Eigen::Vector2i> &GetInnerPixels() const {
    return m_pixels;
  }
  [[nodiscard]] const std::vector<int> &GetNeighbors() const {
    return m_neighbors;
  }
  [[nodiscard]] const Eigen::Vector2d &GetCenterOfMass() const {
    return m_centerOfMass;
  }
  [[nodiscard]] const Eigen::Vector3d &GetFeatureVector() const {
    return m_featureVector;
  }

private:
  // Add coordinates of a pixel to the list of pixels
  void AddPixel(Eigen::Vector2i const &position, SeedType type) {
    m_pixels.push_back(position);
    if (type == SeedType::Background) {
      ++m_numBackgroundPixels;
    } else if (type == SeedType::Foreground) {
      ++m_numForegroundPixels;
    }
  }

  // Add neighbor id if not a neighbor already
  bool AddNeighbor(int id) {
    return std::ranges::find(m_neighbors, id) == m_neighbors.end()
           ? m_neighbors.push_back(id),
           true : false;
  }

  // Compute the center of mass with respect to the pixels
  void UpdateCenterOfMass() {
    if (!m_pixels.empty()) {
      auto sum{std::accumulate(m_pixels.begin(), m_pixels.end(),
                               Eigen::Vector2i(0, 0))};
      m_centerOfMass = sum.cast<double>() / m_pixels.size();
    }
  }

  void SetSeedType(SeedType type) { m_seedType = type; }
  void SetFeatureVector(Eigen::Vector3d const &featureVector) {
    m_featureVector = featureVector;
  }
};

class SLIC {
  QImage const *m_image;
  int const m_width;
  int const m_height;
  int const m_numPixels;
  int const m_superpixelSize;
  double const m_compactness;
  int const m_numIterations;

  int m_numLabels{};                  // Number of superpixels
  std::vector<int> m_pixelLabels;     // Superpixel label for each pixel
  std::vector<double> m_labEdges;     // Magnitude of gradient of Lab image
  std::vector<SLICNode> m_labelNodes; // Superpixel label graph
  int m_minVertexValency{}; // Minimum vertex valency of the superpixel graph
  int m_maxVertexValency{}; // Maximum vertex valency of the superpixel graph

  struct Center {
    int m_index{};
    Eigen::Vector2d m_position{Eigen::Vector2d::Zero()}; // Pixel position
    Eigen::Vector3d m_Lab{Eigen::Vector3d::Zero()};      // Lab color
    Eigen::Vector3d m_normLab{Eigen::Vector3d::Zero()};  // Normalized Lab
  };

  std::vector<Center> m_imageVec;
  std::vector<Center> m_seeds;
  int m_step{};

public:
  explicit SLIC(QImage const *image, int superpixelSize, double compactness,
                int numIterations = 10);

  void ComputeLabels();
  void ComputeGraph(QImage const *seeds, QColor foreground = Qt::red,
                    QColor background = Qt::blue);
  void ComputeFeatures();
  [[nodiscard]] int GetMinVertexValency() const { return m_minVertexValency; }
  [[nodiscard]] int GetMaxVertexValency() const { return m_maxVertexValency; }
  [[nodiscard]] int GetNumLabels() const { return m_numLabels; }
  [[nodiscard]] std::vector<SLICNode> const &GetLabelNodes() const {
    return m_labelNodes;
  }

private:
  void ComputeHexSeeds();
  void DetectLabEdges();
  void ComputeSuperpixels();
  void UpdateConnectivity();
};

} // namespace lc

#endif // SLIC_H
