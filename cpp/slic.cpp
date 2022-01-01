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

#include "slic.h"

#include <QtConcurrent>

#include "util.h"

namespace lc {

SLIC::SLIC(QImage const *image, int superpixelSize, double compactness,
           int numIterations)
    : m_image{image},
      m_width{image->width()},
      m_height{image->height()},
      m_numPixels{m_width * m_height},
      m_superpixelSize{superpixelSize},
      m_compactness{compactness},
      m_numIterations{numIterations} {
  Q_ASSERT(image != nullptr);

  m_imageVec.clear();
  m_imageVec.reserve(m_numPixels);

  auto *st{reinterpret_cast<QRgb const *>(m_image->bits())};
  std::span imageSpan{st, static_cast<std::size_t>(m_numPixels)};
  for (auto const &p : imageSpan) {
    auto const color{QColor{QRgb{p}}};
    m_imageVec.emplace_back(
        0, Eigen::Vector2d{},
        Eigen::Vector3d{color.redF(), color.greenF(), color.blueF()},
        Eigen::Vector3d{});
  }

  auto kernel{[](Center &sc) {
    using Lab_t = std::remove_reference_t<decltype(*sc.m_Lab.data())>;
    RGBToLab(std::span<Lab_t, 3>{sc.m_Lab.data(), 3});

    sc.m_normLab << sc.m_Lab.x() / 100.0, (sc.m_Lab.y() + 100.0) / 200.0,
        (sc.m_Lab.z() + 100.0) / 200.0;
  }};

  auto future{QtConcurrent::map(m_imageVec.begin(), m_imageVec.end(), kernel)};
  future.waitForFinished();
}

void SLIC::ComputeLabels() {
  DetectLabEdges();
  ComputeHexSeeds();
  ComputeSuperpixels();
  UpdateConnectivity();
}

void SLIC::ComputeHexSeeds() {
  std::array const neighborOffsets{
      Eigen::Vector2i{-1, 0}, Eigen::Vector2i{-1, -1}, Eigen::Vector2i{0, -1},
      Eigen::Vector2i{1, -1}, Eigen::Vector2i{1, 0},   Eigen::Vector2i{1, 1},
      Eigen::Vector2i{0, 1},  Eigen::Vector2i{-1, 1}};

  m_step = static_cast<int>(std::round(std::sqrt(m_superpixelSize)));

  auto xStrips{
      static_cast<int>(std::round(static_cast<double>(m_width) / m_step))};
  auto yStrips{
      static_cast<int>(std::round(static_cast<double>(m_height) / m_step))};

  auto xErr{m_width - m_step * xStrips};
  if (xErr < 0) xErr = m_width - m_step * --xStrips;

  auto yErr{m_height - m_step * yStrips};
  if (yErr < 0) yErr = m_height - m_step * --yStrips;

  m_seeds.clear();
  m_seeds.reserve(static_cast<std::size_t>(xStrips) * yStrips);

  auto const xErrPerStrip{static_cast<double>(xErr) / xStrips};
  auto const yErrPerStrip{static_cast<double>(yErr) / yStrips};
  auto const xOffset{m_step / 2};
  auto const yOffset{xOffset};
  auto const wm1{m_width - 1};
  for (auto y{0}; y < yStrips; ++y) {
    auto const ye{static_cast<int>(y * yErrPerStrip)};
    for (auto x{0}; x < xStrips; ++x) {
      auto const xe{static_cast<int>(x * xErrPerStrip)};
      auto const seedx{std::min(
          x * m_step +
              static_cast<int>(static_cast<unsigned int>(xOffset)
                               << (static_cast<unsigned int>(y) & 0x1U)) +
              xe,
          wm1)};
      auto const seedy{y * m_step + yOffset + ye};

      Center sc{m_imageVec[seedy * m_width + seedx]};
      sc.m_index = m_seeds.size();
      sc.m_position = {seedx, seedy};

      m_seeds.push_back(std::move(sc));
    }
  }

  Eigen::Vector2i const imageSize{m_width, m_height};
  for (auto &seed : m_seeds) {
    Eigen::Vector2i const originalPos{seed.m_position.cast<int>()};
    auto const originalIndex{originalPos.y() * m_width + originalPos.x()};

    auto savedIndex{originalIndex};
    for (auto const &offset : neighborOffsets) {
      Eigen::Vector2i newPos{originalPos + offset};

      if (inBounds(newPos, Eigen::Vector2i{0, 0}, imageSize)) {
        auto const newIndex{newPos.y() * m_width + newPos.x()};
        if (m_labEdges[newIndex] < m_labEdges[savedIndex])
          savedIndex = newIndex;
      }
    }
    if (savedIndex != originalIndex) {
      seed.m_Lab = m_imageVec[savedIndex].m_Lab;
      seed.m_position = {savedIndex % m_width, savedIndex / m_width};
    }
  }
}

void SLIC::DetectLabEdges() {
  m_labEdges.clear();
  m_labEdges.resize(static_cast<size_t>(m_numPixels), -1.0);

  auto offset{m_width + 1};
  for (auto i{1}; i < m_height - 1; ++i) {
    for (auto j{1}; j < m_width - 1; ++j) {
      m_labEdges[offset] = offset;
      ++offset;
    }
    offset += 2;
  }

  auto kernel{[this](double &index) {
    if (index < 0.0) return;
    auto const i{static_cast<int>(index)};
    auto const &w{m_width};
    auto const dx{(m_imageVec[i - 1].m_Lab - m_imageVec[i + 1].m_Lab).array()};
    auto const dy{(m_imageVec[i - w].m_Lab - m_imageVec[i + w].m_Lab).array()};
    auto const dx2s{dx.square().sum()};
    auto const dy2s{dy.square().sum()};
    index = dx2s * dx2s + dy2s * dy2s;
  }};

  auto future{QtConcurrent::map(m_labEdges.begin(), m_labEdges.end(), kernel)};
  future.waitForFinished();
}

void SLIC::ComputeSuperpixels() {
  auto const numSeeds{static_cast<int>(m_seeds.size())};
  auto const stepByComp{m_step / m_compactness};
  auto const invWeight{1.0 / (stepByComp * stepByComp)};

  std::vector<double> clusterSize(numSeeds, 0);
  std::vector<double> distVec(m_numPixels, 0);
  std::vector<double> inv(numSeeds, 0);
  std::vector<Center> sigma(numSeeds);

  m_pixelLabels.resize(m_numPixels, 0);

  auto kernel{[&](auto const &seed) {
    auto const &seedPos{seed.m_position};
    auto const x1{std::max(0, static_cast<int>(seedPos.x() - m_step))};
    auto const y1{std::max(0, static_cast<int>(seedPos.y() - m_step))};
    auto const x2{std::min(m_width, static_cast<int>(seedPos.x() + m_step))};
    auto const y2{std::min(m_height, static_cast<int>(seedPos.y() + m_step))};

    for (auto y{y1}; y < y2; ++y) {
      auto const offset{y * m_width};
      for (auto x{x1}; x < x2; ++x) {
        auto const i{offset + x};

        auto const LabDiff{m_imageVec[i].m_Lab - seed.m_Lab};
        auto const LabDist{LabDiff.array().square().sum()};

        Eigen::Vector2d const pos{x, y};
        auto const posDiff{pos - seedPos};
        auto const posDist{posDiff.array().square().sum()};

        auto const dist{LabDist + posDist * invWeight};
        // auto const dist{std::sqrt(LabDist) + std::sqrt(posDist * invWeight)};

        if (dist < distVec[i]) {
          distVec[i] = dist;
          m_pixelLabels[i] = seed.m_index;
        }
      }
    }
  }};

  for (auto iteration{0}; iteration < m_numIterations; ++iteration) {
    distVec.assign(m_numPixels, std::numeric_limits<double>::max());

    auto future{QtConcurrent::map(m_seeds.begin(), m_seeds.end(), kernel)};
    future.waitForFinished();

    sigma.assign(numSeeds, Center{});
    clusterSize.assign(numSeeds, 0);

    auto index{0};
    for (auto i{0}; i < m_height; ++i) {
      for (auto j{0}; j < m_width; ++j) {
        auto const label{m_pixelLabels[index]};
        auto const offset{Eigen::Vector2d(j, i)};

        sigma[label].m_Lab += m_imageVec[index].m_Lab;
        sigma[label].m_position += offset;

        clusterSize[label] += 1.0;
        ++index;
      }
    }

    for (auto k{0}; k < numSeeds; ++k) {
      if (clusterSize[k] <= 0) clusterSize[k] = 1;
      inv[k] = 1.0 / clusterSize[k];
    }

    for (auto k{0}; k < numSeeds; ++k) {
      m_seeds[k].m_Lab = sigma[k].m_Lab * inv[k];
      m_seeds[k].m_position = sigma[k].m_position * inv[k];
    }
  }
}

void SLIC::ComputeGraph(const QImage *seedsImage, QColor foreground,
                        QColor background) {
  std::array const neighborOffsets{
      Eigen::Vector2i{-1, 0}, Eigen::Vector2i{0, -1}, Eigen::Vector2i{1, 0},
      Eigen::Vector2i{0, 1}};

  Eigen::Vector2i const imageSize{m_width, m_height};

  m_labelNodes.clear();
  m_labelNodes.resize(m_numLabels, SLICNode{});

  for (auto i{0}, index{0}; i < m_height; ++i, index += m_width) {
    for (auto j{0}; j < m_width; ++j) {
      auto const label{m_pixelLabels[index + j]};
      SLICNode &node{m_labelNodes[label]};

      // Add current pixel to superpixel node
      QColor const pixelColor{seedsImage->pixelColor(j, i)};
      auto seedType{SeedType::None};
      if (pixelColor == foreground)
        seedType = SeedType::Foreground;
      else if (pixelColor == background)
        seedType = SeedType::Background;
      node.AddPixel(Eigen::Vector2i{j, i}, seedType);

      Eigen::Vector2i const pixel{j, i};
      // Find adjacent label
      for (auto const &offset : neighborOffsets) {
        Eigen::Vector2i pos{pixel + offset};
        if (inBounds(pos, Eigen::Vector2i{0, 0}, imageSize)) {
          auto const neighborLabel{m_pixelLabels[pos.y() * m_width + pos.x()]};
          if (neighborLabel != label) {
            // A neighbor has been found. Update both
            node.AddNeighbor(neighborLabel);
            m_labelNodes[neighborLabel].AddNeighbor(label);
          }
        }
      }
    }
  }

  // Update center of mass and valency
  m_minVertexValency = std::numeric_limits<int>::max();
  m_maxVertexValency = 0;
  for (auto &labelNode : m_labelNodes) {
    labelNode.UpdateCenterOfMass();
    auto const numNeighbors{static_cast<int>(labelNode.GetNeighbors().size())};
    m_maxVertexValency = std::max(m_maxVertexValency, numNeighbors);
    m_minVertexValency = std::min(m_minVertexValency, numNeighbors);
  }

  // Simple voting scheme for deciding whether a superpixel is seeded
  // as foreground or background
  for (auto &labelNode : m_labelNodes) {
    auto const numFgPixels{labelNode.GetNumForegroundPixels()};
    auto const numBgPixels{labelNode.GetNumBackgroundPixels()};

    if (numFgPixels + numBgPixels > 0) {
      // This is a foreground superpixel if the number of internal
      // foreground pixels is greater than the number of background
      // pixels (and vice-versa)
      labelNode.SetSeedType(numFgPixels > numBgPixels ? SeedType::Foreground
                                                      : SeedType::Background);
    }
  }
}

void SLIC::ComputeFeatures() {
  for (auto &labelNode : m_labelNodes) {
    Eigen::Vector3i sum{Eigen::Vector3i::Zero()};
    for (auto const &pt : labelNode.GetInnerPixels()) {
      QRgb const qRGB{m_image->pixel(pt.x(), pt.y())};
      sum += Eigen::Vector3i{qRed(qRGB), qGreen(qRGB), qBlue(qRGB)};
    }
    sum /= static_cast<int>(labelNode.GetInnerPixels().size());
    labelNode.SetFeatureVector(sum.cast<double>() / 255.0);
  }
}

void SLIC::UpdateConnectivity() {
  std::array const neighborOffsets{
      Eigen::Vector2i{-1, 0}, Eigen::Vector2i{0, -1}, Eigen::Vector2i{1, 0},
      Eigen::Vector2i{0, 1}};
  auto const numSuperpixels{
      static_cast<int>(m_numPixels / static_cast<double>(m_step * m_step))};
  auto const segmentSizeLimit{m_numPixels / numSuperpixels /
                              neighborOffsets.size()};
  Eigen::Vector2i const imageSize{m_width, m_height};

  std::vector<int> newLabels;
  newLabels.resize(m_numPixels, -1);

  auto newNumLabels{0};
  auto adjacentLabel{0};
  std::vector<Eigen::Vector2i> posVec(m_numPixels);
  for (auto pixelIndex{0}; pixelIndex < m_numPixels; ++pixelIndex) {
    if (newLabels[pixelIndex] >= 0) continue;

    posVec.front() = {pixelIndex % m_width, pixelIndex / m_width};
    newLabels[pixelIndex] = newNumLabels;

    for (auto const &offset : neighborOffsets) {
      Eigen::Vector2i pos{posVec.front() + offset};

      if (inBounds(pos, Eigen::Vector2i{0, 0}, imageSize)) {
        auto const neighborIndex{pos.y() * m_width + pos.x()};
        if (newLabels[neighborIndex] >= 0)
          adjacentLabel = newLabels[neighborIndex];
      }
    }

    auto labelSize{1};
    for (auto k{0}; k < labelSize; ++k) {
      for (auto const &offset : neighborOffsets) {
        Eigen::Vector2i pos{posVec[k] + offset};

        if (inBounds(pos, Eigen::Vector2i{0, 0}, imageSize)) {
          auto const neighborIndex{pos.y() * m_width + pos.x()};
          if (newLabels[neighborIndex] >= 0) continue;

          if (m_pixelLabels[pixelIndex] == m_pixelLabels[neighborIndex]) {
            posVec[labelSize++] = pos;
            newLabels[neighborIndex] = newNumLabels;
          }
        }
      }
    }

    if (labelSize <= segmentSizeLimit) {
      for (auto const &pos : std::ranges::take_view{posVec, labelSize}) {
        newLabels[pos.y() * m_width + pos.x()] = adjacentLabel;
      }
    } else {
      ++newNumLabels;
    }
  }

  m_numLabels = newNumLabels;
  m_pixelLabels = newLabels;
}

}  // namespace lc
